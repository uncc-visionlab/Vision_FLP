/*
Copyright (C) 2006 Pedro Felzenszwalb

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
 */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
//#include <conio.h>
#include <algorithm>
#include <assert.h>
#include <image.h>
#include <misc.h>
#include <pnmfile.h>
#include <filter.h>
#include <imconv.h>

#include "depth.h"

#define ITER 5       // number of BP iterations at each scale
#define LEVELS 5     // number of scales

#define DISC_K 2.0F         // truncation of discontinuity cost
#define DATA_K 0.5F        // truncation of data cost
#define LAMBDA 0.07F        // weighting of data cost

#define INF 1E20     // large cost
#define VALUES 3    // number of possible disparities
#define SCALE 16     // scaling from disparity to graylevel in output

#define SIGMA 0.7    // amount to smooth the input images

static void dt(float f[VALUES]) {
    for (int q = 1; q < VALUES; q++) {
        float prev = f[q - 1] + 1.0F;
        if (prev < f[q])
            f[q] = prev;
    }
    for (int q = VALUES - 2; q >= 0; q--) {
        float prev = f[q + 1] + 1.0F;
        if (prev < f[q])
            f[q] = prev;
    }
}

// compute message

void msg(float s1[VALUES], float s2[VALUES],
        float s3[VALUES], float s4[VALUES],
        float dst[VALUES]) {
    float val;

    // aggregate and find min
    float minimum = INF;
    for (int value = 0; value < VALUES; value++) {
        dst[value] = s1[value] + s2[value] + s3[value] + s4[value];
        if (dst[value] < minimum)
            minimum = dst[value];
    }

    // dt
    dt(dst);

    // truncate 
    minimum += DISC_K;
    for (int value = 0; value < VALUES; value++)
        if (minimum < dst[value])
            dst[value] = minimum;

    // normalize
    val = 0;
    for (int value = 0; value < VALUES; value++)
        val += dst[value];

    val /= VALUES;
    for (int value = 0; value < VALUES; value++)
        dst[value] -= val;
}

// computation of data costs

image<float[VALUES]> *comp_data(image<float> *imgX, image<float> *imgY,
        image<float> *imgZ, std::vector<Plane3f> planeCoeffs) {
    image<float[VALUES]> *data = new image<float[VALUES]>(imgX->width(), imgX->height());
    int height = imgX->height();
    int width = imgX->width();
    float *xVals, *yVals, *zVals;

    for (int y = 0; y < height; y++) {
        xVals = imPtr(imgX, 0, y);
        yVals = imPtr(imgY, 0, y);
        zVals = imPtr(imgZ, 0, y);
        for (int x = 0; x < width; x++) {
            for (int value = 0; value < VALUES; value++) {
                Point3f p(xVals[x], yVals[x], zVals[x]);
                //Point3f p(imRef(imgX, x, y), imRef(imgY, x, y), imRef(imgZ, x, y));
                float val = abs(planeCoeffs[value].evaluate(p)); // sq intensity diff
                //std::cout << "point = " << p << std::endl;
                //std::cout << "plane = " << planeCoeffs[value] << std::endl;
                //std::cout << "error = "  << val << std::endl;
                imRef(data, x, y)[value] = LAMBDA * std::min(val, DATA_K);
            }
        }
    }

    return data;
}

// generate output from current messages

image<uchar> *output(image<float[VALUES]> *u, image<float[VALUES]> *d,
        image<float[VALUES]> *l, image<float[VALUES]> *r,
        image<float[VALUES]> *data) {
    int width = data->width();
    int height = data->height();
    image<uchar> *out = new image<uchar>(width, height);

    for (int y = 1; y < height - 1; y++) {
        for (int x = 1; x < width - 1; x++) {
            // keep track of best value for current pixel
            int best = 0;
            float best_val = INF;
            for (int value = 0; value < VALUES; value++) {
                float val =
                        imRef(u, x, y + 1)[value] +
                        imRef(d, x, y - 1)[value] +
                        imRef(l, x + 1, y)[value] +
                        imRef(r, x - 1, y)[value] +
                        imRef(data, x, y)[value];
                if (val < best_val) {
                    best_val = val;
                    best = value;
                }
            }
            imRef(out, x, y) = best * SCALE;
        }
    }

    return out;
}

// belief propagation using checkerboard update scheme

void bp_cb(image<float[VALUES]> *u, image<float[VALUES]> *d,
        image<float[VALUES]> *l, image<float[VALUES]> *r,
        image<float[VALUES]> *data,
        int iter) {
    int width = data->width();
    int height = data->height();

    for (int t = 0; t < ITER; t++) {
        //std::cout << "iter " << t << "\n";

        for (int y = 1; y < height - 1; y++) {
            for (int x = ((y + t) % 2) + 1; x < width - 1; x += 2) {

                msg(imRef(u, x, y + 1), imRef(l, x + 1, y), imRef(r, x - 1, y),
                        imRef(data, x, y), imRef(u, x, y));

                msg(imRef(d, x, y - 1), imRef(l, x + 1, y), imRef(r, x - 1, y),
                        imRef(data, x, y), imRef(d, x, y));

                msg(imRef(u, x, y + 1), imRef(d, x, y - 1), imRef(r, x - 1, y),
                        imRef(data, x, y), imRef(r, x, y));

                msg(imRef(u, x, y + 1), imRef(d, x, y - 1), imRef(l, x + 1, y),
                        imRef(data, x, y), imRef(l, x, y));

            }
        }
    }
}

// multiscale belief propagation for image restoration

image<uchar> *plane_ms(image<float> *imgX, image<float> *imgY,
        image<float> *imgZ, std::vector<Plane3f> planecoeffs) {
    image<float[VALUES]> *u[LEVELS];
    image<float[VALUES]> *d[LEVELS];
    image<float[VALUES]> *l[LEVELS];
    image<float[VALUES]> *r[LEVELS];
    image<float[VALUES]> *data[LEVELS];

    // data costs
    data[0] = comp_data(imgX, imgY, imgZ, planecoeffs);

    // data pyramid
    for (int i = 1; i < LEVELS; i++) {
        int old_width = data[i - 1]->width();
        int old_height = data[i - 1]->height();
        int new_width = (int) ceil(old_width / 2.0);
        int new_height = (int) ceil(old_height / 2.0);

        assert(new_width >= 1);
        assert(new_height >= 1);

        data[i] = new image<float[VALUES]>(new_width, new_height);
        for (int y = 0; y < old_height; y++) {
            for (int x = 0; x < old_width; x++) {
                for (int value = 0; value < VALUES; value++) {
                    imRef(data[i], x / 2, y / 2)[value] += imRef(data[i - 1], x, y)[value];
                }
            }
        }
    }

    // run bp from coarse to fine
    for (int i = LEVELS - 1; i >= 0; i--) {
        int width = data[i]->width();
        int height = data[i]->height();

        // allocate & init memory for messages
        if (i == LEVELS - 1) {
            // in the coarsest level messages are initialized to zero
            u[i] = new image<float[VALUES]>(width, height);
            d[i] = new image<float[VALUES]>(width, height);
            l[i] = new image<float[VALUES]>(width, height);
            r[i] = new image<float[VALUES]>(width, height);
        } else {
            // initialize messages from values of previous level
            u[i] = new image<float[VALUES]>(width, height, false);
            d[i] = new image<float[VALUES]>(width, height, false);
            l[i] = new image<float[VALUES]>(width, height, false);
            r[i] = new image<float[VALUES]>(width, height, false);

            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    for (int value = 0; value < VALUES; value++) {
                        imRef(u[i], x, y)[value] = imRef(u[i + 1], x / 2, y / 2)[value];
                        imRef(d[i], x, y)[value] = imRef(d[i + 1], x / 2, y / 2)[value];
                        imRef(l[i], x, y)[value] = imRef(l[i + 1], x / 2, y / 2)[value];
                        imRef(r[i], x, y)[value] = imRef(r[i + 1], x / 2, y / 2)[value];
                    }
                }
            }
            // delete old messages and data
            delete u[i + 1];
            delete d[i + 1];
            delete l[i + 1];
            delete r[i + 1];
            delete data[i + 1];
        }

        // BP
        bp_cb(u[i], d[i], l[i], r[i], data[i], ITER);
    }

    image<uchar> *out = output(u[0], d[0], l[0], r[0], data[0]);

    delete u[0];
    delete d[0];
    delete l[0];
    delete r[0];
    delete data[0];

    return out;
}

class ASCIIMat {
public:
    float size[2];
    float *data;

    ASCIIMat() {
        size[0] = 0;
        size[1] = 0;
    }

    ASCIIMat(const ASCIIMat& other) {
        size[0] = other.size[0];
        size[1] = other.size[1];
        data = new float[ (int) (size[0] * size[1])];
        memcpy(data, other.data, (int) (size[0] * size[1]) * sizeof (float));
    }

    ASCIIMat& operator=(const ASCIIMat& other) {
        if (this != &other) {
            size[0] = other.size[0];
            size[1] = other.size[1];
            data = new float[ (int) (size[0] * size[1])];
            memcpy(data, other.data, (int) (size[0] * size[1]) * sizeof (float));
        }
        return *this;
    }

    ~ASCIIMat() {
        delete [] data;
    }
};

static image<float> *loadDepthImageData(char *name, ASCIIMat &planeCoeffMat,
        ASCIIMat &kMat) {
    std::string sizeline, dataline;
    int height, width, col;
    /* read header */
    std::vector<ASCIIMat> matVec;
    std::ifstream pFile(name);
    if (pFile.is_open()) {
        while (!pFile.eof()) {
            ASCIIMat mat;
            std::getline(pFile, sizeline);
            std::stringstream s1(sizeline);
            col = 0;
            while (s1 >> mat.size[col])
                col++;
            std::cout << col << std::endl;
            mat.data = new float[ (int) (mat.size[0] * mat.size[1])];
            std::getline(pFile, dataline);
            std::stringstream s2(dataline);
            col = 0;
            while (s2 >> mat.data[col])
                col++;
            std::cout << col << std::endl;
            matVec.push_back(mat);
        }
        pFile.close();
    } else {
        std::cout << "Unable to open file";
    }
    image<float> *im;
    if (matVec.size() == 3) {
        kMat = matVec[0];
        height = kMat.size[0];
        width = kMat.size[1];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                std::cout << kMat.data[y * width + x] << " ";
            }
            std::cout << "\n";
        }
        planeCoeffMat = matVec[1];
        height = planeCoeffMat.size[0];
        width = planeCoeffMat.size[1];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                std::cout << planeCoeffMat.data[y * width + x] << " ";
            }
            std::cout << "\n";
        }
        height = matVec[2].size[0];
        width = matVec[2].size[1];
        im = new image<float>(width, height);
        for (int y = 0; y < height; y++) {
            //float *rowPtr = imPtr(im, 0, y);
            for (int x = 0; x < width; x++) {
                //std::cout << "(" << x << ", " << y << ")" << std::endl;
                //rowPtr[x] = matVec[2].data[y * width + x];
                //imRef(im, x, y) = matVec[2].data[x * height + y];
                imRef(im, x, y) = matVec[2].data[y * width + x];
            }
        }
    }
    return im;
}

int main(int argc, char **argv) {
    image<float> *imgX, *imgY, *imgZ;
    image<uchar> *out, *edges, *sc_imgZ;

    if (argc != 3) {
        std::cerr << "usage: " << argv[0] << " depth(dat) out(pgm)\n";
        exit(1);
    }

    // load input
    ASCIIMat kMat;
    ASCIIMat planeCoeffMat;
    imgZ = loadDepthImageData(argv[1], planeCoeffMat, kMat);

    imgX = new image<float>(imgZ->width(), imgZ->height());
    imgY = new image<float>(imgZ->width(), imgZ->height());
    sc_imgZ = new image<uchar>(imgZ->width(), imgZ->height());
    float cx = kMat.data[0 * 2 + 2], cy = kMat.data[1 * 3 + 2];
    float ifx = 1.0f / kMat.data[0 * 3 + 0], ify = 1.0f / kMat.data[1 * 3 + 1];
    for (int y = 0; y < imgZ->height(); y++) {
        for (int x = 0; x < imgZ->width(); x++) {
            float Z = imRef(imgZ, x, y);
            imRef(imgX, x, y) = Z * (x - cx) * ifx;
            imRef(imgY, x, y) = Z * (y - cy) * ify;
            if (Z > 0 && Z < 6) {
                imRef(sc_imgZ, x, y) = (uchar) (Z * 255.0f / 6.0f);
            } else {
                imRef(sc_imgZ, x, y) = 0;
            }
        }
    }

    std::vector<Plane3f> planes;
    for (int idx = 0; idx < planeCoeffMat.size[0]; idx++) {
        Plane3f newplane(planeCoeffMat.data[idx * 4 + 0], planeCoeffMat.data[idx * 4 + 1],
                planeCoeffMat.data[idx * 4 + 2], planeCoeffMat.data[idx * 4 + 3]);
        planes.push_back(newplane);
    }

    for (int i = 0; i < 5; i++) {
        clock_t begin = clock();
        // compute segmentation
        out = plane_ms(imgX, imgY, imgZ, planes);
        clock_t end = clock();
        double elapsed_secs = 1000.0 * double(end - begin) / CLOCKS_PER_SEC;
        std::cout << "Elapsed time = " << elapsed_secs << " ms." << std::endl;
    }
    // save output
    image<uchar> *sc_imgOut = new image<uchar>(out->width(), out->height());
    for (int y = 0; y < out->height(); y++) {
        for (int x = 0; x < out->width(); x++) {
            imRef(sc_imgOut, x, y) = (255 / VALUES) * imRef(out, x, y);
        }
    }

    //savePGM(sc_imgZ, argv[2]);
    savePGM(sc_imgOut, argv[2]);

    delete imgX;
    delete imgY;
    delete imgZ;
    delete sc_imgZ;
    delete out;
    delete sc_imgOut;
    return 0;
}
