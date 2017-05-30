#ifndef __TYPETRUNCATEDQUADRATIC2D_H__
#define __TYPETRUNCATEDQUADRATIC2D_H__

#include <string.h>
#include <assert.h>

// #include <stdio.h>

// This file is a stub, to make TRWS compile.

struct TypeTruncatedQuadratic2D {
    typedef double REAL;

    struct Edge {

//        void DistanceTransformL2(int /*K*/, int /*stride*/, REAL /*alpha*/, REAL* /*source*/, REAL* /*dest*/,
//                int* /*parabolas*/, int* /*intersections*/) {
/*            printf("\n\
+-------------------------------------------------------------+\n\
|   In order to run TRW-S with truncted L2 terms,             |\n\
|   you need to download the implementation from              |\n\
|      http://pub.ist.ac.at/~vnk/papers/TRW-S.html            |\n\
|   and copy file  typeTruncatedQuadratic2D.h                 |\n\
|   to the main directory (i.e. replace the existing file)    |\n\
+-------------------------------------------------------------+\n\
			");
 */
//            exit(1);
//        }

        void DistanceTransformL2(int K, int stride, REAL alpha, REAL* source, REAL* dest, int* parabolas, int* intersections) {
            assert(alpha >= 0);

            if (alpha == 0) {
                REAL* ptr;
                REAL vMin = source[0];

                for (ptr = source + stride; ptr < source + K * stride; ptr += stride) {
                    if (vMin > ptr[0]) {
                        vMin = ptr[0];
                    }
                }
                for (ptr = dest; ptr < dest + K * stride; ptr += stride) {
                    ptr[0] = vMin;
                }
                return;
            }

            int i, j, k, p;
            int r = 0; // = number of parabolas minus 1
            // parabolas[p] will be base of parabola p (0<=p<=r)
            // intersections[p] will be intersection between parabolas p-1 and p (1<=p<=r)
            // intersections[0] will be always 0
            parabolas[0] = 0;
            intersections[0] = 0;

            for (i = 1; i < K; i++) {
                while (1) {
                    j = parabolas[r]; // base of previous rightmost visible parabola
                    // k is intersection of parabolas i and j
                    k = (int) (1 + ((i + j) + (source[stride * j] - source[stride * i]) / (alpha * (j - i))) / 2);
                    if (k >= K) {
                        // i is not visible
                        break;
                    }
                    if (k < 0) {
                        k = 0;
                    }

                    if (k > intersections[r]) {
                        // intersection is rightmost, add it to end
                        r++;
                        parabolas[r] = i;
                        intersections[r] = k;
                        break;
                    }
                    // j is not visible
                    if (r == 0) {
                        parabolas[0] = i;
                        break;
                    }
                    r--;
                }
            }

            intersections[r + 1] = K;

            i = 0;
            for (p = 0; p <= r; p++) {
                j = parabolas[p];
                // i values in [intersections[p], intersections[p+1]) are assigned to j
                for (; i < intersections[p + 1]; i++) {
                    dest[stride * i] = source[stride * j] + alpha * (i - j)*(i - j);
                }
            }
        }

    };
};

#endif
