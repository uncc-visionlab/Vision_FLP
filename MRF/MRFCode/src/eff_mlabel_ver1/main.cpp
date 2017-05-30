
/***  Please see README.txt  ***/

#include <cstdio>
#include "image.h"

int main() {
    image* Image;

    //   b. Create an image object using a ppm file:
    //      image* Image = image("file.ppm", number_of_labels);
    //      pgm files are also supported. Change image::image(char*, int) accordingly.
    int num_labels = 16;
    // Example
    Image = new image("images/tsukuba_l.ppm", num_labels);

    //   c. Solve the energy function
    //Image -> kovtun(); // Compute the partially optimal solution only
    //Image -> kovtun(true); // Compute the partially optimal solution and 
    // solve the remaining nodes using alpha expansion.
    Image -> trw(); // Compute standard TRW/BP solutions.
    //Image -> trw(true); // Compute the partially optimal solution and 
    // solve the remaining nodes using TRW/BP.

    //   d. Write your own method to save the result, eg. image::save_disparity();
    Image ->save_disparity("output.ppm");

    //################################################################################

}

//int main-old()
//{
//  image* Image;
//
//  int num_labels = 4; // number of labels
//
//  // NOTE: Write code to initialize the energy parameters in image::image(char*, int) 
//  Image = new image("file.ppm", num_labels);
//
//  Image -> kovtun(true);
//  //Image -> kovtun();
//  //Image -> trw();
//  //Image -> trw(true);
//
//  return 0;
//}
