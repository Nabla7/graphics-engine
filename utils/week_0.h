//
// Created by Pim Van den Bosch on 2023/02/24.
//

#ifndef ENGINE_WEEK_0_H
#define ENGINE_WEEK_0_H

#include "easy_image.h"
#include <vector>

using namespace std;

img::EasyImage ColorRectangle(int width, int height);

img::EasyImage Blocks(int width,
                      int height,
                      int nrXBlocks,
                      int nrYBlocks,
                      vector<double> colorWhite,
                      vector<double> colorBlack,
                      bool invertColors);

img::EasyImage QuarterCircle(int width,
                             int height,
                             int nrLines,
                             vector<double> backoundcolor,
                             vector<double> linecolor);

#endif //ENGINE_WEEK_0_H




