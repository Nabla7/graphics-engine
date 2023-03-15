//
// Created by Pim Van den Bosch on 2023/02/24.
//

#include "easy_image.h"
#include <vector>

using namespace std;

img::EasyImage ColorRectangle(int width, int height)
{
    img::EasyImage image(width, height);
    for(unsigned int i = 0; i < width; i++)
    {
        for(unsigned int j = 0; j < height; j++)
        {
            image(i,j).red = i%256;
            image(i,j).green = j%256;
            image(i,j).blue = (i+j)%256;
        }
    }
    return image;
}

img::EasyImage Blocks(int width,
                      int height,
                      int nrXBlocks,
                      int nrYBlocks,
                      vector<double> colorWhite,
                      vector<double> colorBlack,
                      bool invertColors)
{
    img::EasyImage image(width, height);
    int w_b = width/nrXBlocks ;
    int h_b = height/nrYBlocks;

    if (invertColors)
    {
        vector<double> cw = colorWhite;
        vector<double> cb = colorBlack;

        colorWhite = cb;
        colorBlack = cw;
    }

    for(unsigned int i = 0; i < width; i++)
    {
        int b_x = i/w_b;
        for(unsigned int j = 0; j < height; j++)
        {
            int b_y = j/h_b;
            int b_xy = b_x + b_y;
            if (b_xy % 2 == 0)
            {
                image(i,j).red = colorWhite[0]*255;
                image(i,j).green = colorWhite[1]*255;
                image(i,j).blue = colorWhite[2]*255;
            }
            else
            {
                image(i,j).red = colorBlack[0]*255;
                image(i,j).green = colorBlack[1]*255;
                image(i,j).blue = colorBlack[2]*255;
            }
        }

    }
    return image;
}

img::EasyImage QuarterCircle(int width,
                             int height,
                             int nrLines,
                             vector<double> backoundcolor,
                             vector<double> linecolor)
{
    img::EasyImage image(width, height);
    int h_s = height/(nrLines-1);
    int w_s = width/(nrLines-1);

    img::Color color = img::Color (255,0,0);

    int x_0 = 0;
    int y_1 = height-1;

    for (int i = 0; i <  nrLines-1; ++i)
    {
        int y_0 = i * h_s;
        int x_1 = i * w_s;
        image.draw_line(x_0, y_0, x_1, y_1, color);
    }
    return image;
}