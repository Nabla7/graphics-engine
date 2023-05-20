//
// Created by Pim Van den Bosch on 2023/02/24.
//


#ifndef ENGINE_REAL_LINES_H
#define ENGINE_REAL_LINES_H

#include "easy_image.h"
#include <list>


using namespace img;
using namespace std;

class Point2D {
public:
    double x;
    double y;
    double z;
};



class Line2D {
public:
    Point2D p1;
    Point2D p2;
    Color color;

    Line2D(Point2D P1, Point2D P2, Color linecolor) {
        p1 = P1;
        p2 = P2;
        color = linecolor;
    }
};

using Lines2D = std::list<Line2D>;

Lines2D L_System2D (const string &inputfile,
                    const Color color);

img::EasyImage draw2DLines(Lines2D &lines,
                           const int &size,
                           Color &backgroundcolor);

double interpolateZ(const Point2D& p1, const Point2D& p2, double t);

img::EasyImage draw2DLinesWithZBuffer(Lines2D &lines,
                                      const int &size,
                                      Color &backgroundcolor,
                                      vector<vector<double>> &zBuffer);

#endif //ENGINE_REAL_LINES_H
