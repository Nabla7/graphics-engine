//
// Created by Pim Van den Bosch on 2023/02/24.
//


#ifndef ENGINE_REAL_LINES_H
#define ENGINE_REAL_LINES_H

#include "easy_image.h"
#include "3d_linedrawings.h"
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

class Triangle2D {
public:
    Point2D p1;
    Point2D p2;
    Point2D p3;
    Color color;
    size_t figureIndex;
    size_t faceIndex;

    Triangle2D(Point2D P1, Point2D P2, Point2D P3, Color triangleColor) {
        p1 = P1;
        p2 = P2;
        p3 = P3;
        color = triangleColor;
        size_t figureIndex;
        size_t faceIndex;
    }
};

typedef vector<Triangle2D> Triangles2D;


double max(std::initializer_list<double> values);
double min(std::initializer_list<double> values);
std::tuple<double, double, double, double> calculateMinMax(Lines2D& lines);
std::pair<double, double> calculateImageDimensions(double current_x_max, double current_x_min, double current_y_max, double current_y_min, double size);
void scaleAndShiftPoints(Lines2D& lines, double current_x_min, double current_y_min, double x_range, double y_range, double image_x, double image_y);

Lines2D L_System2D (const string &inputfile,
                    const Color color);

img::EasyImage draw2DLines(Lines2D &lines,
                           const int &size,
                           Color &backgroundcolor);

double interpolateZ(const Point2D& p1, const Point2D& p2, double t);

img::EasyImage draw2DLinesWithZBuffer(Lines2D &lines,
                                      const int &size,
                                      Color &backgroundcolor);

class Figure3D;
typedef vector<Figure3D> Figures3D;

img::EasyImage draw2DTrianglesWithZBuffer(Triangles2D &triangles,
                                          Figures3D TriangulatedFigures3D,
                                          const int &size,
                                          Color &backgroundcolor);

#endif //ENGINE_REAL_LINES_H
