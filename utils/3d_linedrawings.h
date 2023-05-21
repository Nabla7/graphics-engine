//
// Created by Pim Van den Bosch on 2023/03/24.
//

#ifndef ENGINE_3D_LINEDRAWINGS_H
#define ENGINE_3D_LINEDRAWINGS_H

using namespace std;

#include "vector/vector3d.h"
#include "ini_configuration.h"
#include "real_lines.h"
#include <list>


class Face{
public:
    vector<int> point_indices;
};

class Figure3D{
public:
    vector<Vector3D> points;
    vector<Face> faces;
    img::Color lineColor;
    double angleX;
    double angleY;
    double angleZ;
    double scale;
    vector<double> center;
    vector<double> bgColor;
    int n;
    int m;
    double r;
    double R;
    double height;
};

typedef vector<Figure3D> Figures3D;

img::EasyImage linedrawer3D(ini::Configuration &configuration);
img::EasyImage linedrawer3DWithZBuffer(ini::Configuration &configuration);
img::EasyImage linedrawer3DWithZBufferTriangles(ini::Configuration &configuration);

Figure3D triangulateFigure(Figure3D& figure);

#endif //ENGINE_3D_LINEDRAWINGS_H
