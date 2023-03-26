#include "vector/vector3d.h"
#include "3d_linedrawings.h"
#include "easy_image.h"
#include "real_lines.h"
#include "ini_configuration.h"



Matrix scaleFigure(const double scale) {
    Matrix scaleMatrix;
    scaleMatrix(1, 1) = scale;
    scaleMatrix(2, 2) = scale;
    scaleMatrix(3, 3) = scale;
    return scaleMatrix;
}


Matrix rotateX(const double angle) {
    Matrix rotateXMatrix;
    double cosAngle = cos(angle);
    double sinAngle = sin(angle);
    rotateXMatrix(2, 2) = cosAngle;
    rotateXMatrix(2, 3) = -sinAngle;
    rotateXMatrix(3, 2) = sinAngle;
    rotateXMatrix(3, 3) = cosAngle;
    return rotateXMatrix;
}

Matrix rotateY(const double angle) {
    Matrix rotateYMatrix;
    double cosAngle = cos(angle);
    double sinAngle = sin(angle);
    rotateYMatrix(1, 1) = cosAngle;
    rotateYMatrix(1, 3) = sinAngle;
    rotateYMatrix(3, 1) = -sinAngle;
    rotateYMatrix(3, 3) = cosAngle;
    return rotateYMatrix;
}

Matrix rotateZ(const double angle) {
    Matrix rotateZMatrix;
    double cosAngle = cos(angle);
    double sinAngle = sin(angle);
    rotateZMatrix(1, 1) = cosAngle;
    rotateZMatrix(1, 2) = -sinAngle;
    rotateZMatrix(2, 1) = sinAngle;
    rotateZMatrix(2, 2) = cosAngle;
    return rotateZMatrix;
}

Matrix translate(const vector<double> &vector) {
    Matrix translationMatrix;
    translationMatrix(1, 4) = vector[0];
    translationMatrix(2, 4) = vector[1];
    translationMatrix(3, 4) = vector[2];
    return translationMatrix;
}

void toPolar(Vector3D v, double &r, double &theta, double &phi){
    r = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    theta = atan2(v.y, v.x);
    phi = acos(v.z/r);
}

Matrix eyeMatrix(Vector3D eye){
    double r, theta, phi;
    toPolar(eye, r, theta, phi);

    // Initialize 4x4 identity matrix
    Matrix transformMatrix;
    // Construct eye point transformation matrix
    transformMatrix(1, 1) = -sin(theta);
    transformMatrix(1, 2) = -cos(theta) * cos(phi);
    transformMatrix(1, 3) = cos(theta) * sin(phi);
    transformMatrix(2, 1) = cos(theta);
    transformMatrix(2, 2) = -sin(theta) * cos(phi);
    transformMatrix(2, 3) = sin(theta) * sin(phi);
    transformMatrix(3, 2) = sin(phi);
    transformMatrix(3, 3) = cos(phi);
    transformMatrix(4, 3) = -r;

    return transformMatrix;
}


Figure3D applyTransform(const Figure3D &figure, const Matrix &transformation) {
    Figure3D transformedFigure = figure;
    for (Vector3D &point : transformedFigure.points) {
        point = point*transformation;
    }
    return transformedFigure;
}

Vector3D project(Vector3D point, const Vector3D &eye){
    Matrix projectionMatrix = eyeMatrix(eye); // Assuming eyeMatrix() creates the transformation matrix
    Vector3D projectedPoint = point*projectionMatrix;
    return projectedPoint;
};

#include <iostream>

Lines2D projectFigure(const Figure3D &figure, const Vector3D &eye, double d, Color lineColor) {
    const double epsilon = 1e-6;
    const double maxCoord = 1e6; // Define a maximum allowed coordinate value
    Lines2D lines2D;

    for (const auto& edge : figure.faces) {
        Vector3D startPoint = figure.points[edge.point_indices.front()];
        Vector3D endPoint = figure.points[edge.point_indices.back()];

        startPoint = project(startPoint, eye);
        startPoint.x = startPoint.x * d / (startPoint.z + epsilon);
        startPoint.y = startPoint.y * d / (startPoint.z + epsilon);

        endPoint = project(endPoint, eye);
        endPoint.x = endPoint.x * d / (endPoint.z + epsilon);
        endPoint.y = endPoint.y * d / (endPoint.z + epsilon);

        Point2D startPoint2D = {startPoint.x, startPoint.y};
        Point2D endPoint2D = {endPoint.x, endPoint.y};

        if (std::abs(startPoint2D.x) < maxCoord && std::abs(startPoint2D.y) < maxCoord &&
            std::abs(endPoint2D.x) < maxCoord && std::abs(endPoint2D.y) < maxCoord) {
            lines2D.push_back(Line2D(startPoint2D, endPoint2D, lineColor));
        }
    }

    return lines2D;
}



Figures3D parseiniFigures(ini::Configuration &configuration) {
    Figures3D figures;

    int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
    for (int i = 0; i < nrFigures; i++) {
        std::string figureSection = "Figure" + std::to_string(i);
        Figure3D figure;

        figure.angleX = configuration[figureSection]["rotateX"].as_double_or_die();
        figure.angleY = configuration[figureSection]["rotateY"].as_double_or_die();
        figure.angleZ = configuration[figureSection]["rotateZ"].as_double_or_die();
        figure.scale = configuration[figureSection]["scale"].as_double_or_die();
        figure.center = configuration[figureSection]["center"].as_double_tuple_or_die();
        // Parse the color
        auto colorTuple = configuration[figureSection]["color"].as_double_tuple_or_die();
        figure.lineColor.red = colorTuple[0]*255;
        figure.lineColor.green = colorTuple[1]*255;
        figure.lineColor.blue = colorTuple[2]*255;

        // Parse the points
        int nrPoints = configuration[figureSection]["nrPoints"].as_int_or_die();
        for (int j = 0; j < nrPoints; j++) {
            std::string pointKey = "point" + std::to_string(j);
            auto pointTuple = configuration[figureSection][pointKey].as_double_tuple_or_die();
            Vector3D point = Vector3D::point(pointTuple[0], pointTuple[1], pointTuple[2]);
            figure.points.push_back(point);
        }

        // Parse the lines (faces)
        int nrLines = configuration[figureSection]["nrLines"].as_int_or_die();
        for (int j = 0; j < nrLines; j++) {
            std::string lineKey = "line" + std::to_string(j);
            auto lineTuple = configuration[figureSection][lineKey].as_int_tuple_or_die();
            Face face;
            face.point_indices.push_back(lineTuple[0]);
            face.point_indices.push_back(lineTuple[1]);
            figure.faces.push_back(face);
        }

        figures.push_back(figure);
    }

    return figures;
}

std::tuple<Vector3D, Color, int> parseGeneralSettings(ini::Configuration &configuration) {
    // Parse eye vector
    auto eyeTuple = configuration["General"]["eye"].as_double_tuple_or_die();
    Vector3D eye = Vector3D::vector(eyeTuple[0], eyeTuple[1], eyeTuple[2]);

    // Parse center vector
    vector<double> placeholder {0,0,0};
    auto centerTuple = configuration["General"]["center"].as_double_tuple_or_default(placeholder);
    Vector3D center = Vector3D::vector(centerTuple[0], centerTuple[1], centerTuple[2]);

    // Parse background color
    auto bgColorTuple = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
    Color bgColor(bgColorTuple[0]*255, bgColorTuple[1]*255, bgColorTuple[2]*255);

    // Size
    double size = configuration["General"]["size"].as_int_or_die();


    return {eye, bgColor, size};
}

img::EasyImage linedrawer3D(ini::Configuration &configuration) {

    img::EasyImage image;

    auto [eye,
          bgColor,
          size] = parseGeneralSettings(configuration);

    // Create and initialize your 3D figures
    Figures3D figures = parseiniFigures(configuration);
    double d = 1.0; // Distance parameter for projection

    Lines2D lines2D;

    for (const Figure3D &figure : figures) {
        vector<double> translationVector = figure.center;
        double scale = figure.scale;
        Color lineColor = figure.lineColor;
        Matrix transformation = scaleFigure(scale) *
                rotateX((figure.angleX)) *
                rotateY(figure.angleY) *
                rotateZ(figure.angleZ) *
                translate(translationVector);
        Figure3D transformedFigure = applyTransform(figure, transformation);
        lines2D = projectFigure(transformedFigure, eye, d, lineColor);
    }


    image = draw2DLines(lines2D, size, bgColor);

    return image;
}

