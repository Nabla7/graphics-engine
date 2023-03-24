#include "vector/vector3d.h"
#include "3d_linedrawings.h"
#include "easy_image.h"
#include "real_lines.h"


Matrix scaleFigure(const double scale) {
    Matrix scaleMatrix;
    scaleMatrix(0, 0) = scale;
    scaleMatrix(1, 1) = scale;
    scaleMatrix(2, 2) = scale;
    return scaleMatrix;
}

Matrix rotateX(const double angle) {
    Matrix rotateXMatrix;
    double cosAngle = cos(angle);
    double sinAngle = sin(angle);
    rotateXMatrix(1, 1) = cosAngle;
    rotateXMatrix(1, 2) = -sinAngle;
    rotateXMatrix(2, 1) = sinAngle;
    rotateXMatrix(2, 2) = cosAngle;
    return rotateXMatrix;
}

Matrix rotateY(const double angle) {
    Matrix rotateYMatrix;
    double cosAngle = cos(angle);
    double sinAngle = sin(angle);
    rotateYMatrix(0, 0) = cosAngle;
    rotateYMatrix(0, 2) = sinAngle;
    rotateYMatrix(2, 0) = -sinAngle;
    rotateYMatrix(2, 2) = cosAngle;
    return rotateYMatrix;
}

Matrix rotateZ(const double angle) {
    Matrix rotateZMatrix;
    double cosAngle = cos(angle);
    double sinAngle = sin(angle);
    rotateZMatrix(0, 0) = cosAngle;
    rotateZMatrix(0, 1) = -sinAngle;
    rotateZMatrix(1, 0) = sinAngle;
    rotateZMatrix(1, 1) = cosAngle;
    return rotateZMatrix;
}

Matrix translate(const Vector3D &vector) {
    Matrix translationMatrix;
    translationMatrix(0, 3) = vector.x;
    translationMatrix(1, 3) = vector.y;
    translationMatrix(2, 3) = vector.z;
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
    transformMatrix(0,0) = -sin(theta);
    transformMatrix(0,1) = -cos(theta)*-cos(phi);
    transformMatrix(0,2) = cos(theta)*sin(phi);
    transformMatrix(1,0) = cos(theta);
    transformMatrix(1,1) = -sin(theta)*cos(phi);
    transformMatrix(1,2) = sin(theta)*sin(phi);
    transformMatrix(2,1) = sin(phi);
    transformMatrix(2,2) = cos(phi);
    transformMatrix(3,2) = -r;

    return transformMatrix;
}

Vector3D multiplyMatrixVector(const Matrix &mat, const Vector3D point) {
    Vector3D result;
    result.x = mat(1, 1) * point.x + mat(1, 2) * point.y + mat(1, 3) * point.z;
    result.y = mat(2, 1) * point.x + mat(2, 2) * point.y + mat(2, 3) * point.z;
    result.z = mat(3, 1) * point.x + mat(3, 2) * point.y + mat(3, 3) * point.z;
    return result;
}

Figure3D applyTransform(const Figure3D &figure, const Matrix &transformation) {
    Figure3D transformedFigure = figure;
    for (Vector3D &point : transformedFigure.points) {
        point = multiplyMatrixVector(transformation, point);
    }
    return transformedFigure;
}

void project(const Figure3D &figure, double d, Color lineColor, Lines2D &lines2D) {
    for (const Face &face : figure.faces) {
        Point2D previousPoint;
        bool isFirstPoint = true;
        for (int i = 0; i < face.point_indices.size(); i++) {
            Vector3D projectedPoint = figure.points[face.point_indices[i]];
            projectedPoint.x = projectedPoint.x * d / projectedPoint.z;
            projectedPoint.y = projectedPoint.y * d / projectedPoint.z;

            Point2D currentPoint;
            currentPoint.x = projectedPoint.x;
            currentPoint.y = projectedPoint.y;

            if (isFirstPoint) {
                isFirstPoint = false;
            } else {
                lines2D.push_back(Line2D(previousPoint, currentPoint, lineColor));
            }

            previousPoint = currentPoint;
        }
        // Connect the last point to the first point
        if (!isFirstPoint) {
            Vector3D firstProjectedPoint = figure.points[face.point_indices[0]];
            firstProjectedPoint.x = firstProjectedPoint.x * d / firstProjectedPoint.z;
            firstProjectedPoint.y = firstProjectedPoint.y * d / firstProjectedPoint.z;
            Point2D firstPoint;
            firstPoint.x = firstProjectedPoint.x;
            firstPoint.y = firstProjectedPoint.y;
            lines2D.push_back(Line2D(previousPoint, firstPoint, lineColor));
        }
    }
}


Lines2D doProjection(const Figures3D &figures, const Vector3D &eye, double d){
    Lines2D lines2D;
    Matrix transformMatrix = eyeMatrix(eye); // Assuming eyeMatrix() creates the transformation matrix
    Color lineColor; // Set the desired line color

    for (const Figure3D &figure : figures){
        Figure3D transformedFigure = figure; // Create a copy of the figure to apply the transformation
        applyTransform(transformedFigure, transformMatrix);
    }

    return lines2D;
};


