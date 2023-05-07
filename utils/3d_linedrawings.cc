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
    rotateXMatrix(2, 3) = sinAngle;
    rotateXMatrix(3, 2) = -sinAngle;
    rotateXMatrix(3, 3) = cosAngle;
    return rotateXMatrix;
}

Matrix rotateY(const double angle) {
    Matrix rotateYMatrix;
    double cosAngle = cos(angle);
    double sinAngle = sin(angle);
    rotateYMatrix(1, 1) = cosAngle;
    rotateYMatrix(1, 3) = -sinAngle;
    rotateYMatrix(3, 1) = sinAngle;
    rotateYMatrix(3, 3) = cosAngle;
    return rotateYMatrix;
}

Matrix rotateZ(const double angle) {
    Matrix rotateZMatrix;
    double cosAngle = cos(angle);
    double sinAngle = sin(angle);
    rotateZMatrix(1, 1) = cosAngle;
    rotateZMatrix(1, 2) = sinAngle;
    rotateZMatrix(2, 1) = -sinAngle;
    rotateZMatrix(2, 2) = cosAngle;
    return rotateZMatrix;
}

Matrix translate(const vector<double> &vector) {
    Matrix translationMatrix;
    translationMatrix(4, 1) = vector[0];
    translationMatrix(4, 2) = vector[1];
    translationMatrix(4, 3) = vector[2];
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

Lines2D projectFigure(const Figure3D &figure, const Vector3D &eye, double d, Color lineColor) {
    Lines2D lines2D;

    for (const auto& face : figure.faces) {
        for (size_t i = 0; i < face.point_indices.size(); ++i) {
            Vector3D startPoint = figure.points[face.point_indices[i]];
            Vector3D endPoint = figure.points[face.point_indices[(i + 1) % face.point_indices.size()]];

            startPoint = project(startPoint, eye);
            startPoint.x = -startPoint.x * d / (startPoint.z);
            startPoint.y = -startPoint.y * d / (startPoint.z);

            endPoint = project(endPoint, eye);
            endPoint.x = -endPoint.x * d / (endPoint.z);
            endPoint.y = -endPoint.y * d / (endPoint.z);

            Point2D startPoint2D = {startPoint.x, startPoint.y};
            Point2D endPoint2D = {endPoint.x, endPoint.y};

            lines2D.push_back(Line2D(startPoint2D, endPoint2D, lineColor));
        }
    }

    return lines2D;
}

Figure3D createCube(Figure3D cube) {

    cube.points.push_back(Vector3D::point(1, -1, -1));
    cube.points.push_back(Vector3D::point(-1, 1, -1));
    cube.points.push_back(Vector3D::point(1, 1, 1));
    cube.points.push_back(Vector3D::point(-1, -1, 1));
    cube.points.push_back(Vector3D::point(1, 1, -1));
    cube.points.push_back(Vector3D::point(-1, -1, -1));
    cube.points.push_back(Vector3D::point(1, -1, 1));
    cube.points.push_back(Vector3D::point(-1, 1, 1));

    Face face1;
    face1.point_indices.push_back(0);
    face1.point_indices.push_back(4);
    face1.point_indices.push_back(2);
    face1.point_indices.push_back(6);
    cube.faces.push_back(face1);

    Face face2;
    face2.point_indices.push_back(4);
    face2.point_indices.push_back(1);
    face2.point_indices.push_back(7);
    face2.point_indices.push_back(2);
    cube.faces.push_back(face2);

    Face face3;
    face3.point_indices.push_back(1);
    face3.point_indices.push_back(5);
    face3.point_indices.push_back(3);
    face3.point_indices.push_back(7);
    cube.faces.push_back(face3);

    Face face4;
    face4.point_indices.push_back(5);
    face4.point_indices.push_back(0);
    face4.point_indices.push_back(6);
    face4.point_indices.push_back(3);
    cube.faces.push_back(face4);

    Face face5;
    face5.point_indices.push_back(6);
    face5.point_indices.push_back(2);
    face5.point_indices.push_back(7);
    face5.point_indices.push_back(3);
    cube.faces.push_back(face5);

    Face face6;
    face6.point_indices.push_back(0);
    face6.point_indices.push_back(5);
    face6.point_indices.push_back(1);
    face6.point_indices.push_back(4);
    cube.faces.push_back(face6);

    return cube;
}
Figure3D createTetrahedron(Figure3D tetrahedron) {

    tetrahedron.points.push_back(Vector3D::point(1, -1, -1));
    tetrahedron.points.push_back(Vector3D::point(-1, 1, -1));
    tetrahedron.points.push_back(Vector3D::point(1, 1, 1));
    tetrahedron.points.push_back(Vector3D::point(-1, -1, 1));

    Face face1;
    face1.point_indices.push_back(0);
    face1.point_indices.push_back(1);
    face1.point_indices.push_back(2);
    tetrahedron.faces.push_back(face1);

    Face face2;
    face2.point_indices.push_back(1);
    face2.point_indices.push_back(3);
    face2.point_indices.push_back(2);
    tetrahedron.faces.push_back(face2);

    Face face3;
    face3.point_indices.push_back(0);
    face3.point_indices.push_back(3);
    face3.point_indices.push_back(1);
    tetrahedron.faces.push_back(face3);

    Face face4;
    face4.point_indices.push_back(0);
    face4.point_indices.push_back(2);
    face4.point_indices.push_back(3);
    tetrahedron.faces.push_back(face4);

    return tetrahedron;
}
Figure3D createOctahedron(Figure3D octahedron) {

    octahedron.points.push_back(Vector3D::point(1, 0, 0));
    octahedron.points.push_back(Vector3D::point(0, 1, 0));
    octahedron.points.push_back(Vector3D::point(-1, 0, 0));
    octahedron.points.push_back(Vector3D::point(0, -1, 0));
    octahedron.points.push_back(Vector3D::point(0, 0, -1));
    octahedron.points.push_back(Vector3D::point(0, 0, 1));


    Face face1;
    face1.point_indices.push_back(0);
    face1.point_indices.push_back(1);
    face1.point_indices.push_back(5);
    octahedron.faces.push_back(face1);

    Face face2;
    face2.point_indices.push_back(1);
    face2.point_indices.push_back(2);
    face2.point_indices.push_back(5);
    octahedron.faces.push_back(face2);

    Face face3;
    face3.point_indices.push_back(2);
    face3.point_indices.push_back(3);
    face3.point_indices.push_back(5);
    octahedron.faces.push_back(face3);

    Face face4;
    face4.point_indices.push_back(3);
    face4.point_indices.push_back(0);
    face4.point_indices.push_back(5);
    octahedron.faces.push_back(face4);

    Face face5;
    face5.point_indices.push_back(1);
    face5.point_indices.push_back(0);
    face5.point_indices.push_back(4);
    octahedron.faces.push_back(face5);

    Face face6;
    face6.point_indices.push_back(2);
    face6.point_indices.push_back(1);
    face6.point_indices.push_back(4);
    octahedron.faces.push_back(face6);

    Face face7;
    face7.point_indices.push_back(3);
    face7.point_indices.push_back(2);
    face7.point_indices.push_back(4);
    octahedron.faces.push_back(face7);

    Face face8;
    face8.point_indices.push_back(0);
    face8.point_indices.push_back(3);
    face8.point_indices.push_back(4);
    octahedron.faces.push_back(face8);

    return octahedron;
}
Figure3D createIcosahedron(Figure3D icosahedron) {

    icosahedron.points.push_back(Vector3D::point(0, 0, sqrt(5)/2));

    icosahedron.points.push_back(Vector3D::point( cos((2-2)*2*M_PI/5), sin((2-2)*2*M_PI/5), 0.5));
    icosahedron.points.push_back(Vector3D::point( cos((3-2)*2*M_PI/5), sin((3-2)*2*M_PI/5), 0.5));
    icosahedron.points.push_back(Vector3D::point( cos((4-2)*2*M_PI/5), sin((4-2)*2*M_PI/5), 0.5));
    icosahedron.points.push_back(Vector3D::point( cos((5-2)*2*M_PI/5), sin((5-2)*2*M_PI/5), 0.5));
    icosahedron.points.push_back(Vector3D::point( cos((6-2)*2*M_PI/5), sin((6-2)*2*M_PI/5), 0.5));

    icosahedron.points.push_back(Vector3D::point(cos((M_PI/5)+(7-7)*(2*M_PI)/5), sin((M_PI/5)+(7-7)*(2*M_PI)/5), -0.5));
    icosahedron.points.push_back(Vector3D::point(cos((M_PI/5)+(8-7)*(2*M_PI)/5), sin((M_PI/5)+(8-7)*(2*M_PI)/5), -0.5));
    icosahedron.points.push_back(Vector3D::point(cos((M_PI/5)+(9-7)*(2*M_PI)/5), sin((M_PI/5)+(9-7)*(2*M_PI)/5), -0.5));
    icosahedron.points.push_back(Vector3D::point(cos((M_PI/5)+(10-7)*(2*M_PI)/5), sin((M_PI/5)+(10-7)*(2*M_PI)/5), -0.5));
    icosahedron.points.push_back(Vector3D::point(cos((M_PI/5)+(11-7)*(2*M_PI)/5), sin((M_PI/5)+(11-7)*(2*M_PI)/5), -0.5));

    icosahedron.points.push_back(Vector3D::point(0, 0, -sqrt(5)/2));


    const int face_indices[][3] = {
            {0, 1, 2},
            {0, 2, 3},
            {0, 3, 4},
            {0, 4, 5},
            {0, 5, 1},
            {1, 6, 2},
            {2, 6, 7},
            {2, 7, 3},
            {3, 7, 8},
            {3, 8, 4},
            {4, 8, 9},
            {4, 9, 5},
            {5, 9, 10},
            {5, 10, 1},
            {1, 10, 6},
            {11, 7, 6},
            {11, 8, 7},
            {11, 9, 8},
            {11, 10, 9},
            {11, 6, 10},
    };

    for (int i = 0; i < 20; ++i) {
        Face face;
        face.point_indices.push_back(face_indices[i][0]);
        face.point_indices.push_back(face_indices[i][1]);
        face.point_indices.push_back(face_indices[i][2]);
        icosahedron.faces.push_back(face);
    }

    return icosahedron;
}
Figure3D createDodecahedron(Figure3D icosahedron) {
    Figure3D dodecahedron;
    // Step 0: Ensure that the dodecahedron inherits the properties from icosahedron.
    dodecahedron.center = icosahedron.center;
    dodecahedron.lineColor = icosahedron.lineColor;
    dodecahedron.angleX = icosahedron.angleX;
    dodecahedron.angleY = icosahedron.angleY;
    dodecahedron.angleZ = icosahedron.angleZ;
    dodecahedron.scale = icosahedron.scale;

    // Step 1: Compute the centroid of each face of the icosahedron.
    for (const Face& face : icosahedron.faces) {
        Vector3D centroid = Vector3D::point(0, 0, 0);
        for (int idx : face.point_indices) {
            centroid += icosahedron.points[idx];
        }
        centroid /= face.point_indices.size();

        // Step 2: Normalize the centroid points to form the vertices of the dodecahedron.
        double radius = sqrt(3); // Assuming icosahedron edge length is 2
        centroid = Vector3D::point(centroid.x * radius / centroid.length(),
                                   centroid.y * radius / centroid.length(),
                                   centroid.z * radius / centroid.length());
        dodecahedron.points.push_back(centroid);
    }

    // Step 3: Define the faces of the dodecahedron using the new vertices.
    const int face_indices[][5] = {
            {0,  1,  2,  3,  4},
            {0,  5, 6, 7,  1},
            {1,  7, 8, 9,  2},
            {2, 9, 10, 11,  3},
            {3, 11, 12,  13,  4},
            {4, 13, 14,  5,  0},
            {19, 18, 17,  16,  15},
            {19, 14, 13, 12, 18},
            {18,  12,  11,  10,  17},
            {17,  10,  9, 8, 16},
            {16,  8, 7, 6,  15},
            {15,  6, 5,  14,  19}
    };

    for (int i = 0; i < 12; ++i) {
        Face face;
        for (int j = 0; j < 5; ++j) {
            face.point_indices.push_back(face_indices[i][j]);
        }
        dodecahedron.faces.push_back(face);
    }

    return dodecahedron;
}
Figure3D createCylinder(Figure3D cylinder) {

    int n = cylinder.n;
    double h = cylinder.height;
    // Create top and bottom circle points
    for (int i = 0; i < n; ++i) {
        double angle = 2 * M_PI * i / n;
        double x = cos(angle);
        double z = sin(angle);

        // Bottom circle point
        cylinder.points.push_back(Vector3D::point(x, z, -h / 2));
        // Top circle point
        cylinder.points.push_back(Vector3D::point(x, z, h / 2));
    }

    // Create bottom and top faces
    Face bottom_face;
    Face top_face;
    for (int i = 0; i < n; ++i) {
        bottom_face.point_indices.push_back(2 * i);
        top_face.point_indices.push_back(2 * i + 1);
    }
    cylinder.faces.push_back(bottom_face);
    cylinder.faces.push_back(top_face);

    // Create side faces
    for (int i = 0; i < n; ++i) {
        int i_next = (i + 1) % n;
        Face side_face;
        side_face.point_indices.push_back(2 * i);          // Bottom point
        side_face.point_indices.push_back(2 * i_next);     // Next bottom point
        side_face.point_indices.push_back(2 * i_next + 1); // Next top point
        side_face.point_indices.push_back(2 * i + 1);      // Top point
        cylinder.faces.push_back(side_face);
    }

    return cylinder;
}
Figure3D createCone(Figure3D cone) {

    int n = cone.n;
    double h = cone.height;
    // Add the top point of the cone
    cone.points.push_back(Vector3D::point(0, 0, h));

    // Add the points of the base in counterclockwise order
    for (int i = 0; i < n; ++i) {
        double angle = 2 * M_PI * i / n;
        double x = cos(angle);
        double y = sin(angle);
        cone.points.push_back(Vector3D::point(x, y, 0));
    }

    // Add the triangular faces connecting the top point to the base
    for (int i = 1; i <= n; ++i) {
        Face face;
        face.point_indices.push_back(0); // The top point of the cone
        face.point_indices.push_back(i);
        face.point_indices.push_back(i % n + 1);
        cone.faces.push_back(face);
    }

    // Add the base face
    Face baseFace;
    for (int i = 1; i <= n; ++i) {
        baseFace.point_indices.push_back(i);
    }
    cone.faces.push_back(baseFace);

    return cone;
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
        figure.n = configuration[figureSection]["n"].as_int_or_default(0);
        figure.height = configuration[figureSection]["height"].as_double_or_default(0.0);

        // Parse the color
        auto colorTuple = configuration[figureSection]["color"].as_double_tuple_or_die();
        figure.lineColor.red = colorTuple[0] * 255;
        figure.lineColor.green = colorTuple[1] * 255;
        figure.lineColor.blue = colorTuple[2] * 255;

        std::string figure_type = configuration[figureSection]["type"].as_string_or_die();

        if(figure_type == "LineDrawing"){
            int nrPoints = configuration[figureSection]["nrPoints"].as_int_or_die();
            int nrLines = configuration[figureSection]["nrLines"].as_int_or_die();

            // Parse the points
            for (int j = 0; j < nrPoints; j++) {
                std::string pointKey = "point" + std::to_string(j);
                auto pointTuple = configuration[figureSection][pointKey].as_double_tuple_or_die();
                Vector3D point = Vector3D::point(pointTuple[0], pointTuple[1], pointTuple[2]);
                figure.points.push_back(point);
            }

            // Parse the lines (faces)
            for (int j = 0; j < nrLines; j++) {
                std::string lineKey = "line" + std::to_string(j);
                auto lineTuple = configuration[figureSection][lineKey].as_int_tuple_or_die();
                Face face;
                face.point_indices.push_back(lineTuple[0]);
                face.point_indices.push_back(lineTuple[1]);
                figure.faces.push_back(face);
            }
        }
        else if (figure_type == "Cube"){
            figure = createCube(figure);
        }
        else if (figure_type == "Tetrahedron"){
            figure = createTetrahedron(figure);
        }
        else if (figure_type == "Octahedron"){
            figure = createOctahedron(figure);
        }
        else if (figure_type == "Icosahedron"){
            figure = createIcosahedron(figure);
        }
        else if (figure_type == "Dodecahedron"){
            figure = createIcosahedron(figure);
            figure = createDodecahedron(figure);
        }
        else if (figure_type == "Cylinder"){
            figure = createCylinder(figure);
        }
        else if (figure_type == "Cone"){
            figure = createCone(figure);
        }
        else if (figure_type == "Sphere"){
            //figure = createSphere();
        }
        else if (figure_type == "Torus"){
            //figure = createTorus();
        }
        else if (figure_type == "3DLSystem"){
            //figure = create3DLSystem();
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
    Lines2D allLines2D;

    for (const Figure3D &figure : figures) {
        vector<double> translationVector = figure.center;
        double scale = figure.scale;
        Color lineColor = figure.lineColor;
        Matrix transformation = scaleFigure(scale) *
                                rotateX((figure.angleX * M_PI) / 180) *
                                rotateY((figure.angleY * M_PI) / 180) *
                                rotateZ((figure.angleZ * M_PI) / 180) *
                                translate(translationVector);
        Figure3D transformedFigure = applyTransform(figure, transformation);
        lines2D = projectFigure(transformedFigure, eye, d, lineColor);

        for(const Line2D &line : lines2D){
            allLines2D.push_back(line);
        }
    }

    image = draw2DLines(allLines2D, size, bgColor);

    return image;
}


