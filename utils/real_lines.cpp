//
// Created by Pim Van den Bosch on 2023/02/24.
//

#include <list>
#include "easy_image.h"
#include "real_lines.h"
#include "cmath"
#include "l_parser.h"
#include <fstream>
#include <limits>

using namespace std;
using namespace LParser;

double max(std::initializer_list<double> values) {
    double max_val = *(values.begin()); // Initialize the maximum value
    for (double value : values) {
        if (value > max_val) max_val = value;
    }
    return max_val;
}

double min(std::initializer_list<double> values) {
    double min_val = *(values.begin()); // Initialize the minimum value
    for (double value : values) {
        if (value < min_val) min_val = value;
    }
    return min_val;
}

std::tuple<double, double, double, double> calculateMinMax(Lines2D& lines) {
    double current_x_max = lines.begin()->p1.x;
    double current_x_min = lines.begin()->p1.x;
    double current_y_max = lines.begin()->p1.y;
    double current_y_min = lines.begin()->p1.y;

    for (Line2D &line : lines) {
        current_x_max = max({line.p1.x, line.p2.x, current_x_max});
        current_x_min = min({line.p1.x, line.p2.x, current_x_min});
        current_y_max = max({line.p1.y, line.p2.y, current_y_max});
        current_y_min = min({line.p1.y, line.p2.y, current_y_min});
    }

    return std::make_tuple(current_x_max, current_x_min, current_y_max, current_y_min);
}

std::tuple<double, double, double, double> calculateMinMaxTriangle(Triangles2D& triangles) {
    double current_x_max = triangles.begin()->p1.x;
    double current_x_min = triangles.begin()->p1.x;
    double current_y_max = triangles.begin()->p1.y;
    double current_y_min = triangles.begin()->p1.y;

    for (Triangle2D &triangle : triangles) {
        current_x_max = max({triangle.p1.x, triangle.p2.x, triangle.p3.x, current_x_max});
        current_x_min = min({triangle.p1.x, triangle.p2.x, triangle.p3.x, current_x_min});
        current_y_max = max({triangle.p1.y, triangle.p2.y, triangle.p3.y, current_y_max});
        current_y_min = min({triangle.p1.y, triangle.p2.y, triangle.p3.y, current_y_min});
    }

    return std::make_tuple(current_x_max, current_x_min, current_y_max, current_y_min);
}


std::pair<double, double> calculateImageDimensions(double current_x_max, double current_x_min, double current_y_max, double current_y_min, double size) {
    double x_range = current_x_max - current_x_min;
    double y_range = current_y_max - current_y_min;

    double image_x = size*(x_range)/max(x_range, y_range);
    double image_y = size*(y_range)/max(x_range, y_range);

    return std::make_pair(image_x, image_y);
}

void scaleAndShiftPoints(Lines2D& lines, double current_x_min, double current_y_min, double x_range, double y_range, double image_x, double image_y) {
    double d = 0.95 * image_x/x_range;

    for (Line2D &line : lines) {
        line.p1.x = d * line.p1.x;
        line.p1.y = d * line.p1.y;
        line.p2.x = d * line.p2.x;
        line.p2.y = d * line.p2.y;

        line.p1.x = line.p1.x - (current_x_min * d) + ((image_x - x_range * d)/2);
        line.p1.y = line.p1.y - (current_y_min * d) + ((image_y - y_range * d)/2);
        line.p2.x = line.p2.x - (current_x_min * d) + ((image_x - x_range * d)/2);
        line.p2.y = line.p2.y - (current_y_min * d) + ((image_y - y_range * d)/2);

        line.p1.x = lround(line.p1.x);
        line.p1.y = lround(line.p1.y);
        line.p2.x = lround(line.p2.x);
        line.p2.y = lround(line.p2.y);
    }
}

void scaleAndShiftPointsTriangles(Triangles2D& triangles, double current_x_min, double current_y_min, double x_range, double y_range, double image_x, double image_y) {
    double d = 0.95 * image_x/x_range;

    for (Triangle2D &triangle : triangles) {
        for(Point2D* p : {&triangle.p1, &triangle.p2, &triangle.p3}){
            p->x = d * p->x;
            p->y = d * p->y;

            p->x = p->x - (current_x_min * d) + ((image_x - x_range * d)/2);
            p->y = p->y - (current_y_min * d) + ((image_y - y_range * d)/2);

            //p->x = lround(p->x);
            //p->y = lround(p->y);
        }
    }
}

// This function takes an input file and a color, reads an L-system from the file,
// and returns a list of 2D lines that represent the L-system after iterating
// it a specified number of times and applying certain drawing rules.
Lines2D L_System2D (const string &inputfile, const Color color)
{
    // Open the input file
    ifstream file(inputfile);

    // Read the L-system from the file
    LSystem2D system;
    file >> system;

    // Get the alphabet, angle, initiator, replacement rules, and other parameters
    set<char> alphabet = system.get_alphabet();
    double angle = system.get_angle(); // The angle in degrees to turn the turtle
    string initiator = system.get_initiator(); // The initial string to start the L-system from
    string replacement;
    double starting_angle = system.get_starting_angle(); // The starting angle of the turtle
    int nr_iterations = system.get_nr_iterations(); // The number of iterations to apply the L-system

    // Iterate the L-system to produce the final string
    for (int i = 0; i < nr_iterations; ++i)
    {
        string next_string = "";

        // Iterate through the characters in the initiator string.
        for (char c: initiator)
        {
            // If the character is in the alphabet and should be drawn, replace it with the replacement string.
            if (alphabet.count(c) && system.draw(c) )
            {
                replacement = system.get_replacement(c);
                next_string += replacement;
            }
            else{
                // Otherwise, keep the character in the string.
                next_string += c;
            }
        }
        // Update the initiator string for the next iteration.
        initiator = next_string;
    }

    // The final string after all iterations
    string final_string = initiator;

    // The starting position of the turtle
    Point2D start_pos = {0,0};

    // The ending position of each line segment drawn by the turtle
    Point2D end_pos;

    // The list of lines representing the L-system
    Lines2D line_list;

    // The current angle of the turtle
    double current_angle = starting_angle;
    // Lists to keep track of points and directions associated with a bracket
    list<Point2D> bracket_points;
    list<double> bracket_angle;

    // Iterate the final string and draw lines or turn the turtle
    for (char c: final_string)
    {
        // If the character is in the alphabet and should be drawn
        if ( alphabet.count(c) && system.draw(c) )
        {
            // Calculate the ending position of the line segment
            end_pos = {start_pos.x + cos(current_angle * M_PI / 180),
                       start_pos.y + sin(current_angle * M_PI / 180)};

            // Create the line segment and add it to the list
            Line2D line(start_pos, end_pos, color);
            line_list.push_back(line);

            // Update the starting position for the next line segment
            start_pos = end_pos;
        }
            // If the character is a turn command
        else
        {
            if (c == '+'){
                // Turn to the left
                current_angle += angle;
            }
            else if (c == '-'){
                // Turn to the right
                current_angle -= angle;
            }
            else if (c == '('){
                // Keep track of the current point and angle
                bracket_points.push_front(start_pos);
                bracket_angle.push_front(current_angle);
            }
            else if (c == ')'){
                // Set current angle and starting position to previously saved bracket
                start_pos = bracket_points.front();
                current_angle = bracket_angle.front();

                bracket_points.pop_front();
                bracket_angle.pop_front();
            }
            else{
                // Move forward without drawing a line
                start_pos = {start_pos.x + cos(current_angle * M_PI / 180),
                             start_pos.y + sin(current_angle * M_PI / 180)};
            }
        }
    }

    // Close the file
    file.close();

    // Return the list of lines representing the L-system
    return line_list;
}


// A function that takes a list of 2D lines and a size parameter
// and returns an EasyImage of the lines drawn on it
img::EasyImage draw2DLines(Lines2D &lines,
                           const int &size,
                           Color &backgroundcolor) {

    // Calculate the maximum and minimum x and y values among all the points in the lines
    auto [current_x_max, current_x_min, current_y_max, current_y_min] = calculateMinMax(lines);

    // Calculate the range of x and y values
    double x_range = current_x_max - current_x_min;
    double y_range = current_y_max - current_y_min;

    // Calculate the dimensions of the image based on the size parameter
    auto [image_x, image_y] = calculateImageDimensions(current_x_max, current_x_min, current_y_max, current_y_min, size);

    // Scale and shift all the points
    scaleAndShiftPoints(lines, current_x_min, current_y_min, x_range, y_range, image_x, image_y);

    // Create a new EasyImage of the appropriate size
    img::EasyImage image(lround(image_x), lround(image_y), backgroundcolor);

    // Draw all the lines on the image
    for (Line2D line: lines){
        img::Color line_color = line.color;
        image.draw_line(line.p1.x, line.p1.y, line.p2.x, line.p2.y, line.color);
    }

    // Return the image
    return image;
}


img::EasyImage draw2DLinesWithZBuffer(Lines2D &lines,
                                      const int &size,
                                      Color &backgroundcolor) {
    // Calculate the maximum and minimum x and y values among all the points in the lines
    auto [current_x_max, current_x_min, current_y_max, current_y_min] = calculateMinMax(lines);

    // Calculate the range of x and y values
    double x_range = current_x_max - current_x_min;
    double y_range = current_y_max - current_y_min;

    // Calculate the dimensions of the image based on the size parameter
    auto [image_x, image_y] = calculateImageDimensions(current_x_max, current_x_min, current_y_max, current_y_min, size);

    // Create a z-buffer with the same dimensions as the image
    vector<vector<double>> zBuffer(image_x, vector<double>(image_y, numeric_limits<double>::infinity()));

    // Scale and shift all the points
    scaleAndShiftPoints(lines, current_x_min, current_y_min, x_range, y_range, image_x, image_y);

    // Create a new EasyImage of the appropriate size
    img::EasyImage image(image_x, image_y, backgroundcolor);

    // Draw all the lines on the image
    for (Line2D line : lines) {
        img::Color line_color = line.color;
        // Custom line drawing function that updates the z-buffer and checks depth values
        image.draw_zbuff_line(line.p1.x, line.p1.y, line.p2.x, line.p2.y, line.color,
                              [&](int x, int y, double t) {
                                  double z = interpolateZ(line.p1, line.p2, t);
                                  if (z < zBuffer[x][y]) {
                                      zBuffer[x][y] = z;
                                      return true;
                                  }
                                  return false;
                              });

    }

    // Return the image
    return image;
}

double interpolateZ(const Point2D& p1, const Point2D& p2, double t) {
    return   t/p1.z + (1 - t)/p2.z ;
}

img::EasyImage draw2DTrianglesWithZBuffer(Triangles2D &triangles,
                                          Figures3D triangulatedFigures3D,
                                          const int &size,
                                          Color &backgroundcolor) {
    // Calculate the maximum and minimum x and y values among all the points in the triangles
    auto [current_x_max, current_x_min, current_y_max, current_y_min] = calculateMinMaxTriangle(triangles);

    // Calculate the range of x and y values
    double x_range = current_x_max - current_x_min;
    double y_range = current_y_max - current_y_min;

    // Calculate the dimensions of the image based on the size parameter
    auto [image_x, image_y] = calculateImageDimensions(current_x_max, current_x_min, current_y_max, current_y_min, size);

    // Create a z-buffer with the same dimensions as the image
    vector<vector<double>> zBuffer(image_x, vector<double>(image_y, -numeric_limits<double>::infinity()));

    // Scale and shift all the points
    scaleAndShiftPointsTriangles(triangles, current_x_min, current_y_min, x_range, y_range, image_x, image_y);

    // Create a new EasyImage of the appropriate size
    img::EasyImage image(image_x, image_y, backgroundcolor);

    // Sort all the triangles based on average z value (depth)
    std::sort(triangles.begin(), triangles.end(), [&](const Triangle2D& a, const Triangle2D& b) {
        // Get corresponding Figure3D objects
        Figure3D& figA = triangulatedFigures3D[a.figureIndex];
        Figure3D& figB = triangulatedFigures3D[b.figureIndex];

        // Compute average depth for each triangle
        double zA = (figA.points[figA.faces[a.faceIndex].point_indices[0]].z +
                     figA.points[figA.faces[a.faceIndex].point_indices[1]].z +
                     figA.points[figA.faces[a.faceIndex].point_indices[2]].z) / 3.0;

        double zB = (figB.points[figB.faces[b.faceIndex].point_indices[0]].z +
                     figB.points[figB.faces[b.faceIndex].point_indices[1]].z +
                     figB.points[figB.faces[b.faceIndex].point_indices[2]].z) / 3.0;

        // Compare
        return zA > zB;  // Change to < for different coordinate system
    });

    // Draw all the triangles on the image
    for (Triangle2D &triangle : triangles) {
        // use the indices stored in the triangle to look up the corresponding Figure3D and face
        Figure3D &figure = triangulatedFigures3D[triangle.figureIndex];
        Face &face = figure.faces[triangle.faceIndex];

        // The unprojected z-values for the vertices of the triangle
        double z1 = figure.points[face.point_indices[0]].z;
        double z2 = figure.points[face.point_indices[1]].z;
        double z3 = figure.points[face.point_indices[2]].z;

        // Custom triangle drawing function that updates the z-buffer and checks depth values
        image.draw_zbuff_triangle(triangle.p1.x, triangle.p1.y, z1,
                                  triangle.p2.x, triangle.p2.y, z2,
                                  triangle.p3.x, triangle.p3.y, z3,
                                  triangle.color,
                                  [&](int x, int y, double z) {
                                      if (z > zBuffer[x][y]) {
                                          zBuffer[x][y] = z;
                                          return true;
                                      }
                                      return false;
                                  });
    }



    // Return the image
    return image;
}








