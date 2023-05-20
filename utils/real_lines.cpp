//
// Created by Pim Van den Bosch on 2023/02/24.
//

#include <list>
#include "easy_image.h"
#include "real_lines.h"
#include "cmath"
#include "l_parser.h"
#include <fstream>

using namespace std;
using namespace LParser;

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

    // Find the maximum and minimum x and y values among all the points in the lines
    double current_x_max = lines.begin()->p1.x;
    double current_x_min = lines.begin()->p1.x;
    double current_y_max = lines.begin()->p1.y;
    double current_y_min = lines.begin()->p1.y;

    // Loop over all the lines to find the maximum and minimum x and y values
    for (Line2D &line : lines) {
        current_x_max = max({line.p1.x, line.p2.x, current_x_max});
        current_x_min = min({line.p1.x, line.p2.x, current_x_min});
        current_y_max = max({line.p1.y, line.p2.y, current_y_max});
        current_y_min = min({line.p1.y, line.p2.y, current_y_min});
    }

    // Find the range of x and y values
    double x_range = current_x_max - current_x_min;
    double y_range = current_y_max - current_y_min;

    // Calculate the dimensions of the image based on the size parameter
    double image_x = size*(x_range)/max(x_range, y_range);
    double image_y = size*(y_range)/max(x_range, y_range);

    // Calculate the scaling factor
    double d = 0.95 * image_x/x_range;

    // Scale and shift all the points
    for (Line2D &line : lines) {
        line.p1.x = d * line.p1.x;
        line.p1.y = d * line.p1.y;
        line.p2.x = d * line.p2.x;
        line.p2.y = d * line.p2.y;

        line.p1.x = line.p1.x - (current_x_min * d) + ((image_x - x_range * d)/2);
        line.p1.y = line.p1.y - (current_y_min * d) + ((image_y - y_range * d)/2);
        line.p2.x = line.p2.x - (current_x_min * d) + ((image_x - x_range * d)/2);
        line.p2.y = line.p2.y - (current_y_min * d) + ((image_y - y_range * d)/2);

        // Round the coordinates to the nearest integer
        line.p1.x = lround(line.p1.x);
        line.p1.y = lround(line.p1.y);
        line.p2.x = lround(line.p2.x);
        line.p2.y = lround(line.p2.y);
    }

    // Create a new EasyImage of the appropriate size
    img::EasyImage image(image_x, image_y, backgroundcolor);

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
                                      Color &backgroundcolor,
                                      vector<vector<double>> &zBuffer) {
    // Find the maximum and minimum x and y values among all the points in the lines
    double current_x_max = lines.begin()->p1.x;
    double current_x_min = lines.begin()->p1.x;
    double current_y_max = lines.begin()->p1.y;
    double current_y_min = lines.begin()->p1.y;

    // Loop over all the lines to find the maximum and minimum x and y values
    for (Line2D &line : lines) {
        current_x_max = max({line.p1.x, line.p2.x, current_x_max});
        current_x_min = min({line.p1.x, line.p2.x, current_x_min});
        current_y_max = max({line.p1.y, line.p2.y, current_y_max});
        current_y_min = min({line.p1.y, line.p2.y, current_y_min});
    }

    // Find the range of x and y values
    double x_range = current_x_max - current_x_min;
    double y_range = current_y_max - current_y_min;

    // Calculate the dimensions of the image based on the size parameter
    double image_x = size*(x_range)/max(x_range, y_range);
    double image_y = size*(y_range)/max(x_range, y_range);

    // Calculate the scaling factor
    double d = 0.95 * image_x/x_range;

    // Scale and shift all the points
    for (Line2D &line : lines) {
        line.p1.x = d * line.p1.x;
        line.p1.y = d * line.p1.y;
        line.p2.x = d * line.p2.x;
        line.p2.y = d * line.p2.y;

        line.p1.x = line.p1.x - (current_x_min * d) + ((image_x - x_range * d)/2);
        line.p1.y = line.p1.y - (current_y_min * d) + ((image_y - y_range * d)/2);
        line.p2.x = line.p2.x - (current_x_min * d) + ((image_x - x_range * d)/2);
        line.p2.y = line.p2.y - (current_y_min * d) + ((image_y - y_range * d)/2);


        // Round the coordinates to the nearest integer
        line.p1.x = lround(line.p1.x);
        line.p1.y = lround(line.p1.y);
        line.p2.x = lround(line.p2.x);
        line.p2.y = lround(line.p2.y);
    }


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


