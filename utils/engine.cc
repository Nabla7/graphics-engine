#include "easy_image.h"
#include "ini_configuration.h"
#include "week_0.h"
#include "real_lines.h"
#include "3d_linedrawings.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace std;


img::EasyImage generate_image(ini::Configuration &configuration)
{
    img::EasyImage image;

    string type;
    int size;
    int width;
    int height;
    Color backgroundcolor;
    vector<double> bgc_vector;
    Color linecolor;
    vector<double> lc_vector;
    string inputfile;

    try {
        type = configuration["General"]["type"].as_string_or_die();
        size = configuration["General"]["size"].as_int_or_die();
        width = configuration["ImageProperties"]["width"].as_int_or_default(0);
        height = configuration["ImageProperties"]["height"].as_int_or_default(0);
        bgc_vector = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        lc_vector = configuration["2DLSystem"]["color"].as_double_tuple_or_die();
        inputfile = configuration["2DLSystem"]["inputfile"].as_string_or_die();

        backgroundcolor.red = bgc_vector[0]*255;
        backgroundcolor.green= bgc_vector[1]*255;
        backgroundcolor.blue = bgc_vector[2]*255;

        linecolor.red = lc_vector[0]*255;
        linecolor.green = lc_vector[1]*255;
        linecolor.blue = lc_vector[2]*255;

    }
    catch (ini::NonexistentEntry){

    }

    if (type == "ZBuffering"){
        image = linedrawer3DWithZBufferTriangles(configuration);
    }

    if (type == "ZBufferedWireframe"){
        image = linedrawer3DWithZBuffer(configuration);
    }

    if (type == "Wireframe"){
        image = linedrawer3D(configuration);
    }

    if (type == "2DLSystem")
    {
        Lines2D line_list = L_System2D(inputfile, linecolor);
        image = draw2DLines(line_list, size, backgroundcolor);
    }

    if (type == "IntroColorRectangle")
    {
        image = ColorRectangle(width, height);
    }

    if (type == "IntroBlocks")
    {
        int nrXBlocks = configuration["BlockProperties"]["nrXBlocks"].as_int_or_die();
        int nrYBlocks = configuration["BlockProperties"]["nrYBlocks"].as_int_or_die();
        vector<double> colorWhite = configuration["BlockProperties"]["colorWhite"].as_double_tuple_or_die();
        vector<double> colorBlack = configuration["BlockProperties"]["colorBlack"].as_double_tuple_or_die();
        bool invertColors = configuration["BlockProperties"]["invertColors"].as_bool_or_die();

        image = Blocks(width, height, nrXBlocks, nrYBlocks, colorWhite, colorBlack, invertColors);
    }

    if (type == "IntroLines")
    {
        int nrLines = configuration["LineProperties"]["nrLines"].as_int_or_die();
        vector<double> backgroundcolor = configuration["LineProperties"]["backgroundcolor"].as_double_tuple_or_die();
        vector<double> lineColor = configuration["LineProperties"]["lineColor"].as_double_tuple_or_die();
        image = QuarterCircle(width, height, nrLines, backgroundcolor, lineColor);
    }

    return image;
}

int main(int argc, char const* argv[])
{
        int retVal = 0;
        try
        {
                std::vector<std::string> args = std::vector<std::string>(argv+1, argv+argc);
                if (args.empty()) {
                        std::ifstream fileIn("filelist");
                        std::string filelistName;
                        while (std::getline(fileIn, filelistName)) {
                                args.push_back(filelistName);
                        }
                }
                for(std::string fileName : args)
                {
                        ini::Configuration conf;
                        try
                        {
                                std::ifstream fin(fileName);
                                if (fin.peek() == std::istream::traits_type::eof()) {
                                    std::cout << "Ini file appears empty. Does '" <<
                                    fileName << "' exist?" << std::endl;
                                    continue;
                                }
                                fin >> conf;
                                fin.close();
                        }
                        catch(ini::ParseException& ex)
                        {
                                std::cerr << "Error parsing file: " << fileName << ": " << ex.what() << std::endl;
                                retVal = 1;
                                continue;
                        }

                        img::EasyImage image = generate_image(conf);
                        if(image.get_height() > 0 && image.get_width() > 0)
                        {
                                std::string::size_type pos = fileName.rfind('.');
                                if(pos == std::string::npos)
                                {
                                        //filename does not contain a '.' --> append a '.bmp' suffix
                                        fileName += ".bmp";
                                }
                                else
                                {
                                        fileName = fileName.substr(0,pos) + ".bmp";
                                }
                                try
                                {
                                        std::ofstream f_out(fileName.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                                        f_out << image;

                                }
                                catch(std::exception& ex)
                                {
                                        std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                                        retVal = 1;
                                }
                        }
                        else
                        {
                                std::cout << "Could not generate image for " << fileName << std::endl;
                        }
                }
        }
        catch(const std::bad_alloc &exception)
        {
    		//When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
    		//Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
    		//(Unless of course you are already consuming the maximum allowed amount of memory)
    		//If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
		//mark the test as failed while in reality it just needed a bit more memory
                std::cerr << "Error: insufficient memory" << std::endl;
                retVal = 100;
        }
        return retVal;
}
