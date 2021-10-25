#include "project3/file.hpp"

std::ofstream get_file(std::string filename) {
    std::ofstream ofile;
    ofile.open("data/" + filename);
    ofile << std::scientific;
    return ofile;
}

std::ofstream get_solution_file(std::string filename) {
    std::ofstream ofile = get_file(filename);
    ofile << "x,y,z" << std::endl;
    return ofile;
}

std::ofstream get_solution_file_double_x_y_z(std::string filename) {
    std::ofstream ofile = get_file(filename);
    ofile << "x1,y1,z1,x2,y2,z2" << std::endl;
    return ofile;
}
