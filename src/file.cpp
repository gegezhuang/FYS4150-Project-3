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
