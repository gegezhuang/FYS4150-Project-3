#pragma once

#include <armadillo>

// std::vector<arma::vec> solve_analytically(arma::vec t, double q, double m,
std::vector<arma::vec> solve_analytically(arma::vec t, double q, double m,
                                          double x0, double z0, double v0,
                                          double B0, double V0, double d);


// std::vector<arma::vec> solve_analytically(arma::vec t, double q, double m,
//                                           double x0, double z0, double v0,
//                                           double B0, double V0, double d);
    // Computes the the analytical solution of a penning trap with one charged 
    // particle in it, assuming y0 = 0, and the velocity in x and z directions 
    // are also 0

    // TODO: Consider changing the format so that indexing is the same here as 
    // for the numerical solution
