#include "project3/analytical.hpp"

// std::vector<arma::vec> solve_analytically(arma::vec t, double q, double m, double x0, double z0, double v0, double B0, double V0, double d) {
// int solve_analytically(arma::vec t, double q, double m, double x0, double z0, double v0, double B0, double V0, double d) {
//     std::cout << "I:JFEOI" << std::endl;
//
//     // std::vector<arma::vec> solution = {1, 2, 3};
//     // std::vector<arma::vec> solution;
//     //
//     // return solution;
//     return 69;
// }

std::vector<arma::vec> solve_analytically(arma::vec t, double q, double m,
                                          double x0, double z0, double v0,
                                          double B0, double V0, double d) {
    // Computes the the analytical solution of a penning trap with one qd
    // particle in it, assuming y0 = 0, and the velocity in x and z directions
    // are also 0

    arma::vec p = { x0, 0.0, z0 };
    arma::vec v = { 0.0, v0, 0 };

    double omega0 = q * B0 / m;
    double omegaz2 = 2 * q * V0 / (m * pow(d, 2));
    double omegaz = sqrt(omegaz2);

    double omegap = (omega0 + sqrt(pow(omega0, 2) - 2*omegaz2))/2;
    double omegam = (omega0 - sqrt(pow(omega0, 2) - 2*omegaz2))/2;

    double Ap = (v0 + omegam * x0)/(omegam - omegap);
    double Am = - (v0 + omegap * x0)/(omegam - omegap);

    // int N = 100; double dt = 1e-4;

    arma::vec x(t.size()), y(t.size()), z(t.size());
    for (int i = 0; i < t.size(); i++) {
        x[i] = Ap * cos(- omegap * t[i]) + Am * cos(- omegam * t[i]);
        y[i] = Ap * sin(- omegap * t[i]) + Am * sin(- omegam * t[i]);
        z[i] = z0 * cos(omegaz * t[i]);
    }

    std::vector<arma::vec> solution = {x, y, z};

    return solution;
}
