// #include <assert.h>
//
// #include "project3/Particle.hpp"
// #include "project3/PenningTrap.hpp"
// #include "project3/analytical.hpp"
//
// int main() {
//     double q = 1, m = 1;
//     double x0 = 0.5, z0 = 1.0, v0 = 1.0;
//     arma::vec p = { x0, 0.0, z0 };
//     arma::vec v = { 0.0, v0, 0 };
//
//     Particle particle(q, m, p, v);
//
//     double B0 = 4, V0 = 2.5, d = 3;
//     PenningTrap pt(B0, V0, d);
//     pt.add_particle(particle);
//
//     int N = 100; double dt = 1e-4;
//     pt.solve_RK4(N, dt);
//     std::vector<arma::mat> numerical = pt.get_solution();
//     std::vector<double> t = pt.get_time();
//
//     std::vector<arma::vec> analytical = solve_analytically(t, q,  m, x0, z0, v0, B0, V0, d);
//
//     double x_err, y_err;
//     for (int i = 1; i < N; i++) {
//         x_err = fabs((analytical[0][i] - numerical[i][0])/analytical[0][i]);
//         y_err = fabs((analytical[1][i] - numerical[i][1])/analytical[1][i]);
//         // TODO: Add z_err when that is computed analytically
//
//         // Asserts that the relative error stays below 10^-3
//         assert (x_err < 1e-3);
//         assert (y_err < 1e-3);
//     }
//     return 0;
// }
