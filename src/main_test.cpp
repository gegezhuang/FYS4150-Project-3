#include <assert.h>

#include "project3/Particle.hpp"
#include "project3/PenningTrap.hpp"

int main() {
    double charge = 1, mass = 1;
    double x0 = 0.5, z0 = 1.0, v0 = 1.0;
    arma::vec pos = { x0, 0.0, z0 };
    arma::vec vel = { 0.0, v0, 0 };

    Particle particle(charge, mass, pos, vel);

    double B0 = 4, V0 = 2.5, d = 3;
    PenningTrap pt(B0, V0, d);
    pt.add_particle(particle);

    double omega0 = charge * B0 / mass;
    double omegaz2 = 2 * charge * V0 / (mass * pow(d, 2));

    double omegap = (omega0 + sqrt(pow(omega0, 2) - 2*omegaz2))/2;
    double omegam = (omega0 - sqrt(pow(omega0, 2) - 2*omegaz2))/2;

    double Ap = (v0 + omegam * x0)/(omegam - omegap);
    double Am = - (v0 + omegap * x0)/(omegam - omegap);

    int N = 100; double dt = 1e-4;
    pt.solve_RK4(N, dt);
    std::vector<arma::mat> solution = pt.get_solution();
    std::vector<double> t = pt.get_time();
    double x_analytic, y_analytic, x_err, y_err;
    for (int i = 1; i < N; i++) {
        x_analytic = Ap * cos(- omegap * t[i]) + Am * cos(- omegam * t[i]);
        y_analytic = Ap * sin(- omegap * t[i]) + Am * sin(- omegam * t[i]);

        x_err = fabs((x_analytic - solution[i][0])/x_analytic);
        y_err = fabs((y_analytic - solution[i][1])/y_analytic);

        // Asserts that the relative error stays below 10^-3
        assert(x_err < 1e-3);
        assert(y_err < 1e-3);
    }
    return 0;
}
