#include <fstream>

#include <iostream>
using namespace std; // TODO: Remove

#include <assert.h>
#include <armadillo>
#include "project3/file.hpp"
#include "project3/config.hpp"
#include "project3/Particle.hpp"
#include "project3/analytical.hpp"
#include "project3/PenningTrap.hpp"
#include "tests/test_runge_kutta_forward_euler.hpp"

PenningTrap initialize_penning_trap_one_particle() {
    arma::vec pos = { 1.0, 0.0, 0.0 };
    arma::vec vel = { 0.0, 1.0, 0.0 };

    Particle a(charge, mass, pos, vel);
    // Particle b(2*charge, mass/2, -pos, -vel);

    PenningTrap pt(B_0, V_0, d);
    pt.add_particle(a);
    // pt.add_particle(b);

    return pt;
}

void test_runge_kutta() {
    arma::vec pos = { 1.0, 0.0, 0.0 };
    arma::vec vel = { 0.0, 1.0, 0.0 };

    PenningTrap pt = initialize_penning_trap_one_particle();

    pt.solve_RK4(N, step_size);

    std::vector<double> time_vec = pt.get_time();

    for (long unsigned int i = 0; i < time_vec.size(); i++) {
        double t = time_vec[i];
        double t_i = ((double) i) * step_size;
        assert (1e-8 > fabs(t - t_i));
    }

    std::vector<arma::mat> sol = pt.get_solution();

    ofstream ofile = get_solution_file("runge_kutta_positons_one_particle.csv");
    for (int i = 0; i < N+1; i++){
        ofile << sol[i](0) << "," << sol[i](1) << "," << sol[i](2) << endl;
    }
    ofile.close();
}

/*
ofstream get_file(string name) {

    return ofile;
};

*/

void test_forward_euler() {
    PenningTrap pt = initialize_penning_trap_one_particle();

    pt.solve_forward_Euler(N, step_size);

    std::vector<double> time_vec = pt.get_time();

    std::vector<arma::mat> sol = pt.get_solution();

    ofstream ofile = get_solution_file("forward_euler_one_particle.csv");
    for (int i = 0; i < N+1; i++){
        ofile << sol[i](0) << "," << sol[i](1) << "," << sol[i](2) << endl;
    }
    ofile.close();
}

void test_analytical() {
    // TODO: FIX THE POSITIONS
    arma::vec pos = { 1.0, 0.0, 0.0 };
    arma::vec vel = { 0.0, 1.0, 0.0 };

    // arma::vec t = { 0.0, 0.3,  0.6, 1.0 };
    arma::vec t(N+1);
    for (auto i = 0; i < N+1; i++) {
        t(i) = ((double) i) * step_size;
    }

    std::vector<arma::vec> sol = solve_analytically(t, charge, mass, pos(0), pos(2), vel(1), B_0, V_0, d);

    ofstream ofile = get_solution_file("analytical_positons.csv");
    for (int i = 0; i < N+1; i++) {
        ofile << sol[0][i] << "," << sol[1][i] << "," << sol[2][i] << endl;
    }
    ofile.close();
}

void test_runge_kutta_forward_euler() {
    test_runge_kutta();
    test_forward_euler();
    test_analytical();
}
