#include <iostream>
using namespace std; // TODO: Remove

#include <assert.h>
#include <armadillo>
#include "project3/analytical.hpp"
#include "project3/Particle.hpp"
#include "project3/PenningTrap.hpp"
#include "tests/test_runge_kutta_forward_euler.hpp"
#include "tests/config.hpp"

PenningTrap initialize_penning_trap() {
    Particle a(charge, mass, pos, vel);
    Particle b(2*charge, mass/2, -pos, -vel);

    PenningTrap pt(B_0, V_0, d);
    pt.add_particle(a);
    pt.add_particle(b);

    return pt;
}

void test_runge_kutta() {
    arma::vec pos = { 1.0, 0.0, 0.0 };
    arma::vec vel = { 0.0, 1.0, 0.0 };

    PenningTrap pt = initialize_penning_trap();

    pt.solve_RK4(N, step_size);

    std::vector<double> time_vec = pt.get_time();

    for (long unsigned int i = 0; i < time_vec.size(); i++) {
        double t = time_vec[i];
        double t_i = ((double) i) * step_size;
        assert (1e-10 > abs(t - t_i));
    }


    std::vector<arma::mat> sol = pt.get_solution();

    for (int i = 0; i < 3; i++){
        sol[i].print();
        std::cout << std::endl;
    }

    for (int i = N-3; i < N; i++){
        sol[i].print();
        std::cout << std::endl;
    }
}

void test_forward_euler() {
    PenningTrap pt = initialize_penning_trap();

    pt.solve_forward_Euler(N, step_size);

    std::vector<double> time_vec = pt.get_time();
}

void test_analytical() {
    // TODO: FIX THE POSITIONS
    arma::vec pos = { 1.0, 0.0, 0.0 };
    arma::vec vel = { 0.0, 1.0, 0.0 };

    arma::vec t = { 0.0, 0.3,  0.6, 1.0 };

    std::vector<arma::vec> sol = solve_analytically(t, charge, mass, pos(0), pos(2), vel(1), B_0, V_0, d);
    arma::mat solution_matrix  = std::move(arma::join_rows( sol[0], sol[1], sol[2] ));
    cout << "        x        y        z" << endl;
    solution_matrix.print();
}

void test_runge_kutta_forward_euler() {
    test_runge_kutta();
    test_forward_euler();
    test_analytical();
}
