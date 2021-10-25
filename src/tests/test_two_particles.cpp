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
#include "tests/test_two_particles.hpp"

PenningTrap initialize_penning_trap_two_particle(bool coulomb_interactions) {
    PenningTrap pt(B_0, V_0, d, coulomb_interactions);
    pt.add_random_particles(charge, mass, 2);
    return pt;
}

void test_runge_kutta_two_particles(bool coulomb_interactions) {
    PenningTrap pt = initialize_penning_trap_two_particle(coulomb_interactions);
    pt.solve_RK4(N, step_size);

    std::vector<double> time_vec = pt.get_time();
    for (long unsigned int i = 0; i < time_vec.size(); i++) {
        double t = time_vec[i];
        double t_i = ((double) i) * step_size;
        assert (1e-8 > fabs(t - t_i));
    }

    std::vector<arma::mat> sol = pt.get_solution();

    ofstream ofile;
    if (coulomb_interactions) {
        ofile = get_solution_file_double_x_y_z("runge_kutta_positons_two_particles_with_interactions.csv");
    } else {
        ofile = get_solution_file_double_x_y_z("runge_kutta_positons_two_particles_without_interactions.csv");
    }
    for (int i = 0; i < N+1; i++){
        ofile << sol[i](0) << "," << sol[i](1) << "," << sol[i](2) << "," << sol[i](3) << "," << sol[i](4) << "," << sol[i](5) << endl;
    }
    ofile.close();
}

// void test_forward_euler_two_particles(bool coulomb_interactions) {
//     PenningTrap pt = initialize_penning_trap_two_particle(coulomb_interactions);
//     pt.solve_forward_Euler(N, step_size);
//
//     std::vector<double> time_vec = pt.get_time();
//     std::vector<arma::mat> sol = pt.get_solution();
//     ofstream ofile;
//     if (coulomb_interactions) {
//         ofile = get_solution_file_double_x_y_z("forward_euler_two_particles_with_interactions.csv");
//     } else {
//         ofile = get_solution_file_double_x_y_z("forward_euler_two_particles_without_interactions.csv");
//     }
//     for (int i = 0; i < N+1; i++){
//         ofile << sol[i](0) << "," << sol[i](1) << "," << sol[i](2) << "," << sol[i](3) << "," << sol[i](4) << "," << sol[i](5) << endl;
//     }
//     ofile.close();
// }

// void test_analytical() {
//     arma::vec t(N+1);
//     for (auto i = 0; i < N+1; i++) {
//         t(i) = ((double) i) * step_size;
//     }
//
//     std::vector<arma::vec> sol = solve_analytically(t, charge, mass, pos(0), pos(2), vel(1),
//                                                     B_0, V_0, d);
//     ofstream ofile = get_solution_file("analytical_positons.csv");
//     for (int i = 0; i < N+1; i++) {
//         ofile << sol[0][i] << "," << sol[1][i] << "," << sol[2][i] << endl;
//     }
//     ofile.close();
// }

void test_two_particles() {
    arma::arma_rng::set_seed(42);
    test_runge_kutta_two_particles(false);
    arma::arma_rng::set_seed(42);
    test_runge_kutta_two_particles(true);
    // arma::arma_rng::set_seed(42);
    // test_forward_euler_two_particles(false);
    // arma::arma_rng::set_seed(42);
    // test_forward_euler_two_particles(true);
    // PenningTrap pt = initialize_penning_trap_two_particle();
    
}
