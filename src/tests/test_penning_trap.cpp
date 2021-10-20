#include <iostream>
#include <armadillo>
#include "project3/Particle.hpp"
#include "project3/PenningTrap.hpp"
#include "tests/test_penning_trap.hpp"

PenningTrap test_penning_trap_creation() {
    double charge = 1, mass = 2;
    arma::vec pos = { 3.0, 0.0, 0.0 };
    arma::vec vel = { 0.0, 1.0, 1.0 };
    // PenningTrap pn (1, 2, 3);

    Particle p1 (charge, mass, pos, vel);
    Particle p2 (2*charge, mass/2, -pos, 0.5*vel);

    PenningTrap pt1(1.0, 2.0, 3.0);
    PenningTrap pt2(1.0, 2.0, 3.0, 2.0, 6.0);

    return pt1;
}

void test_penning_trap() {
    std::cout << "B" << std::endl;
    PenningTrap pt = test_penning_trap_creation();
}
