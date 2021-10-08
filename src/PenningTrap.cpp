#include <iostream>
#include <armadillo>
#include "project3/Particle.hpp"
#include "project3/PenningTrap.hpp"

PenningTrap::PenningTrap(double magnetic_field_strength, double applied_potential, double characteristic_dimension) {
    B_0 = magnetic_field_strength;
    V_0 = applied_potential;
    d = characteristic_dimension;
}

void PenningTrap::add_particle(Particle p){
    particles.push_back(p);
}

void PenningTrap::run() {
    std::cout << particles[0].q << std::endl;
}
