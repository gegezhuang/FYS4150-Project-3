#include <iostream>
#include <armadillo>
#include "Particle.hpp"
#include "PenningTrap.hpp"

PenningTrap::PenningTrap(double magnetic_field_strength, double applied_potential,
        double characteristic_dimension) {
            B_0 = magnetic_field_strength;
            V_0 = applied_potential;
            d = characteristic_dimension;
        }

void PenningTrap::add_particle(Particle p){
    particles.push_back(p);
}

