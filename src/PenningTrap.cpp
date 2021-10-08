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

void PenningTrap::evolve_forward_Euler(double dt){
    int n = particles.size();
    for (int i = 0; i < n; i++)
    {
        Particle particle = particles[i];
        arma::vec F_ext = particle.q * (external_E_field(particle.r) 
                                      + cross(particle.v, external_B_field(particle.r));
        arma::vec F_particles(3, arma::fill::zeroes);
        for (int j = 0; )
        
    }
}
