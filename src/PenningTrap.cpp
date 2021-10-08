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

arma::vec PenningTrap::external_E_field(arma::vec r){
    arma::vec E(3);
    E(0) = - 2 * r(0);
    E(1) = - 2 * r(1);
    E(2) = 4 * r(2);
    return E;
}

arma::vec PenningTrap::external_B_field(arma::vec r){
    arma::vec B(3, arma::fill::zeros);
    B(2) = B_0;
    return B;
}

arma::vec PenningTrap::force_particle(int i, int j){
    double coloumb_constant = 1.38935333e5;
    Particle p_i = particles[i];
    Particle p_j = particles[j];
    arma::vec dist = p_i.r - p_j.r;
    double abs_dist = arma::norm(dist);
    return coloumb_constant * p_i.q * p_j.q * (dist) / (abs_dist * abs_dist * abs_dist);
}
