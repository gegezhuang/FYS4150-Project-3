#include <iostream>
#include <armadillo>
#include "project3/Particle.hpp"

/**
 * Particle class keeping state of particle
 * 
 * @param charge Charge of the particle (q) [unit: the elementary charge - e]
 * @param mass Mass of the particle (m) [unit: Atomic mass unit - u]
 * @param position Initial postion of the particle in 3d space (r) [unit: micrometer - &micro;m]
 * @param velocity Initial velocity of the particle in 3d space (v) [unit: metetes per second - m/s]
 */
Particle::Particle(double charge, double mass, arma::vec position, arma::vec velocity) {
    if (position.size() != 3 || velocity.size() != 3)
        throw std::invalid_argument("Position and velocity needs to be 3-dimensional vectors");
    q = charge;
    m = mass;
    r = position;
    v = velocity;
}
