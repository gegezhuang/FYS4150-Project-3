#include <iostream>
#include <armadillo>
#include "project3/Particle.hpp"


Particle::Particle(double charge, double mass, arma::vec position, arma::vec velocity) {
    if (position.size() != 3 || velocity.size() != 3)
        throw std::invalid_argument("Position and velocity needs to be 3-dimensional vectors");
    q = charge;
    m = mass;
    r = position;
    v = velocity;
}
