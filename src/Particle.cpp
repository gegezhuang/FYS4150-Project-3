#include <iostream>
#include <armadillo>
#include "Particle.hpp"


Particle::Particle(double charge, double mass, arma::vec position, arma::vec velocity) {
    // require r anv v to be size 3
    q = charge;
    m = mass;
    r = position;
    v = velocity;
}