#pragma once

#include <armadillo>

class Particle {
    public:
    /**
    * Particle class keeping state of particle
    * 
    * @param charge Charge of the particle (q) [unit: the elementary charge - e]
    * @param mass Mass of the particle (m) [unit: Atomic mass unit - u]
    * @param position Initial postion of the particle in 3d space (r) [unit: micrometer - &micro;m]
    * @param velocity Initial velocity of the particle in 3d space (v) [unit: metetes per second - m/s]
    */
        Particle(double charge, double mass, arma::vec position, arma::vec velocity);
    private:
        double q;
        double m;
        arma::vec r;
        arma::vec v;

        friend class PenningTrap;
};
