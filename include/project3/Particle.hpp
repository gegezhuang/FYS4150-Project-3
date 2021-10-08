#pragma once

#include <armadillo>

class Particle {
    public:
        Particle(double charge, double mass, arma::vec position, arma::vec velocity);
    private:
        double q;
        double m;
        arma::vec r;
        arma::vec v;

        friend class PenningTrap;
};
