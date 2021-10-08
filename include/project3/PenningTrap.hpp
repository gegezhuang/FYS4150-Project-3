#pragma once

#include <armadillo>
#include "project3/Particle.hpp"

class PenningTrap {
    public:
        PenningTrap(double magnetic_field_strength, double applied_potential,
        double characteristic_dimension);
        void add_particle(Particle p);
        void run();
        arma::vec external_E_field(arma::vec r);
        arma::vec external_B_field(arma::vec r); // TODO: What is this?
    private:
        double B_0;
        double V_0;
        double d;
        std::vector<Particle> particles;
};

