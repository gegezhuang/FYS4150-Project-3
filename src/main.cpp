#include "project3/Particle.hpp"
#include "project3/PenningTrap.hpp"

int main() {
    double charge, mass; charge = 1.0; mass = 2.0;
    arma::vec pos;
    arma::vec vel;

    Particle a(charge, mass, pos, vel);
    Particle b(2*charge, mass/2, pos, vel);

    PenningTrap pt(1.0, 2.0, 3.0);
    pt.add_particle(a);
    pt.add_particle(b);

    pt.run();

    return 0;
}
