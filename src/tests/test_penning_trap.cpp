#include <assert.h>
#include <armadillo>
#include "project3/Particle.hpp"
#include "project3/PenningTrap.hpp"
#include "tests/test_penning_trap.hpp"

PenningTrap test_penning_trap_creation() {
    PenningTrap pt(1.0, 2.0, 3.0);

    return pt;
}

PenningTrap test_penning_trap_creations_all_parameters() {
    PenningTrap pt(1.0, 2.0, 3.0, 2.0, 6.0);
    return pt;
}

void test_add_particles(PenningTrap pt) {

    double charge = 1, mass = 2;

    arma::vec pos = { 1.0, 0.0, 0.0 };
    arma::vec vel = { 0.0, 1.0, 1.0 };

    assert (pt.count_particles_in_region() == 0);

    Particle p1 (charge, mass, pos, vel);
    Particle p2 (2*charge, mass/2, -pos, 0.5*vel);
    pt.add_particle(p1);
    assert (pt.count_particles_in_region() == 1);
    pt.add_particle(p2);
    assert (pt.count_particles_in_region() == 2);

    arma::vec pos3 = { 100.0, 0.0, 0.0 };
    Particle p3 (charge, mass, pos3, vel);
    pt.add_particle(p3);

    assert (pt.count_particles_in_region() == 2);
}

void test_add_random_particles(PenningTrap pt) {
    assert (pt.count_particles_in_region() == 0);
    pt.add_random_particles(5);
    assert (pt.count_particles_in_region() == 5);
    pt.add_random_particles(3);
    assert (pt.count_particles_in_region() == 8);
}

void test_penning_trap() {
    PenningTrap pt1 = test_penning_trap_creation();
    test_add_particles(pt1);

    PenningTrap pt2 = test_penning_trap_creations_all_parameters();
    test_add_random_particles(pt2);
}
