#ifndef PARTICLE_H
#define PARTICLE_H

class Particle {
    public:
        double q;
        double m;
        arma::vec r;
        arma::vec v;
        Particle(double charge, double mass, arma::vec position, arma::vec velocity);
};

#endif