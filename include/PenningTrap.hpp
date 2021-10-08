#ifndef PENNINGTRAP_H
#define PENNINGTRAP_H

class PenningTrap {
    private:
        double B_0;
        double V_0;
        double d;
        std::vector<Particle> particles;

    public:
        PenningTrap(double magnetic_field_strength, double applied_potential,
         double characteristic_dimension);
        void add_particle(Particle p);
        arma::vec external_E_field(arma::vec r);
        arma::vec external_B_field(arma::vec r); // TODO: What is this?
        void evolve_RK4(double dt);
        void evolve_forward_Euler(double dt);
};

#endif