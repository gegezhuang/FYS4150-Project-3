#pragma once

#include <armadillo>
#include <vector>
#include "project3/Particle.hpp"

class PenningTrap {
    public:
        PenningTrap(double magnetic_field_strength, double applied_potential,
        double characteristic_dimension);
        void add_particle(Particle p);
        arma::vec external_E_field(arma::vec r);
        arma::vec external_B_field(arma::vec r);
        arma::vec force_particle(int i, int j);
        arma::vec total_force_external(int i);
        arma::vec total_force_particles(int i);
        arma::vec total_force(int i);
        //Solves ODE with forward Euler. N is number of points in solution.
        void solve_forward_Euler(int N, double dt);
        //Solves ODE with Runge Kutta 4. N is number of points in solution.
        void solve_RK4(int N, double dt);
        std::vector<arma::mat> get_solution();
        std::vector<double> get_time();
        void add_random_particles(int n);
    private:
        double B_0;
        double V_0;
        double d;
        double t;
        std::vector<Particle> particles;
        //solution[k] contains matrix with the solution in t_k.
        //Column i in the matrix is the position of particle i in t_k.
        std::vector<arma::mat> solution;
        std::vector<double> t_sol;
        //Calculate position of particles in t_(k+1) and add to solution
        void evolve_RK4(double dt);
        void evolve_forward_Euler(double dt);
        int count_particles_in_region();
};
