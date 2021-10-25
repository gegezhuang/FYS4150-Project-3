#pragma once

#include <armadillo>
#include <vector>

class PenningTrap {
    public:
        
        /**
         * Class holding the state of the PenningTrap and methods for simulating the system
         * 
         * @param magnetic_field_strength Magnetic field strength (B_0) [unit: u / (&micro;s * e)]
         * @param applied_potential The applied potential to the elocrodes (V_0) [unit: u * &micro;m^2 / (&micro;s^2 * e) ]
         * @param characteristic_dimension The length scale for the region between the electrodes (d) [unit: micrometer - &micro;m ]
         * @param coloumb_interactions If false, the force due to coloumb interaction between particle is ignored. Default value is true
         */
        PenningTrap(double magnetic_field_strength, double applied_potential,
        double characteristic_dimension, bool coloumb_interactions = true);
        /**
        * Class holding the state of the PenningTrap and methods for simulating the system
        * 
        * @param magnetic_field_strength Magnetic field strength (B_0) [unit: u / (&micro;s * e)]
        * @param applied_potential The applied potential to the elocrodes (V_0) [unit: u * &micro;m^2 / (&micro;s^2 * e) ]
        * @param characteristic_dimension The length scale for the region between the electrodes (d) [unit: micrometer - &micro;m ]
        * @param amplitude Amplitude of the time-dependent potential term
        * @param angular_frequency Frequency of the time-dependent potential term [unit: MHz]
        * @param coloumb_interactions If False, the force due to coloumb interaction between particle is ignored. Default value is true
        */
        PenningTrap(double magnetic_field_strength, double applied_potential,
        double characteristic_dimension, double amplitude, double angular_frequency, bool coloumb_interactions = true);
        /**
        * Adds particle PenningTrap
        */
        void add_particle(Particle p);
        /**
        * Calculates strength of electrical field at specified postion in 3d space. (E(r)) 
        * 
        * 
        * @param r Positon in 3d space (x, y, z) [unit: micrometer - &micro;m]
        * @return (-2x, -2y, 4z). 
        * Or zero vecotor if position is outside PenningTrap (i.e if |r| > d) [unit: TODO: what's the unit?]
        */
        arma::vec external_E_field(arma::vec& r);
        /**
        * Calculates strength of magnetic field at specified postion in 3d space. B(r)
        * 
        * @param r Positon in 3d space (x, y, z) [unit: micrometer - &micro;m]
        * @return Calculated strength of the magneticfield field at positon (0, 0, B_0). 
        * Or zero vecotor if position is outside PenningTrap (i.e if |r| > d) [unit: u / (&micro;s * e)]
        */
        arma::vec external_B_field(arma::vec& r);
        /**
        * Force on particle i from particle j (agrees with force on particle j from particle i)
        * 
        * @param i Index of the particle - particles are stored in the order they were added
        * @param j Index of the particle - particles are stored in the order they were added
        * @return Vector form of coloum's law (&kappa; * q_i * q_j * (r_i - r_j) / |r_i - r_j|^3) (unit: N TODO: check)
        */
        arma::vec force_particle(int i, int j);
        /**
        * Calculates the total external force on the i-th particle
        * 
        * @param i Index of the particle - particles are stored in the order they were added
        * @return q_i * external_E_field(r_i) + v_i x external_B_field(r_i)
        */
        arma::vec total_force_external(int i);
        /**
        * Calculates the total force from other particles on the i-th particle. 
        * Returns zero if Coloumb interactions is to be ignored
        * 
        * @param i Index of the particle - particles are stored in the order they were added
        * @return for (j &ne; i) &sum; force_particle(i, j) or 0 if Coloumb interactions is to be ignored 
        */
        arma::vec total_force_particles(int i);
        /**
        * Calculates the total force on the i-th particle
        * 
        * @param i Index of the particle - particles are stored in the order they were added
        * @return sum of external force and force from other particles
        */
        arma::vec total_force(int i);
        /**
        * Clears previous solutions, reserves memory and calculates the solution for 
        * a system by evolving using forward euler repeatedly
        *
        * @param N The number of iterations to evolve the system
        * @param dt The length for each time step
        */
        void solve_forward_Euler(int N, double dt);
        /**
        * Clears previous solutions, reserves memory and calculates the solution for 
        * a system by evolving using RK4 repeatedly
        *
        * @param N The number of iterations to evolve the system
        * @param dt The length for each time step
        */
        void solve_RK4(int N, double dt);
        std::vector<arma::mat> get_solution();
        std::vector<double> get_time();
        /**
        * Adds n particles with random starting position and velocity
        * 
        * Random values are sampled using Gaussian distribution, and scaled with d.
        * 
        * @param charge the charge of the particles [unit: u]
        * @param mass the mass of the particle [unit: e]
        */
        void add_random_particles(double charge, double mass, int n);
        int count_particles_in_region();
    private:
        double B_0;
        double d;
        double V_0_by_d_squared;
        double t;
        double f;
        double omega_V;
        bool ci;
        std::vector<Particle> particles;
        //solution[k] contains matrix with the solution in t_k.
        //Column i in the matrix is the position of particle i in t_k.
        std::vector<arma::mat> solution;
        std::vector<double> t_sol;
        /**
        * Calculates the system one time step forward in time using Forward Euler
        *
        * @param dt The length for the time step
        */
        void evolve_RK4(double dt);
        /**
        * Calculates the system one time step forward in time using RK4
        *
        * @param dt The length for the time step
        */
        void evolve_forward_Euler(double dt);
};
