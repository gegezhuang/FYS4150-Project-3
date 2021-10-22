#include <iostream>
#include <armadillo>
#include "project3/Particle.hpp"
#include "project3/PenningTrap.hpp"


PenningTrap::PenningTrap(double magnetic_field_strength, double applied_potential,
 double characteristic_dimension, double amplitude, double angular_frequency) {
    B_0 = magnetic_field_strength;
    V_0 = applied_potential;
    d = characteristic_dimension;
    f = amplitude;
    omega_V = angular_frequency;
}

PenningTrap::PenningTrap(double magnetic_field_strength, double applied_potential,
 double characteristic_dimension) : PenningTrap(magnetic_field_strength, applied_potential, characteristic_dimension, 0, 0){
     // TODO: initialize f and omega_V in hpp?
 }


void PenningTrap::add_particle(Particle p){
    particles.push_back(p);
}

arma::vec PenningTrap::external_E_field(arma::vec& r){
    arma::vec E(3, arma::fill::zeros);
    // r.print();
    if (arma::norm(r) > d) return E; // return 0 if particle is outside PenningTrap
    E(0) = + r(0);
    E(1) = + r(1);
    E(2) = - 2.0 * r(2);
    E *= (V_0 * (1 + f * cos(omega_V * t))) / (d * d);
    // E.print();
    return E;
}


arma::vec PenningTrap::external_B_field(arma::vec& r){
    arma::vec B(3, arma::fill::zeros);
    if (arma::norm(r) > d) return B; // return 0 if particle is outside PenningTrap
    B(2) = B_0;
    return B;
}


arma::vec PenningTrap::force_particle(int i, int j){
    double coloumb_constant = 1.38935333e5;
    Particle p_i = particles[i];
    Particle p_j = particles[j];
    arma::vec dist = p_i.r - p_j.r;
    double abs_dist = arma::norm(dist);
    return coloumb_constant * p_i.q * p_j.q * (dist) / (abs_dist * abs_dist * abs_dist);
}


arma::vec PenningTrap::total_force_external(int i){
    return particles[i].q * (external_E_field(particles[i].r) 
                        + arma::cross(particles[i].v, external_B_field(particles[i].r)));
}



arma::vec PenningTrap::total_force_particles(int i){
    arma::vec force_particles(3, arma::fill::zeros);
    int n = particles.size();
    for (int j = 0; j < n; j++){
        if (i != j){
            force_particles += force_particle(i, j);
        }
    }
    return force_particles;
}


arma::vec PenningTrap::total_force(int i){
    return total_force_particles(i) + total_force_external(i);
}


void PenningTrap::evolve_forward_Euler(double dt){
    int n = particles.size();
    arma::vec a_i;
    arma::mat positions(3, n);
    arma::mat velocities(3, n);
    //Calculate v_(i+1), r_(i+1) using forward Euler
    for (int i = 0; i < n; i++){
        a_i = total_force(i) / particles[i].m;
        velocities.col(i) = particles[i].v + dt*a_i;
        positions.col(i) = particles[i].r + dt*particles[i].v;
    }
    //add to solution
    t += dt;
    t_sol.push_back(t);
    solution.push_back(positions);
    //update position in Particle objects 
    for (int i=0; i<n; i++){
        particles[i].v = velocities.col(i);
        particles[i].r = positions.col(i);
    }
}


void PenningTrap::solve_forward_Euler(int N, double dt){
    //remove previous solution and reserve memory
    solution.clear();
    t_sol.clear();
    solution.reserve(N);
    t_sol.reserve(N);
    //add initial condition to solution
    t = 0;
    evolve_forward_Euler(0.0);
    //solve ODE
    for (int i=0; i<N-1; i++){
        evolve_forward_Euler(dt);
    }
}

void PenningTrap::evolve_RK4(double dt){
    int n = particles.size();
    Particle particle = particles[0];
    arma::mat positions(3, n);
    arma::mat velocities(3, n);
    arma::mat original_positions(3, n);
    arma::mat original_velocities(3, n);
    arma::mat k1_r(3, n); arma::mat k1_v(3, n); 
    arma::mat k2_r(3, n); arma::mat k2_v(3, n); 
    arma::mat k3_r(3, n); arma::mat k3_v(3, n); 
    arma::mat k4_r(3, n); arma::mat k4_v(3, n); 
    for (int i=0; i<n; i++){
        original_positions.col(i) = particles[i].r;
        original_velocities.col(i) = particles[i].v;
    }
    //calculate k1
    for (int i=0; i<n; i++){
        k1_v.col(i) = total_force(i)/particles[i].m;
    }
    k1_r = original_velocities;
    t += dt / 2;
    //Set new values for v, r to calculate k2
    for (int i=0; i<n; i++){
        particles[i].v = original_velocities.col(i) + dt*k1_v.col(i)/2;
        particles[i].r = original_positions.col(i) + dt*k1_r.col(i)/2;
    }
    for (int i=0; i<n; i++){
        k2_v.col(i) = total_force(i)/particles[i].m;
        k2_r.col(i) = particles[i].v;
    }
    //Set new values for v, r to calculate k3
    for (int i=0; i<n; i++){
        particles[i].v = original_velocities.col(i) + dt*k2_v.col(i)/2;
        particles[i].r = original_positions.col(i) + dt*k2_r.col(i)/2;
    }
    for (int i=0; i<n; i++){
        k3_v.col(i) = total_force(i)/particles[i].m;
        k3_r.col(i) = particles[i].v;
    }
    //set new values for v, r to calculate k4
    for (int i=0; i<n; i++){
        particles[i].v = original_velocities.col(i) + dt*k3_v.col(i);
        particles[i].r = original_positions.col(i) + dt*k3_r.col(i);
    }
    t += dt / 2;
    for (int i=0; i<n; i++){
        k4_v.col(i) = total_force(i)/particles[i].m;
        k4_r.col(i) = particles[i].v;
    }
    //Compute v_(i+1) and r_(i+1) and add to solution-vector
    velocities = original_velocities + dt*(k1_v + 2*k2_v + 2*k3_v + k4_v)/6;
    positions = original_positions + dt*(k1_r + 2*k2_r + 2*k3_r + k4_r)/6;
    t_sol.push_back(t);
    solution.push_back(positions);
    //update positions + velocities in Particle objects
    for (int i=0; i<n; i++){
        particles[i].v = velocities.col(i);
        particles[i].r = positions.col(i);
    }
}

void PenningTrap::solve_RK4(int N, double dt){
    //remove previous solution and reserve memory
    solution.clear();
    t_sol.clear();
    solution.reserve(N);
    t_sol.reserve(N);
    //add initial condition to solution
    t = 0;
    evolve_RK4(0.0);
    //solve ODE
    for (int i=0; i<N-1; i++){
        evolve_RK4(dt);
    }
}

std::vector<arma::mat> PenningTrap::get_solution(){
    return solution;
}

std::vector<double> PenningTrap::get_time(){
    return t_sol;
}

/**
 * TODO: This is not used anywhere?
 * @return The number of particles inside the PenningTrap (i.e particles with |r| < d)
 */
int PenningTrap::count_particles_in_region(){
    int count = 0;
    for (int i = 0; i < particles.size(); i++)
        if (arma::norm(particles[i].r) < d) 
            count++;
    return count;
}


void PenningTrap::add_random_particles(float charge, float mass, int n){
    for (int i = 0; i < n; i++){
        arma::vec r = arma::vec(3).randn() * 0.1 * d;  // random initial position
        arma::vec v = arma::vec(3).randn() * 0.1 * d;  // random initial velocity
        Particle p(charge, mass, r, v);
        add_particle(p);
    }
}
