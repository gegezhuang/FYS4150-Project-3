#pragma GCC optimize("O2,unroll-loops")
#pragma GCC target("avx,avx2,sse,sse2,ssse3,sse4,mmx,popcnt,fma")
#include <armadillo>
#include <iostream>
#include "tests/testing.hpp"
#include "project3/config.hpp"
#include "project3/Particle.hpp"
#include "project3/PenningTrap.hpp"

using namespace std;
// TODO: V_0 and d only appears as V_0 / (d*d)

bool has_flag(const std::string& option, char** begin, char** end){
    return std::find(begin, end, option) != end;
}

void print_help_message() {
    cout << "Usage" << endl;
    cout << "\t./runner [flags]" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "\t-h\tShow this help message" << endl;
    cout << "\t-t\tRun all tests" << endl;
    cout << "\t-p10\tSolves problem 10 in the problems set" << endl;
}



void simulate_resonance(double f){
    double mass = 40.78; // mass of CA-ion
    double charge = 2; // charge of CA-ion

    double omega_V_start = .2;
    double omega_V_end = 2.5;

    double omega_V_stepsize = .02;
    double N = (omega_V_end - omega_V_start) / omega_V_stepsize;
    char buffer[50];
    sprintf(buffer, "output/amplitude%f.csv", f);
    string filename = buffer;
    ofstream outfile;
    outfile.open(filename);
    outfile << "omega_V,particles left" << endl;

    for (int i = 0; i < N; i++){
        double omega_V = omega_V_start + omega_V_stepsize * i; 
        PenningTrap pt(1 * 9.65e1, .0025 * 9.65e7, 500., f, omega_V, false);
        pt.add_random_particles(charge, mass, 100);
        pt.solve_RK4(500 / 1e-2, 1e-2);
        outfile << omega_V <<"," << pt.count_particles_in_region() << endl;
    }
    outfile.close();
}

void problem10(){ //todo probably give different name
    simulate_resonance(0.1);
    simulate_resonance(0.4);
    simulate_resonance(0.7);
}

int main(int argc, char * argv[]) {
    arma::arma_rng::set_seed(42);

    if (has_flag("-h", argv, argv+argc)) print_help_message();
    if (has_flag("-t", argv, argv+argc)) run_testing();
    if (has_flag("-p10",argv, argv+argc)) problem10();
    
    return 0;
}
