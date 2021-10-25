#pragma GCC optimize("O2,unroll-loops")
#pragma GCC target("avx,avx2,sse,sse2,ssse3,sse4,mmx,popcnt,fma")
#include <armadillo>
#include <iostream>
#include "tests/testing.hpp"
#include "project3/config.hpp"
#include "project3/Particle.hpp"
#include "project3/PenningTrap.hpp"
#include "project3/file.hpp"
#include "project3/analytical.hpp"

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

void estimate_error() {
    double hs[5] = {1e-2, 5e-2, 1e-1, 5e-1, 1.};
    string hstr[5] = {"1e-2", "5e-2", "1e-1", "5e-1", "1"};
    

        for (int i = 0; i < 5; i++){
        Particle aRK(charge, mass, pos, vel);
        Particle aFE(charge, mass, pos, vel);
        PenningTrap ptRK(B_0, V_0, d);
        PenningTrap ptFE(B_0, V_0, d);
        ptRK.add_particle(aRK);
        ptFE.add_particle(aFE);
        string filenameRK = "rk4h=" + hstr[i] +  ".csv";
        string filenameFE = "feh=" + hstr[i] + ".csv";
        string fileNameAnalytical = "analyticalh=" + hstr[i] + ".csv";
        ofstream outfileRK = get_solution_file(filenameRK);
        ofstream outfileFE = get_solution_file(filenameFE);
        ofstream outfileAnalytical = get_solution_file(fileNameAnalytical);
        double h = hs[i];
        int N = 100 / h;
        ptRK.solve_RK4(N, h);
        ptFE.solve_RK4(N, h);
        arma::vec t(N+1);
        for (auto i = 0; i < N+1; i++) {
            t(i) = ((double) i) * h;
        }
        std::vector<arma::vec> solAnalytical = solve_analytically(t, charge, mass, pos(0), pos(2), vel(1), 
                                                    B_0, V_0, d);
        vector<arma::mat> solRK = ptRK.get_solution();
        vector<arma::mat> solFE = ptFE.get_solution();
        for (int j = 0; j < N + 1; j++){
            outfileRK << solRK[j](0) << "," << solRK[j](1) << "," << solRK[j](2) << endl;
            outfileFE << solFE[j](0) << "," << solFE[j](1) << "," << solFE[j](2) << endl;
            outfileAnalytical << solAnalytical[0][j] << "," << solAnalytical[1][j] << "," << solAnalytical[2][j] << endl;
        }
        outfileRK.close();
        outfileFE.close();
        outfileAnalytical.close();
    }
}



void simulate_resonance(double f){
    double omega_V_start = .2;
    double omega_V_end = 2.5;
    double omega_V_stepsize = .02;
    int N_omega_V = (omega_V_end - omega_V_start) / omega_V_stepsize;

    char buffer[50];
    sprintf(buffer, "output/amplitude%f.csv", f);
    string filename = buffer;
    ofstream outfile;
    outfile.open(filename);
    outfile << "omega_V,particles left" << endl;

    for (int i = 0; i < N_omega_V+1; i++){
        double omega_V = omega_V_start + omega_V_stepsize * i; 
        PenningTrap pt(1 * 9.65e1, .0025 * 9.65e7, 500., f, omega_V, false);
        pt.add_random_particles(charge, mass, 100);
        pt.solve_RK4(500 / 1e-2, 1e-2);
        outfile << omega_V << "," << pt.count_particles_in_region() << endl;
    }
    outfile.close();
}

void simulate_resonance_fine_grained(){

    double omega_V_start = .29;
    double omega_V_end = .31;
    double omega_V_stepsize = .002;
    double N_omega = (omega_V_end - omega_V_start) / omega_V_stepsize;
    double f = 0.1; //note that we are only running for one amplitude f=0.1


    double T = 500;
    double h = 1e-2;
    int N = T/h;

    ofstream outfile1;
    outfile1.open("output/fine_grained_no_coulomb_interactions.csv");
    outfile1 << "omega_V,particles left" << endl;

    for (int i=0; i < N_omega+1; i++){
        double omega_V = omega_V_start + omega_V_stepsize * i; 
        PenningTrap pt(1 * 9.65e1, .0025 * 9.65e7, 500., f, omega_V, false);
        pt.add_random_particles(charge, mass, 100);
        pt.solve_RK4(N, h);
        outfile1 << omega_V << "," << pt.count_particles_in_region() << endl;
    }
    outfile1.close();

    ofstream outfile2;
    outfile2.open("output/fine_grained_with_coulomb_interactions.csv");
    outfile2 << "omega_V,particles left" << endl;

    for (int i=0; i<N_omega+1; i++){
        double omega_V = omega_V_start + omega_V_stepsize * i; 
        PenningTrap pt(1 * 9.65e1, .0025 * 9.65e7, 500., f, omega_V, true);
        pt.add_random_particles(charge, mass, 100);
        pt.solve_RK4(N, h);
        outfile2 << omega_V <<"," << pt.count_particles_in_region() << endl;
    }
    outfile2.close();
}

void problem10(){
    simulate_resonance(0.1);
    simulate_resonance(0.4);
    simulate_resonance(0.7);
}

int main(int argc, char * argv[]) {
    arma::arma_rng::set_seed(42);

    if (has_flag("-h", argv, argv+argc)) print_help_message();
    if (has_flag("-t", argv, argv+argc)) run_testing();
    if (has_flag("-p9", argv, argv+argc)) estimate_error();
    if (has_flag("-p10",argv, argv+argc)) problem10();
    if (has_flag("-fg", argv, argv+argc)) simulate_resonance_fine_grained();
    
    return 0;
}
