#include <vector>
#include <array>
#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <cmath>
#include <armadillo>
#include <algorithm>
#include <omp.h>

#include "timer.h"

using namespace std;

// basic atom information
struct atom{
    array<double, 3> coords;
    int element;
};

struct basis{
    array<double, 3> a_k;
    array<double, 3> d_k;
};

struct basis_function{
    array<double, 3> center;
    array<double, 3> quantum_nums;

    array<double, 3> a_k;
    array<double, 3> d_k;
    array<double, 3> N_k;

    int atom_index;
    int orbital_type;
};

class CNDO2{
    private:
        vector<atom> atoms; // atoms in the molecule
        int N, electrons, num_atoms, charge; // molecule properties
                                             // N is the number of molecular orbitals
        int p, q; // alpha and beta electrons
        bool outside_tolerance;
        int iterations;

        vector<int> s_basises;
        map<int, vector<basis>> basises; // basis parameters for each atom type

        vector<basis_function> basis_functions; // basis functions for each orbital in order

        vector<double> total_density; // total density calculated from density matrices

        arma::mat overlap_M, gamma_ab, Pa, Pb, Pa_old, Pb_old, fock_A, fock_B, Ca, Cb;
        arma::mat core_hamiltonian;
        arma::vec Ea, Eb;

        double total_energy;

    public:
        CNDO2(string file);

        void read_input(string file);

        void read_basis();
        void construct_basis_functions();
        basis_function populate_basis_function(int atom_num, array<double, 3> quantum_nums);
        
        void calculate_N();
        void initialize_matrices();

        void calc_overlap_matrix();
        void calc_gamma();

        void calc_total_density();

        void calculate_fock_all();

        void execute_SCF();
        void calc_new_P();

        void calc_energies();
        void calc_hamiltonian();

        void print();
        void print_SCF_step();
        void print_fock();
};