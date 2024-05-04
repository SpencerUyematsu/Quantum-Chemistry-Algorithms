#include "CNDO2.h"

// --------------------------------------------------------------------------------------------------------------------------------------
// PROGRAM PARAMETERS
// --------------------------------------------------------------------------------------------------------------------------------------
map<char, int> elements = {
    {'H', 1},
    {'C', 6},
    {'N', 7},
    {'O', 8},
    {'F', 9},
    {'P', 15}
};

map<int, int> valence_electrons = {
    {1, 1},
    {6, 4},
    {7, 5},
    {8, 6},
    {9, 7},
    {15, 5}
};

// Empirical parameters 
map<int, string> basis_function_files = {
    {1, "basis/H_STO3G.txt"},
    {6, "basis/C_STO3G.txt"},
    {7, "basis/N_STO3G.txt"},
    {8, "basis/O_STO3G.txt"},
    {9, "basis/F_STO3G.txt"},
    {15, "basis/P_STO3G.txt"}
};

map<int, double> I_A_s = {
    {1, -7.176},
    {6, -14.051},
    {7, -19.316},
    {8, -25.390},
    {9, -32.272},
    {15, -14.033}
};

map<int, double> I_A_p = {
    {6, -5.572},
    {7, -7.275},
    {8, -9.111},
    {9, -11.08},
    {15, -5.464}
};

map<int, double> I_A_d = {
    {15, 0.5}
};

map<int, double> Beta = {
    {1, -9},
    {6, -21},
    {7, -25},
    {8, -31},
    {9, -39},
    {15, -15.07}
};
// --------------------------------------------------------------------------------------------------------------------------
// OVERLAP CALCULATIONS
//     This section contains the functions used to calculate the 3D overlap integral of primitive Gaussians
// --------------------------------------------------------------------------------------------------------------------------
int double_factorial(int n){
    int result;
    result = 1;
    for(int i = n; i > 0; i -= 2){
        result *= i;
    }
    return result;
}

int single_factorial(int n){
    int result;
    result = 1;
    for(int i = n; i > 0; i--){
        result *= i;
    }
    return result;
}

double integer_prefactor(int m, int n){
    return (single_factorial(m)/(single_factorial(n)*(single_factorial(m-n))));
}

// Center of the product
std::vector<double> center(basis_function& funct1, basis_function& funct2, int k, int l){
    std::vector<double> Rp_values;
    double Rp;
    double divisor = funct1.a_k[k] + funct2.a_k[l];
    for(int i = 0; i < 3; i++){
        Rp = funct1.a_k[k] * funct1.center[i] + funct2.a_k[l] * funct2.center[i];
        Rp_values.push_back(Rp / divisor);
    }
    return Rp_values;
}

// Exponential prefactor specific to the integration.
double exponential_prefactor(basis_function& funct1, basis_function& funct2, int k, int l, int i){
    double result;
    double exponent_intermediate;
    double root_intermediate;

    exponent_intermediate = -1 * (funct1.a_k[k] * funct2.a_k[l]);
    exponent_intermediate *= pow((funct1.center[i] - funct2.center[i]), 2);
    exponent_intermediate /= (funct1.a_k[k] + funct2.a_k[l]);
    exponent_intermediate = exp(exponent_intermediate);

    root_intermediate = sqrt(M_PI / (funct1.a_k[k] + funct2.a_k[l]));
    
    result = exponent_intermediate * root_intermediate;

    return result;
}

// Overall 3D integral calculation
double integrate_analytical(basis_function& funct1, basis_function& funct2, int k, int l){
    double e_prefactor;
    double summation;

    std::vector<double> Rp = center(funct1, funct2, k, l);
    double i_prefactor, numerator, denom;
    double result = 1;

    for(int dim = 0; dim < 3; dim++){
        summation = 0;
        e_prefactor = exponential_prefactor(funct1, funct2, k, l, dim);
        for(int i = 0; i <= funct1.quantum_nums[dim]; i++){
            for(int j = 0; j <= funct2.quantum_nums[dim]; j++){
                if((i + j) % 2 == 0){
                    i_prefactor = (integer_prefactor(funct1.quantum_nums[dim], i) * integer_prefactor(funct2.quantum_nums[dim], j));
                    numerator = double_factorial(i + j -1) * pow((Rp[dim] - funct1.center[dim]), (funct1.quantum_nums[dim] - i)) * 
                                pow((Rp[dim] - funct2.center[dim]), (funct2.quantum_nums[dim] - j));
                    denom = pow(2*(funct1.a_k[k] + funct2.a_k[l]), (i + j)/2);
                    summation += (i_prefactor * numerator / denom);
                }
            }
        }
        result *= summation * e_prefactor;
    }
    return result;
}

// Euclidean distance between two particles
double distance(array<double, 3> coords1, array<double, 3> coords2){
    double total_distance = 0;

    for(int i = 0; i < 3; i++){
        total_distance += (coords1[i] - coords2[i]) * (coords1[i] - coords2[i]);
    }
    return sqrt(total_distance);
}

// --------------------------------------------------------------------------------------------------------------------------------------
// CNDO2
// --------------------------------------------------------------------------------------------------------------------------------------

// Constructor (facilitates the algorithm's calculations)
CNDO2::CNDO2(string file){
    // reads and sets up molecule from an input file
    read_input(file);

    // read and construct basis functions to define Gaussians for each orbital 
    read_basis();
    construct_basis_functions();

    // normalization constants
    calculate_N();

    initialize_matrices();

    calc_overlap_matrix();

    // calculate gamma matrix representing electron repulsion
    calc_gamma();

    // Atomic Orbital basis total density matrix. Initial guess is 0.
    calc_total_density();

    // print CNDO/2 calculation setup
    print();

    // initialize variables associated with the SCF algorithm
    outside_tolerance = true;
    iterations = 0;

    // set number of OpenMP threads to use
    omp_set_num_threads(8);

    // iterative portion
    #pragma omp parallel
    {
        while(outside_tolerance){
            calculate_fock_all();
            execute_SCF();            
        }
    }

    // extract final results
    calc_energies();
}

// Reads molecules from input file
void CNDO2::read_input(string file){
    N = 0;
    electrons = 0;

    ifstream input(file);

    if(!input.is_open()){
        throw std::invalid_argument("File cannot be found.");
    }

    input >> num_atoms;
    atoms.resize(num_atoms);
    input >> charge;

    for(int i = 0; i < num_atoms; i++){
        input >> atoms[i].element;

        electrons += valence_electrons[atoms[i].element];

        if(atoms[i].element == 1) N += 1;
        else N += 4;

        input >> atoms[i].coords[0] >> atoms[i].coords[1] >> atoms[i].coords[2];
    }

    electrons -= charge;

    p = (electrons / 2) + (electrons % 2);
    q = electrons / 2;
}

// Reads parameters from basis files
void CNDO2::read_basis(){
    double exponent;
    for(auto it: basis_function_files){
        basises[it.first];
        if(it.first == 1) basises[it.first].resize(1);
        else basises[it.first].resize(2);

        ifstream basis_file(it.second);

        for(int i = 0; i < 3; i++){
            basis_file >> exponent;
            for(int j = 0; j < basises[it.first].size(); j++){
                basises[it.first][j].a_k[i] = exponent;
                basis_file >> basises[it.first][j].d_k[i];
            }
        }
    }
}

// constructs the basis functions corresponding to each atom in the molecule
void CNDO2::construct_basis_functions(){
    array<double, 3> temp_quantum_nums = {0,0,0};

    for(int i = 0; i < atoms.size(); i++){
        basis_functions.push_back(populate_basis_function(i, temp_quantum_nums));
        s_basises.push_back(basis_functions.size() - 1);

        // orbitals with Quantum number > 1
        if(atoms[i].element != 1){
            for(int j = 0; j < 3; j++){
                temp_quantum_nums[j] = 1;
                basis_functions.push_back(populate_basis_function(i, temp_quantum_nums));
                temp_quantum_nums[j] = 0;
            }
        }
    }
}

// subroutine used to add the appropriate parameters to each basis function
basis_function CNDO2::populate_basis_function(int atom_num, array<double, 3> quantum_nums){
    basis_function temp_basis_function;
    int orbital;

    temp_basis_function.center = atoms[atom_num].coords;
    temp_basis_function.quantum_nums = quantum_nums;

    if(quantum_nums[0] == 0 && 
       quantum_nums[1] == 0 && 
       quantum_nums[2] == 0) {
       orbital = 0;
       temp_basis_function.orbital_type = 0;
    }
    else {
        orbital = 1;
        temp_basis_function.orbital_type = 1;
    }

    temp_basis_function.a_k = basises[atoms[atom_num].element][orbital].a_k;
    temp_basis_function.d_k = basises[atoms[atom_num].element][orbital].d_k;
    temp_basis_function.atom_index = atom_num;

    return temp_basis_function;
}

// calculate normalization constants corresponding to each basis function
void CNDO2::calculate_N(){
    double S;      

    for(int j = 0; j < basis_functions.size(); j++){      
        for(int i = 0; i < 3; i++){   
            S = integrate_analytical(basis_functions[j], basis_functions[j], i, i);   
            basis_functions[j].N_k[i] = sqrt(1 / S);
        }     
    }     
}

// overlap matrix between all orbitals
void CNDO2::calc_overlap_matrix(){
    double overlap_value;     

    for(int i = 0; i < N; i++){   
        for(int j = i; j < N; j++){   
            overlap_value = 0; 

            // all combinations of dimensions
            for(int k = 0; k < 3; k++){   
                for(int l = 0; l < 3; l++){   
                    overlap_value += basis_functions[i].d_k[k] * basis_functions[j].d_k[l] *  
                                     basis_functions[i].N_k[k] * basis_functions[j].N_k[l] * 
                                     integrate_analytical(basis_functions[i], basis_functions[j], k, l);     
                }     
            }     
            overlap_M(i,j) = overlap_value;   

            if(i != j){   
                overlap_M(j,i) = overlap_value;   
            }     
        }     
    }     
}    

// calculate electron-electron repulsions
void CNDO2::calc_gamma(){
    double sigma_A, sigma_B, V_2, U_A, U_B, T, gamma_value, boys_value, coefficients;
    basis_function A_basis, B_basis;

    for(int A = 0; A < num_atoms; A++){
        A_basis = basis_functions[s_basises[A]];
        for(int B = 0; B < num_atoms; B++){
            B_basis = basis_functions[s_basises[B]];

            gamma_value = 0;

            for(int k = 0; k < 3; k++){
                for(int k_ = 0; k_ < 3; k_++){
                    for(int l = 0; l < 3; l++){
                        for(int l_ = 0; l_ < 3; l_++){
                            sigma_A = 1 / (A_basis.a_k[k] + A_basis.a_k[k_]);
                            sigma_B = 1 / (B_basis.a_k[l] + B_basis.a_k[l_]);

                            V_2 = 1 / (sigma_A + sigma_B);

                            U_A = pow((M_PI * sigma_A), 1.5);
                            U_B = pow((M_PI * sigma_B), 1.5);

                            T = V_2 * distance(atoms[A].coords, atoms[B].coords) * distance(atoms[A].coords, atoms[B].coords);

                            if(A == B){
                                boys_value = U_A * U_B * sqrt(2 * V_2) * sqrt(2/ M_PI);
                            } else{
                                boys_value = U_A * U_B * sqrt(1 / pow((distance(atoms[A].coords, atoms[B].coords)), 2)) * erf(sqrt(T));
                            }

                            coefficients = A_basis.d_k[k] * A_basis.d_k[k_] * B_basis.d_k[l] * B_basis.d_k[l_];
                            coefficients *= A_basis.N_k[k] * A_basis.N_k[k_] * B_basis.N_k[l] * B_basis.N_k[l_];

                            gamma_value +=  27.211324570237 * coefficients * boys_value; 
                        }
                    }
                }
            }
            gamma_ab(A, B) = gamma_value;
        }
    }
}

// set size for required matrices
void CNDO2::initialize_matrices(){   
    overlap_M.set_size(N, N);   
    gamma_ab.set_size(num_atoms, num_atoms);
    
    Pa.set_size(N, N);
    Pa.zeros();
    Pb.set_size(N, N);
    Pb.zeros();
    Pa_old.set_size(N, N);
    Pb_old.set_size(N, N);

    total_density.resize(num_atoms);

    fock_A.set_size(N, N);
    fock_B.set_size(N, N);
}

void CNDO2::calc_total_density(){
    fill(total_density.begin(), total_density.end(), 0);
    for(int i = 0; i < N; i++){
        int index = basis_functions[i].atom_index;
        total_density[index] += Pa(i, i) + Pb(i, i);
    }
}

// fock matrices calculation with OpenMP parallelization.
void CNDO2::calculate_fock_all(){
    double I_A, onsite_A, onsite_B, offsite;
    int A_element, B_element, A, B;
    
    #pragma omp for schedule(dynamic, 18)
    for(int v = 0; v < N; v++){
        B = basis_functions[v].atom_index;
        B_element = atoms[B].element;
        for(int u = v; u < N; u++){
            if(v == u){
                if(basis_functions[u].orbital_type) I_A = I_A_p[B_element];
                else I_A = I_A_s[B_element];
                
                onsite_A = ((total_density[B] - valence_electrons[B_element]) - (Pa(u, v) - 0.5)) * gamma_ab(B, B);
                onsite_B = ((total_density[B] - valence_electrons[B_element]) - (Pb(u, v) - 0.5)) * gamma_ab(B, B);

                offsite = 0;

                for(A = 0; A < num_atoms; A++){
                    if(A != B){
                        offsite += (total_density[A] - valence_electrons[atoms[A].element]) * gamma_ab(A, B);
                    }
                }
                fock_A(u, v) = I_A + onsite_A + offsite;
                fock_B(u, v) = I_A + onsite_B + offsite;
            } else{
                A = basis_functions[u].atom_index;
                A_element = atoms[A].element;

                double temp = 0.5 * (Beta[A_element] + Beta[B_element]) * overlap_M(u, v);

                fock_A(u, v) = fock_A(v, u) = temp - Pa(u, v) * gamma_ab(A, B); 
                fock_B(u, v) = fock_B(v, u) = temp - Pb(u, v) * gamma_ab(A, B); 
            }
        }
    }
}

// execute remainder of SCF step
void CNDO2::execute_SCF(){
    #pragma omp single
    {
        if(iterations == 0){
            print_fock();
        }
        arma::eig_sym(Ea, Ca, fock_A);
        arma::eig_sym(Eb, Cb, fock_B);
        
        iterations++;
    
        calc_new_P();

        print_SCF_step();
    }
}

// new density matrices
void CNDO2::calc_new_P(){
    Pa_old = Pa;
    Pb_old = Pb;

    Pa.zeros();
    Pb.zeros();

    for(int u = 0; u < N; u++){
        for(int v = 0; v< N; v++){
            for(int i = 0; i < p; i++){
                Pa(u, v) += Ca(u, i) * Ca(v, i);
            }
            for(int i = 0; i < q; i++){
                Pb(u, v) += Cb(u, i) * Cb(v, i);
            }
        }
    }

    // check new density matrix against old matrix to determine if the algorithm has converged
    if(arma::approx_equal(Pa, Pa_old, "absdiff", 10e-6) && arma::approx_equal(Pb, Pb_old, "absdiff", 10e-6)){
        outside_tolerance = false;
    }

    calc_total_density();
}

// 
void CNDO2::calc_energies(){
    calc_hamiltonian();

    double a_value, b_value, electronic_energy, nuclear_repulsion;
    a_value = b_value = 0;

    // Calculate electron repulsion energy from calculated core hamiltonian
    for(int u = 0; u < N; u++){
        for(int v = 0; v< N; v++){
            a_value += Pa(u, v) * (core_hamiltonian(u, v) + fock_A(u, v));
            b_value += Pb(u, v) * (core_hamiltonian(u, v) + fock_B(u, v));
        }
    }
    electronic_energy = 0.5 * (a_value + b_value);

    cout << "_____________________________________________________________" << endl;
    cout << "Electron Energy: " << electronic_energy << " eV" << endl;

    nuclear_repulsion = 0;
    for(int A = 0; A < num_atoms; A++){
        for(int B = 0; B < A; B++){
            nuclear_repulsion += (valence_electrons[atoms[A].element] * valence_electrons[atoms[B].element]) / (distance(atoms[A].coords, atoms[B].coords));
        }
    }
    nuclear_repulsion *= 27.211324570237;
    cout << "Nuclear Repulsion Energy: " << nuclear_repulsion << " eV" << endl;

    total_energy = electronic_energy + nuclear_repulsion;
    cout << "Total energy for the molecule: " << total_energy << " eV" << endl;
}

// calculates core hamiltonian to extract the molecule's energy
void CNDO2::calc_hamiltonian(){
    core_hamiltonian.set_size(N, N);
    double I_A, onsite, offsite;
    int A_element, B_element, A, B;
    
    #pragma omp parallel for schedule(dynamic, 18)
    for(int v = 0; v < N; v++){
        for(int u = 0; u < N; u++){
            if(v == u){
                A = basis_functions[u].atom_index;
                A_element = atoms[A].element;

                if(basis_functions[u].orbital_type){
                    I_A = I_A_p[A_element];
                } else {
                    I_A = I_A_s[A_element];
                }
                onsite = (valence_electrons[A_element] - 0.5) * gamma_ab(A, A);

                offsite = 0;

                for(B = 0; B < num_atoms; B++){
                    if(A != B){
                        offsite += valence_electrons[atoms[B].element] * gamma_ab(A, B);
                    }
                }
                core_hamiltonian(u, v) = I_A - onsite - offsite;
            } else{
                A = basis_functions[u].atom_index;
                B = basis_functions[v].atom_index;
                A_element = atoms[A].element;
                B_element = atoms[B].element;

                core_hamiltonian(u, v) = 0.5 * (Beta[A_element] + Beta[B_element]) * overlap_M(u, v); 
            }
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// Print results and progress
// -----------------------------------------------------------------------------------------------------------------------------------------

void CNDO2::print(){
    std::cout << "Atoms in this molecule:" << std::endl;
    for(int i = 0; i < num_atoms; i++){   
        std::cout << atoms[i].element << "\t" << atoms[i].coords[0] << "\t" << atoms[i].coords[1] << "\t" << atoms[i].coords[2] << std::endl;    
    }     

    std::cout << std::endl << "Number of basis functions: " << N << std::endl;    
    std::cout << "Number of electrons: " << electrons << std::endl;   
    std::cout << "Number of electron pairs: " << electrons / 2 << std::endl << std::endl;
    std::cout << "p = " << p << endl;
    std::cout << "q = " << q << endl;

    if(num_atoms < 10){
        for(int i = 0; i < basis_functions.size(); i++){      
            std::cout << "Basis Function: " << i << std::endl;    
            std::cout << "Center: " << basis_functions[i].center[0] << ", " << basis_functions[i].center[1] << ", " << basis_functions[i].center[2] << std::endl;     
            std::cout << "Quantum Numbers: " << basis_functions[i].quantum_nums[0] << ", " << basis_functions[i].quantum_nums[1] << ", " << basis_functions[i].quantum_nums[2] << std::endl;      
            std::cout << "Exponents: " << basis_functions[i].a_k[0] << ", " << basis_functions[i].a_k[1] << ", " << basis_functions[i].a_k[2] << std::endl;   
            std::cout << "Contraction Coefficients: " << basis_functions[i].d_k[0] << ", " << basis_functions[i].d_k[1] << ", " << basis_functions[i].d_k[2] << std::endl;    
            std::cout << "N values: " << basis_functions[i].N_k[0] << ", " << basis_functions[i].N_k[1] << ", " << basis_functions[i].N_k[2] << std::endl;    
            std::cout << std::endl;   
        }
        std::cout << "Overlap Matrix: " << std::endl;     
        std::cout << overlap_M << std::endl;      

        std::cout << "Gamma: " << std::endl;     
        std::cout << gamma_ab << std::endl; 
    }
}

void CNDO2::print_fock(){
    cout << "Fa: " << endl;
    cout << fock_A << endl;

    cout << "Fb: " << endl;
    cout << fock_B << endl;
}

void CNDO2::print_SCF_step(){
    if(num_atoms < 10){
        cout << endl << "_____________________________________________________________" << endl;
        cout << "Iteration: " << iterations << endl;

        print_fock();

        cout << "After solving eigen equation: " << endl;
        cout << "Ca: " << endl;
        cout << Ca << endl;

        cout << "Cb: " << endl;
        cout << Cb << endl;

        cout << "Pa new: " << endl;
        cout << Pa << endl;

        cout << "Pb new: " << endl;
        cout << Pb << endl;

        cout << "P total: " << endl;
        for(int i = 0; i < total_density.size(); i++){
            cout << total_density[i] << endl;
        }
    }
}
