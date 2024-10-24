#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>  
#include <tuple>
#include <armadillo>
using namespace std;

struct Atom{

    int _atomic_number;
    double _x, _y, _z;

};


// Struct for a Primitive Gaussian Function
struct PrimitiveGaussian {
    double exponent;  // Alpha, exponent of the Gaussian
    double coefficient;  // dk, contraction coefficient
    double normalization;  // Normalization constant
};


// Struct for a Basis Function
struct BasisFunction {
    Atom center;   // Center of the basis function (coordinates of the atom)
    int l, m, n;   // Angular momentum quantum numbers
    vector<PrimitiveGaussian> primitives; // Collection of primitive gaussians
};


// Primitive Gaussian exponents and coefficients for H (STO-3G basis set)
vector<PrimitiveGaussian> hydrogenBasis = {
    {3.42525091, 0.15432897, 0.0},
    {0.62391373, 0.53532814, 0.0},
    {0.16885540, 0.44463454, 0.0}
};


// Primitive Gaussian exponents and coefficients for C 2s (STO-3G basis set)
vector<PrimitiveGaussian> carbon2sBasis = {
    {2.94124940, -0.09996723, 0.0},
    {0.68348310, 0.39951283, 0.0},
    {0.22228990, 0.70011547, 0.0}
};


// Primitive Gaussian exponents and coefficients for C 2p (STO-3G basis set)
vector<PrimitiveGaussian> carbon2pBasis = {
    {2.94124940, 0.15591627, 0.0},
    {0.68348310, 0.60768372, 0.0},
    {0.22228990, 0.39195739, 0.0}
};


// Compute normalization constant for a single primitive Gaussian
double computeNormalization(const PrimitiveGaussian& pg, int l, int m, int n) {
    // Normalization constant formula for Gaussian functions
    double alpha = pg.exponent;
    double L = l + m + n;
    double norm_const = pow(2 * alpha / M_PI, 0.75) * pow(4 * alpha, L / 2.0);
    return norm_const;
}


// Build basis functions for a given atom with normalization
vector<BasisFunction> buildBasisFunctions(const Atom &atom) {
    vector<BasisFunction> basisFunctions;

    // Check the atom type (currently supports only H and C)
    if (atom._atomic_number == 1) {
        // Hydrogen has 1 s orbital (l=0, m=0, n=0)
        BasisFunction h1s;
        h1s.center = atom;
        h1s.l = h1s.m = h1s.n = 0;  // s orbital

        // Set up primitives and their normalization constants
        for (auto& pg : hydrogenBasis) {
            pg.normalization = computeNormalization(pg, 0, 0, 0);
            h1s.primitives.push_back(pg);
        }
        basisFunctions.push_back(h1s);
    } else if (atom._atomic_number == 6) {
        // Carbon has 1 s orbital and 3 p orbitals
        // 1. Add the 2s orbital (l=0, m=0, n=0)
        BasisFunction c2s;
        c2s.center = atom;
        c2s.l = c2s.m = c2s.n = 0;  // s orbital
        for (auto& pg : carbon2sBasis) {
            pg.normalization = computeNormalization(pg, 0, 0, 0);
            c2s.primitives.push_back(pg);
        }
        basisFunctions.push_back(c2s);

        // 2. Add the 2px, 2py, 2pz orbitals
        BasisFunction c2px, c2py, c2pz;

        // 2px (l=1, m=0, n=0)
        c2px.center = atom;
        c2px.l = 1; c2px.m = 0; c2px.n = 0;
        for (auto& pg : carbon2pBasis) {
            pg.normalization = computeNormalization(pg, 1, 0, 0);
            c2px.primitives.push_back(pg);
        }
        basisFunctions.push_back(c2px);

        // 2py (l=0, m=1, n=0)
        c2py.center = atom;
        c2py.l = 0; c2py.m = 1; c2py.n = 0;
        for (auto& pg : carbon2pBasis) {
            pg.normalization = computeNormalization(pg, 0, 1, 0);
            c2py.primitives.push_back(pg);
        }
        basisFunctions.push_back(c2py);

        // 2pz (l=0, m=0, n=1)
        c2pz.center = atom;
        c2pz.l = 0; c2pz.m = 0; c2pz.n = 1;
        for (auto& pg : carbon2pBasis) {
            pg.normalization = computeNormalization(pg, 0, 0, 1);
            c2pz.primitives.push_back(pg);
        }
        basisFunctions.push_back(c2pz);
    }

    return basisFunctions;
}


// Primitive Gaussian overlap calculation
double primitive_overlap(const PrimitiveGaussian& pgA, const PrimitiveGaussian& pgB, const Atom& centerA, const Atom& centerB) {
    // Gaussian exponents (alpha and beta)
    double alpha = pgA.exponent;
    double beta = pgB.exponent;

    // Distance squared between the centers A and B
    double dx = centerA._x - centerB._x;
    double dy = centerA._y - centerB._y;
    double dz = centerA._z - centerB._z;
    double distance_squared = dx * dx + dy * dy + dz * dz;

    // Calculate the Gaussian overlap integral in 3D
    double prefactor = pow(M_PI / (alpha + beta), 1.5);  // (π/(α+β))^3/2
    double exponent = exp(((-alpha * beta / (alpha + beta))) * distance_squared);  // e^[-(αβ)/(α+β) * |A-B|^2]

    return prefactor * exponent;
}


// Compute overlap between two contracted basis functions
double contracted_overlap(const BasisFunction& bfA, const BasisFunction& bfB) {
    double S_mu_nu = 0.0;

    // Loop over all primitives in both basis functions
    for (const auto& pgA : bfA.primitives) {
        for (const auto& pgB : bfB.primitives) {
            // Compute unnormalized primitive overlap Skl
            double Skl = primitive_overlap(pgA, pgB, bfA.center, bfB.center);

            // Compute normalized overlap S_mu_nu using the contraction coefficients and normalization constants
            double N_A = pgA.normalization;
            double N_B = pgB.normalization;
            S_mu_nu += (pgA.coefficient * pgB.coefficient * N_A * N_B * Skl);
        }
    }

    return S_mu_nu;
}


vector<vector<double>> buildOverlapMatrix(const vector<BasisFunction>& basisFunctions) {
    int N = basisFunctions.size();
    vector<vector<double>> overlapMatrix(N, vector<double>(N, 0.0));

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) { 
            // Compute overlap between basis function i and basis function j
            double S_ij = contracted_overlap(basisFunctions[i], basisFunctions[j]);
            overlapMatrix[i][j] = S_ij;
            if (i != j) {
                overlapMatrix[j][i] = S_ij;  // Symmetry S_ij = S_ji
            }
        }
    }

    return overlapMatrix;
}


void Evaluate_Basis_and_Electrons(vector<Atom> atoms){
    int Carbon_count = 0;
    int Hydrogen_count = 0;
    int total_electrons_pairs = 0;
    int total_electrons = 0;
    int N = 0;

    for(auto atom : atoms){
        if(atom._atomic_number == 1) Hydrogen_count++;
        if(atom._atomic_number == 6) Carbon_count++;
    }

    total_electrons = (Carbon_count * 6) + (Hydrogen_count * 1);
    total_electrons_pairs = total_electrons / 2;

    if(total_electrons_pairs != int(total_electrons_pairs)){
        throw runtime_error("Error: Number of electron paris isn't an integer!");
    }

    N = (4 * Carbon_count) + Hydrogen_count;

    total_electrons = (Carbon_count * 4) + (Hydrogen_count * 1); // Valence electrons

}


// Function to compute the kinetic energy integral (simplified)
double computeKineticEnergyIntegral(const BasisFunction& bfA, const BasisFunction& bfB) {
    double K_mu_nu = 0.0;
    // This is a simplified placeholder; the actual computation would involve detailed integrals
    // For now, return a simple placeholder value
    for (const auto& pgA : bfA.primitives) {
        for (const auto& pgB : bfB.primitives) {
            double Kkl = 0.0; // Calculate the kinetic energy integral for pgA and pgB
            K_mu_nu += (pgA.coefficient * pgB.coefficient * Kkl);
        }
    }
    return K_mu_nu;
}


// Function that builds the Hamiltonian Matrix
vector<vector<double>> buildHamilMatrix(const vector<BasisFunction> & basis_funcs){
    int N = basis_funcs.size(); // Hamiltonian Matrix has the same dimensions as the Basis functions

    vector<vector<double>> Hamils(N, vector<double>(N, 0.0)); // Initialize all the values to 0 for now 

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            double overlap = contracted_overlap(basis_funcs[i], basis_funcs[j]);
            double kinetic_E = computeKineticEnergyIntegral(basis_funcs[i], basis_funcs[j]);

            Hamils[i][j] = kinetic_E;
        }
    }

    return Hamils;

}


// Function that will read in the file & create a vector of 'Atoms'
vector<Atom> Read_Atom_Para(string & filepath){

    ifstream file(filepath);
    if(!file){
        throw runtime_error("Unable to open file!"); // Program should let user know if file cant be accessed
    }

    vector<Atom> atoms; // Vector of class instances will be used throughout the program
    string line;

    getline(file, line);

    while(getline(file, line)){
        stringstream ss(line);
        
        Atom atom;
        
        ss >> atom._atomic_number >> atom._x >> atom._y >> atom._z; // Reading in all the needed info from the file
        atoms.push_back(atom);
    }

    return atoms;

}


// Function that converts standard C++ vector into an armadillo vector
arma::mat convertToArmaMatrix(const vector<vector<double>>& vec) {
    // Get the dimensions of the input vector
    size_t rows = vec.size();
    size_t cols = rows > 0 ? vec[0].size() : 0;

    // Create an Armadillo matrix with the same dimensions
    arma::mat mat(rows, cols);

    // Fill the Armadillo matrix with values from the std::vector
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            mat(i, j) = vec[i][j];
        }
    }

    return mat;
}


// Function to construct the Hamiltonian matrix for the molecule- used for H2 primarily
arma::mat constructHamiltonian( arma::mat& S, vector<tuple<string, string>>& orbitals) {
    int n = S.n_rows;
    arma::mat H(n, n, arma::fill::zeros);

    // Define the diagonal values for each atomic orbital in eV
    map<string, double> ionizationEnergies = {
        {"H", -13.6}, // Hydrogen atom
        {"C1", -21.4}, // Carbon 2s
        {"C2", -11.4}
    };

    // Empirical constant K for the off-diagonal elements
    double K = 1.75;

    // Fill in the diagonal elements
    for (int i = 0; i < n; i++) {
        // std::string atom_orbital = orbitals[i].first + "_" + orbitals[i].second;
        string atom_orbital = get<0>(orbitals[i]) + "_" + get<1>(orbitals[i]);

    // Check if the orbital is known
    if (ionizationEnergies.find(get<0>(orbitals[i])) == ionizationEnergies.end()) {
        throw runtime_error("Unknown orbital: " + get<0>(orbitals[i]));
    }

        H(i, i) = ionizationEnergies[get<0>(orbitals[i])];
    }

    // Fill in the off-diagonal elements
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            // Get the average of the ionization energies
            double avg_ionization_energy = 0.5 * (ionizationEnergies[get<0>(orbitals[i])] + ionizationEnergies[get<0>(orbitals[j])]);

            // Compute the off-diagonal Hamiltonian element
            H(i, j) = K * avg_ionization_energy * S(i, j);
            H(j, i) = H(i, j); // Since H is symmetric
        }
    }

    return H;
}


// Function to construct Hamiltonian matrix- for any molecule more complex than H2
arma::mat constructHamiltonian_complex(arma::mat& S, vector<tuple<string, string>>& orbitals) {
    int n = S.n_rows;
     if (n != orbitals.size()) {
        // This is the make sure the overlap matrix & amount of orbitals are the same
        throw runtime_error("Size mismatch: Overlap matrix size (" + to_string(n) + 
                            ") does not match number of orbitals (" + to_string(orbitals.size()) + ")");
    }

    arma::mat H(n, n, arma::fill::zeros);

    // Define the diagonal values for each atomic orbital in eV
    map<string, double> ionizationEnergies = {
        {"C1", -21.4}, // Carbon 2s
        {"C2", -11.4},
        {"H", -13.6}
    };


    // Empirical constant K for the off-diagonal elements
    double K = 1.75;

    // Fill in the diagonal elements
    for (int i = 0; i < n; i++) {
    string atom_symbol = get<0>(orbitals[i]);

    if (atom_symbol.empty()) {
        throw runtime_error("Atom symbol is empty for index: " + to_string(i));
    }

    // Check if the orbital is known in the map
    if (ionizationEnergies.find(atom_symbol) == ionizationEnergies.end()) {
        throw runtime_error("Unknown orbital: '" + atom_symbol + "'"); // Program should throw error if an uknown/incorrect orbital is found 
    }

    // Assign diagonal elements in Hamiltonian matrix H
    H(i, i) = ionizationEnergies[atom_symbol];
}

    // Fill in the off-diagonal elements
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            // Get the average of the ionization energies
            double avg_ionization_energy = 0.5 * (ionizationEnergies[get<0>(orbitals[i])] + ionizationEnergies[get<0>(orbitals[j])]);

            // Compute the off-diagonal Hamiltonian element
            H(i, j) = K * avg_ionization_energy * S(i, j);
            H(j, i) = H(i, j); // Since H is symmetric
        }
    }

    return H;
}


// Function to solve the generalized eigenvalue problem: HC = SCΛ
void solveGeneralizedEigenProblem(const arma::mat& H, const arma::mat& S, arma::vec& eigenvalues, arma::mat& MO_coefficients, arma::mat& X) {
    // Compute S^(-1/2)
    arma::vec S_eigenvalues;
    arma::mat S_eigenvectors;
    arma::eig_sym(S_eigenvalues, S_eigenvectors, S); // Eigen decomposition of S

    // Form the inverse square root of the overlap matrix
    arma::mat S_inv_sqrt = S_eigenvectors * arma::diagmat(1.0 / arma::sqrt(S_eigenvalues)) * S_eigenvectors.t();
    X = S_inv_sqrt;

    // Transform Hamiltonian
    arma::mat H_orthog = S_inv_sqrt.t() * H * S_inv_sqrt;

    // Diagonalize the transformed Hamiltonian
    arma::eig_sym(eigenvalues, MO_coefficients, H_orthog);

    // Compute molecular orbital coefficients in the original basis
    MO_coefficients = S_inv_sqrt * MO_coefficients;
}


// Function that will return the X, MO overlap, MO Coefficients vector
void solveGeneralizedEigenProblem_complex(const arma::mat& H, const arma::mat& S, arma::vec& eigenvalues, arma::mat& MO_coefficients, arma::mat& X) {
    // Ensure symmetry of S and H
    arma::mat H_sym = 0.5 * (H + H.t());
    arma::mat S_sym = 0.5 * (S + S.t());

    // Compute S^(-1/2)
    arma::vec S_eigenvalues;
    arma::mat S_eigenvectors;
    arma::eig_sym(S_eigenvalues, S_eigenvectors, S_sym); // Eigen decomposition of S

    // Handle small eigenvalues to avoid numerical issues
    arma::vec S_eigenvalues_reg = S_eigenvalues;
    for (size_t i = 0; i < S_eigenvalues_reg.n_elem; ++i) {
        if (S_eigenvalues_reg(i) < 1e-6) {
            S_eigenvalues_reg(i) = 1e-6; // Regularize small eigenvalues
        }
    }

    // Form the inverse square root of the overlap matrix
    arma::mat S_inv_sqrt = S_eigenvectors * arma::diagmat(1.0 / arma::sqrt(S_eigenvalues_reg)) * S_eigenvectors.t();
    X = S_inv_sqrt;

    // Transform Hamiltonian
    arma::mat H_orthog = S_inv_sqrt.t() * H_sym * S_inv_sqrt;

    // Diagonalize the transformed Hamiltonian
    arma::eig_sym(eigenvalues, MO_coefficients, H_orthog);

    // Compute molecular orbital coefficients in the original basis
    MO_coefficients = S_inv_sqrt * MO_coefficients;
}

// Function to compute the total energy based on the lowest n eigenvalues
double computeTotalEnergy(const arma::vec& eigenvalues, int num_electrons) {
    double total_energy = 0.0;

    // Sum the lowest n/2 eigenvalues (each orbital contributes 2 electrons)
    for (int i = 0; i < num_electrons / 2; ++i) {
        total_energy += 2 * eigenvalues(i);
    }

    return total_energy;
}