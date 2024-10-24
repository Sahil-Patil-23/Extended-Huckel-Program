#include "Basis_Hamil_Code.hpp"

int main(){

    // Using 'H2.txt'
    string filepath2 = "../sample_input/H2.txt";
    vector<Atom> h2 = Read_Atom_Para(filepath2);


    Evaluate_Basis_and_Electrons(h2);

    vector<BasisFunction> basisFunctions_H2;
    for (const auto& atom : h2) {
        vector<BasisFunction> atomBasis = buildBasisFunctions(atom);
        basisFunctions_H2.insert(basisFunctions_H2.end(), atomBasis.begin(), atomBasis.end());
    }

    vector<vector<double>> overlapMatrix_H2 = buildOverlapMatrix(basisFunctions_H2);

    // Doing arma calcs
    arma::mat overlap_H2 = convertToArmaMatrix(overlapMatrix_H2);

    // Atomic orbitals for H2 (two H atoms, each with a 1s orbital)
    vector<tuple<string, string>> orbitals = {{"H", "1s"}, {"H", "1s"}};

    // Step 1: Construct the Hamiltonian matrix
    arma::mat H = constructHamiltonian(overlap_H2, orbitals);

    // Step 2: Solve the generalized eigenvalue problem
    arma::vec eigenvalues;
    arma::mat MO_coefficients, X;
    solveGeneralizedEigenProblem(H, overlap_H2, eigenvalues, MO_coefficients, X);

    // Step 3: Compute the total energy (H2 has 2 electrons)
    int num_electrons = 2; // H2 has 2 electrons
    double total_energy = computeTotalEnergy(eigenvalues, num_electrons);

    // Print results in desired format
    cout << "Overlap matrix:" << endl;
    overlap_H2.print();
    cout << "Hamiltonian matrix:" << endl;
    H.print();
    cout << "X_mat:" << endl;
    X.print();
    cout << "MO coefficients (C matrix):" << endl;
    MO_coefficients.print();
    cout << "MO overlap matrix:" << endl;
    arma::mat MO_overlap = MO_coefficients.t() * overlap_H2 * MO_coefficients;
    MO_overlap.print();

    // Output total energy
    cout << "The molecule in file " << filepath2 << " has energy " << setprecision(8) << total_energy << " eV" << endl;

    cout << endl << endl;


    // Using 'C2H2.txt'
    string filepath_c2h2 = "../sample_input/C2H2.txt";
    vector<Atom> c2h2 = Read_Atom_Para(filepath_c2h2);

    Evaluate_Basis_and_Electrons(c2h2);

    vector<BasisFunction> basisFunctions_C2H2;
    for (const auto& atom : c2h2) {
        vector<BasisFunction> atomBasis = buildBasisFunctions(atom);
        basisFunctions_C2H2.insert(basisFunctions_C2H2.end(), atomBasis.begin(), atomBasis.end());
    }

    vector<vector<double>> overlapMatrix_C2H2 = buildOverlapMatrix(basisFunctions_C2H2);

    arma::mat overlap_C2H2 = convertToArmaMatrix(overlapMatrix_C2H2);

    // vector<tuple<string, string>> orbs = {{"H", "1s"}, {"C1", "2s"}, {"C2", "2p"}, {"H", "1s"}};
    vector<tuple<string, string>> orbs = {
    {"C1", "2s"}, {"C1", "2px"}, {"C1", "2py"}, {"C1", "2pz"},
    {"C2", "2s"}, {"C2", "2px"}, {"C2", "2py"}, {"C2", "2pz"},
    {"H", "1s"}, {"H", "1s"}
    };

    arma::mat Hamil = constructHamiltonian_complex(overlap_C2H2, orbs);

    arma::vec eigenval_C2H2;
    arma::mat MO_C2H2, X_C2H2;
    solveGeneralizedEigenProblem_complex(Hamil, overlap_C2H2, eigenval_C2H2, MO_C2H2, X_C2H2);

    int num_elec = 5;
    double total_E = computeTotalEnergy(eigenval_C2H2, num_elec);


    cout << "Overlap matrix:" << endl;
    overlap_C2H2.print();
    cout << "Hamiltonian matrix:" << endl;
    Hamil.print();
    cout << "X_mat:" << endl;
    X_C2H2.print();
    cout << "MO coefficients (C Matrix):" << endl;
    MO_C2H2.print();
    cout << "MO overlap matrix:" << endl;
    arma::mat MO_overlap_C2H2 = MO_C2H2.t() * overlap_C2H2 * MO_C2H2;
    MO_overlap_C2H2.print();

    cout << "The molecule in file " << filepath_c2h2 << " has energy " << setprecision(8) << total_E << " eV" << endl;


    cout << endl << endl;


    // Using 'C2H4.txt'
    string filepath_c2h4 = "../sample_input/C2H4.txt";
    vector<Atom> c2h4 = Read_Atom_Para(filepath_c2h4);

    Evaluate_Basis_and_Electrons(c2h4);

    vector<BasisFunction> basisFunctions_C2H4;
    for (const auto& atom : c2h4) {
        vector<BasisFunction> atomBasis = buildBasisFunctions(atom);
        basisFunctions_C2H4.insert(basisFunctions_C2H4.end(), atomBasis.begin(), atomBasis.end());
    }

    vector<vector<double>> overlapMatrix_C2H4 = buildOverlapMatrix(basisFunctions_C2H4);

    arma::mat overlap_C2H4 = convertToArmaMatrix(overlapMatrix_C2H4);

    vector<tuple<string, string>> orbitals_c2h4 = {
    {"C1", "2s"}, {"C1", "2px"}, {"C1", "2py"}, {"C1", "2pz"},
    {"C2", "2s"}, {"C2", "2px"}, {"C2", "2py"}, {"C2", "2pz"},
    {"H", "1s"}, {"H", "1s"},  {"H", "1s"},  {"H", "1s"}
    };

    arma::mat Ham = constructHamiltonian_complex(overlap_C2H4, orbitals_c2h4);

    arma::vec eigenval_C2H4;
    arma::mat MO_C2H4, X_C2H4;
    solveGeneralizedEigenProblem_complex(Ham, overlap_C2H4, eigenval_C2H4, MO_C2H4, X_C2H4);

    int num_E = 12;
    double total_EN = computeTotalEnergy(eigenval_C2H4, num_E);

    cout << "Overlap matrix:" << endl;
    overlap_C2H4.print();
    cout << "Hamiltonian matrix:" << endl;
    Ham.print();
    cout << "X_mat:" << endl;
    X_C2H4.print();
    cout << "MO coefficients (C Matrix):" << endl;
    MO_C2H4.print();
    cout << "MO overlap matrix:" << endl;
    arma::mat MO_overlap_C2H4 = MO_C2H4.t() * overlap_C2H4 * MO_C2H4;
    MO_overlap_C2H4.print();    

    cout << "The molecule in file " << filepath_c2h4 << " has energy " << total_EN << " eV" << endl;

    return 0;
}