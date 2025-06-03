#include "hamiltonian.hpp"


void Hamiltonian::eigensystem(const Vec& k, Eigen::VectorXd& evals, Eigen::MatrixXcd& evecs) const {
    // Compute the Hamiltonian matrix at wavevector k
    Mat Hk_matrix = Hk(k);
    
    // Use Eigen's SelfAdjointEigenSolver for Hermitian matrices
    Eigen::SelfAdjointEigenSolver<Mat> solver(Hk_matrix);
    
    if (solver.info() != 0) {
        throw std::runtime_error("Eigensystem computation failed");
    }
    
    evals = solver.eigenvalues();
    evecs = solver.eigenvectors();
}