#include "kubo.hpp"
#include <cmath>
#include <complex>

KuboSolver::KuboSolver(const Hamiltonian& H, const Mesh& mesh, double eta,
                       bool temperature_in_kelvin, double energy_scale)
    : H(H), mesh(mesh), eta(eta),
      temperature_in_kelvin(temperature_in_kelvin),
      energy_scale(energy_scale) {}

static inline double fermi(double e, double Ef, double beta) {
    double x = beta * (e - Ef);
    return 1.0 / (std::exp(x) + 1.0);
}

Eigen::Matrix3d KuboSolver::conductivity(double Ef, double T) {
    using namespace Eigen;
    using std::complex;

    const double kB = 8.617333262e-5; // eV/K
    double beta = 1.0 / (temperature_in_kelvin ? (kB * T) : T); // Unitless if T not in K

    Matrix3d sigma = Matrix3d::Zero();
    VectorXd evals;
    MatrixXcd evecs;
    const double hbar = 1.0; // Unitless
    const double prefactor = 2.0 * M_PI; // 2πe²/hbar omitted for unitless

    for (const auto& k : mesh.getKPoints()) {
        H.eigensystem(k, evals, evecs);
        int N = evals.size();

        // Compute velocity operator matrices numerically
        std::array<VectorXcd, 3> vel;
        for (int i = 0; i < 3; ++i) {
            Vector3d dk = Vector3d::Zero();
            dk(i) = 1e-5;
            MatrixXcd H_plus = H.Hk(k + dk);
            MatrixXcd H_minus = H.Hk(k - dk);
            MatrixXcd dH = (H_plus - H_minus) / (2.0 * dk(i));

            vel[i] = VectorXcd(N);
            for (int m = 0; m < N; ++m) {
                vel[i](m) = (evecs.col(m).adjoint() * dH * evecs.col(m))(0, 0);
            }
        }

        // Double sum over bands
        for (int n = 0; n < N; ++n) {
            for (int m = 0; m < N; ++m) {
                if (n == m) continue;
                double deltaE = evals[n] - evals[m];
                double f_diff = fermi(evals[n], Ef, beta) - fermi(evals[m], Ef, beta);
                double denom = deltaE * deltaE + eta * eta;

                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        complex<double> vnm_i = (evecs.col(n).adjoint() * (evecs.col(m) * vel[i](m).real()))(0, 0);
                        complex<double> vmn_j = (evecs.col(m).adjoint() * (evecs.col(n) * vel[j](n).real()))(0, 0);
                        sigma(i, j) += (f_diff * std::imag(vnm_i * vmn_j)) / denom;
                    }
                }
            }
        }
    }

    sigma *= (prefactor / mesh.size()) * energy_scale * energy_scale;
    return sigma;
}
