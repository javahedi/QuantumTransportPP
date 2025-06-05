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
    double beta = 1.0 / (temperature_in_kelvin ? (kB * T) : T);

    Matrix3d sigma = Matrix3d::Zero();
    VectorXd evals;
    MatrixXcd evecs;
    const double prefactor = 2.0 * M_PI; // Unitless version of prefactor

    for (const auto& k : mesh.getKPoints()) {
        H.eigensystem(k, evals, evecs);
        int N = evals.size();

        // Compute full velocity operator matrices v_i(n, m) = <n|∂H/∂k_i|m>
        std::array<MatrixXcd, 3> v;
        for (int i = 0; i < 3; ++i) {
            Vector3d dk = Vector3d::Zero();
            dk(i) = 1e-5;
            MatrixXcd dH = (H.Hk(k + dk) - H.Hk(k - dk)) / (2.0 * dk(i));
            v[i] = evecs.adjoint() * dH * evecs;  // transformed to eigenbasis
        }

        // Kubo formula sum over bands
        for (int n = 0; n < N; ++n) {
            for (int m = 0; m < N; ++m) {
                if (n == m) continue;

                double deltaE = evals[n] - evals[m];
                double f_diff = fermi(evals[n], Ef, beta) - fermi(evals[m], Ef, beta);
                double denom = deltaE * deltaE + eta * eta;

                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        complex<double> vnm_i = v[i](n, m);
                        complex<double> vmn_j = v[j](m, n);
                        sigma(i, j) += (f_diff * std::imag(vnm_i * vmn_j)) / denom;
                    }
                }
            }
        }
    }

    sigma *= (prefactor / mesh.size()) * energy_scale * energy_scale;

    // To get conductance in units of e^2/h , you could multiply by:
    const double e2_over_h = 1.0 / (2 * M_PI);  // ≈ 0.1592
    sigma *= e2_over_h;

    return sigma;
}
