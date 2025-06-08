#include "kubo.hpp"
#include <cmath>
#include <complex>

/// @brief  KuboSolver constructor
/// @details Initializes the KuboSolver with a Hamiltonian, mesh, and parameters.
/// @param H 
/// @param mesh 
/// @param eta 
/// @param temperature_in_kelvin 
/// @param energy_scale 
KuboSolver::KuboSolver(const Hamiltonian& H, const Mesh& mesh, double eta,
                       bool temperature_in_kelvin, double energy_scale)
    : H(H), mesh(mesh), eta(eta),
      temperature_in_kelvin(temperature_in_kelvin),
      energy_scale(energy_scale) {}

static inline double fermi(double e, double Ef, double beta) {
    const double x = beta * (e - Ef);
    return 1.0 / (std::exp(x) + 1.0);
}

std::tuple<Eigen::Matrix3d, Eigen::Matrix3d, Eigen::Matrix3d> KuboSolver::computeTransportTensors(double Ef, double T) {
    using namespace Eigen;
    using std::complex;

    constexpr double kB = 8.617333262e-5; // eV/K
    const double beta = 1.0 / (temperature_in_kelvin ? (kB * T) : T);
    const double norm_factor = 1.0/ mesh.size(); // Normalization factor for averaging over k-points
    constexpr double e2_over_h = 1.0 / (2 * M_PI); // e^2/h in unitless form


    Matrix3d L0 = Matrix3d::Zero(); // Electrical conductivity
    Matrix3d L1 = Matrix3d::Zero(); // Thermoelectric
    Matrix3d L2 = Matrix3d::Zero(); // Thermal

    //VectorXd evals;
    //MatrixXcd evecs;

    for (const auto& k : mesh.getKPoints()) {
        H.eigensystem(k, evals, evecs);
        const int N = evals.size();

        // Compute full velocity operator matrices v_i(n, m) = <n|∂H/∂k_i|m>
        std::array<MatrixXcd, 3> v;
        for (int i = 0; i < 3; ++i) {
            Vector3d dk = Vector3d::Zero();
            dk(i) = 1e-5;
            // Reuse dH_buffer for finite difference calculation
            dH_buffer = (H.Hk(k + dk) - H.Hk(k - dk)) / (2.0 * dk(i));
            velocity_matrices[i] = evecs.adjoint() * dH_buffer * evecs;  // transformed to eigenbasis
        }

        // Kubo formula sum over bands
        for (int n = 0; n < N; ++n) {
            const double E_n = evals[n];
            const double f_n = fermi(E_n, Ef, beta);

            for (int m = 0; m < N; ++m) {
                if (n == m) continue;

                const double E_m = evals[m];
                const double deltaE = E_n - E_m;
                const double f_diff =f_n - fermi(E_m, Ef, beta);
                const double denom = deltaE * deltaE + eta * eta;
                const double factor = f_diff / denom;
                const double omega = (E_n + E_m) / 2.0 - Ef;

                for (int i = 0; i < 3; ++i) {
                    const complex<double> vnm_i = velocity_matrices[i](n, m);
                    for (int j = 0; j < 3; ++j) {
                        const complex<double> vmn_j = velocity_matrices[j](m, n);

                        double imagPart = std::imag(vnm_i * vmn_j);

                        L0(i, j) += factor * imagPart;
                        L1(i, j) += factor * imagPart * omega;
                        L2(i, j) += factor * imagPart * omega * omega;
                    }
                }
            }
        }
    }

    const double scaling = 2.0 * M_PI * norm_factor * energy_scale * energy_scale * e2_over_h;

    L0 *= scaling;
    L1 *= scaling;
    L2 *= scaling;

    Matrix3d sigma = L0;
    Matrix3d alpha = L1 / (T);
    Matrix3d kappa = (L2 - (L1 * L1.transpose()).cwiseQuotient(L0)) / T;


    return {sigma, alpha, kappa};

}
