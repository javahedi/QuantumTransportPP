#include "boltzmann.hpp"
#include <cmath>
#include <iostream>


/// @brief Fermi-Dirac distribution function
/// @param e: energy level
/// @param Ef: Fermi energy
/// @param T: temperature (in eV if temperature_in_kelvin is false)
/// @param temperature_in_kelvin: if true, T is interpreted in Kelvin
/// @param scale_energy: energy unit (e.g., t1), used when converting k_B T 
inline double fermi_derivative(double e, double Ef, double T,
                               bool temperature_in_kelvin, double scale_energy) {
    if (temperature_in_kelvin) {
        const double kB = 8.617333262e-5; // eV/K
        T = kB * T / scale_energy;
    }

   
    double x = (e - Ef) / T;
    double sech2 = 1.0 / std::cosh(x / 2.0);
    return -0.25  * sech2 * sech2 / T;
};


Eigen::Matrix3d BoltzmannSolver::secondDerivatives(const Eigen::Vector3d& k, 
                                                    int band, double dk) const {

    Eigen::Matrix3d hessian = Eigen::Matrix3d::Zero();
    Eigen::Vector3d dki, dkj;

    for (int i = 0; i < 3; ++i) {
        for (int j = i; j < 3; ++j) {
            dki = Eigen::Vector3d::Zero();
            dkj = Eigen::Vector3d::Zero();
            dki(i) = dk;
            dkj(j) = dk;

            double epp, epm, emp, emm;

            H.eigensystem(k + dki + dkj, evals_plus, evecs_dummy);
            epp = evals_plus(band);
            H.eigensystem(k + dki - dkj, evals_plus, evecs_dummy);
            epm = evals_plus(band);
            H.eigensystem(k - dki + dkj, evals_plus, evecs_dummy);
            emp = evals_plus(band);
            H.eigensystem(k - dki - dkj, evals_plus, evecs_dummy);
            emm = evals_plus(band);

            double second_derivative = (epp - epm - emp + emm) / (4.0 * dk * dk);
            hessian(i, j) = second_derivative;
            hessian(j, i) = second_derivative;
        }
    }

    return hessian;
};



/// @brief Calculate the group velocity and anomalous velocity at a given k-point
/// @param k k-point in reciprocal space
/// @param band band index (0 for lowest band) 
/// @param E electric field vector
/// @param B magnetic field vector
/// @param dk small displacement in k-space for numerical differentiation
/// @return    VelocityResult containing group velocity and phase-space factor
VelocityResult BoltzmannSolver::velocity(
                    double energy, double Ef, double T,
                    const Eigen::Vector3d& k, int band,
                    const Eigen::Vector3d& gradT,
                    const Eigen::Vector3d& E, 
                    const Eigen::Vector3d& B,  
                    double dk) const {
    
    Eigen::Vector3d v_group;

    for (int dim = 0; dim < 3; ++dim) {
        Eigen::Vector3d dk_vec = Eigen::Vector3d::Zero();
        dk_vec(dim) = dk; // Perturb in the current dimension

        //Eigen::VectorXd evals_plus, evals_minus;
        //Eigen::MatrixXcd evecs_dummy;

        H.eigensystem(k + dk_vec, evals_plus, evecs_dummy);
        H.eigensystem(k - dk_vec, evals_minus, evecs_dummy);

        v_group(dim) = (evals_plus(band) - evals_minus(band)) / (2.0 * dk);
    }

    // Get Berry curvature for band
    const double omega_z = berryCurvatureFHS(H, k, 1e-3, band); 
    const Eigen::Vector3d omega(0, 0, omega_z);

    // 3. Phase-space factor (D_n = 1 + B·Ω_n)
    double D_n = 1.0 + B.dot(omega);
    if (D_n <= 0.0) D_n = 1.0; // Fallback if unphysical

    // Anomalous velocity: -e/ħ (E × Ω)
    const Eigen::Vector3d v_anomalous = -E.cross(omega) / D_n;

    // Lorentz force correction: -e/ħ (v × B) × Ω
    const Eigen::Vector3d v_lorentz = -v_group.cross(B).cross(omega) / D_n;

    // Temperature gradient anomalous velocity
    const Eigen::Vector3d v_gradT = ((energy - Ef) / T) * gradT.cross(omega);


    VelocityResult result;
    result.velocity = (v_group + v_anomalous + v_lorentz + v_gradT) / D_n;
    result.phaseSpaceFactor = D_n;
    return result;
};



/// @brief Compute transport tensors for conductivity and thermopower
/// @param Ef fermi energy
/// @param T temperature (in Kelvin if temperature_in_kelvin is true)
/// @param Efield 
/// @param Bfield 
/// @return Pair of conductivity tensor and thermopower tensor
std::tuple<Eigen::Matrix3d, Eigen::Matrix3d> 
BoltzmannSolver::computeTransportTensors(double Ef, double T, 
                                        const Eigen::Vector3d& gradT,
                                        const Eigen::Vector3d& Efield,
                                        const Eigen::Vector3d& Bfield) const {
    
    Eigen::Matrix3d sigma = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d alpha = Eigen::Matrix3d::Zero();
    const double norm = 1.0 / mesh.size();
    const double e = 1.0;           // unit charge (in units of e)
    const double hbar = 1.0;        // assume ℏ = 1 in natural units or scale appropriately
    const double dk = 1e-2;  // step size in momentum space for second derivatives


    for (const auto& k : mesh.getKPoints()) {
        //Eigen::VectorXd evals;
        //Eigen::MatrixXcd evecs;
        H.eigensystem(k, evals, evecs);
        
        for (int band = 0; band < evals.size(); ++band) {
            const double energy       = evals(band);
            const double dfde         = -fermi_derivative(energy, Ef, T, 
                                                    temperature_in_kelvin, 
                                                    energy_scale);

            VelocityResult vres = velocity(energy, Ef, T, k, band, gradT, Efield, Bfield); // <-- returns velocity & D_n
            Eigen::Vector3d v = vres.velocity;
            double D_n = vres.phaseSpaceFactor;

            Eigen::Matrix3d vvT = v * v.transpose(); // Reuse for both tensors

        
           // σ_ij = e²τ ∑ D_n (-∂f/∂E) v_i v_j
            sigma += tau * D_n * dfde  * vvT;

            // α_ij = -eτ ∑ D_n (-∂f/∂E) (E-Ef)/T v_i v_j
            alpha += tau * D_n * dfde  * ((energy - Ef) / T) * vvT;

            /*
            // 2. B·correction terms (second order in τ)
            if (Bfield.norm() > 1e-12) {
                Eigen::Matrix3d hess = secondDerivatives(k, band, dk);

                // Assuming B = (0, 0, Bz) only
                double Bz = Bfield.z();

                Eigen::Vector3d C;
                C.x() = -v.x() * hess(0,1) + v.y() * hess(0,0);  // -vx ∂xy ε + vy ∂xx ε
                C.y() = -v.x() * hess(0,1) + v.y() * hess(0,0);  // same as x in 2D symmetry
                C.z() = 0.0;

                double prefactor = -std::pow(e, 3) * std::pow(tau, 2) * Bz / std::pow(hbar, 3);
                sigma += prefactor * dfde * (v * C.transpose());

                double energyShift = (energy - Ef) / T;
                double alphaPrefactor = prefactor * energyShift;
                alpha += alphaPrefactor * dfde * (v * C.transpose());
            }
            */
        }
    }

    
    return {sigma * norm, alpha * norm};
}