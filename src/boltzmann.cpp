#include "boltzmann.hpp"
#include <cmath>

/*
Linear response theory for Boltzmann transport

\sigma_{ij} = e^2 \tau \sim_{n,k} v_i^n(k) v_j^n(k) \left( -\frac{\partial f}{\partial E} \right)

where:
- \sigma_{ij} is the conductivity tensor
- e is the electron charge
- \tau is the relaxation time
- v_i^n(k) is the group velocity of band n at wavevector k in direction i
- f is the Fermi-Dirac distribution function
- E is the energy of band n at wavevector k
- \frac{\partial f}{\partial E} is the derivative of the Fermi-Dirac distribution with respect to energy

- Include anomalous velocity from Berry curvature in velocity:
- the group velocity 
- $\widetilde{\mathbf{v}}_n = \mathbf{v}_n +  \mathbf{v}^{\rm an}_n = \left(\partial_\mathbf{k}{\varepsilon_n} - e \mathbf{E} \times \mathbf{\Omega}_n \right)/\hbar$ for the \( n \)-th band is the sum of
ordinary velocity term and anomalous velocity associated with the Berry curvature $\mathbf{\Omega}_n$
where:
*/

// Fermi-Dirac distribution derivative
inline double fermi_derivative(double e, double Ef, double T,
                               bool temperature_in_kelvin, double scale_energy) {
    if (temperature_in_kelvin) {
        const double kB = 8.617333262e-5; // eV/K
        T = kB * T / scale_energy;
    }

    double beta = 1.0 / T;
    double x = beta * (e - Ef);
    double sech2 = 1.0 / std::cosh(x / 2.0);
    return -0.25 * beta * sech2 * sech2;
}



Eigen::Vector3d BoltzmannSolver::velocity(const Eigen::Vector3d& k, int band, const Eigen::Vector3d& E, double dk) const {
    
    Eigen::Vector3d v_group;

    for (int dim = 0; dim < 3; ++dim) {
        Eigen::Vector3d dk_vec = Eigen::Vector3d::Zero();
        dk_vec(dim) = dk; // Perturb in the current dimension

        Eigen::VectorXd evals_plus, evals_minus;
        Eigen::MatrixXcd evecs_dummy;

        H.eigensystem(k + dk_vec, evals_plus, evecs_dummy);
        H.eigensystem(k - dk_vec, evals_minus, evecs_dummy);

        v_group(dim) = (evals_plus(band) - evals_minus(band)) / (2.0 * dk);
    }

    // Get Berry curvature for band
    double omega_z = berryCurvatureDifferential(H, k); // You could generalize to 3D later
    Eigen::Vector3d omega(0, 0, omega_z);

    // Anomalous velocity: -E × Ω
    Eigen::Vector3d v_anomalous = -E.cross(omega);

    return v_group + v_anomalous;


}

Eigen::Matrix3d BoltzmannSolver::conductivity(double Ef, double T, const Eigen::Vector3d& Efield) {
    using Mat3 = Eigen::Matrix3d;

    
    
    
    // Initialize conductivity tensor
    Mat3 sigma = Eigen::Matrix3d::Zero();

    

    for (const auto& k : mesh.getKPoints()) {
        Eigen::VectorXd evals;
        Eigen::MatrixXcd evecs;
        H.eigensystem(k, evals, evecs);
        

        for (int band = 0; band < evals.size(); ++band) {
            double energy = evals(band);
            
            double dfde = -fermi_derivative(energy, Ef, T,
                               temperature_in_kelvin,
                               energy_scale);


            Eigen::Vector3d v = velocity(k, band, Efield);

            // Add v_i v_j weighted by -df/de
            sigma += tau * dfde * (v * v.transpose());
        }
    }

    sigma *= (1.0 / mesh.size()); // Normalize by # of k-points

    return sigma ; // Scale by relaxation time
}