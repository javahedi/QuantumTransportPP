#pragma once
#include "hamiltonian.hpp"
#include <complex>
#include <cmath>

class HaldaneModel : public Hamiltonian {
public:
    double t1;     // NN hopping
    double t2;     // NNN hopping
    double phi;    // NNN complex phase
    double M;      // sublattice mass

    // Constructor with default parameters
    HaldaneModel(double t1 = 1.0, double t2 = 0.1, double phi = M_PI / 2.0, double M = 0.2)
        : t1(t1), t2(t2), phi(phi), M(M) {}

    Mat Hk(const Vec& k) const override {
        using namespace std::complex_literals;
        using std::cos, std::sin;
        const std::complex<double> I(0,1);


        double kx = k(0), ky = k(1);
        std::complex<double> f = t1 * (1.0 + std::exp(-I * kx) + std::exp(-I * (kx / 2 + std::sqrt(3) * ky / 2)));
        double d_x = f.real();
        double d_y = f.imag();
        double d_z = M - 2.0 * t2 * sin(phi) * (sin(kx) - sin(kx / 2 + std::sqrt(3) * ky / 2));

        Mat H(2, 2);
        H << d_z, d_x - I * d_y,
             d_x + I * d_y, -d_z;
        return H;
    }
};
