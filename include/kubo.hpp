/*
Kubo-Greenwood 
*/
#pragma once
#include "hamiltonian.hpp"
#include "mesh.hpp"
#include <Eigen/Dense>

class KuboSolver {
public:
    KuboSolver(const Hamiltonian& H, const Mesh& mesh, double eta = 1e-3,
               bool temperature_in_kelvin = false, double energy_scale = 1.0);

    // Compute real part of conductivity tensor σ_ij
    Eigen::Matrix3d conductivity(double Ef, double T);

private:
    const Hamiltonian& H;
    const Mesh& mesh;
    double eta;
    bool temperature_in_kelvin;
    double energy_scale;
};
