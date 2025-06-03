#include "mesh.hpp"
#include "tightbinding_k.hpp"
#include <iostream>

int main() {
    TightBindingSquare tb(1.0); // Create a tight-binding model with t = 1.0
    
    Mesh mesh(5, 5); // Create a 2D mesh with 5x5 k-points
    const auto& kpoints = mesh.getKPoints();

    std::cout << "Generated " << mesh.size() << " k-points:\n";
    for (const auto& kp : kpoints) {
        double energy = tb.energy(kp); // Calculate energy for each k-point
        Eigen::Vector3d velocity = tb.groupVelocity(kp); // Calculate group velocity for each k-point
        std::cout << "k: " << kp.transpose() 
                  << " | Energy: " << energy 
                  << " | Velocity: " << velocity.transpose() << "\n";
    }
    return 0;
}
