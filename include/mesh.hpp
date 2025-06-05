#pragma once // Include guard to prevent multiple inclusions
#include <vector>
#include <Eigen/Dense>

/// @brief Represents a k-point mesh for Brillouin zone sampling
class Mesh {
public:
    using Vec = Eigen::Vector3d;
    // int is a signed integer, which can be negative
    // size_t is an unsigned integer, which cannot be negative

    /// @brief Construct a mesh with given dimensions
    Mesh(size_t nx, size_t ny, size_t nz = 1, double kmax = 2 * M_PI) ;
        
    /// @brief Return all k-points in the mesh
    const std::vector<Vec>& getKPoints() const; 

    /// @brief Return the number of k-points
    size_t size() const;


private:
    size_t nx_, ny_, nz_;
    double kmax_;
    std::vector<Vec> kpoints_;
    // Uniform grid: k = 2 * pi * (i, j, [,k]) / N

    /// @brief Generate the mesh grid
    void generateMesh(); // Generate the k-point mesh based on nx, ny, nz, and kmax

};