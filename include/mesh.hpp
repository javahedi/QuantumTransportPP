#pragma once // Include guard to prevent multiple inclusions
// #ifndef BANDSTRUCTURE_HPP  
// #define BANDSTRUCTURE_HPP
//  ....
//#endif // BANDSTRUCTURE_HPP

#include <vector>
#include <Eigen/Dense>

class Mesh {
public:
    using Vec = Eigen::Vector3d;
    // int is a signed integer, which can be negative
    // size_t is an unsigned integer, which cannot be negative
    Mesh(size_t nx, size_t ny, size_t nz = 1, double kmax = 2 * M_PI) ;
        

    const std::vector<Vec>& getKPoints() const; 
    size_t size() const;


private:
    size_t nx_, ny_, nz_;
    double kmax_;
    std::vector<Vec> kpoints_;
    // Uniform grid: k = 2 * pi * (i, j, [,k]) / N

    void generateMesh(); // Generate the k-point mesh based on nx, ny, nz, and kmax

};