#include "mesh.hpp"

Mesh::Mesh(size_t nx, size_t ny, size_t nz, double kmax)
    : nx_(nx), ny_(ny), nz_(nz), kmax_(kmax) {
    generateMesh();
}

const std::vector<Mesh::Vec>& Mesh::getKPoints() const {
    return kpoints_;
}

size_t Mesh::size() const {
    return kpoints_.size();
}

void Mesh::generateMesh() {
    kpoints_.clear();
    for (size_t i = 0; i < nx_; ++i) {
        for (size_t j = 0; j < ny_; ++j) {
            for (size_t k = 0; k < nz_; ++k) {
                Vec kp;
                kp(0) = -kmax_ + 2 * kmax_ * i / (nx_ - 1);
                kp(1) = -kmax_ + 2 * kmax_ * j / (ny_ - 1);
                kp(2) = -kmax_ + 2 * kmax_ * k / (nz_ - 1);
                kpoints_.emplace_back(kp);
            }
        }
    }
}
