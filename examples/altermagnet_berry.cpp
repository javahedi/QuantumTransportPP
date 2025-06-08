#include "mesh.hpp"
#include "altermagnet.hpp"
#include "geometry.hpp" 
#include "logInfo.hpp"
#include <fstream>


using namespace QuantumTransport;

int main() {
    AltermagnetModel altermagnet(1.0, 0.5, 0.1); // Example parameters
    size_t N = 1000; // Grid resolution
    Mesh mesh(N, N, 1, M_PI); // 2D mesh


    logInfo("Starting altermagnet Berry curvature calculation...");
    std::ofstream out("altermagnet_berry.csv");
    out << "kx,ky,Berry_curvature_FHS, berryCurvatureDiff\n";

    for (const auto& k : mesh.getKPoints()) {
        double berry_curvature_fhs  = berryCurvatureFHS(altermagnet, k, 1e-2, 1); // Band index 1
        double berry_curvature_diff = berryCurvatureDifferential(altermagnet, k, 1e-2, 1); // Band index 1
        out << k(0) << "," << k(1) << "," << berry_curvature_fhs << "," << berry_curvature_diff << "\n";
    }
    out.close();
    logSuccess("Altermagnet Berry curvature calculation completed successfully. Results saved to altermagnet_berry.csv");
    return 0;
}



