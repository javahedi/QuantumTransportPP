# Quantum Transport in Altermagnetic Materials 🧲⚡

![Project Banner](https://via.placeholder.com/800x200/1a237e/ffffff?text=Quantum+Transport+in+Altermagnets)

A modern C++ toolkit for simulating quantum transport phenomena in altermagnetic materials using both Kubo and Boltzmann formalisms.

## 🌟 Key Features

- **Advanced Transport Calculations**: Compute conductivity (σ) and thermoelectric (α) tensors
- **Dual Solver Support**: Both Kubo and Boltzmann transport formalisms implemented
- **Parallel Computing**: OpenMP-accelerated calculations
- **Comprehensive Analysis**: Berry curvature, DOS, and band structure capabilities
- **Modern C++**: Clean, templated code with Eigen integration
- **Visualization Ready**: Jupyter notebooks for post-processing included

## 📦 Installation & Setup

### Prerequisites

- **Compiler**: GCC ≥ 9.0 or Clang ≥ 10.0
- **Dependencies**:
  - Eigen3 (≥ 3.3.7)
  - OpenMP
  - CMake (≥ 3.12)
  - Python (for post-processing)

### 🚀 Quick Start (Linux/macOS)

```bash
# Clone the repository
git clone https://github.com/yourusername/quantum-transport-altermagnets.git
cd quantum-transport-altermagnets

# Create build directory
mkdir build && cd build

# Configure and build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4

# Run an example
./examples/altermagnet_conductivity
```

### 🏗️ Build Options

| Option | Description | Default |
|--------|-------------|---------|
| `-DUSE_OPENMP=ON` | Enable parallel computation | ON |
| `-DBUILD_TESTS=ON` | Build test suite | OFF |
| `-DBUILD_DOCS=ON` | Build Doxygen documentation | OFF |

## 🧮 Running Simulations

### Basic Usage

```bash
# Run conductivity calculations (4 threads)
./examples/altermagnet_conductivity 4

# Compute Berry curvature
./examples/altermagnet_berry
```

### Example Parameter File

Create `params.ini`:
```ini
[System]
J = 0.0:1.0:0.025      # J values (start:end:step)
lambda = 0.0,0.1,0.5   # Spin-orbit coupling values
mesh_points = 250      # k-point resolution

[Transport]
Ef = 0.5               # Fermi energy
temperature = 0.02     # Reduced temperature
eta = 1e-2             # Kubo broadening
tau = 100.0            # Relaxation time
```

## 📊 Post-Processing

Jupyter notebooks are provided for visualization:

```bash
jupyter notebook postprocessing/plot_conductivity.ipynb
```

![Berry Curvature Visualization](postprocessing/berry_curvature_pyplot.png)
## 🧩 Code Structure

```
quantum-transport-altermagnets/
├── include/            # Header files
├── src/               # Core implementation
├── examples/          # Ready-to-run examples
├── data/              # Output data
├── postprocessing/    # Visualization scripts
├── tests/             # Unit tests
└── docs/              # Documentation
```

## 🤝 Contributing

We welcome contributions! Please see our [Contribution Guidelines](CONTRIBUTING.md).

## 📜 License

MIT License - See [LICENSE](LICENSE) for details.

---

**Get Started Now** - Unlock the secrets of quantum transport in altermagnets! ✨

```bash
git clone https://github.com/yourusername/quantum-transport-altermagnets.git
```