# QuantumTransport++

[![docs](https://img.shields.io/badge/docs-Doxygen-blue.svg)](https://javahedi.github.io/QuantumTransportPP/)


🚀 **QuantumTransport++** is a high-performance C++/Python hybrid framework for studying quantum transport phenomena in tight-binding models. It enables accurate and efficient simulations of **semiclassical dynamics**, **Berry curvature effects**, **Kubo linear response**, and **Boltzmann transport** — with quantum geometry built-in.

---

## 🔬 Key Features

- ⚛️ **Tight-binding models in momentum space** — general N-band Hamiltonians
- 🌀 **Berry curvature**, **quantum metric**, and **Berry dipole** computations
- 📈 **Kubo conductivity** solver (linear response)
- 🚗 **Semiclassical Boltzmann solver** with Berry corrections
- 🔁 (Planned) **Floquet field support**: time-periodic `A(t)` fields
- 🧠 Automatic **Brillouin zone meshing** (Monkhorst-Pack, hexagonal)
- 🐍 Python bindings via **pybind11**
- 📊 Example Jupyter notebooks for reproducible plots

---
## Tight-binding models

| Model               | Type                 | Key Features                             |
| ------------------- | -------------------- | ---------------------------------------- |
| **Kane-Mele**       | 2D Topological       | SOC + sublattice potential               |
| **BHZ Model**       | 2D Quantum Spin Hall | 4×4 Dirac-like                           |
| **Weyl Model**      | 3D Topological       | Linear crossings, chirality              |
| **Graphene**        | 2D Semimetal         | Simple honeycomb TB                      |
| **Chern Insulator** | 2-band               | Trivial/Topological phase tuning         |
| **Haldane**         | 2D Chern Insulator   | Complex hopping, breaks TRS, topological |


## 📦 Installation

### Prerequisites

- C++17 compiler
- CMake ≥ 3.15
- Python ≥ 3.7
- [pybind11](https://github.com/pybind/pybind11)

### Build Instructions

```bash
git clone https://github.com/yourusername/QuantumTransportPP.git
cd QuantumTransportPP
mkdir build && cd build
cmake ..
make
````

This builds the Python module `qtpp`.

---

## 🐍 Python Usage Example

```python
import qtpp

model = qtpp.TightBindingModel(2)  # e.g. 2-band model
res = qtpp.compute_boltzmann_transport(model, Ef=0.1, Nk=50)

print("σ_xx:", res.conductivity[0])
print("σ_xy (Berry curvature):", res.hall_conductivity)
print("Berry dipole:", res.berry_dipole_magnitude)
```

See `python/example_kubo.ipynb` and `example_boltzmann.ipynb` for full tutorials.

---

## 📂 Project Structure

```
QuantumTransport++/
├── include/           # C++ headers
├── src/               # C++ source files
├── python/            # Example notebooks
├── bindings.cpp       # pybind11 interface
├── CMakeLists.txt     # Build system
├── README.md
```

---

## 📈 Output Formats

Simulation results (band structure, Berry curvature, conductivities) are exported as CSV or NumPy arrays for easy visualization using:

* `matplotlib`, `seaborn` (Python)
* `gnuplot`, `paraview` (optional)

---

## 📌 Roadmap

* [x] Band structure + eigen solver in k-space
* [x] Berry curvature, quantum metric
* [x] Kubo and Boltzmann transport
* [ ] Time-dependent fields + Floquet formalism
* [ ] Python package (`pip install quantumtransportpp`)
* [ ] GUI / interactive viewer (optional)

---

## 🤝 Contributing

Pull requests are welcome! For major changes, please open an issue first to discuss your ideas.

We welcome:

* New tight-binding models
* Additional transport quantities
* Numerical stability improvements
* Parallelization or GPU backends


---

## 📜 License

MIT License. See `LICENSE` file.

---

## 🧠 Acknowledgments

Inspired by the quantum geometry formalism, semiclassical theory, and nonlinear response frameworks.

---

🧪 Built for researchers. Optimized for speed. Designed for insight.
