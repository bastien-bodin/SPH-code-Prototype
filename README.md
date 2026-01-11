# SPH-Prototype v2.0: High-Performance Cryolava Solver

<!-- [![Fortran](https://img.shields.io/badge/Language-Fortran%202003%2F2008-734f96.svg)](https://gcc.gnu.org/fortran/)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE) -->

**SPH-Prototype** is a Weakly Compressible Smoothed Particle Hydrodynamics (WCSPH) solver designed for modeling cryolava flows on planetary surfaces. While its primary focus is the low-gravity environment of **Europa** ($g \approx 1.31 \, \text{m/s}^2$), v2.0 introduces robust support for terrestrial validation cases.

---

## Previous version

First prototype for making SPH Simulations. Test with a Dam Break problem.
based on the work done by Goffin, Louis. "Development of a didactic SPH model." (2013).

It served as a playground to learn the bases of the SPH method.
For the current model, please go see the [CryoPy repository](https://github.com/bastien-bodin/CryoPy)

---

## 🎓 How to Cite

If you use **SPH-Prototype** in your research, papers, or reports, please cite it as follows:

> Bodin, B. (2026). *SPH-Prototype: High-Performance WCSPH Solver for Cryolava Flows* (Version 2.0). GitHub. https://github.com/bastien-bodin/SPH-Europa

Alternatively, you can use the **"Cite this repository"** button provided by GitHub in the right sidebar to export the citation in BibTeX or APA format.

---

## 🚀 Key Features in v2.0

* **Optimized Architecture:** Fully refactored into **Structure of Arrays (SoA)** to improve cache locality and enable SIMD vectorization.
* **Hydrostatic Initialization:** Particles are initialized with consistent pressure/density gradients to eliminate initial numerical shocks.
* **Dynamic Boundary Conditions (DBC):** Implementation of staggered double-layer boundaries following **Cabrera-Crespo (2007)**.
* **Advanced Integrators:** Choice between 2nd-order Runge-Kutta (**RK22**) and **Symplectic Velocity Verlet**.
* **Parallel Computing:** Multi-threaded execution via **OpenMP** with a serial toggle for debugging.
* **Flexible Physics:** Automated gravity vector $\vec{g}$ scaling based on problem dimensionality ($nDim$).

---

## 🔬 Physical Model

### 1. Governing Equations
The fluid motion is governed by the Lagrangian form of the Navier-Stokes equations:

* **Continuity Equation (Mass Conservation):**
    $$\frac{d\rho_a}{dt} = \sum_b m_b (\vec{v}_a - \vec{v}_b) \cdot \vec{\nabla} W_{ab}$$

* **Momentum Equation:**
    $$\frac{d\vec{v}_a}{dt} = -\sum_b m_b \left( \frac{P_a}{\rho_a^2} + \frac{P_b}{\rho_b^2} + \Pi_{ab} \right) \vec{\nabla} W_{ab} + \vec{g}$$



### 2. Equation of State (Tait EoS)
To ensure weakly compressible behavior, we use the Tait equation:
$$P = B \left[ \left( \frac{\rho}{\rho_0} \right)^\gamma - 1 \right] \quad \text{where} \quad B = \frac{\rho_0 c_0^2}{\gamma}$$

---

## 💻 Code Structure (SoA Layout)

The SoA layout ensures that coordinates, velocities, and properties are stored in contiguous memory blocks, significantly boosting performance in the force calculation loop.

| Module | Purpose |
| :--- | :--- |
| `parameters.f90` | Global constants, gravity vectors, and simulation flags. |
| `particles.f90` | Definition of the `ParticleSystem` type and SoA memory management. |
| `geometries.f90` | Generators for tanks, fluid blocks, and staggered DBC layers. |
| `forces.f90` | SPH kernels, neighbor search interface, and force evaluation. |
| `integrator.f90` | Time-stepping logic (RK22, Verlet, Euler) and data shifting. |
| `application.f90` | Simulation orchestration and hydrostatic setup. |


---

## 🛠 Compilation & Installation

### Requirements
* **Compiler:** `gfortran` (GCC 9+ recommended).
* **Parallelism:** OpenMP 4.5+ support.

### Building the Project
```bash
# Compile parallel version (Default)
make clean && make

# Compile serial version (For debugging/profiling)
make clean && make OMP=0
```

---

## 📊 Running & Visualization

### Executing a Simulation
By default, the solver is configured for the **4m x 3m Terrestrial Dam-Break** validation case. 
To run the simulation:
```bash
./sph_europa
```
The solver will generate a sequence of `.csv` files (e.g., `output_0000.csv`, `output_0050.csv`, etc.) containing particle coordinates, velocities, density, and pressure.

_Note: Ensure you have compiled the code using `make` before running._

### Visualization workflow (Paraview)
To visualize the fluid dynamics and the impact splash, follow this optimized ParaView workflow:

1. Load Data: Open ParaView and select the group of output_*.csv files. Click Apply.
2. Geometry Conversion: Apply the Table To Points filter:
    - Set `X Column` to `x`.
    - Set `Y Column` to `y`.
    - Check `2D Points`.
    - Click `Apply`.
3. Fluid Rendering: To achieve a realistic fluid look instead of raw points:
    - Change the Representation (top toolbar) from Surface to Point Gaussian.
    - In the Properties panel, set the Gaussian Radius to your simulation spacing (e.g., 0.04).
    - Use the Shader Preset "Sphere" for better depth.
4. Data Analysis:
    - Pressure Field: Color by p to observe the hydrostatic gradient and the pressure spikes during the wall impact.
    - Velocity Field: Color by v_mag to track the wavefront progression and the "splash" height.
    - Boundary Check: Use a Threshold filter on the mobile field (set to 0) to isolate and inspect the Dynamic Boundary Particles (DBC).

---
## 📜 Academic Reference

This version implements the staggered double-layer DBC methodology as described in:

> Cabrera Crespo, A. J., Gómez Gesteira, R., & Dalrymple, R. A. (2007). Boundary conditions generated by dynamic particles in SPH methods. Computers, Materials, & Continua.

---
**Author**: Bastien Bodin, PhD in Physics

**Contact**: [bastien.bodin@proton.me](mailto:bastien.bodin@proton.me)