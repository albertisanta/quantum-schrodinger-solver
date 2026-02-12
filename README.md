# Time-Dependent Schrödinger Equation Solver

A Fortran-based numerical solver for the 1D and 2D Time-Dependent Schrödinger Equation (TDSE). This project was developed as part of the Computational Physics course at the University of Barcelona.

## Features
* **1D Solver:** Implements the **Crank-Nicolson scheme**, ensuring a unitary evolution and conservation of the wavefunction norm.
* **Numerical Method:** Utilizes the **Thomas Algorithm** (Tridiagonal Matrix Algorithm) for efficient computation.
* **2D Extension:** Solves the 2D TDSE using the **Alternating Direction Implicit (ADI)** method.
* **Physics Applications:** Includes simulations of wavepacket scattering against potential barriers, with analysis of reflection and transmission coefficients.

## Physics Context
The solver handles the evolution of:
$$\psi(x, t+\Delta t) = \frac{1 - \frac{i\Delta t}{2\hbar}\hat{H}}{1 + \frac{i\Delta t}{2\hbar}\hat{H}}\psi(x, t)$$

## Requirements
* Fortran compiler (gfortran recommended).
* Gnuplot or Python (for visualization of the exported data).
