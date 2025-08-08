Quantum Burgers' Equation Solver

Project Overview
This repository contains the full submission for a project that implements a quantum algorithm to solve the 1-D viscous Burgers' equation. The project uses the Hamiltonian Simulation using an Expander (HSE) framework, leveraging the Cole-Hopf transformation to map the non-linear problem to a linear diffusion equation solvable on a quantum computer.

The submission includes:

A detailed Algorithm Design Brief describing the theoretical framework, gate decomposition, and resource estimates.

Prototype Code for both a classical finite-difference solver and a quantum solver implemented with Qiskit.

A Validation and Benchmark report comparing the quantum solver's performance against the classical reference.

A Resource and Noise Analysis that includes a Zero Noise Extrapolation (ZNE) strategy for error mitigation.

A Scalability Study and an Algorithm Comparison between QTN and HSE.
Getting Started
Prerequisites
To run the code, you'll need Python 3.8+ and the following libraries. You can install them using pip:

pip install -r requirements.txt

The requirements.txt file should contain:

numpy
qiskit
qiskit-aer
matplotlib
Code Description
classical_solver.py: Implements a high-fidelity finite-difference solver for the 1-D Burgers' equation.

quantum_solver.py: Contains the core quantum algorithm, including the Cole-Hopf transformation, the Trotterized time-evolution circuit, and the Zero Noise Extrapolation functions.

main.py: The entry point for the project, which orchestrates the entire simulation, benchmarking, and visualization process.