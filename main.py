import numpy as np
import matplotlib.pyplot as plt
import time
# Assuming classical_solver.py and quantum_solver.py are in the same directory
from classical_solver import classical_burgers_solver
from quantum_solver import (
    get_initial_phi_statevector,
    create_time_evolution_circuit,
    post_process_results,
    scale_circuit_noise,
    extrapolate_zne_results,
)
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import StatePreparation

# --- Main Execution and Benchmarking ---
if __name__ == "__main__":
    # --- Simulation Parameters ---
    nu = 0.01  # Viscosity coefficient
    L = 1.0    # Domain length
    nx = 16    # Number of spatial grid points (must be a power of 2)
    dx = L / (nx - 1)
    
    total_time = 0.05
    dt_classical = 0.0001
    num_time_steps_classical = int(total_time / dt_classical)
    
    dt_quantum_step = 0.005
    num_quantum_steps = int(total_time / dt_quantum_step)
    
    if num_quantum_steps < 3:
        num_quantum_steps = 3
        dt_quantum_step = total_time / num_quantum_steps

    num_trotter_steps_per_quantum_step = 10
    num_qubits = int(np.log2(nx))
    
    if 2**num_qubits != nx:
        raise ValueError("Number of spatial points (nx) must be a power of 2.")

    print(f"--- Simulation Setup ---")
    print(f"Viscosity (nu): {nu}")
    print(f"Domain: x in [0, {L}] with {nx} points (dx={dx:.4f})")
    print(f"Total Simulation Time: {total_time}")
    print(f"Quantum simulation with {num_qubits} qubits, {num_quantum_steps} quantum steps (dt_q={dt_quantum_step}), {num_trotter_steps_per_quantum_step} Trotter steps/quantum step.")
    print(f"Classical simulation with {num_time_steps_classical} steps (dt_c={dt_classical}).")

    # --- Classical Pre-processing: Initial Condition ---
    x = np.linspace(0, L, nx)
    # Riemann step initial condition: u(x,0) = 1 for x <= 0.5, 0 otherwise
    u0_initial = np.where(x <= 0.5, 1.0, 0.0)

    # --- Run Classical Solver for Reference ---
    print("\n--- Running Classical Solver ---")
    start_time_classical = time.perf_counter()
    u_classical_reference = classical_burgers_solver(u0_initial, nu, dx, dt_classical, num_time_steps_classical)
    end_time_classical = time.perf_counter()
    classical_wall_time = end_time_classical - start_time_classical
    print(f"Classical solver wall-clock time: {classical_wall_time:.6f} seconds")

    # --- Quantum Solver Workflow ---
    print("\n--- Running Quantum Solver ---")
    
    # Initialize phi with the initial u0
    current_phi_vector = get_initial_phi_statevector(u0_initial, nu, dx)

    # Store benchmark results
    benchmark_results = []
    
    # Define noise model for noisy simulations
    noise_model = NoiseModel()
    error_1_qubit = depolarizing_error(0.001, 1)
    noise_model.add_all_qubit_quantum_error(error_1_qubit, ['rz', 'rx', 'ry'])
    error_2_qubit = depolarizing_error(0.01, 2)
    noise_model.add_all_qubit_quantum_error(error_2_qubit, ['cx', 'rzz'])
    
    # Define noise scaling factors for ZNE
    noise_factors = [1.0, 2.0, 3.0]
    
    for q_step in range(num_quantum_steps):
        current_sim_time = (q_step + 1) * dt_quantum_step
        print(f"\n  Quantum Step {q_step + 1}/{num_quantum_steps} (Simulated Time: {current_sim_time:.4f})")
        
        current_phi_vector = current_phi_vector / np.linalg.norm(current_phi_vector)
        initial_phi_circuit = StatePreparation(current_phi_vector) 
        
        base_qc = QuantumCircuit(num_qubits)
        base_qc.append(initial_phi_circuit, base_qc.qubits)
        evolution_circuit = create_time_evolution_circuit(num_qubits, nu, dt_quantum_step, num_trotter_steps_per_quantum_step, dx)
        base_qc.compose(evolution_circuit, inplace=True)
        base_qc.save_statevector()
        
        simulator_noiseless = AerSimulator()
        transpiled_qc_noiseless = transpile(base_qc, simulator_noiseless)
        job_noiseless = simulator_noiseless.run(transpiled_qc_noiseless, shots=1)
        result_noiseless = job_noiseless.result()
        final_phi_vector_noiseless = np.real(result_noiseless.get_statevector(base_qc).data)
        u_quantum_noiseless = post_process_results(final_phi_vector_noiseless, dx, nu)
        
        u_vectors_noisy_for_zne = []
        for factor in noise_factors:
            scaled_qc = scale_circuit_noise(base_qc, factor)
            simulator_noisy = AerSimulator(noise_model=noise_model)
            transpiled_qc_noisy = transpile(scaled_qc, simulator_noisy)
            job_noisy = simulator_noisy.run(transpiled_qc_noisy, shots=1) 
            result_noisy = job_noisy.result()
            final_phi_vector_noisy = np.real(result_noisy.get_statevector(scaled_qc).data)
            u_noisy_result = post_process_results(final_phi_vector_noisy, dx, nu)
            u_vectors_noisy_for_zne.append(u_noisy_result)
        
        u_quantum_zne = extrapolate_zne_results(u_vectors_noisy_for_zne, noise_factors)
        current_phi_vector = final_phi_vector_noiseless

        u_classical_step_reference = classical_burgers_solver(u0_initial, nu, dx, dt_classical, int(current_sim_time / dt_classical))
        
        l2_error_noiseless = np.linalg.norm(u_quantum_noiseless - u_classical_step_reference) / np.linalg.norm(u_classical_step_reference)
        l2_error_zne = np.linalg.norm(u_quantum_zne - u_classical_step_reference) / np.linalg.norm(u_classical_step_reference)
        
        benchmark_results.append({
            'time_step': q_step + 1,
            'sim_time': current_sim_time,
            'l2_error_noiseless': l2_error_noiseless,
            'l2_error_zne': l2_error_zne,
            'u_quantum_noiseless': u_quantum_noiseless,
            'u_quantum_zne': u_quantum_zne,
            'u_classical_reference': u_classical_step_reference
        })
        
        print(f"    L2-error (Noiseless): {l2_error_noiseless:.4f}")
        print(f"    L2-error (ZNE): {l2_error_zne:.4f}")

    # --- Final Post-processing and Reporting ---
    print("\n--- Final Benchmark Report ---")
    print(f"Classical Solver Total Wall-Clock Time: {classical_wall_time:.6f} seconds")
    print("\n| Time Step | Sim. Time (s) | L2-Error (Noiseless) | L2-Error (ZNE) |")
    print("|-----------|---------------|----------------------|----------------|")
    for res in benchmark_results:
        print(f"| {res['time_step']:<9} | {res['sim_time']:<13.4f} | {res['l2_error_noiseless']:<20.4f} | {res['l2_error_zne']:<14.4f} |")

    # --- Visualization of Final Results ---
    try:
        plt.figure(figsize=(12, 8))
        plt.plot(x, u_classical_reference, 'o-', label='Classical Solver (Reference)')
        plt.plot(x, benchmark_results[-1]['u_quantum_noiseless'], 'x--', label='Quantum Solver (HSE - Noiseless)')
        plt.plot(x, benchmark_results[-1]['u_quantum_zne'], 's--', label='Quantum Solver (HSE - ZNE Mitigated)')
        plt.title(f'1D Burgers\' Equation Shock Tube at t={total_time}')
        plt.xlabel('x')
        plt.ylabel('u(x,t)')
        plt.legend()
        plt.grid(True)
        plt.show()
    except ImportError:
        print("\nMatplotlib not installed. Cannot plot results.")
    except Exception as e:
        print(f"\nAn error occurred during plotting: {e}")

    # --- Quantum Circuit Drawing (for the last quantum step) ---
    print("\n--- Example Quantum Circuit for one time step ---")
    print(base_qc.draw(output='text', idle_wires=False))