import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit.circuit.library import StatePreparation
from qiskit.quantum_info import Statevector
import matplotlib.pyplot as plt
import time
from qiskit_aer.noise import NoiseModel, depolarizing_error

# This is a placeholder for the classical solver to be used for validation.
# It is assumed to be in a separate file (classical_solver.py).
def classical_burgers_solver(u0, nu, dx, dt, num_time_steps):
    """
    Placeholder for the classical solver function. 
    It is used here for benchmarking and validation.
    """
    # This is a basic forward-Euler finite difference solver.
    u = u0.copy()
    for _ in range(num_time_steps):
        u_new = u.copy()
        for i in range(1, len(u) - 1):
            convection = u[i] * (u[i] - u[i-1]) / dx if u[i] >= 0 else u[i] * (u[i+1] - u[i]) / dx
            diffusion = nu * (u[i+1] - 2 * u[i] + u[i-1]) / dx**2
            u_new[i] = u[i] - dt * convection + dt * diffusion
        u = u_new
    return u

# --- Quantum Solver Components (HSE Framework) ---

def get_initial_phi_statevector(u0, nu, dx):
    """
    Applies the Cole-Hopf transformation to the initial velocity field u0
    to obtain the initial phi(x,0) statevector for the quantum algorithm.
    
    The transformation is: phi(x, 0) = A * exp(-1/(2*nu) * integral(u(xi, 0) dxi))
    
    Args:
        u0 (np.array): Initial velocity field.
        nu (float): Viscosity coefficient.
        dx (float): Spatial step size.
        
    Returns:
        np.array: Normalized statevector representing the initial phi.
    """
    nx = len(u0)
    phi = np.zeros(nx)
    
    # Numerically integrate u0 using cumulative sum (a simple trapezoidal rule).
    integral_u = np.cumsum(u0) * dx
    
    # Apply the exponential transformation.
    phi = np.exp(-integral_u / (2 * nu))
    
    # Normalize the statevector for quantum amplitude encoding.
    phi = phi / np.linalg.norm(phi)
    
    return phi

def create_time_evolution_circuit(num_qubits, nu, dt, num_trotter_steps, dx):
    """
    Creates a quantum circuit for time evolution of phi using a Trotter-Suzuki decomposition.
    The Hamiltonian is proportional to the discretized Laplacian operator.
    
    Args:
        num_qubits (int): Number of qubits (log2 of grid size).
        nu (float): Viscosity coefficient.
        dt (float): Total time step size for evolution.
        num_trotter_steps (int): Number of Trotter steps.
        dx (float): Spatial step size for Laplacian discretization.
    
    Returns:
        QuantumCircuit: Circuit for time evolution.
    """
    qc = QuantumCircuit(num_qubits)
    trotter_dt = dt / num_trotter_steps
    coupling_strength = nu * trotter_dt / (dx**2)
    
    for _ in range(num_trotter_steps):
        # Simulating nearest-neighbor coupling from the Laplacian with RZZ gates.
        for i in range(num_qubits - 1):
            qc.rzz(2 * coupling_strength, i, i+1)
        
        # Simulating diagonal terms from the Laplacian with RZ gates.
        for i in range(num_qubits):
            qc.rz(2 * coupling_strength * 2, i)
            
    return qc

def post_process_results(final_phi, dx, nu):
    """
    Applies the inverse Cole-Hopf transformation to get u(x,t) from the phi statevector.
    
    u(x, t) = -2*nu * (1/phi(x,t)) * d(phi)/dx
    
    Args:
        final_phi (np.array): Final phi statevector (amplitudes).
        dx (float): Spatial step size.
        nu (float): Viscosity coefficient.
    
    Returns:
        np.array: Final velocity field u(x,t).
    """
    dphi_dx = np.zeros_like(final_phi, dtype=float)
    
    # Use central difference for interior points.
    dphi_dx[1:-1] = (final_phi[2:] - final_phi[:-2]) / (2 * dx)
    # Use forward difference for the first point.
    dphi_dx[0] = (final_phi[1] - final_phi[0]) / dx
    # Use backward difference for the last point.
    dphi_dx[-1] = (final_phi[-1] - final_phi[-2]) / dx
    
    with np.errstate(divide='ignore', invalid='ignore'):
        u = -2 * nu * (dphi_dx / final_phi)
        u = np.nan_to_num(u, nan=0.0, posinf=0.0, neginf=0.0) 
    
    return u
    
def scale_circuit_noise(circuit, scale_factor):
    """
    Scales the noise in a circuit for Zero Noise Extrapolation (ZNE) by repeating two-qubit gates.
    
    Args:
        circuit (QuantumCircuit): The original quantum circuit.
        scale_factor (float): The factor by which to scale the noise (e.g., 1.0, 2.0, 3.0).
    
    Returns:
        QuantumCircuit: The scaled circuit.
    """
    if scale_factor == 1.0:
        return circuit.copy()
        
    scaled_qc = QuantumCircuit(circuit.num_qubits)
    num_repetitions = int(scale_factor)
    
    for instruction, qargs, cargs in circuit.data:
        # We assume two-qubit gates are the main source of error.
        if len(qargs) == 2 and instruction.name in ['rzz', 'cx']:
            for _ in range(num_repetitions):
                scaled_qc.append(instruction, qargs, cargs)
        else:
            scaled_qc.append(instruction, qargs, cargs)
            
    return scaled_qc

def extrapolate_zne_results(u_vectors_noisy, noise_factors):
    """
    Performs linear extrapolation on a set of noisy u-vectors to estimate the zero-noise result.
    
    Args:
        u_vectors_noisy (list of np.array): A list of u-vectors corresponding to each noise factor.
        noise_factors (list of float): The list of noise scale factors.
    
    Returns:
        np.array: The extrapolated u-vector at the zero-noise limit.
    """
    num_points = u_vectors_noisy[0].shape[0]
    u_zne_extrapolated = np.zeros(num_points)
    
    # Extrapolate each grid point's value independently using a linear fit.
    for i in range(num_points):
        y_values = [u[i] for u in u_vectors_noisy]
        
        try:
            poly = np.polyfit(noise_factors, y_values, 1)
            extrapolated_value = poly[1]
            u_zne_extrapolated[i] = extrapolated_value
        except Exception:
            # Fallback if the fitting fails.
            u_zne_extrapolated[i] = u_vectors_noisy[0][i]

    return u_zne_extrapolated


# --- Main Execution and Benchmarking ---
if __name__ == "__main__":
    # --- Simulation Parameters ---
    nu = 0.01
    L = 1.0
    nx = 16
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

    print("--- Simulation Setup ---")
    print(f"Viscosity (nu): {nu}")
    print(f"Quantum simulation with {num_qubits} qubits, {num_quantum_steps} quantum steps.")
    
    # --- Classical Pre-processing: Initial Condition ---
    x = np.linspace(0, L, nx)
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
    current_phi_vector = get_initial_phi_statevector(u0_initial, nu, dx)
    benchmark_results = []
    
    noise_model = NoiseModel()
    error_1_qubit = depolarizing_error(0.001, 1)
    noise_model.add_all_qubit_quantum_error(error_1_qubit, ['rz', 'rx', 'ry'])
    error_2_qubit = depolarizing_error(0.01, 2)
    noise_model.add_all_qubit_quantum_error(error_2_qubit, ['cx', 'rzz'])
    
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