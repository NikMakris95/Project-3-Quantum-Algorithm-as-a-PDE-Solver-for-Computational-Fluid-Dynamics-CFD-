import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit.circuit.library import StatePreparation
from qiskit.quantum_info import Statevector
import matplotlib.pyplot as plt
import time # Import time module for wall-clock time measurement
from qiskit_aer.noise import NoiseModel, depolarizing_error # For building custom noise models for noisy simulation

# Import for IBM Quantum services (uncomment and configure if using real QPU)
# from qiskit_ibm_provider import IBMProvider 
# from qiskit_ibm_provider.job import IBMQJob

# --- Classical Solver for High-Fidelity Reference ---
def classical_burgers_solver(u0, nu, dx, dt, num_time_steps):
    """
    Classical finite-difference solver for 1D Burgers' equation.
    Uses an explicit upwind scheme for convection and central difference for diffusion.
    
    Args:
        u0 (np.array): Initial velocity field.
        nu (float): Viscosity coefficient.
        dx (float): Spatial step size.
        dt (float): Time step size.
        num_time_steps (int): Number of time steps to evolve.
        
    Returns:
        np.array: Final velocity field after num_time_steps.
    """
    u = u0.copy()
    nx = len(u)
    
    # Check CFL condition for stability of the explicit scheme
    max_u_abs = np.max(np.abs(u))
    cfl_convection = max_u_abs * dt / dx
    cfl_diffusion = nu * dt / dx**2
    
    if cfl_convection > 1.0 or cfl_diffusion > 0.5:
        print(f"Warning: CFL condition might be violated. Convection CFL: {cfl_convection:.2f}, Diffusion CFL: {cfl_diffusion:.2f}")
        print("Results may be unstable or inaccurate. Consider reducing dt or dx.")
        
    for _ in range(num_time_steps):
        u_new = u.copy()
        for i in range(1, nx - 1):
            # Upwind scheme for convection: depends on the sign of u[i]
            if u[i] >= 0: # Use forward difference if u[i] is positive or zero
                convection_term = u[i] * (u[i] - u[i-1]) / dx
            else: # Use backward difference if u[i] is negative
                convection_term = u[i] * (u[i+1] - u[i]) / dx
            
            # Central difference for diffusion
            diffusion_term = nu * (u[i+1] - 2 * u[i] + u[i-1]) / dx**2
            
            u_new[i] = u[i] - dt * convection_term + dt * diffusion_term
        
        # Apply Dirichlet boundary conditions (fixed at ends)
        u_new[0] = u0[0] # Left boundary fixed to initial left state
        u_new[nx-1] = u0[nx-1] # Right boundary fixed to initial right state
        u = u_new
        
    return u

# --- Quantum Solver Components (HSE Framework) ---
def get_initial_phi_statevector(u0, nu, dx):
    """
    Applies the Cole-Hopf transformation to the initial velocity field u0
    to obtain the initial phi(x,0) statevector.
    
    phi(x, 0) = A * exp(-1/(2*nu) * integral(u(xi, 0) dxi))
    
    Args:
        u0 (np.array): Initial velocity field.
        nu (float): Viscosity coefficient.
        dx (float): Spatial step size.
        
    Returns:
        np.array: Normalized statevector representing the initial phi.
    """
    nx = len(u0)
    phi = np.zeros(nx)
    
    # Numerically integrate u0 using cumulative sum
    # We choose x0 = 0 for the integration start point.
    integral_u = np.cumsum(u0) * dx
    
    # Apply the exponential transformation
    phi = np.exp(-integral_u / (2 * nu))
    
    # Normalize the statevector for quantum amplitude encoding
    phi = phi / np.linalg.norm(phi)
    
    return phi

def create_time_evolution_circuit(num_qubits, nu, dt, num_trotter_steps, dx):
    """
    Creates a quantum circuit for time evolution of phi using a Trotter-Suzuki decomposition.
    The Hamiltonian for the diffusion equation is H = -i * nu * Laplacian.
    
    Args:
        num_qubits (int): Number of qubits (log2 of grid size).
        nu (float): Viscosity coefficient.
        dt (float): Total time step size for evolution.
        num_trotter_steps (int): Number of Trotter steps.
        dx (float): Spatial step size for Laplacian discretization.
    
    Returns:
        QuantumCircuit: Circuit for time evolution.
    """
    num_qubits=4
    qc = QuantumCircuit(num_qubits)
    
    # The Hamiltonian for the diffusion equation (after mapping to quantum system) is
    # proportional to the discretized Laplacian operator.
    # For a 1D system with central differences, the discretized Laplacian is a tridiagonal matrix
    # with -2 on the diagonal and 1 on the off-diagonals (scaled by 1/dx^2).
    # H_eff = -nu * (1/dx^2) * L_matrix
    # where L_matrix has -2 on diag, 1 on off-diag.
    #
    # This Hamiltonian can be decomposed into a sum of Pauli terms.
    # For a 4-qubit system (16 points), a common decomposition strategy for a 1D chain
    # involves terms like Z_i Z_{i+1} and single-ququbit Z_i terms.
    #
    # Below is a more structured placeholder for Trotter steps that mimics
    # the application of exponentials of Pauli strings, rather than a random sequence.
    # This is NOT the exact decomposition of the Laplacian, but demonstrates the
    # *type* of gates used in such a decomposition.
    # A full implementation would derive the exact Pauli decomposition of the
    # discretized Laplacian and then apply the corresponding R_P(theta) gates.
    
    trotter_dt = dt / num_trotter_steps
    
    # Example coefficients for illustrative Pauli terms
    # These coefficients would be derived from the actual Laplacian decomposition
    # The factor of 2*pi is often included for Rz gates in Qiskit if angle is in cycles,
    # but here we use radians, so 2*theta is for exp(-i*theta*Z)
    
    # Example: A simple 1D chain interaction
    coupling_strength = nu * trotter_dt / (dx**2) # Scale factor for the Hamiltonian terms
    
    for _ in range(num_trotter_steps):
        # Apply ZZ interactions (simulating nearest-neighbor coupling from Laplacian)
        for i in range(num_qubits - 1):
            qc.rzz(2 * coupling_strength, i, i+1) # RZZ gate for exp(-i*theta*Z_i*Z_{i+1})
        
        # Apply Z rotations (simulating diagonal terms from Laplacian)
        for i in range(num_qubits):
            qc.rz(2 * coupling_strength * 2, i) # RZ gate for exp(-i*theta*Z_i)
            
    return qc

def post_process_results(final_phi, dx, nu):
    """
    Applies the inverse Cole-Hopf transformation to get u(x,t) from phi.
    
    u(x, t) = -2*nu * (1/phi(x,t)) * d(phi)/dx
    
    Args:
        final_phi (np.array): Final phi statevector (amplitudes).
        dx (float): Spatial step size.
        nu (float): Viscosity coefficient.
    
    Returns:
        np.array: Final velocity field u(x,t).
    """
    # Numerically differentiate phi using central difference for interior points
    # For boundary points, a forward/backward difference might be needed or handled by BCs.
    dphi_dx = np.zeros_like(final_phi, dtype=float)
    
    # Central difference for interior points
    dphi_dx[1:-1] = (final_phi[2:] - final_phi[:-2]) / (2 * dx)
    
    # Forward difference for the first point (x=0)
    dphi_dx[0] = (final_phi[1] - final_phi[0]) / dx
    
    # Backward difference for the last point (x=L)
    dphi_dx[-1] = (final_phi[-1] - final_phi[-2]) / dx
    
    # Apply the inverse transformation formula
    # Handle potential division by zero or very small phi values
    with np.errstate(divide='ignore', invalid='ignore'):
        u = -2 * nu * (dphi_dx / final_phi)
        # Replace NaN (from 0/0) or inf (from X/0) with a reasonable default, e.g., 0
        # In a real scenario, these points might need special handling based on physics.
        u = np.nan_to_num(u, nan=0.0, posinf=0.0, neginf=0.0) 
    
    return u
    
def scale_circuit_noise(circuit, scale_factor):
    """
    Scales the noise in a circuit by repeating two-qubit gates.
    A scale factor of 1.0 returns the original circuit.
    A scale factor of 2.0 would repeat each two-qubit gate once.
    This is a common method for ZNE on noisy simulators.
    
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
        # Only scale two-qubit gates for simplicity, as they are the primary source of error
        if len(qargs) == 2 and instruction.name in ['rzz', 'cx']:
            for _ in range(num_repetitions):
                scaled_qc.append(instruction, qargs, cargs)
        else:
            scaled_qc.append(instruction, qargs, cargs)
            
    return scaled_qc

def extrapolate_zne_results(u_vectors_noisy, noise_factors):
    """
    Performs linear extrapolation on a set of noisy u-vectors.
    
    Args:
        u_vectors_noisy (list of np.array): A list of u-vectors corresponding to
                                            each noise factor.
        noise_factors (list of float): The list of noise scale factors.
    
    Returns:
        np.array: The extrapolated u-vector at the zero-noise limit.
    """
    num_points = u_vectors_noisy[0].shape[0]
    u_zne_extrapolated = np.zeros(num_points)
    
    # Extrapolate each grid point's value independently
    for i in range(num_points):
        y_values = [u[i] for u in u_vectors_noisy]
        
        # Fit a linear polynomial: u = a*lambda + b, where lambda is the noise factor
        # We want to find the value at lambda=0, which is 'b'.
        try:
            poly = np.polyfit(noise_factors, y_values, 1) # 1st degree polynomial
            extrapolated_value = poly[1] # 'b' is the y-intercept
            u_zne_extrapolated[i] = extrapolated_value
        except Exception as e:
            # Fallback if fitting fails (e.g., not enough data points)
            print(f"Warning: Failed to extrapolate for point {i}. Using the least noisy result.")
            u_zne_extrapolated[i] = u_vectors_noisy[0][i]

    return u_zne_extrapolated


# --- Main Execution and Benchmarking ---
if __name__ == "__main__":
    # --- Simulation Parameters ---
    nu = 0.01  # Viscosity coefficient (for diffusion)
    L = 1.0    # Domain length
    nx = 16    # Number of spatial grid points (must be a power of 2 for amplitude encoding)
    dx = L / (nx - 1) # Spatial step size
    
    # Time parameters
    total_time = 0.05 # Total time to simulate
    dt_classical = 0.0001 # Time step for classical solver (smaller for accuracy/stability)
    num_time_steps_classical = int(total_time / dt_classical)
    
    dt_quantum_step = 0.005 # Quantum time step per evolution circuit
    num_quantum_steps = int(total_time / dt_quantum_step) # Number of quantum evolution steps
    
    # Ensure at least 3 quantum time steps for reporting as per challenge
    if num_quantum_steps < 3:
        num_quantum_steps = 3
        dt_quantum_step = total_time / num_quantum_steps
        print(f"Adjusted quantum steps to {num_quantum_steps} for reporting, new dt_quantum_step: {dt_quantum_step:.4f}")

    num_trotter_steps_per_quantum_step = 10 # Number of Trotter steps within each quantum dt
    num_qubits = int(np.log2(nx)) # Number of qubits required for amplitude encoding
    
    if 2**num_qubits != nx:
        raise ValueError("Number of spatial points (nx) must be a power of 2 for amplitude encoding.")

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
    print(f"Classical solution at t={total_time}: {u_classical_reference}")
    print(f"Classical solver wall-clock time: {classical_wall_time:.6f} seconds")

    # --- Quantum Solver Workflow ---
    print("\n--- Running Quantum Solver ---")
    
    # Initialize phi with the initial u0
    current_phi_vector = get_initial_phi_statevector(u0_initial, nu, dx)

    # Store benchmark results
    benchmark_results = []
    
    # Define noise model for noisy simulations
    noise_model = NoiseModel()
    error_1_qubit = depolarizing_error(0.001, 1) # 0.1% error rate on single-qubit gates
    noise_model.add_all_qubit_quantum_error(error_1_qubit, ['rz', 'rx', 'ry'])
    error_2_qubit = depolarizing_error(0.01, 2) # 1% error rate on two-qubit gates
    noise_model.add_all_qubit_quantum_error(error_2_qubit, ['cx', 'rzz'])
    
    # Define noise scaling factors for ZNE
    noise_factors = [1.0, 2.0, 3.0]
    
    # Loop through multiple quantum time steps
    for q_step in range(num_quantum_steps):
        current_sim_time = (q_step + 1) * dt_quantum_step
        print(f"\n  Quantum Step {q_step + 1}/{num_quantum_steps} (Simulated Time: {current_sim_time:.4f})")
        
        # Re-normalize the statevector before StatePreparation
        current_phi_vector = current_phi_vector / np.linalg.norm(current_phi_vector)
        
        # 1. Prepare initial state for phi from the current phi_vector
        initial_phi_circuit = StatePreparation(current_phi_vector) 
        
        # 2. Create the base quantum circuit for one quantum time step
        base_qc = QuantumCircuit(num_qubits)
        base_qc.append(initial_phi_circuit, base_qc.qubits)
        evolution_circuit = create_time_evolution_circuit(num_qubits, nu, dt_quantum_step, num_trotter_steps_per_quantum_step, dx)
        base_qc.compose(evolution_circuit, inplace=True)
        base_qc.save_statevector()
        
        # --- Run noiseless simulation for a perfect reference ---
        simulator_noiseless = AerSimulator()
        transpiled_qc_noiseless = transpile(base_qc, simulator_noiseless)
        job_noiseless = simulator_noiseless.run(transpiled_qc_noiseless, shots=1)
        result_noiseless = job_noiseless.result()
        final_phi_vector_noiseless = np.real(result_noiseless.get_statevector(base_qc).data)
        u_quantum_noiseless = post_process_results(final_phi_vector_noiseless, dx, nu)
        
        # --- Perform Zero Noise Extrapolation ---
        print(f"    Running Zero Noise Extrapolation with factors: {noise_factors}")
        
        u_vectors_noisy_for_zne = []
        for factor in noise_factors:
            print(f"      Running simulation for noise factor: {factor}")
            # Scale the circuit by repeating two-qubit gates
            scaled_qc = scale_circuit_noise(base_qc, factor)
            
            # Use the noisy simulator
            simulator_noisy = AerSimulator(noise_model=noise_model)
            transpiled_qc_noisy = transpile(scaled_qc, simulator_noisy)

            # NOTE: For ZNE with statevector, we can still use shots=1.
            # For real hardware, you would need to run with many shots and
            # reconstruct the state from counts.
            job_noisy = simulator_noisy.run(transpiled_qc_noisy, shots=1) 
            result_noisy = job_noisy.result()
            
            # Get the final statevector from the noisy simulation
            final_phi_vector_noisy = np.real(result_noisy.get_statevector(scaled_qc).data)
            u_noisy_result = post_process_results(final_phi_vector_noisy, dx, nu)
            u_vectors_noisy_for_zne.append(u_noisy_result)
        
        # Extrapolate to zero noise
        u_quantum_zne = extrapolate_zne_results(u_vectors_noisy_for_zne, noise_factors)
        
        # Update current_phi_vector for the next quantum step using the noiseless result
        current_phi_vector = final_phi_vector_noiseless

        # --- Benchmark Metrics for this step ---
        # Get classical reference for this specific time step
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
