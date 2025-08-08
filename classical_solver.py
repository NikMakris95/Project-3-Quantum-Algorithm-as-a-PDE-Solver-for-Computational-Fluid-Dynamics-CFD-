import numpy as np

def classical_burgers_solver(u0, nu, dx, dt, num_time_steps):
    """
    Classical finite-difference solver for 1D Burgers' equation.
    
    This solver uses an explicit upwind scheme for the convection term
    and a central difference scheme for the diffusion term. It is a
    reliable reference for comparison with the quantum solution.
    
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
    
    # --- Stability Check (CFL condition) ---
    # The Courant-Friedrichs-Lewy (CFL) condition is crucial for the stability
    # of explicit numerical methods. We check it here to warn if the time step
    # might be too large, which could lead to an unstable solution.
    max_u_abs = np.max(np.abs(u))
    cfl_convection = max_u_abs * dt / dx
    cfl_diffusion = nu * dt / dx**2
    
    if cfl_convection > 1.0 or cfl_diffusion > 0.5:
        print(f"Warning: CFL condition might be violated. Convection CFL: {cfl_convection:.2f}, Diffusion CFL: {cfl_diffusion:.2f}")
        print("Results may be unstable or inaccurate. Consider reducing dt or dx.")
        
    # --- Time Stepping Loop ---
    for _ in range(num_time_steps):
        u_new = u.copy()
        
        # We iterate over the interior points, as boundary points are fixed.
        for i in range(1, nx - 1):
            # --- Convection Term: u * du/dx ---
            # Using an upwind scheme based on the sign of u[i] to ensure stability.
            if u[i] >= 0:
                # Use a forward difference for the spatial derivative if velocity is positive
                convection_term = u[i] * (u[i] - u[i-1]) / dx
            else:
                # Use a backward difference for the spatial derivative if velocity is negative
                convection_term = u[i] * (u[i+1] - u[i]) / dx
            
            # --- Diffusion Term: nu * d^2u/dx^2 ---
            # Using a central difference scheme, which is typically stable for this term.
            diffusion_term = nu * (u[i+1] - 2 * u[i] + u[i-1]) / dx**2
            
            # --- Update Equation ---
            # This is the explicit forward Euler update for the value at the next time step.
            u_new[i] = u[i] - dt * convection_term + dt * diffusion_term
        
        # --- Boundary Conditions ---
        # We apply Dirichlet boundary conditions by keeping the ends fixed.
        u_new[0] = u0[0]
        u_new[nx-1] = u0[nx-1]
        
        # Update the velocity field for the next time step
        u = u_new
        
    return u

# --- Example Usage (Optional) ---
if __name__ == '__main__':
    # Define parameters for a sample shock tube problem
    nu = 0.01
    L = 1.0
    nx = 16
    dx = L / (nx - 1)
    dt = 0.0001
    total_time = 0.05
    num_time_steps = int(total_time / dt)

    x = np.linspace(0, L, nx)
    u0_initial = np.where(x <= 0.5, 1.0, 0.0) # Riemann step initial condition

    print("Running classical solver for 1D Burgers' equation...")
    u_final = classical_burgers_solver(u0_initial, nu, dx, dt, num_time_steps)
    print("Final velocity field:", u_final)
