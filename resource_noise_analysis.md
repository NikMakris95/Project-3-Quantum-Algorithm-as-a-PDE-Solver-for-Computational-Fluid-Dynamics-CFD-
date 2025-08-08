Resource & Noise Analysis
1. Resource Table
Metric

Grid Size N=16

Two-qubit-gate depth

T-count

Value

4 qubits

~30 per time step

0

2. Noise Mitigation Strategy
The primary error mitigation strategy employed is Zero Noise Extrapolation (ZNE). This technique scales the noise in the quantum circuit by repeating the primary error-prone gates (in this case, the RZZ gates). We run the scaled circuits on a noisy simulator to obtain results at different noise levels. These results are then used to fit a linear model, and the value is extrapolated back to the theoretical zero-noise limit to estimate the ideal, noiseless outcome. The prototype code implements this by running the circuit with noise scale factors of 1.0, 2.0, and 3.0.