Quantum Hardware Run
1. Hardware Execution Details
Quantum Processing Unit (QPU): [e.g., IBM Quantum Eagle, Quantinuum H1-2]

Date of Execution: [Date]

Backend Version: [e.g., ibm_washington version 1.2.3]

Mitigation Applied: [e.g., Zero Noise Extrapolation (ZNE)]

2. Raw and Error-Mitigated Results
To demonstrate the impact of noise, the circuit for one time step was executed on the physical QPU. The raw results from the noisy device show a significant deviation from the classical solution. After applying ZNE, the results show a marked improvement, bringing the quantum solution closer to the reference.

3. Diagnostics
Runtime: The wall-clock time for the QPU execution was [Time in minutes/hours].

Noise Diagnostics: The QPU backend reported average gate fidelities of [e.g., 99.7% for CNOT, 99.9% for single-qubit gates]. These metrics highlight the importance of error mitigation for obtaining meaningful results on current-generation hardware.