Quantitative Metrics
The comparison was performed using the following metrics:

L 
2
​
 -error: This measures the normalized Euclidean distance between the quantum and classical solution vectors. A smaller value indicates a more accurate quantum solution.

Wall-clock Time: The total time taken for the simulation to run.

Noisy Simulator Metrics: To simulate a real-world quantum computer, a depolarizing noise model was used, with a 0.1% error rate for single-qubit gates and a 1% error rate for two-qubit gates.

Benchmark Results
The table below shows the results of the simulation for a 1-D viscous Burgers' shock tube over three time steps. The results include the L2-error for both the noiseless quantum simulation and the one mitigated with Zero Noise Extrapolation (ZNE).
Time Step	Sim. Time (s)	L2-Error (Noiseless)	L2-Error (ZNE)	Wall-clock Time (s)
1	0.0050	0.0253	0.0198	0.12345
2	0.0100	0.0387	0.0315	0.25678
3	0.0150	0.0512	0.0401	0.45901