# TEBD
Time-evolving block decimation algorithm from quantum systems based on the G. Vidal paper 'Efficient Classical Simulation of Slightly Entangled Quantum Computations', PRL 91, 147902 (2003).
Currently the program is set to solve the 1D Ising model with transverse field with open boundary conditions. 

## How it works?
In main_tebd.m we define the parameters.We create a MPS with random variables with initial_maps.m. We start evolving the MPS with the sweep programs and at the end of two sweeps we calculate the expectation value of the energy for the actual MPS. The convergence of the energy is calculated. 
We finish by plotting the results and comparing them with the theoretical values for the OBC (only valid for g=1) and the PBC.
