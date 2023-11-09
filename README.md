# HD-QPSK-CVQKD
Hello! Here are MATLAB-based code packages with functions for calculating the secret key rate when excuting the homodyne detection quadrature phase shift keying continuous-variable quantum key distribution protocol.

The preprint that introduces this protocol is given on PRX Quantum 2, 040334 (arXiv with Num. 2104.11152).

Code 1 and Code 2 are two code packages based on this protocol. They are alternative choices if you want to simulate this protocol.

# Introduction of Code 1
There is one controller and nine functions:

-> Controller - afin_control.m

This is a function with the key rate as the final output. The inputs are experimental parameters.

-> Function - S_1.m

This is a function to generate all Hermitian matrices in the constraints.

-> Function - S_2.m

This is a function to calculate the value of each constraint.

-> Function - S_3.m

This is a function to generate a primal density matrix \rho_{AB} under the constraints.

-> Function - S_4.m

This is a function to generate the matrices of interval operators in the photon number representation.

-> Function - S_5.m

This is a function to calculate \Del\rho, which is the solution of the primal problem in the step 1.

-> Function - S_6.m

This is a function to search the parameter \lambda that makes the key rate minimum.

-> Function - S_7.m

This is a function to solve the dual problem in step 2.

-> Function - S_8.m

This is a function to calculate the error correction for each quadrature.

-> Function - abstractep.m

This is a function to gain the parameter ep2 used in dual problem.

# Introduction of Code 2
There is one controller and eleven functions:

-> Controller - protocol_simulation.m

This is a script with the key rate as the final result. It sets the parameters and controls other functions.

-> Function - K_pq_new.m

This is a function to generate Kraus operators in the key rate formula.

-> Function - Gamma.m

This is a function to calculate the constraint matrices.

-> Function - gamma_exp.m

This is a function to calculate the value of each constraint.

-> Function - search_rho0.m

This is a function to generate a primal density matrix \rho_{AB} under the constraints.

-> Function - algorithm1.m

This is a function to do the step 1 (primal problem).

-> Function - abstractep.m

This is a function to gain the parameter eps used in dual problem.

-> Function - G_perturb_deg.m

This function is used to calculte the hermitian conjugate of the perturbed G.

-> Function - G_perturb.m

This function is used to calculte the perturbed G.

-> Function - max_gamma_y.m

This is a function to solve the dual problem in step 2.

-> Function - rel_ent.m

This is a function to calcualte the relative entropy of the inputs.

-> Function - error_correction.m

This is a function to calculate the error correction for each quadrature.
