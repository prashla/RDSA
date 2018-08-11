# RDSA
Random direction stochastic approximation for simulation optimization
=====================================================================
											August, 2018
											-------------
Contents
--------
1 Introduction                                                                                                  
2 Notes on usage 														    
3 References
                                                       
1 Introduction
--------------
This software package provides a class of algorithms for simulation optimization using random direction stochastic approximation (RDSA). The package is in two part: a matlab part that implements RDSA with random and deterministic perturbations schemes on five synthetic functions and a traffic simulation part that incorporates RDSA schemes for adaptively tuning the thresholds of a traffic light control (TLC) algorithm.

The implementation includes first-order (gradient) as well as second-order (Newton) schemes. The RDSA variants with random perturbations incorporate both continuous-valued as well as discrete-valued perturbations. The former are chosen to be independent and identically distributed (i.i.d.) symmetric uniformly distributed random variables (r.v.), while the latter are i.i.d., asymmetric Bernoulli r.v.s. See [3] for a detailed description. While the RDSA variants with deterministic perturbations uses two choices for perturbations - a semi-lexicographic sequence and a perturbation matrix-based sequence. See [4] for a detailed description.

2a Matlab part
--------------
The main files in the distribution are:

First-order schemes:
-------------------
i) onespsa.m --> this contains the implementation of first-order SPSA algorithm [1] and is a simplified version of that available on J.C. Spall's ISSO book website http://www.jhuapl.edu/ISSO

ii) onerdsa_unif.m --> An RDSA variant of the 1SPSA code from J.C. Spall. The primary difference is in the generation of perturbation r.v.s. In this case, the latter are sampled from an uniform [-1,1] distribution.

iii) onerdsa_asymber.m --> Similar to 1RDSA-Unif, except that the perturbation r.v.s. follow an asymmetric Bernoulli distribution.

iv) onerdsa_lex_dp.m --> An 1RDSA variant with deterministic perturbation. The primary difference is in the generation of perturbation r.v.s. Deterministic Perturbations are generated from lexicographic sequence.

v) onerdsa_perm_dp.m --> Similar to 1RDSA-Lex-DP, except that the deterministic perturbations r.v.s. are generated from permuatation matrix.

vi) onerdsa_kw_dp.m --> Similar to 1RDSA-Perm-DP, except that in this case unkike Perm-DP, KW-DP algorithm independently updates the individual coordinates.


Second-order schemes:
---------------------
i) twospsa.m --> this contains the implementation of second-order SPSA algorithm [2] and is a simplified version of that available on J.C. Spall's ISSO book website 

ii) twordsa_unif.m --> An RDSA variant of the 2SPSA code from J.C. Spall. The primary difference is in the generation of perturbation r.v.s. In this case, the latter are sampled from an uniform [-1,1] distribution.

iii) twordsa_asymber.m --> Similar to 2RDSA-Unif, except that the perturbation r.v.s. follow an asymmetric Bernoulli distribution.

iv) twordsa_lex_dp.m --> An 2RDSA variant with deterministic perturbation. The primary difference is in the generation of perturbation r.v.s. Deterministic Perturbations are generated from lexicographic sequence.

v) twordsa_perm_dp.m --> Similar to 2RDSA-Lex-DP, except that the deterministic perturbations r.v.s. are generated from permuatation matrix.



On input parameters: Most of the algorithms above take as input the following:

p -> dimension of the problem

sigma -> noise parameter. Noise is (p+1)-dimensional Gaussian with variance sigma^2

type -> 1 for quadratic, 2 for fourth-order loss, 3 = Powell singular function, 4 = Rosenbrock function, 5 = Rastrigin function (see loss_myexample.m)

numSimulation -> this is the simulation budget that impacts the number of 2SPSA iterations

replications -> number of independent simulations

theta_0 -> initial point

2b Traffic simulation part
-------------------------

This part is in Java and provides the implementation of adaptive threshold tuning algorithms based on RDSA. 
In practice, obtaining exact queue length information is difficult, but one can obtain coarse estimates along the lanes of the road network, for instance, by placing magnetic sensor loops at some distance from the junction. The challenge is to choose the optimal locations for placing sensor loops to infer congestion information for any lane in the road network considered and the threshold tuning algorithms cater to this need.

This package is based on the source code of the Green Light District (GLD) traffic simulator [2]. GLD codebase is modified to include the RL based TLC algorithm. The files relevant to RDSA schemes in the distribution are: 

i) src.gld.tt.RDSARunner.java --> Wrapper for running the 2RDSA based threshold tuning algorithm

ii) gld.tt.TwoRDSA_OuterLoop.java --> Implementation of the  threshold tuning outer loop that is based on 2RDSA 

iii) gld.algo.tlc.PTLCL12345T123 --> Simulates traffic for a given threshold parameter with a simple priority based traffic light control scheme

On input parameters for 2RDSA based schemes (to be set in the main function of RDSARunner.java): 

i) trainingBudget -> the number of function evaluations in the training phase 

ii) numReplicationsForTesting -> After the training phase, the thresholds are fixed and then a number of independent simulations (that this variable is set to) are run and the empirical average cost from each simulation is recorded 

iii) trajectoryLengthTrainingPhase -> length of each simulated trajectory during training phase

iv) trajectoryLengthTestingPhase --> length of each simulated trajectory during testing phase

See Section IV-E of [3] for a detailed description.

3 References
------------
[1] J. C. Spall, "Multivariate stochastic approximation using a simultaneous perturbation gradient approximation", IEEE Trans. Auto. Cont., vol. 37, no. 3, pp. 332-341, 1992.

[2] J. C. Spall, "Adaptive stochastic approximation by the simultaneous perturbation method", IEEE Trans. Autom. Contr., vol. 45, pp. 1839-1853, 2000.

[3] Prashanth L.A., S. Bhatnagar, Michael Fu and Steve Marcus, "Adaptive system optimization using (simultaneous) random directions stochastic approximation", arXiv:1502.05577, 2015.

[4] Prashanth L A, Shalabh Bhatnagar, Nirav Bhavsar, Michael Fu and Steven I. Marcus, "Random directions stochastic approximation with deterministic perturbations", 	arXiv:1808.02871, 2018.

4 For more information, whom do I contact?
------------------------------------------

    Prashanth L.A. email: prashla@umd.edu
    Nirav Bhavsar email: cs17s016@smail.iitm.ac.in
