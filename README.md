# RDSA
Random direction stochastic approximation for simulation optimization
=====================================================================
											August, 2015
											-------------
Contents
--------
1 Introduction                                                                                                  
2 Notes on usage	
3 References
                                                       
1 Introduction
--------------
This software package (in Matlab) provides a class of algorithms for simulation optimization using random direction stochastic approximation (RDSA). 
These include first-order (gradient) as well as second-order (Newton) schemes. The algorithms incorporate both continuous-valued as well as discrete-valued perturbations into both our algorithms. The former are chosen to be independent and identically distributed (i.i.d.) symmetric uniformly distributed random variables (r.v.), while the latter are i.i.d., asymmetric Bernoulli r.v.s. See [3] for a detailed description.

2 Notes on Usage
----------------
The main files in the distribution are:

First-order schemes:
-------------------
i) onespsa.m --> this contains the implementation of first-order SPSA algorithm [1] and is a simplified version of that available on J.C. Spall's ISSO book website http://www.jhuapl.edu/ISSO

ii) onerdsa_unif.m --> An RDSA variant of the 1SPSA code from J.C. Spall. The primary difference is in the generation of perturbation r.v.s. In this case, the latter are sampled from an uniform [-1,1] distribution.

iii) onerdsa_asymber.m --> Similar to 1RDSA-Unif, except that the perturbation r.v.s. follow an asymmetric Bernoulli distribution.

Second-order schemes:
---------------------
i) twospsa.m --> this contains the implementation of second-order SPSA algorithm [2] and is a simplified version of that available on J.C. Spall's ISSO book website 

ii) twordsa_unif.m --> An RDSA variant of the 2SPSA code from J.C. Spall. The primary difference is in the generation of perturbation r.v.s. In this case, the latter are sampled from an uniform [-1,1] distribution.

iii) twordsa_asymber.m --> Similar to 2RDSA-Unif, except that the perturbation r.v.s. follow an asymmetric Bernoulli distribution.

On input parameters: Most of the algorithms above take as input the following:
p -> dimension of the problem
sigma -> noise parameter. Noise is (p+1)-dimensional Gaussian with variance sigma^2
type -> 1 for quadratic, 2 for fourth-order loss (see loss_myexample.m)
numSimulation -> this is the simulation budget that impacts the number of 2SPSA iterations
replications -> number of independent simulations
theta_0 -> initial point

4 References
------------
[1] J. C. Spall, "Multivariate stochastic approximation using a simultaneous perturbation gradient approximation", IEEE Trans. Auto. Cont., vol. 37, no. 3, pp. 332-341, 1992.

[2] J. C. Spall, "Adaptive stochastic approximation by the simultaneous perturbation method", IEEE Trans. Autom. Contr., vol. 45, pp. 1839-1853, 2000.

[3] Prashanth L.A., S. Bhatnagar, Michael Fu and Steve Marcus, "Adaptive system optimization using (simultaneous) random directions stochastic approximation", arXiv:1502.05577, 2015.
