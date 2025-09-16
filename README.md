# Parameter Reconstruction in CD-SAXS with Full-Wave Maxwell Solver and Bayesian Optimization

## Author: Ligang Sun

---

## Overview

We couple a rigorous full-wave finite-element Maxwell forward solver with Bayesian optimization algorithms to efficiently reconstruct grating geometry from critical dimension small-angle X-ray scattering (CD-SAXS), then quantify parameter uncertainties with Markov chain Monte Carlo (MCMC) algorithms. The approach makes high-fidelity inversion computationally viable by using Gaussian-process surrogates, so only a few dozen expensive forward solves are required.

---

## What’s included
- a forward model based on JCMsuite to compute intensity in the reciprocal space.
- scripts to run Bayesian optimization (BO) on a full-wave forward model and to run MCMC on the trained surrogate.
- example that reconstructs 5 grating parameters (CD, height, SWA, r_top, r_bot) for a tungsten grating.
- short notebook / plotting utilities to visualize convergence and posterior samples.

---

## Requirements
- JCMsuite — commercial FEM Maxwell solver (required). Download and documentation: https://jcmwave.com/  
- JCMoptimizer — optimization modules (required) for Bayesian optimization and MCMC analysis: https://optimizer.jcmwave.com/ 
- Python 3.8+

---


## Remarks

Please acknowledge this work by citing the following article:

*Sun, L., Fukuto, M., Yager, K. G., 2025, September. Parameter Reconstruction in CD-SAXS with Full-Wave Maxwell Solver and Bayesian Optimization.
In Photomask Technology + Extreme Ultraviolet Lithography, SPIE.*

Comments and questions: lsun1@bnl.gov

