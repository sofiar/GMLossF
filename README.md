# Gompertz Model with observation Error 

## Context 
Stochastic population dynamics models are a crucial component of applied and theoretical ecology. Statistical inference for these models can be challenging, especially when both process error and sampling error are present. Ignoring sampling variability can lead to biased estimations and erroneous conclusions about system behavior. The Gompertz model is widely used to describe the growth of animals, plants, or cells, and it can be adapted to account for sampling variability.

This project aims to develop an inference method for estimating the parameters of the model when sampling error follows a Poisson distribution.

## Description

This repository serves two main purposes: it contains the `gse` library for applying the methods described in the manuscript *Likelihood-based inference for the Gompertz model with Poisson errors*, as well as the code to reproduce all analyses presented there. 


## Contents 
* 📁 `gse`: R library. 
* 📁 `Real_data`: Scripts to analyze the American Redstart counts dataset.
* 📁 `Simulation_studies`: Script to run the simulation analysis presented in the paper 


## Overview of the `gse` library 
- [Installation](#installation)
- [Quick Start](#quick-start)

### **Installation**

```r
devtools::install_github("sofiar/GMLossF/gse")
```

### **Quick Start**
```r
library(gse)

# simulate model 
data = GompPois_rng(T = 100, theta1 = 1.92, theta2= 0.22, b =-0.24)

# fit using methods of moments
GompPois_MoM(data)

# Get posterior simulation using elliptical slice sampling within Gibbs
# under flat and mixture of normal-inverse gamma priors
GompPois_flatMixNIG_mcmc(nsim = 1e+4, Nstar = data, kappa = 1, psi1 = 0.5, psi2 =0.5 , 
                        nu = 2 , zeta1 = 0, zeta2 = 1, starter = NULL, burn = 1,
                        thin = 1, verbose = +Inf)

# Get posterior simulation using elliptical slice sampling within Gibbs 
# under flat and normal-inverse gamma priors
GompPois_flatNIG_mcmc(nsim = 1e+4, Nstar =data , phi1 = 0.5, phi2 = 0.5, 
                      eta1 = 0, eta2 = 1, starter = NULL, burn = 1, thin = 1,
                      verbose = +Inf)


```