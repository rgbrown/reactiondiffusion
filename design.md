# Code Design

We want to be able to solve 1D reaction diffusion systems with the
following features

* Any number of the variables diffusing
* Spatially variable diffusion coefficients
* Spatially variable parameters
* Periodic, zero-flux, or Dirichlet boundary conditions

We also want the user to be able to choose between finite-difference and
Spectral methods. Further if using a spectral method, it should
automatically use the Fourier method if we are using periodic boundary
conditions, and Chebyshev method otherwise.

We will impose a [-1, 1] computational domain for each method

## Problem specification
User specifies one of the following boundary condition types. For each
variable that has diffusion, we need to specify boundary conditions. These
will be of the following

* *Periodic* This encompasses both ends of the domain, obviously
* *Dirichlet* a time-varying Dirichlet boundary condition. For each end:
    + Value ul, ur, or, 
    + function of time ul = fl(t), ur = fr(t)
    + A mixture of the two is fine
* *Neumann* For now, simply zero-flux

If Dirichlet, there is a 



