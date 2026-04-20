# Discontinuous Galerkin Methods for Convection-Dominated Problems

## Overview

This project implements a discontinuous Galerkin (DG) method for solving one-dimensional convection-dominated transport equations, with a particular focus on systems exhibiting shock formation and propagation.

Conservation laws underpin many physical systems, leading naturally to transport equations governing the evolution of quantities such as mass, momentum, or energy. In convection-dominated regimes, analytical solutions are often unavailable—particularly in the presence of non-linear fluxes or discontinuities—necessitating robust numerical methods.

The discontinuous Galerkin method provides a flexible framework that combines high-order accuracy in smooth regions with stability in the presence of shocks. This project explores its implementation, behaviour, and limitations across a range of test problems.

---

## Mathematical Model

We consider conservation laws of the form:

u_t + f(u)_x = 0

with particular focus on:

- Linear advection: f(u) = u  
- Inviscid Burgers’ equation: f(u) = 0.5 u²  

These systems exhibit fundamentally different behaviour:
- Linear transport preserves smoothness
- Nonlinear transport leads to shock formation

---

## Numerical Method

The implementation uses a discontinuous Galerkin formulation with:

- Piecewise linear polynomial approximation (2 DOFs per element)
- Weak formulation with numerical flux coupling between elements
- Two-point Gaussian quadrature for flux integration
- Periodic boundary conditions

### Fluxes

- Lax–Friedrichs-type numerical flux
- Flux derivative used to control stability

### Time Integration

- Forward Euler (baseline)
- Fourth-order Runge–Kutta (RK4)

### Shock Handling

- Generalised minmod limiter
- Suppresses oscillations near discontinuities
- Tunable parameter balancing stability and numerical diffusion

---

## Project Structure

```
.
├── src/
│   ├── Galerkin.cpp
│   ├── mvector.h
│   ├── mmatrix.h
│
├── results/
│   ├── burgersin.pdf
│   ├── burgersquare.pdf
│   ├── MCLimiter.pdf
│   ├── Erroragainsthandt.jpeg
│
├── Discontinuous_Galerkin_Methods_for_Convection_Dominated_Problems
└── README.md
```

---

## How to Run

### Requirements

- C++ compiler (e.g. g++)

### Compile

```bash
g++ -O2 main.cpp -o dg_solver
```

### Run

```bash
./dg_solver
```

The `main()` function contains a collection of test cases (commented/uncommented as required), including:

- Linear advection with sinusoidal initial condition  
- Square wave advection (with/without limiter)  
- Burgers’ equation (shock formation and propagation)  
- Error analysis experiments  

To reproduce specific experiments, uncomment the corresponding blocks in `main()`.

---

## Results

### Shock Formation (Burgers’ Equation)

A smooth sinusoidal initial condition evolves into a discontinuity due to nonlinear steepening. The method captures the onset of shock formation and its subsequent evolution.

### Shock Propagation

Square wave initial conditions demonstrate how discontinuities propagate:

- Downward jumps remain sharp  
- Upward jumps diffuse over time  

### Limiter Effects

Without a limiter, the scheme exhibits Gibbs-type oscillations near discontinuities.  
The minmod limiter suppresses these oscillations at the cost of slight smoothing.

### Error Analysis

Error behaviour is analysed as a function of:

- Spatial resolution (h)  
- Time step (Δt)  

Results confirm:

- Stability governed by the CFL condition  
- Coupling between spatial and temporal discretisation orders  

---

## Key Observations

- Increasing polynomial order does not improve accuracy near discontinuities  
- Limiters are essential for physically meaningful solutions  
- Nonlinear fluxes naturally generate shocks even from smooth initial data  
- Numerical diffusion can be controlled via both limiter and flux parameters  

---

## Notes

- Code is structured around an element-based design with neighbour coupling  
- Implementation is easily extendable to higher-order polynomials or alternative fluxes  
- Suitable for parallelisation due to local element interactions  

---

## Author

Jacek Schmidt
jacek14schmidt@gmail.com

---
