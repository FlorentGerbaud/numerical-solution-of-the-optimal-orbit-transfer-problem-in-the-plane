# Numerical Solution of the Optimal Orbit Transfer Problem in the Plane

## Project Overview

This project aims to solve the optimal orbit transfer problem in the plane using numerical methods. The objective is to transfer a satellite from an initial orbit to a final orbit with minimal fuel consumption by controlling the thrust direction and magnitude.

## Table of Contents

- [Project Overview](#project-overview)
- [Getting Started](#getting-started)
- [Prerequisites](#prerequisites)
- [Implementation Details](#implementation-details)
  - [Theory](#theory)
  - [Implementation](#implementation)
- [Results](#results)

## Getting Started

To run this project, ensure that you have the necessary software and dependencies installed as listed in the [Prerequisites](#prerequisites) section. The project is designed using Scilab for solving differential equations and optimization problems.

## Prerequisites

- Scilab (version 6.0 or later)
- A basic understanding of numerical optimization and orbit dynamics
- Mathematical libraries for solving ordinary differential equations (ODEs)

## Implementation Details

### Theory

The problem involves determining the optimal control law for thrust that minimizes the fuel usage while transferring a satellite from one orbit to another. The Hamiltonian dynamics of the system are modeled, and boundary conditions at the initial and final times are defined to form a two-point boundary value problem.

Key equations used:

- Hamiltonian formulation and adjoint state dynamics.
- Control equation for the thrust angle.
- System of equations to solve for the unknown variables at different time points.

### Implementation

The project is implemented in several steps:

1. **Problem Data Declaration**: Defines initial and final conditions such as initial orbit, mass, velocity, and thrust parameters. The data is normalized for accuracy.
2. **System Dynamics**: The dynamics of the spacecraft are modeled using a set of differential equations that account for thrust, gravitational forces, and changes in velocity and position over time.
3. **Numerical Resolution**: Uses numerical methods (ODE solver and optimization) to compute the optimal trajectory and control history for thrust values ranging from 0.1N to 0.6N.

Key Functions:

- `dynpol`: Calculates the dynamics of the spacecraft at each time step.
- `gnultmin`: Solves the optimization problem by minimizing the Hamiltonian and enforcing boundary conditions.
- **Control History Calculation**: Computes the thrust control angle and plots the control history over time.

## Results

The results of the simulation are plotted for different thrust values ranging from 0.1N to 0.6N. The control history shows the thrust angle evolution over time, and the transfer duration decreases with higher thrust values. Example figures show the thrust profiles for various thrust levels.

