# IISPH Fluid Solver Implementation

This IISPH fluid solver code was developed as part of my Master's thesis "Various Boundary
Handling Techniques for Incompressible SPH".

Due to the limited time frame, the focus was on numerical correctness and performance, not software architecture. For this repository, I improved the structure and documentation.

## Demo

![Demo](docs/demo.gif)

## Features

This code implements:
* a particle-based [IISPH](https://cg.informatik.uni-freiburg.de/publications/2013_TVCG_IISPH.pdf) fluid solver (2D)
* neighborhood search using an Index Sort variant of [Compressed Neighbor Lists for SPH (Band et al. 2018)](https://cg.informatik.uni-freiburg.de/publications/2019_CGF_CompressedNeighbors.pdf)
* boundary handling: [pressure mirroring (Akinci et al. 2012)](https://cg.informatik.uni-freiburg.de/publications/2012_SIGGRAPH_rigidFluidCoupling.pdf)
* semi-implicit time integration



The solver
* numerically solves the Navier-Stokes equations
* strongly enforces incompressibility by solving a discretized Poisson pressure equation using a relaxed Jacobi method
* uses the [SFML](https://www.sfml-dev.org/) library for visualization


## Build

```bash
mkdir build
cd build
cmake ..
make
