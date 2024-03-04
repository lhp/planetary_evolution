This repository contains code which simulates the evolution of a simplified planet. The planet has a mantle and a core. The mantle is a rectangle with an isothermal cold surface. The mantle may contain multiple layers: each layer may be volumetrically heated and may be mobile or not. Volumetric heating is uniform within a layer and decays exponentially over time. Within the mobile layer(s), temperature perturbations may result in the development of motion (thermal convection with a Boussinesq approximation). The numerical solutions of the coupled momentum and heat transport problems are implemented using Lattice-Boltzmann methods. The core is isothermal and is coupled to the mantle at the core-mantle boundary by continuity of temperature and heat flux. 

A description of files.

1. Setup and execution
- run_directly.m: setup of simulation, entry point.
- setup_model.m: creates objects used to model the planet.
- run_guts.m: manages execution of simulation, including visualization and data output.
- Scaling.m, PhysicalProperties.m, ScaledPhysicalProperties.m: manages translation between lattice properties and another reference frame (e.g., real-world properties).

2. Implementation
- PlanetModel.m: couples mantle transport solvers and core.
- Core.m: implementation of the core
- AdvectionDiffusion.m: implementation of heat transport problem.
- NavierStokes.m: implementation of momentum transport problem.
- settle_velocities.m: algorithm for finding steady-state velocities, including managing scaling of forces.

3. Output
- ModelHistory.m: tracks results of simulation.
- Visualizer.m: visualizes simulation during execution.

The coupled momentum and heat transport solvers are capable of modeling thermal convection; we verified the transport properties of the modeled system against the benchmark values of Blankenbach et al. (1989). Our implementation provides an excellent match, illustrated below.

![Ra_Nu](https://github.com/lhp/planetary_evolution/assets/1531596/35f630d2-711e-4295-8eda-69f1d2130053)
