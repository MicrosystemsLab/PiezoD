# PiezoD

Modeling and optimization of piezoresistive and piezoelectric sensors and actuators.

## Overview

PiezoD is a tool for modeling the performance and optimizing the design of piezoresistive and piezoelectric sensors and actuators. It is designed to make it fast and easy to design high performance microdevices. While originally developed for the design of cantilever beams, the code is modular and can be applied to a variety of problems.

## Implementations

| Implementation | Status | Notes |
|----------------|--------|-------|
| [MATLAB](matlab/) | Complete | Primary implementation, full feature set |
| [Python](python/) | In progress | Port of MATLAB implementation |

Future development will focus on the Python version.

## Features

### Materials
- p-type and n-type single crystal silicon piezoresistive cantilevers
- epitaxial, diffused and ion implanted piezoresistors
- aluminum nitride (AlN) and lead zirconate titanate (PZT) piezoelectric cantilevers
- polycrystalline silicon and thin metal film cantilevers

### Model components
- noise: 1/f (Hooge), Johnson, amplifier, thermomechanical
- mechanics: Bernoulli beam bending with accurate stiffness and resonant frequency for multilayer, segmented structures
- concentration and temperature dependent effects: thermal conductivity, carrier density and mobility, piezoresistance factor
- fluid damping: resonant frequency and quality factor
- temperature profile: finite difference, lumped, and linearized circuit models for self-heating

### Constraints
- simple bounds: all design parameters, e.g. cantilever length, bias voltage
- general nonlinear constraints: power dissipation, temperature, resonant frequency, resistance, spring constant, etc

### Optimization goal functions
- integrated noise (force, displacement, surface stress)
- noise power spectral density (force, displacement, surface stress, acceleration, etc)

## Documentation

Model details and example usage are available at [microsystems.stanford.edu/piezod](http://microsystems.stanford.edu/piezod).

## Attribution

If you use PiezoD in your research, we only ask that you cite the first paper that we wrote on the subject:

Joseph C. Doll, Sung-Jin Park and Beth L. Pruitt
Design optimization of piezoresistive cantilevers for force sensing in air and water
Journal of Applied Physics 106.6 (2009): 064310-064310.

## License

piezoD is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

piezoD is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
