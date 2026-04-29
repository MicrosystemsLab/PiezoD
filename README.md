# PiezoD

Modeling and optimization of piezoresistive and piezoelectric sensors and actuators.

## Overview

PiezoD is a tool for modeling the performance and optimizing the design of piezoresistive and piezoelectric sensors and actuators. It is designed to make it fast and easy to design high performance microdevices. While originally developed for the design of cantilever beams, the code is modular and can be applied to a variety of problems.

## Implementations

**Python** is the primary implementation and the only version under active development.

```
pip install piezod
```

See [python/README.md](python/README.md) and [python/examples/quickstart.py](python/examples/quickstart.py) to get started.

**MATLAB** ([matlab/](matlab/)) is the original implementation, kept for historical reference. It is no longer maintained -- bug reports and feature requests should target the Python version.

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

### Design optimization constraints
- simple bounds: all design parameters, e.g. cantilever length, bias voltage
- general nonlinear constraints: power dissipation, temperature, resonant frequency, resistance, spring constant, etc

### Optimization goal functions
- integrated noise (force, displacement, surface stress)
- noise power spectral density (force, displacement, surface stress, acceleration, etc)

## Documentation

- [Tutorial](docs/tutorial.md) - Step-by-step Python walkthrough
- [MATLAB tutorial](matlab/tutorial.md) - Legacy walkthrough for the MATLAB implementation

## Contributing

Contributions are welcome. Please report issues and feature requests on [GitHub Issues](https://github.com/MicrosystemsLab/PiezoD/issues).

## Attribution

If you use PiezoD in your research, we only ask that you cite the first paper that we wrote on the subject:

Joseph C. Doll, Sung-Jin Park and Beth L. Pruitt
Design optimization of piezoresistive cantilevers for force sensing in air and water
Journal of Applied Physics 106.6 (2009): 064310-064310.

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or
  http://www.apache.org/licenses/LICENSE-2.0)
- MIT license ([LICENSE-MIT](LICENSE-MIT) or
  http://opensource.org/licenses/MIT)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in the work by you, as defined in the Apache-2.0
license, shall be dual licensed as above, without any additional terms or
conditions.
