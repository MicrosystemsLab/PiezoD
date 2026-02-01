# piezod

Modeling and optimization of piezoresistive and piezoelectric sensors and actuators.

## Installation

```bash
pip install piezod
```

For development:

```bash
git clone https://github.com/MicrosystemsLab/PiezoD.git
cd PiezoD/python
uv sync
```

## Usage

```python
from piezod import CantileverEpitaxy

# Create cantilever and set geometry
c = CantileverEpitaxy()
c.l = 300e-6  # length: 300 um
c.w = 44e-6   # width: 44 um
c.t = 89e-9   # thickness: 89 nm
c.fluid = "water"

# Calculate properties
print(f"Stiffness: {c.stiffness() * 1e3:.3f} mN/m")
print(f"Resonant frequency: {c.omega_vacuum_hz() / 1e3:.1f} kHz")
freq_hz, Q = c.omega_damped_hz_and_Q()
print(f"Damped frequency: {freq_hz / 1e3:.1f} kHz, Q = {Q:.1f}")
```

Run the example:

```bash
uv run python examples/quickstart.py
```

## Development

Install dependencies:

```bash
uv sync
```

Run tests:

```bash
uv run pytest
```

Lint and format:

```bash
uvx ruff check src tests --fix && uvx ruff format src tests
```

Type check:

```bash
uvx ty check src
```

## Documentation

See the main project documentation at [microsystems.stanford.edu/piezod](http://microsystems.stanford.edu/piezod).

## Citation

If you use piezod in your research, please cite:

> Joseph C. Doll, Sung-Jin Park and Beth L. Pruitt
> Design optimization of piezoresistive cantilevers for force sensing in air and water
> Journal of Applied Physics 106.6 (2009): 064310-064310.

## License

GPL-3.0-or-later
