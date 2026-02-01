# piezod

Modeling and optimization of piezoresistive and piezoelectric sensors and actuators.

## Installation

```bash
pip install piezod
```

For development:

```bash
git clone https://github.com/certus/PiezoD.git
cd PiezoD/python
uv sync
```

## Usage

```python
from piezod import Cantilever, CantileverEpitaxy

c = CantileverEpitaxy()
c.doping_profile()
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
uvx ruff check src tests --fix --extend-select I,B,SIM,C4,ISC,PIE && uvx ruff format src tests
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
