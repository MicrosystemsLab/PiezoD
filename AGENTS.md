# Overview

Piezoresistor design and optimization tool.

## Repository layout

```
PiezoD/
├── matlab/                      # MATLAB implementation (primary, complete)
│   ├── PiezoD/                  # Core library
│   │   ├── cantilever.m         # Base class
│   │   ├── cantileverDiffusion.m
│   │   ├── cantileverEpitaxy.m
│   │   ├── cantileverImplantation.m
│   │   ├── cantileverPiezoelectric.m
│   │   ├── cantileverPoly.m
│   │   ├── ionImplantLookupTable.mat  # Ion implant lookup data
│   │   └── sensorSimulation.mdl # Simulink model
│   ├── tests/                   # MATLAB test suite (4 cantilever*Test.m files)
│   ├── sampleCode.m             # Usage example
│   └── README.md
├── python/                      # Python port (in progress)
│   ├── src/piezod/              # Core library
│   │   ├── cantilever.py        # Base class
│   │   ├── cantilever_diffusion.py
│   │   ├── cantilever_epitaxy.py
│   │   ├── cantilever_implantation.py
│   │   ├── cantilever_piezoelectric.py
│   │   ├── cantilever_poly.py
│   │   ├── default.mplstyle
│   │   └── data/                # Bundled data files
│   │       ├── ionImplantLookupTable_tsuprem.h5
│   │       └── ionImplantLookupTable_dopedealer.h5
│   ├── tests/                   # pytest suite
│   ├── examples/                # Usage examples (quickstart.py, etc.)
│   ├── archive/                 # Legacy/experimental scripts
│   ├── pyproject.toml
│   └── README.md
├── lookupTableGeneration/       # FLOOXS TCAD evaluation (archived, not viable
│                                # as a TSUPREM-4 replacement). See its
│                                # README.md and docs/status.md for details.
├── docs/
│   └── tutorial.md              # User-facing tutorial
└── .github/workflows/           # CI and release automation
    ├── ci.yml                   # Tests + lint on push/PR
    ├── release.yml              # Bump version, create GitHub Release (manual)
    └── publish.yml              # Build wheel, upload to PyPI (on release)
```

## Tech Stack

- Python
	- Environment: `uv`
	- Lint: `ruff`
	- Test: `pytest`
	- Type checking: `ty`
	- Build: `hatchling`
- Matlab
	- Lint: `checkcode`
	- Test: `runtests`

## Commands

- Python
	- Install: `uv sync` from pyproject.toml
	- Lint: `uvx ruff check [path] --fix --extend-select I,B,SIM,C4,ISC,PIE && uvx ruff format [path]`
	- Test: `uv run pytest`
	- Type check: `uvx ty check [path]`
- Matlab
	- Test: `matlab -batch "runtests('tests')"

## Safety

- Never run `git checkout` (destructive)
- Never use unicode in code (breaks Python on Windows)

## Workflow

- Plan first, share plan, get explicit approval before implementing
- Never assume a good plan means you should implement it
- If unsure or stuck, ask the user
- If something fails, do not silently move on - ask for clarification
- If `uv` fails due to permissions, consult user

## Testing

- Run tests after changes; all tests must pass before task is complete
- Prefer tests over manual verification

## Release

- Use `.claude/skills/release/skill.md` for every release.
- Do not substitute individual workflow or manual commands for the release skill.
- A release is not complete until the GitHub Release has user-facing release notes,
  the publish workflow succeeds, and PyPI shows the new version.
- Never publish locally with `uv publish` or `twine upload`; PyPI uses CI trusted publishing.

## Code Quality

- Use type hints and docstrings for all functions
- Keep files under 300 LOC; split if longer
- Fix all linter issues

## Parameters

- Single source of truth for defaults (e.g. config.py)
- Use enums instead of bare strings
- Never use `params.get('key', default)` pattern
- Never use `if config then value else default` pattern

## Refactoring

- Delete old interfaces; no legacy wrappers or thin compatibility layers

## Documentation

- No emojis in code or docs
- No excessive bold in markdown; use styling selectively
- No "last updated" dates or authorship
- Commands on single line, not split across lines

## Plotting

- For python, include `addcopyfighandler` for easy matplotlib copy/paste
