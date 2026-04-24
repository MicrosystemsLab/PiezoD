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
│   └── sampleCode.m             # Usage example
├── python/                      # Python port (in progress)
│   ├── src/piezod/              # Core library
│   │   ├── cantilever.py        # Base class
│   │   ├── cantilever_diffusion.py
│   │   ├── cantilever_epitaxy.py
│   │   ├── cantilever_implantation.py
│   │   ├── cantilever_piezoelectric.py
│   │   ├── cantilever_poly.py
│   │   └── data/                # Bundled data files
│   │       ├── ionImplantLookupTable_tsuprem.h5
│   │       └── ionImplantLookupTable_dopedealer.h5
│   ├── tests/                   # Test suite
│   ├── examples/                # Usage examples
│   └── archive/                 # Legacy/experimental scripts
├── lookupTableGeneration/       # FLOOXS TCAD evaluation (archived, not viable)
│   ├── Dockerfile               # FLOOXS build instructions
│   ├── docker-compose.yml       # Container service definition
│   ├── templates/               # FLOOXS simulation templates
│   ├── simulations/             # Input/output directory
│   ├── legacy/                  # TSUPREM-4 reference files
│   └── lookupTable.mat          # Generated source data
├── Docs/                        # Website (jemdoc)
│   ├── html/                    # Generated HTML
│   └── cgi-bin/                 # Download scripts
└── Releases/                    # Version archives
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
- KLayout
	- pcell: use python not ruby

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

- PyPI releases go through CI, never from a local machine. There are no
  credentials on local machines; PyPI uses trusted publishing (OIDC) from
  GitHub Actions.
- Plan of record: invoke the `/release` skill. It bumps the version,
  creates a GitHub Release, and the release event triggers the
  `publish.yml` workflow which builds and uploads to PyPI.
- Equivalent manual trigger: `gh workflow run release.yml -f version=X.Y.Z`
  from the repo root.
- Do not run `uv publish` or `twine upload` locally. There is no local
  auth path; it will fail, and the intent is to keep the release pipeline
  auditable and reproducible via CI only.
- Workflow files: `.github/workflows/release.yml` (create release + tag),
  `.github/workflows/publish.yml` (build + upload to PyPI on release event).

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

## Hardware

- Read project-specific docs and datasheets before making suggestions
- Use exact part numbers, not generic equivalents
- Prioritize safety over convenience; err on side of caution

## Plotting

- For python, include `addcopyfighandler` for easy matplotlib copy/paste
