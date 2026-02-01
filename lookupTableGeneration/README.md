# Lookup Table Generation

Generate dopant profile lookup tables for PiezoD using FLOOXS TCAD process simulation.

## Prerequisites

- Docker (with Compose)
- FLOOXS source code (free for academic/research use, requires registration)

Request access: http://www.flooxs.ece.ufl.edu/index.php/Download

## Setup

1. Download FLOOXS source code after registration
2. Extract to `FLOOXS_2026/` in this folder (gitignored, not committed)
3. Build the Docker image (one time):

```bash
docker compose build
```

## Usage

Run a simulation:

```bash
docker compose run --rm flooxs input.tcl
```

Input/output files go in `simulations/` which is mounted into the container.

Interactive shell (for debugging):

```bash
docker compose run --rm flooxs /bin/bash
```

## Workflow

1. Generate inputs - Python script creates FLOOXS input files in `simulations/`
2. Run simulations - `docker compose run` executes FLOOXS on each input
3. Post-process - Python/MATLAB extracts profiles and builds lookup table

## Parameter Space

| Parameter | Values |
|-----------|--------|
| Dopants | Boron, Phosphorus, Arsenic |
| Doses | 2e14, 2e15, 2e16 cm^-2 |
| Energies | 20, 50, 80 keV |
| Anneal times | 15, 30, 45, 60, 75, 90, 105, 120 min |
| Anneal temps | 900, 1000, 1100 C |
| Oxidation | With/without passivation oxide |

Total: 3 x 3 x 3 x 8 x 3 x 2 = 1296 simulations

## Files

| File | Description |
|------|-------------|
| `docker-compose.yml` | Container service definition |
| `Dockerfile` | FLOOXS build instructions |
| `simulations/` | Input/output directory (mounted in container) |
| `postProcessTables.m` | MATLAB post-processor (legacy) |
| `simulation.template` | TSUPREM-4 template (legacy, needs conversion) |
| `simulationControl.py` | Python 2 runner (legacy, needs rewrite) |
