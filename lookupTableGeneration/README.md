# Lookup Table Generation

Generate dopant profile lookup tables for PiezoD using FLOOXS TCAD process simulation.

## Prerequisites

### Docker

Install Docker Desktop: https://www.docker.com/products/docker-desktop/

Free for personal and educational use.

### FLOOXS

FLOOXS is free for academic/research use but requires registration.

Request access: http://www.flooxs.ece.ufl.edu/index.php/Download

### FLOOXS Documentation

Download the FLOOXS manual for offline reference:

```powershell
scoop install wget
wget -r -l 5 -E -k -np -nH --cut-dirs=1 -e robots=off --reject-regex "(load|api)\.php|action=|Special:|oldid=|title=|User:|Talk:|File:" -P flooxs_docs http://www.flooxs.ece.ufl.edu/index.php/Main_Page
```

## Setup

1. Install Docker Desktop and verify it's running
2. Download FLOOXS source code after registration
3. Extract to `FLOOXS_2026/` in this folder (gitignored, not committed)
4. Build the Docker image (one time):

```bash
docker compose build
```

## Usage

Run a simulation:

```bash
docker compose run --rm flooxs input.tcl
```

Input/output files go in `simulations/` which is mounted into the container.

Interactive FLOOXS prompt:

```bash
docker compose run --rm flooxs
```

Interactive bash shell (for debugging):

```bash
docker compose run --rm --entrypoint /bin/bash flooxs
```

## Workflow

1. Generate inputs - Python script creates FLOOXS input files in `simulations/`
2. Run simulations - `docker compose run` executes FLOOXS on each input
3. Post-process - Python extracts profiles and builds lookup table

## Structure

```
lookupTableGeneration/
├── Dockerfile
├── docker-compose.yml
├── FLOOXS_2026/             # FLOOXS source (gitignored)
├── templates/
│   └── ion_implant.tcl      # FLOOXS simulation template
├── simulations/             # Input/output directory (mounted in container)
├── simulation.template      # TSUPREM-4 reference template
├── simulationControl.py     # TSUPREM-4 batch runner
├── postProcessTables.m      # MATLAB post-processing
└── lookupTable.mat          # Generated lookup table data
```

## Parameter Space

| Parameter | Values |
|-----------|--------|
| Dopants | Boron, Phosphorus, Arsenic |
| Doses | 2e14, 2e15, 2e16 cm^-2 |
| Energies | 20, 50, 80 keV |
| Anneal times | 15, 30, 45, 60, 75, 90, 105, 120, 180, 240, ..., 900 min |
| Anneal temps | 900, 1000, 1100 C |
| Oxidation | With/without passivation oxide |

Total: 3 x 3 x 3 x 21 x 3 x 2 = 3402 simulations
