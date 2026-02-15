# Lookup Table Generation

Generate dopant profile lookup tables for PiezoD using FLOOXS TCAD process simulation.

## Prerequisites

### Docker

Install Docker Desktop: https://www.docker.com/products/docker-desktop/

### FLOOXS

FLOOXS is free for academic/research use but requires registration.

Request access: http://www.flooxs.ece.ufl.edu/index.php/Download

Two versions are supported with separate templates:

| Version | Template | Docker Files |
|---------|----------|--------------|
| FLOOXS_2026 | `templates/ion_implant.tcl` | `Dockerfile`, `docker-compose.yml` |
| FLOOXS_2024 | `templates/ion_implant_2024.tcl` | `Dockerfile.2024`, `docker-compose.2024.yml` |

### FLOOXS Documentation

Local copy in `www.flooxs.ece.ufl.edu/` (gitignored). To download:

```powershell
./grab_flooxs_docs.ps1
```

Or manually:

```powershell
scoop install wget
wget -r -l 5 -E -k -np -e robots=off --reject-regex "(load|api)\.php|action=|Special:|oldid=|title=|User:|Talk:|File:" http://www.flooxs.ece.ufl.edu/index.php/Main_Page
```

## Setup

1. Install Docker Desktop and verify it's running
2. Download FLOOXS source code after registration
3. Extract to `FLOOXS_2024/` or `FLOOXS_2026/` (gitignored)
4. Build the Docker image:

```bash
docker compose build
docker compose -f docker-compose.2024.yml build
```

## Usage

Run a simulation:

```bash
docker compose run --rm flooxs input.tcl
docker compose -f docker-compose.2024.yml run --rm flooxs input.tcl
```

Input/output files go in `simulations/` which is mounted into the container.

## Structure

```
lookupTableGeneration/
├── Dockerfile, docker-compose.yml           # FLOOXS_2026
├── Dockerfile.2024, docker-compose.2024.yml # FLOOXS_2024
├── FLOOXS_2024/, FLOOXS_2026/               # Source (gitignored)
├── templates/                               # Parameterized .tcl files (${dopant}, ${dose}, etc.)
│   ├── ion_implant.tcl                      # FLOOXS_2026
│   └── ion_implant_2024.tcl                 # FLOOXS_2024
├── scripts/                                 # Python: parameter substitution, run Docker, parse output
├── simulations/                             # Working directory mounted into Docker container
│   └── reference_iv_recomb.tcl              # Reference: I-V recombination model
├── simulation.template                      # TSUPREM-4 reference (legacy)
├── simulationControl.py                     # TSUPREM-4 runner (legacy)
└── lookupTable.mat                          # Generated data
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
