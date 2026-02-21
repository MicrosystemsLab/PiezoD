# Lookup Table Generation

Generate dopant profile lookup tables for PiezoD using FLOOXS TCAD process simulation.

## Prerequisites

### Docker

Install Docker Desktop: https://www.docker.com/products/docker-desktop/

### FLOOXS

FLOOXS is free for academic/research use but requires registration.

Request access: http://www.flooxs.ece.ufl.edu/index.php/Download

Two versions are supported with separate templates:

| Version | Docker Files |
|---------|--------------|
| FLOOXS_2026 | `Dockerfile`, `docker-compose.yml` |
| FLOOXS_2024 | `Dockerfile.2024`, `docker-compose.2024.yml` |

No patches to the FLOOXS source code are required. All templates use stock FLOOXS features.

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

Generate a test script from a template and run it:

```bash
cd lookupTableGeneration
bash scripts/gen_test.sh boron 2e14 80 1000 30 templates/ion_implant_5pd.tcl > simulations/test.tcl
MSYS_NO_PATHCONV=1 docker compose run --rm flooxs /work/test.tcl
```

Arguments: `<dopant> <dose> <energy> <temp> <time> [template]`

The template defaults to `ion_implant_fermi.tcl` if not specified.

Input/output files go in `simulations/` which is mounted into the container.

## Templates

| Template | Model | Status |
|----------|-------|--------|
| `ion_implant_5pd.tcl` | 5-stream pair diffusion + {311} clustering (all dopants) | Primary, used for lookup tables |
| `ion_implant_ted.tcl` | 5-stream pair diffusion, no {311} | Reference (under-predicts Xj) |
| `ion_implant_ted_qss.tcl` | QSS/Fermi hybrid | Superseded by 5pd |
| `ion_implant_fermi.tcl` | Fermi-level dependent diffusion, no TED | Baseline (no TED physics) |
| `ion_implant_2024.tcl` | FLOOXS_2024 | Legacy |

See [docs/implementation.md](docs/implementation.md) for physics details. See [comparison.md](comparison.md) for validation against TSUPREM-4.

## Structure

```
lookupTableGeneration/
├── Dockerfile, docker-compose.yml           # FLOOXS_2026
├── Dockerfile.2024, docker-compose.2024.yml # FLOOXS_2024
├── FLOOXS_2024/, FLOOXS_2026/               # Source (gitignored)
├── templates/                               # Parameterized .tcl files (${dopant}, ${dose}, etc.)
│   ├── ion_implant_5pd.tcl                  # 5-stream + {311} clustering (primary)
│   ├── ion_implant_ted.tcl                  # 5-stream pair diffusion (no {311})
│   ├── ion_implant_ted_qss.tcl              # QSS/Fermi hybrid (superseded)
│   ├── ion_implant_fermi.tcl                # Fermi-level diffusion (no TED)
│   └── ion_implant_2024.tcl                 # FLOOXS_2024 (legacy)
├── scripts/                                 # Parameter substitution, run Docker, parse output
│   └── gen_test.sh                          # Generate .tcl from template via sed
├── simulations/                             # Working directory mounted into Docker (gitignored)
├── docs/                                    # Physics details, theory, implementation notes
├── simulation.template                      # TSUPREM-4 input template (legacy)
├── simulationControl.py                     # TSUPREM-4 batch runner (legacy)
└── lookupTable.mat                          # Generated lookup table data
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
