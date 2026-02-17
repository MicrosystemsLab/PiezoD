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
| FLOOXS_2026 | `templates/ion_implant_fermi.tcl`, `templates/ion_implant_react.tcl` | `Dockerfile`, `docker-compose.yml` |
| FLOOXS_2024 | `templates/ion_implant_2024.tcl` | `Dockerfile.2024`, `docker-compose.2024.yml` |

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
bash scripts/gen_test.sh boron 2e14 80 1000 30 > simulations/test.tcl
docker compose run --rm flooxs test.tcl
```

Arguments: `<dopant> <dose> <energy> <temp> <time>`

Input/output files go in `simulations/` which is mounted into the container.

## Templates

### Fermi model (`ion_implant_fermi.tcl`)

Fermi-level dependent diffusion with equilibrium defects. Single solved variable (dopant concentration) with algebraic charge neutrality. Works for all three dopants.

Physics:
- Fermi-level dependent diffusivity via Noni (n/ni) or Poni (p/ni)
- Combined interstitial + vacancy diffusion paths
- 250A protection oxide (matches TSUPREM-4 pre-implant oxidation)
- Background doping: n-type (P) for boron, p-type (B) for phosphorus/arsenic
- 45 min ramp from 800C, dwell, 45 min ramp down (matches TSUPREM-4)

Does NOT model: transient enhanced diffusion (TED), explicit point defect transport, implant damage effects, oxidation-enhanced diffusion.

Validation against TSUPREM-4 `pd.5str` reference (1000C, 30 min, inert):

| Dopant | Dose/Energy | FLOOXS (um) | TSUPREM-4 (um) | Error |
|--------|-------------|-------------|----------------|-------|
| B | 2e14/80keV | 0.652 | 1.230 | -47% |
| B | 2e15/20keV | 0.639 | 0.880 | -27% |
| P | 2e14/80keV | 0.547 | 0.710 | -23% |
| P | 2e15/50keV | 0.832 | 0.800 | +4% |
| As | 2e14/80keV | 0.237 | 0.520 | -54% |
| As | 2e15/50keV | 0.317 | 0.400 | -21% |

The systematic under-prediction is expected: TSUPREM-4 reference uses `pd.5str` (full pair-diffusion with TED from implant damage), while the fermi model uses equilibrium defects. The gap is largest at low dose where TED dominates, and smallest at high dose where Fermi-level enhancement dominates.

### React model (`ion_implant_react.tcl`)

Full pair-diffusion model with explicit defect transport and BIC cluster chemistry. Boron only. Work in progress.

## Structure

```
lookupTableGeneration/
├── Dockerfile, docker-compose.yml           # FLOOXS_2026
├── Dockerfile.2024, docker-compose.2024.yml # FLOOXS_2024
├── FLOOXS_2024/, FLOOXS_2026/               # Source (gitignored)
├── templates/                               # Parameterized .tcl files (${dopant}, ${dose}, etc.)
│   ├── ion_implant_fermi.tcl                # Fermi-level diffusion (all dopants)
│   ├── ion_implant_react.tcl                # Pair-diffusion + BIC clusters (boron, WIP)
│   └── ion_implant_2024.tcl                 # FLOOXS_2024 (legacy)
├── scripts/                                 # Parameter substitution, run Docker, parse output
│   └── gen_test.sh                          # Generate .tcl from template via sed
├── simulations/                             # Working directory mounted into Docker (gitignored)
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
