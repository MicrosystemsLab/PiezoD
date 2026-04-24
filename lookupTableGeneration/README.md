# Lookup Table Generation

PiezoD uses a TSUPREM-4-generated lookup table for dopant profiles (`lookupTable.mat` -> `ionImplantLookupTable_tsuprem.h5`). This directory contains an attempt to replace TSUPREM-4 with FLOOXS (open-source TCAD). The attempt failed due to fundamental parameter and solver issues. See [docs/status.md](docs/status.md) for the full evaluation.

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

Three C++ bugs and three Tcl bugs were found and patched. See [docs/flooxs_bugs.md](docs/flooxs_bugs.md) for details.

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
bash scripts/gen_test.sh boron 2e14 80 1000 30 templates/ion_implant_5pd_v5.tcl > simulations/test.tcl
MSYS_NO_PATHCONV=1 docker compose run --rm flooxs /work/test.tcl
```

Arguments: `<dopant> <dose> <energy> <temp> <time> [template]`

The template defaults to `ion_implant_fermi.tcl` if not specified.

Input/output files go in `simulations/` which is mounted into the container.

## Templates

| Template | Model | Status |
|----------|-------|--------|
| `ion_implant_eff_v4.tcl` | Effective D, FLOOXS pair D, ScaleInter cap, 10um mesh | Best stability (6/6 B) |
| `ion_implant_eff_v3.tcl` | Effective D, Table A-3 D, IV relaxation | Best low-dose B accuracy |
| `ion_implant_eff_v2.tcl` | Effective D, FLOOXS pair D, IV relaxation | Better high-dose stability |
| `ion_implant_eff.tcl` | Effective D, Table A-3 D, no IV sink | First eff template |
| `ion_implant_5pd_v5.tcl` | 5-stream + {311}, FLOOXS params, V_eq init | Best 1000C accuracy |
| `ion_implant_5pd_v5_noVac.tcl` | 5-stream + {311}, no Vac PDE for B | Worse than v5 |
| `ion_implant_5pd_v4.tcl` | 5-stream + {311}, FLOOXS params | Superseded by v5 |
| `ion_implant_5pd.tcl` | 5-stream + {311}, TSUPREM-4 Appendix A params | Failed (pair fraction) |
| `ion_implant_ted.tcl` | 5-stream pair diffusion, no {311} | Under-predicts Xj |
| `ion_implant_ted_qss.tcl` | QSS/Fermi hybrid | Superseded by 5pd |
| `ion_implant_builtin_fermi.tcl` | FLOOXS built-in model procs | Fermi only (no TED) |
| `ion_implant_fermi.tcl` | Fermi-level dependent diffusion, no TED | Baseline |
| `ion_implant_2024.tcl` | FLOOXS_2024 | Legacy |

See [docs/implementation.md](docs/implementation.md) for physics details, [comparison.md](comparison.md) for validation against TSUPREM-4, and [docs/status.md](docs/status.md) for the evaluation conclusion.

## Structure

```
lookupTableGeneration/
├── Dockerfile, docker-compose.yml           # FLOOXS_2026
├── Dockerfile.2024, docker-compose.2024.yml # FLOOXS_2024
├── FLOOXS_2024/, FLOOXS_2026/               # Source (gitignored)
├── templates/                               # Parameterized .tcl files (${dopant}, ${dose}, etc.)
│   ├── ion_implant_eff_v4.tcl               # Effective D, best stability
│   ├── ion_implant_eff_v3.tcl               # Effective D, best low-dose B
│   ├── ion_implant_eff_v2.tcl               # Effective D, FLOOXS pair D
│   ├── ion_implant_eff.tcl                  # Effective D, Table A-3 D
│   ├── ion_implant_5pd_v5.tcl               # 5-stream + {311}, best 1000C
│   ├── ion_implant_5pd_v5_noVac.tcl         # 5-stream, no Vac (worse than v5)
│   ├── ion_implant_5pd_v4.tcl               # 5-stream + {311} (superseded by v5)
│   ├── ion_implant_5pd.tcl                  # 5-stream, TSUPREM-4 params (failed)
│   ├── ion_implant_ted.tcl                  # 5-stream pair diffusion (no {311})
│   ├── ion_implant_ted_qss.tcl              # QSS/Fermi hybrid (superseded)
│   ├── ion_implant_builtin_fermi.tcl        # FLOOXS built-in procs (Fermi only)
│   ├── ion_implant_fermi.tcl                # Fermi-level diffusion (baseline)
│   ├── ion_implant_2024.tcl                 # FLOOXS_2024 (legacy)
│   └── Dopant_patched.tcl                   # Patched FLOOXS Dopant.tcl
├── scripts/                                 # Parameter substitution, run Docker, parse output
│   └── gen_test.sh                          # Generate .tcl from template via sed
├── simulations/                             # Working directory mounted into Docker (gitignored)
├── docs/                                    # Physics details, evaluation status, bug reports
│   ├── status.md                            # Evaluation conclusion
│   ├── flooxs_bugs.md                       # 6 bugs found and patched
│   ├── implementation.md                    # Physics and equation details
│   ├── eff_template_plan.md                 # Effective D derivation
│   └── diffusion_theory.md                  # Background theory
├── comparison.md                            # Validation data vs TSUPREM-4
├── simulation.template                      # TSUPREM-4 input template (legacy)
├── simulationControl.py                     # TSUPREM-4 batch runner (legacy)
├── lookupTable.mat                          # Generated lookup table data
└── postProcessTables.m                      # MATLAB: compute metrics from TSUPREM-4 profiles
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

## Lookup Table Format (`lookupTable.mat`)

Generated by `postProcessTables.m` from TSUPREM-4 simulation output files.

### Variables

| Variable | Shape | Units | Description |
|----------|-------|-------|-------------|
| `z` | (501,) | um | Depth grid, 0 to 5 um at 10 nm spacing |
| `n` | (501, 3, 3, 3, 3, 8, 2) | cm^-3 | Concentration profiles |
| `Xj` | (3, 3, 3, 3, 8, 2) | m | Junction depth |
| `Rs` | (3, 3, 3, 3, 8, 2) | ohm/sq | Sheet resistance |
| `Beta1` | (3, 3, 3, 3, 8, 2) | - | Piezoresistive coefficient (uniform term) |
| `Beta2` | (3, 3, 3, 3, 8, 2) | um | Piezoresistive coefficient (depth term) |
| `Nz` | (3, 3, 3, 3, 8, 2) | m^-2 | Mobility-weighted effective carrier density |
| `Nz_total` | (3, 3, 3, 3, 8, 2) | m^-2 | Total effective carrier density |

### Array indices

Dimensions are (dopant, dose, energy, temp, time, oxidation):

| Dimension | Index 0 | Index 1 | Index 2 | Index 3-7 |
|-----------|---------|---------|---------|-----------|
| dopant | B | P | As | |
| dose | 2e14 | 2e15 | 2e16 | |
| energy | 20 keV | 50 keV | 80 keV | |
| temp | 900 C | 1000 C | 1100 C | |
| time | 15 min | 30 min | 45 min | 60, 75, 90, 105, 120 min |
| oxidation | inert | oxide | | |

### Metric computation

All metrics (Rs, Beta1, Beta2, Nz) are computed by integrating the concentration
profile from the surface to the junction depth only. Below the junction the
substrate is the opposite carrier type and does not contribute to the
piezoresistor. The beta decomposition is:

    beta*(t) = Beta1 - 2 * Beta2 / t

where `t` is cantilever thickness in um.
