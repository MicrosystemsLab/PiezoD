# Open-Source TCAD Process Simulator Comparison

The lookup table generation requires a 1D semiconductor process simulator to compute dopant concentration profiles for ion-implanted piezoresistors.

Original lookup table (generated with TSUPREM-4):
- Dopants: Boron, Phosphorus, Arsenic
- Doses: 2e14, 2e15, 2e16 cm^-2
- Energies: 20, 50, 80 keV
- Anneal times: 15, 30, 45, 60, 75, 90, 105, 120 minutes
- Anneal temperatures: 900, 1000, 1100 C

This is insufficient for new design work and we want to extend it. The original design optimization paper showed data to 900 minute anneal, but current lookup ta bles only go to 120 minutes.

Extended lookup table (planned):
- Dopants: Boron, Phosphorus, Arsenic
- Doses: 1e13 to 8e16 cm^-2 (1, 2, 4, 8 per decade = 16 values)
- Energies: 20, 50, 80 keV
- Anneal times: 5, 10, 15, 30, 45, 60, 75, 90, 105, 120, 180, 240, ..., 900 minutes
- Anneal temperatures: 900 to 1100 C in 50 C steps

TSUPREM-4 requires an expensive commercial license. This document describes two open-source alternatives that can be deployed via Docker. We choose to use Docker for a portable, reproducible build environment.

---

## Option 1: SUPREM-IV.GS

### Overview

SUPREM-IV.GS is Stanford's 1993 open-source process simulator, the academic ancestor of commercial TSUPREM-4. Cogenda ported it to modern Linux; rafael1193 maintains a build-ready fork.

| Attribute | Value |
|-----------|-------|
| Source | https://github.com/rafael1193/suprem4gs |
| License | Stanford/public domain |
| Last active development | 1993 (original), 2016 (port) |
| Dimensions | 1D |
| Dopants | B, P, As (Pearson-IV implant models) |
| TED support | No |

### Limitations

- No TED (Transient Enhanced Diffusion)
    - Dopant profiles will be shallower than reality at high doses.
    - Acceptable for relative comparisons but not absolute accuracy.
- 1993 physics
    - Missing modern point-defect coupling models.
- Syntax compatibility
    - May require minor input file adjustments from TSUPREM-4 format.

---

## Option 2: FLOOXS/FLOOPS

### Overview

FLOOXS (Florida Object-Oriented Process/Device Simulator) is actively maintained by University of Florida. FLOOPS is the process simulation subset.

| Attribute | Value |
|-----------|-------|
| Source | http://www.flooxs.org/ |
| License | Free for academic/research (registration required) |
| Last release | 2026 |
| Dimensions | 1D, 2D, (3D experimental) |
| Dopants | B, P, As, BF2, Sb |
| TED support | Yes |
| Scripting | Tcl + Alagator |

### Obtaining FLOOXS

FLOOXS requires registration:

1. Complete license form at http://www.flooxs.ece.ufl.edu/index.php/Download
2. Email Prof. Mark Law
3. Receive archive password
4. Extract and use source code for Docker build

**IMPORTANT:** We do not have authorization to distribute the source code and so it is explicitly ignored and not mixed with any of the piezod code. New developers need to acquire their own access to the code.

### FLOOPS Syntax

FLOOPS uses Tcl scripting, not TSUPREM-4 syntax. Template conversion is required.

However, the changes are fairly minor and can be updated quickly and easily. The docs for FLOOXS are excellent.

### Advantages Over SUPREM-IV.GS

- TED modeling
    - Critical for accurate B/P/As profiles at doses >1e14 cm^-2
- Point defect coupling
    - Interstitial/vacancy models
- Active maintenance
    - Bug fixes, modern Linux support, has a future
- Better documentation
    - Comprehensive wiki with tutorials

The only downside is the registration flow, but that is an extremely minor inconvenience and most users only depend on the distilled lookup table.

---

## Summary

| Feature | SUPREM-IV.GS | FLOOXS/FLOOPS |
|---------|--------------|---------------|
| Porting effort | Low | Low |
| Physics accuracy | 1993 vintage | Modern (2026) |
| Docker complexity | Simple | Simple |
| Template rewrite | Minor tweaks | Full rewrite to Tcl |
| Maintenance | Dormant | Active |

The choice is clear - we will use FLOOXS to update the lookup tables.

---

## References

- SUPREM-IV.GS GitHub: https://github.com/rafael1193/suprem4gs
- FLOOXS main page: http://www.flooxs.ece.ufl.edu/index.php/Main_Page
- FLOOXS process tutorial: http://www.flooxs.ece.ufl.edu/index.php/Process_Tutorial
- FLOOXS diffuse command: http://www.flooxs.ece.ufl.edu/index.php/Diffuse_Command
- Stanford TSUPREM-4 notes: https://web.stanford.edu/class/ee410/TSUPREM4_Notes.pdf
