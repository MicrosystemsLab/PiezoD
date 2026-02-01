# Open-Source TCAD Process Simulators for PiezoD

## Problem Statement

PiezoD is a MATLAB-based tool for optimizing piezoresistive cantilever designs. Its `lookupTableGeneration` module requires a 1D semiconductor process simulator to compute dopant concentration profiles for ion-implanted piezoresistors. The module sweeps 432 parameter combinations:

- **Dopants**: Boron, Phosphorus, Arsenic
- **Doses**: 2×10¹⁴, 2×10¹⁵, 2×10¹⁶ cm⁻²
- **Energies**: 20, 50, 80 keV
- **Anneal times**: 15, 30, 45, 60, 75, 90, 105, 120 minutes
- **Anneal temperatures**: 900, 1000, 1100°C

The original code calls Synopsys TSUPREM-4, which requires an expensive commercial license. This document describes two open-source alternatives that can be deployed via Docker for cross-platform use.

---

## Option 1: SUPREM-IV.GS (Minimal Porting Effort)

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

### Dockerfile

```dockerfile
FROM ubuntu:22.04

RUN apt-get update && apt-get install -y \
    git \
    gcc \
    gfortran \
    make \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt
RUN git clone https://github.com/rafael1193/suprem4gs.git

WORKDIR /opt/suprem4gs
RUN make depend && make install

# Add to PATH
ENV PATH="/opt/suprem4gs:${PATH}"

WORKDIR /work
ENTRYPOINT ["suprem4gs"]
```

### Build and Run

```bash
# Build image
docker build -t suprem4gs .

# Run simulation (mount input/output directory)
docker run --rm -v $(pwd)/simulations:/work suprem4gs input_file.inp

# Interactive shell
docker run --rm -it -v $(pwd)/simulations:/work --entrypoint /bin/bash suprem4gs
```

### PiezoD Integration

Modify `simulationControl.py`:

```python
# Original
subprocess.call('tsuprem4 ' + simulationFileName, shell=True)

# Docker version
subprocess.call('docker run --rm -v ' + simDir + ':/work suprem4gs ' + simulationFileName, shell=True)
```

### Limitations

- **No TED (Transient Enhanced Diffusion)**: Dopant profiles will be shallower than reality at high doses. Acceptable for relative comparisons but not absolute accuracy.
- **1993 physics**: Missing modern point-defect coupling models.
- **Syntax compatibility**: May require minor input file adjustments from TSUPREM-4 format.

---

## Option 2: FLOOXS/FLOOPS (Better Physics, More Porting Effort)

### Overview

FLOOXS (Florida Object-Oriented Process/Device Simulator) is actively maintained by University of Florida. FLOOPS is the process simulation subset.

| Attribute | Value |
|-----------|-------|
| Source | http://www.flooxs.org/ |
| License | Free for academic/research (registration required) |
| Last release | 2024.0.0 |
| Dimensions | 1D, 2D, (3D experimental) |
| Dopants | B, P, As, BF₂, Sb |
| TED support | Yes |
| Scripting | Tcl + Alagator |

### Obtaining FLOOXS

FLOOXS requires registration:

1. Complete license form at http://www.flooxs.ece.ufl.edu/index.php/Download
2. Email Prof. Mark Law (law@ece.ufl.edu)
3. Receive download credentials
4. Download appropriate .deb package

### Dockerfile

```dockerfile
FROM ubuntu:24.04

RUN apt-get update && apt-get install -y \
    tcl \
    tk \
    libblas3 \
    liblapack3 \
    libsuitesparse-dev \
    libplplot-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy the .deb (must be obtained via registration)
COPY ubuntu_24.04-flooxs_2024.0.0_amd64.deb /tmp/
RUN dpkg -i /tmp/ubuntu_24.04-flooxs_2024.0.0_amd64.deb || apt-get install -f -y

ENV FLXSHOME=/usr/local/flooxs
WORKDIR /work

ENTRYPOINT ["flooxs"]
```

### Build and Run

```bash
# Build image (after placing .deb in build context)
docker build -t flooxs .

# Run simulation
docker run --rm -v $(pwd)/simulations:/work flooxs input_file.tcl

# Interactive shell
docker run --rm -it -v $(pwd)/simulations:/work --entrypoint /bin/bash flooxs
```

### FLOOPS Syntax Examples

FLOOPS uses Tcl scripting, not TSUPREM-4 syntax. Template conversion is required.

**1D mesh setup:**
```tcl
line x loc=0.0  spac=0.001 tag=Top
line x loc=5.0  spac=0.01  tag=Bottom

region Silicon xlo=Top xhi=Bottom
init
```

**Ion implantation:**
```tcl
implant boron dose=1.0e15 energy=20
implant phosphorus dose=2.0e15 energy=50
implant arsenic dose=2.0e16 energy=80
```

**Diffusion/anneal:**
```tcl
# Time in minutes, temp in Celsius
diffuse time=30 temp=1000
diffuse time=60 temp=1100 wetO2  # with oxidation
```

**Extract profile:**
```tcl
sel z=Boron
print.1d
```

### PiezoD Template Conversion

TSUPREM-4 input:
```
INITIALIZE <100> SILICON
IMPLANT BORON DOSE=1E15 ENERGY=20
DIFFUSION TIME=30 TEMP=1000
PRINT LAYERS
```

Equivalent FLOOPS:
```tcl
line x loc=0.0 spac=0.001 tag=Top
line x loc=5.0 spac=0.01  tag=Bottom
region Silicon xlo=Top xhi=Bottom
init

implant boron dose=1.0e15 energy=20

solution add name=Boron pde solve !negative
pdbSet Si Boron Equation "ddt(Boron) - D*grad(Boron)"
diffuse time=30 temp=1000

sel z=Boron
print.1d
```

### Advantages Over SUPREM-IV.GS

- **TED modeling**: Critical for accurate B/P/As profiles at doses >10¹⁴ cm⁻²
- **Point defect coupling**: Interstitial/vacancy models
- **Active maintenance**: Bug fixes, modern Linux support
- **Better documentation**: Comprehensive wiki with tutorials

### Limitations

- **Registration required**: Cannot freely redistribute .deb
- **Syntax incompatibility**: Requires rewriting all PiezoD templates
- **Learning curve**: Alagator PDE language for custom models

---

## Comparison Summary

| Feature | SUPREM-IV.GS | FLOOXS/FLOOPS |
|---------|--------------|---------------|
| Porting effort | Low | Moderate |
| Physics accuracy | 1993 vintage | Modern (2024) |
| TED support | No | Yes |
| Installation | `git clone` + `make` | Registration + .deb |
| Docker complexity | Simple | Simple (once you have .deb) |
| Template rewrite | Minor tweaks | Full rewrite to Tcl |
| Maintenance | Dormant | Active |

## Recommendation

- **Quick prototype / educational use**: SUPREM-IV.GS
- **Accurate profiles / research publication**: FLOOXS/FLOOPS

For PiezoD's use case (relative optimization across parameter space), SUPREM-IV.GS is adequate. If absolute profile accuracy matters (e.g., validating against SIMS data), invest the effort in FLOOXS.

---

## References

- SUPREM-IV.GS GitHub: https://github.com/rafael1193/suprem4gs
- FLOOXS main page: http://www.flooxs.ece.ufl.edu/index.php/Main_Page
- FLOOXS process tutorial: http://www.flooxs.ece.ufl.edu/index.php/Process_Tutorial
- FLOOXS diffuse command: http://www.flooxs.ece.ufl.edu/index.php/Diffuse_Command
- PiezoD repository: https://github.com/MicrosystemsLab/PiezoD
- Stanford TSUPREM-4 notes: https://web.stanford.edu/class/ee410/TSUPREM4_Notes.pdf
