# Ion implantation and diffusion simulation for PiezoD lookup tables
# FLOOXS template with I-V recombination - parameters substituted by runner script
#
# Parameters:
#   ${dopant}      - boron, phosphorus, or arsenic
#   ${dose}        - implant dose (cm^-2)
#   ${energy}      - implant energy (keV)
#   ${temp}        - anneal temperature (C)
#   ${time}        - anneal time (minutes)
#
# Physics:
#   - Interstitial and vacancy point defect dynamics
#   - I-V bulk recombination (diffusion-limited)
#   - Surface recombination at oxide/silicon interface
#   - Dopant diffusion via both interstitial and vacancy mechanisms
#   - Temperature ramp matching TSUPREM-4
#
# All parameters from FLOOXS_2026/Params/Silicon/

math diffuse dim=1 umf none col scale

# 1D mesh (5um depth, fine near surface for implant resolution)
line x loc=-0.001 spac=0.001 tag=TopOx
line x loc=0.0 spac=0.001 tag=Top
line x loc=0.1 spac=0.001
line x loc=0.5 spac=0.005
line x loc=2.0 spac=0.02
line x loc=5.0 spac=0.05 tag=Bottom

region Oxide xlo=TopOx xhi=Top
region Silicon xlo=Top xhi=Bottom
init

# Ion implantation
implant ${dopant} dose=${dose} energy=${energy} tilt=7

# Initialize interstitials from implant damage (+1 model)
sel z=${dopant} name=Inter

# Set up solutions
options !constdatafields storenodes
solution name=Temp !negative add const val=800

# =============================================================================
# POINT DEFECT PARAMETERS (temperature-dependent via Arrhenius)
# Source: FLOOXS_2026/Params/Silicon/Interstitial, Vacancy, Info
# =============================================================================

set lattice 2.714417617e-8

# Equilibrium concentrations: Cstar = prefactor * exp(-Ea/kT)
# exp(11.2) = 7.2967e4, exp(9.0) = 8.1031e3
pdbSetDouble Silicon Inter Cstar {[Arrhenius 3.6484e27 3.7]}
pdbSetDouble Silicon Vac Cstar {[Arrhenius 4.0515e26 3.97]}
set cis [pdbDelayDouble Silicon Inter Cstar]
set cvs [pdbDelayDouble Silicon Vac Cstar]

# Defect diffusivities
pdbSetDouble Silicon Inter Diff {[Arrhenius 0.138 1.37]}
pdbSetDouble Silicon Vac Diff {[Arrhenius 1.18e-4 0.1]}
set DiffI [pdbDelayDouble Silicon Inter Diff]
set DiffV [pdbDelayDouble Silicon Vac Diff]

# I-V bulk recombination rate: kIV = 4*pi*(D_I+D_V)*lattice
# Cannot express sum of Arrhenius as single Arrhenius, so inline it
set kIV "(4.0*3.14159*($DiffI+$DiffV)*$lattice)"

# =============================================================================
# DOPANT DIFFUSION PARAMETERS
# Source: FLOOXS_2026/Params/Silicon/{Dopant}/Interstitial, Vacancy
# =============================================================================

switch ${dopant} {
    boron {
        pdbSetDouble Silicon Boron DiffI {[Arrhenius 1.36 3.56]}
        pdbSetDouble Silicon Boron DiffV {[Arrhenius 0.34 3.56]}
    }
    phosphorus {
        pdbSetDouble Silicon Phosphorus DiffI {[Arrhenius 5.6 3.71]}
        pdbSetDouble Silicon Phosphorus DiffV 0.0
    }
    arsenic {
        pdbSetDouble Silicon Arsenic DiffI {[Arrhenius 0.0666 3.45]}
        pdbSetDouble Silicon Arsenic DiffV {[Arrhenius 12.8 4.05]}
    }
}
set DiffDopantI [pdbDelayDouble Silicon ${dopant} DiffI]
set DiffDopantV [pdbDelayDouble Silicon ${dopant} DiffV]

# =============================================================================
# SURFACE RECOMBINATION PARAMETERS
# Source: FLOOXS_2026/Params/Oxide_Silicon/Interstitial, Vacancy
# =============================================================================

# Interstitial surface recombination: Ksurf = pi * D_I * lattice * KinkSite
# KinkSite ~ 1.3e15 (saturated at anneal temperatures)
set KsurfI "(3.14159*$DiffI*$lattice*1.3e15)"

# Vacancy surface recombination: KinkSite = 1.0e5 (constant)
set KsurfV "(3.14159*$DiffV*$lattice*1.0e5)"

# =============================================================================
# SOLUTION DEFINITIONS
# =============================================================================

solution name=Inter solve !negative add
solution name=Vac solve !negative add
solution name=${dopant} solve !negative

# Initialize vacancies from implant damage (Frenkel pair model)
# Each implanted atom creates one I-V pair, so initial V = initial I
sel z=${dopant} name=Vac

# =============================================================================
# DIFFUSION EQUATIONS
# =============================================================================

# Interstitial diffusion with I-V bulk recombination
# dI/dt = D_I * laplacian(I) + k_IV * (I*V - I*V*)
# Bulk recombination sign: positive drives ddt negative (FLOOXS convention)
pdbSetString Silicon Inter Equation "ddt(Inter) - $DiffI*grad(Inter) + $kIV*(Inter*Vac - ($cis)*($cvs))"

# Vacancy diffusion with I-V bulk recombination
pdbSetString Silicon Vac Equation "ddt(Vac) - $DiffV*grad(Vac) + $kIV*(Inter*Vac - ($cis)*($cvs))"

# Surface recombination at oxide/silicon interface
# Interface recombination sign: negative (FLOOXS convention)
pdbSetString Oxide_Silicon Inter Silicon Equation "-$KsurfI*(Inter(Silicon)-($cis))"
pdbSetString Oxide_Silicon Vac Silicon Equation "-$KsurfV*(Vac(Silicon)-($cvs))"

# Dopant diffusion with interstitial AND vacancy enhancement
# dC/dt = div(D_I*(I/I*) * grad(C)) + div(D_V*(V/V*) * grad(C))
pdbSetString Silicon ${dopant} Equation "ddt(${dopant}) - ($DiffDopantI) * (Inter/($cis)) * grad(${dopant}) - ($DiffDopantV) * (Vac/($cvs)) * grad(${dopant})"

# =============================================================================
# OUTPUT: PRE-ANNEAL PROFILE
# =============================================================================

puts "=== PRE-ANNEAL ==="
sel z=${dopant}
foreach v [print.1d] {
    lassign $v d c
    if {$c == "Value"} continue
    if {$d < 0} continue
    puts "$d $c"
}

# =============================================================================
# ANNEAL SEQUENCE (matches TSUPREM-4)
# =============================================================================

# Calculate ramp rate: (temp - 800) / (45 min * 60 s/min) in C/s
set ramprate [expr {(${temp} - 800.0) / 2700.0}]

puts "=== RAMP UP: 45 min, 800C -> ${temp}C ==="
diffuse time=45 temp=800 ramprate=$ramprate

puts "=== DWELL: ${time} min at ${temp}C ==="
diffuse time=${time} temp=${temp}

puts "=== RAMP DOWN: 45 min, ${temp}C -> 800C ==="
diffuse time=45 temp=${temp} ramprate=-$ramprate

# =============================================================================
# OUTPUT: POST-ANNEAL PROFILES
# =============================================================================

puts "=== POST-ANNEAL ==="
sel z=${dopant}
foreach v [print.1d] {
    lassign $v d c
    if {$c == "Value"} continue
    if {$d < 0} continue
    puts "$d $c"
}

puts "=== POST-ANNEAL INTERSTITIALS ==="
sel z=Inter
foreach v [print.1d] {
    lassign $v d c
    if {$c == "Value"} continue
    if {$d < 0} continue
    if {$d > 0.5} break
    puts "$d $c"
}

puts "=== POST-ANNEAL VACANCIES ==="
sel z=Vac
foreach v [print.1d] {
    lassign $v d c
    if {$c == "Value"} continue
    if {$d < 0} continue
    if {$d > 0.5} break
    puts "$d $c"
}

# Junction depth at 1e15 cm^-3 (typical substrate background)
sel z=${dopant}-1e15
set xj [interpolate silicon val=0.0]
puts "=== RESULTS ==="
puts "Xj: $xj um"

exit
