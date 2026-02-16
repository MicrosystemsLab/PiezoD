# Ion implantation and diffusion simulation for PiezoD lookup tables
# Uses FLOOXS built-in model procs (PotentialEqns, DefectBulk, DopantBulk)
#
# NOTE: This template requires a FLOOXS installation with working pdb
# auto-loading and built-in TclLib procs. It does NOT work in our Docker
# environment due to missing support for:
#   - "term list" command (used by DefectBulk for Pressure term check)
#   - "term name=X print" command (used by DopantBulk for Charge check)
#   - "DiffLimit" function (used at solve time for I-V recombination)
#   - pdb auto-loading from Params/ directory
# Use ion_implant_fermi.tcl or ion_implant_react.tcl for actual simulations.
#
# Parameters:
#   ${dopant}      - boron, phosphorus, or arsenic
#   ${dose}        - implant dose (cm^-2)
#   ${energy}      - implant energy (keV)
#   ${temp}        - anneal temperature (C)
#   ${time}        - anneal time (minutes)
#
# Physics:
#   - Fermi-level dependent dopant diffusivity (PotentialEqns: charge neutrality)
#   - Detailed-balance defect diffusion (DefectBulk)
#   - Dopant Fermi-enhanced diffusion (DopantBulk, DiffModel=Fermi)
#   - I-V bulk recombination
#   - Surface recombination at oxide/silicon interface (DefectBound)
#   - Background doping (phosphorus, 10 ohm-cm n-type)
#   - Temperature ramp matching TSUPREM-4
#
# Built-in proc call sequence (from TclLib/Models/):
#   PotentialEqns  -> creates Potential, Charge, Noni, Poni
#   DefectInit     -> clears defect equation strings
#   DefectBulk     -> creates EqInt/EqVac, ScaleInt/ScaleVac, IVRecomb,
#                     PressureInt/PressureVac; sets defect transport equations
#   InitDopant     -> creates boronInt/boronVac nosolve solutions
#   DopantBulk     -> dispatches to DopantFermi (DiffModel=1);
#                     creates boronActive, Diffboron; adds charge to Charge term
#   DefectBound    -> surface recombination at Oxide_Silicon interface
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

# Background doping: phosphorus at ~1.4e15 cm^-3 (10 ohm-cm n-type)
sel z=1.4e15 name=Phosphorus

# Ion implantation
implant ${dopant} dose=${dose} energy=${energy} tilt=7

# Initialize interstitials from implant damage (+1 model)
sel z=${dopant} name=Int

# Set up solutions
options !constdatafields storenodes
solution name = Temp !negative add const val = 800

solution name = Int solve !negative add
solution name = Vac solve !negative add
solution name = ${dopant} solve !negative

# Initialize vacancies from implant damage (Frenkel pair model)
sel z=${dopant} name=Vac

# =============================================================================
# PDB PARAMETERS
# =============================================================================
# In a standard FLOOXS installation these load automatically from
# Params/Silicon/{Interstitial,Vacancy,Potential,Boron,...}.
# Set manually here for environments where auto-loading is unavailable.

# --- Potential / Fermi level ---
pdbSetDouble Silicon Potential ni {3.87e16 * sqrt((Temp+273.0)*(Temp+273.0)*(Temp+273.0)) * exp(-7014.0/(Temp+273.0))}
pdbSetDouble Silicon Potential Permittivity 11.7

# --- Interstitial parameters ---
pdbSetDouble Silicon Int Cstar {[Arrhenius 3.6484e27 3.7]}
pdbSetDouble Silicon Int neutral 1.0
pdbSetDouble Silicon Int negative {[Arrhenius 5.68 0.48]}
pdbSetDouble Silicon Int positive {[Arrhenius 5.68 0.42]}
pdbSetDouble Silicon Int dnegative 0.0
pdbSetDouble Silicon Int dpositive 0.0
pdbSetDouble Silicon Int D0 {[Arrhenius 0.138 1.37]}
pdbSetDouble Silicon Int Volume 0.0
pdbSetString Silicon Int RecombDefect Vac
pdbSetSwitch Silicon Int DiffModel 1

# --- Vacancy parameters ---
pdbSetDouble Silicon Vac Cstar {[Arrhenius 4.0515e26 3.97]}
pdbSetDouble Silicon Vac neutral 1.0
pdbSetDouble Silicon Vac negative {[Arrhenius 5.68 0.145]}
pdbSetDouble Silicon Vac positive {[Arrhenius 5.68 0.455]}
pdbSetDouble Silicon Vac dnegative {[Arrhenius 32.47 0.62]}
pdbSetDouble Silicon Vac dpositive 0.0
pdbSetDouble Silicon Vac D0 {[Arrhenius 1.18e-4 0.1]}
pdbSetDouble Silicon Vac Volume 0.0
pdbSetString Silicon Vac RecombDefect Int
pdbSetSwitch Silicon Vac DiffModel 1

# --- Surface recombination (Oxide_Silicon interface) ---
pdbSetDouble Oxide_Silicon Int Ksurf {3.14159*[Arrhenius 0.138 1.37]*2.714417617e-8*1.3e15}
pdbSetDouble Oxide_Silicon Int Krat 0.0
pdbSetDouble Oxide_Silicon Int Scale 1.0
pdbSetDouble Oxide_Silicon Int theta 0.0
pdbSetDouble Oxide_Silicon Int Ktrap 0.0
pdbSetDouble Oxide_Silicon Vac Ksurf {3.14159*[Arrhenius 1.18e-4 0.1]*2.714417617e-8*1.0e5}
pdbSetDouble Oxide_Silicon Vac Krat 0.0
pdbSetDouble Oxide_Silicon Vac Scale 1.0
pdbSetDouble Oxide_Silicon Vac theta 0.0
pdbSetDouble Oxide_Silicon Vac Ktrap 0.0
pdbSetDouble Silicon LatticeDensity 5.0e22

# --- Dopant parameters (per-dopant) ---
switch ${dopant} {
    boron {
        pdbSetSwitch Silicon ${dopant} DiffModel 1
        pdbSetSwitch Silicon ${dopant} Charge 1
        pdbSetSwitch Silicon ${dopant} ActiveModel 0
        pdbSetDouble Silicon ${dopant} Abs.Error 1.0e10
        pdbSetString Silicon ${dopant} Defects "Int Vac"
        # Combined I+V Fermi diffusivity: D = D0 + Dp*Poni
        pdbSetDouble Silicon ${dopant} D0 {[Arrhenius 0.743 3.56] + [Arrhenius 0.186 3.56]}
        pdbSetDouble Silicon ${dopant} Dp {[Arrhenius 0.617 3.56] + [Arrhenius 0.154 3.56]}
    }
    phosphorus {
        pdbSetSwitch Silicon ${dopant} DiffModel 1
        pdbSetSwitch Silicon ${dopant} Charge 2
        pdbSetSwitch Silicon ${dopant} ActiveModel 0
        pdbSetDouble Silicon ${dopant} Abs.Error 1.0e10
        pdbSetString Silicon ${dopant} Defects "Int"
        # Fermi diffusivity: D = D0 + Dn*Noni + Dnn*Noni^2
        pdbSetDouble Silicon ${dopant} D0 {[Arrhenius 5.6 3.71]}
        pdbSetDouble Silicon ${dopant} Dn {[Arrhenius 6.38 4.05]}
        pdbSetDouble Silicon ${dopant} Dnn {[Arrhenius 2.45e-2 3.23]}
    }
    arsenic {
        pdbSetSwitch Silicon ${dopant} DiffModel 1
        pdbSetSwitch Silicon ${dopant} Charge 2
        pdbSetSwitch Silicon ${dopant} ActiveModel 0
        pdbSetDouble Silicon ${dopant} Abs.Error 1.0e10
        pdbSetString Silicon ${dopant} Defects "Int Vac"
        # Combined I+V Fermi diffusivity: D = D0 + Dn*Noni
        pdbSetDouble Silicon ${dopant} D0 {[Arrhenius 0.0666 3.45]}
        pdbSetDouble Silicon ${dopant} Dn {[Arrhenius 12.8 4.05]}
    }
}

# =============================================================================
# BUILT-IN MODEL SETUP
# =============================================================================

# 1. Fermi level: creates Potential, Charge, Noni, Poni solutions
#    Charge neutrality: Potential/Vt = log(n/ni)
PotentialEqns Silicon Potential

# 2. Point defects: creates EqInt/EqVac, ScaleInt/ScaleVac, IVRecomb,
#    PressureInt/PressureVac; sets detailed-balance transport equations
DefectInit Silicon Int
DefectInit Silicon Vac
DefectBulk Silicon Int
DefectBulk Silicon Vac

# 3. Dopant diffusion with Fermi enhancement
#    DopantFermi: D = (D0 + Dp*Poni)/Poni, flux = D*grad(Active*Poni)
#    Adds dopant charge to Charge term (closing self-consistent loop)
InitDopant Silicon ${dopant}
DopantBulk Silicon ${dopant}

# 4. Surface recombination boundary conditions
DefectBound Oxide_Silicon Int
DefectBound Oxide_Silicon Vac

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
sel z=Int
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

# Junction depth where implanted dopant equals background (1.4e15 cm^-3)
sel z=${dopant}-1.4e15
set xj [interpolate silicon val=0.0]
puts "=== RESULTS ==="
puts "Xj: $xj um"

exit
