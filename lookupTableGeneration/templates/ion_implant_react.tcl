# Ion implantation and diffusion simulation for PiezoD lookup tables
# FLOOXS template with dopant-defect pair reactions (DiffModel=React)
#
# Parameters:
#   ${dose}        - implant dose (cm^-2)
#   ${energy}      - implant energy (keV)
#   ${temp}        - anneal temperature (C)
#   ${time}        - anneal time (minutes)
#
# Physics (extends ion_implant_fermi.tcl with explicit pair kinetics):
#   - Boron-interstitial (BI) and boron-vacancy (BV) pair formation/dissolution
#   - Pair diffusion with Fermi-level dependent mobility
#   - Solid solubility limiting (boronActive)
#   - Substitutional boron tracking (boronSub = boronActive - boronInt - boronVac)
#   - Reaction terms couple back to interstitial/vacancy equations
#   - All other physics same as ion_implant_fermi.tcl
#
# Boron-only: BIC pair kinetics are specific to boron. For phosphorus/arsenic,
# use ion_implant_fermi.tcl.
#
# All parameters from FLOOXS_2026/Params/Silicon/

math diffuse dim=1 umf none col scale

# 1D mesh (5um depth, fine near surface for implant resolution)
# 250A (25nm) protection oxide matches TSUPREM-4 pre-implant oxidation
line x loc=-0.025 spac=0.005 tag=TopOx
line x loc=0.0 spac=0.0002 tag=Top
line x loc=0.05 spac=0.0002
line x loc=0.15 spac=0.001
line x loc=0.5 spac=0.005
line x loc=2.0 spac=0.02
line x loc=5.0 spac=0.05 tag=Bottom

region Oxide xlo=TopOx xhi=Top
region Silicon xlo=Top xhi=Bottom
init

# Background doping: phosphorus at ~1.4e15 cm^-3 (10 ohm-cm n-type)
sel z=1.4e15 name=Phosphorus

# Ion implantation
implant boron dose=${dose} energy=${energy} tilt=7

# Initialize interstitials from implant damage (+1 model)
sel z=boron name=Inter

# Set up solutions
options !constdatafields storenodes
solution name=Temp !negative add const val=800

# Solver settings
pdbSetDouble Math iterLimit 500
pdbSetDouble Math updateLimit 1e-8
pdbSetDouble Math rhsLimit 1e-20
pdbSetDouble Math rhsMin 1e-2

# =============================================================================
# INTRINSIC CARRIER CONCENTRATION
# =============================================================================

# ni from Law et al. - accurate at process temperatures (900-1100C)
set ni "(3.87e16 * sqrt((Temp+273.0)*(Temp+273.0)*(Temp+273.0)) * exp(-7014.0/(Temp+273.0)))"

# =============================================================================
# POINT DEFECT PARAMETERS
# Source: FLOOXS_2026/Params/Silicon/Interstitial, Vacancy, Info
# =============================================================================

set lattice 2.714417617e-8

# Equilibrium concentrations
pdbSetDouble Silicon Inter Cstar {[Arrhenius 3.6484e27 3.7]}
pdbSetDouble Silicon Vac Cstar {[Arrhenius 4.0515e26 3.97]}
set cis [pdbDelayDouble Silicon Inter Cstar]
set cvs [pdbDelayDouble Silicon Vac Cstar]

# Interstitial charge states
pdbSetDouble Silicon Inter neutral 1.0
pdbSetDouble Silicon Inter negative {[Arrhenius 5.68 0.48]}
pdbSetDouble Silicon Inter positive {[Arrhenius 5.68 0.42]}
set Ineu [pdbDelayDouble Silicon Inter neutral]
set Ineg [pdbDelayDouble Silicon Inter negative]
set Ipos [pdbDelayDouble Silicon Inter positive]
set IdenI "($Ineu + $Ineg + $Ipos)"

# Vacancy charge states
pdbSetDouble Silicon Vac neutral 1.0
pdbSetDouble Silicon Vac negative {[Arrhenius 5.68 0.145]}
pdbSetDouble Silicon Vac positive {[Arrhenius 5.68 0.455]}
pdbSetDouble Silicon Vac dnegative {[Arrhenius 32.47 0.62]}
set Vneu [pdbDelayDouble Silicon Vac neutral]
set Vneg [pdbDelayDouble Silicon Vac negative]
set Vpos [pdbDelayDouble Silicon Vac positive]
set Vdng [pdbDelayDouble Silicon Vac dnegative]
set IdenV "($Vneu + $Vneg + $Vdng + $Vpos)"

# Defect diffusivities
pdbSetDouble Silicon Inter Diff {[Arrhenius 0.138 1.37]}
pdbSetDouble Silicon Vac Diff {[Arrhenius 1.18e-4 0.1]}
set DiffI [pdbDelayDouble Silicon Inter Diff]
set DiffV [pdbDelayDouble Silicon Vac Diff]

# I-V bulk recombination rate: kIV = 4*pi*(D_I+D_V)*lattice
set kIV "(4.0*3.14159*($DiffI+$DiffV)*$lattice)"

# =============================================================================
# BORON-DEFECT PAIR PARAMETERS
# Source: FLOOXS_2026/Params/Silicon/Boron/Interstitial, Vacancy
# =============================================================================

# Boron solubility (from Params/Silicon/Boron/Info)
pdbSetDouble Silicon boron Solubility {[Arrhenius 7.68e22 0.7086]}
set Css [pdbDelayDouble Silicon boron Solubility]

# BI pair: Binding, diffusivity components (D0 + Dp*Poni), reaction rate
pdbSetDouble Silicon boron BI_Binding {[Arrhenius 8.0e-23 -1.0]}
pdbSetDouble Silicon boron BI_D0 {[Arrhenius 0.743 3.56]}
pdbSetDouble Silicon boron BI_Dp {[Arrhenius 0.617 3.56]}
set Bind_BI [pdbDelayDouble Silicon boron BI_Binding]
set D0_BI [pdbDelayDouble Silicon boron BI_D0]
set Dp_BI [pdbDelayDouble Silicon boron BI_Dp]

# BV pair: Binding, diffusivity components (D0 + Dp*Poni), reaction rate
pdbSetDouble Silicon boron BV_Binding {[Arrhenius 8.0e-23 -0.5]}
pdbSetDouble Silicon boron BV_D0 {[Arrhenius 0.186 3.56]}
pdbSetDouble Silicon boron BV_Dp {[Arrhenius 0.154 3.56]}
set Bind_BV [pdbDelayDouble Silicon boron BV_Binding]
set D0_BV [pdbDelayDouble Silicon boron BV_D0]
set Dp_BV [pdbDelayDouble Silicon boron BV_Dp]

# BIC cluster: effective two-body B+I <-> BIClust with 2.5 eV binding energy
# Lumps higher-order clusters (B3I, B4I) into single immobile species
pdbSetDouble Silicon boron Clust_Bind {[Arrhenius 2.0e-23 -2.5]}
set Bind_Cl [pdbDelayDouble Silicon boron Clust_Bind]

# Diffusion-limited reaction rates: Krate = 4*pi*D_defect*lattice
# DiffLimit with barrier=0 reduces to just the prefactor
set KrateBI "(4.0*3.14159*$DiffI*$lattice)"
set KrateBV "(4.0*3.14159*$DiffV*$lattice)"

# =============================================================================
# SURFACE RECOMBINATION PARAMETERS
# Source: FLOOXS_2026/Params/Oxide_Silicon/Interstitial, Vacancy
# =============================================================================

set KsurfI "(3.14159*$DiffI*$lattice*1.3e15)"
set KsurfV "(3.14159*$DiffV*$lattice*1.0e5)"

# =============================================================================
# SOLUTION DEFINITIONS
# =============================================================================

solution name=Inter solve !negative add
solution name=Vac solve !negative add
solution name=boron solve !negative
solution name=boronInt solve !negative add
solution name=boronVac solve !negative add
solution name=BIClust solve !negative add

# Fermi level coupling
solution name = Charge add const val = 0.0 Silicon
solution name = Noni add const val = "0.5*(Charge + sqrt(Charge*Charge + 4.0*$ni*$ni)) / $ni" Silicon
solution name = Poni add const val = "0.5*(-Charge + sqrt(Charge*Charge + 4.0*$ni*$ni)) / $ni" Silicon

# Initialize vacancies near zero (physical: +1 model creates only interstitials)
sel z=1.0 name=Vac

# Cluster equilibrium initialization at 800C (ramp start temperature)
# K_cl >> K_BI (~2800x at 800C), so clusters dominate the B/I partition
# Quadratic: BIClust = K_cl * B_free * I_free with B_free = I_free = B_total - BIClust
# Solution: x = [(2KB+1) - sqrt(4KB+1)] / (2K)
set kT_init [expr {8.617e-5 * (800.0 + 273.15)}]
set K_cl_init [expr {2.0e-23 * exp(2.5 / $kT_init)}]
sel z=boron name=Inter
sel z=(2.0*${K_cl_init}*boron+1.0-sqrt(4.0*${K_cl_init}*boron+1.0))/(2.0*${K_cl_init}) name=BIClust
sel z=boron-BIClust name=Inter
sel z=Inter name=boron
sel z=1.0 name=boronInt
sel z=1.0 name=boronVac

# =============================================================================
# EQUILIBRIUM DEFECT CONCENTRATIONS (Fermi-level dependent)
# =============================================================================

set InumI "($Ineu + Noni*$Ineg + Poni*$Ipos)"
set InumV "($Vneu + Noni*($Vneg + Noni*$Vdng) + Poni*$Vpos)"

term name = EqInter add eqn = "$cis * $InumI / $IdenI" Silicon
term name = EqVac add eqn = "$cvs * $InumV / $IdenV" Silicon

term name = ScaleInter add eqn = "Inter/EqInter" Silicon
term name = ScaleVac add eqn = "Vac/EqVac" Silicon

# =============================================================================
# DOPANT-DEFECT PAIR TERMS (DopantReact model)
# Source: FLOOXS_2026/TclLib/Models/Dopant.tcl DopantReact, DopantDefectReact
# =============================================================================

# Solubility-limited active boron
term name = boronActive add eqn = "($Css) * boron / (($Css) + boron)" Silicon

# Substitutional boron: active minus paired (boron already excludes clustered B)
term name = boronSub add eqn = "boronActive - boronInt - boronVac" Silicon

# Charge neutrality: background donor - substitutional acceptor
term name = Charge add eqn = "1.4e15 - boronSub" Silicon

# Reaction terms: forward (pairing) - reverse (dissociation)
# React = Krate * (boronSub * Defect - boronDefect / Binding)
term name = ReactBI add eqn = "$KrateBI * (boronSub * Inter - boronInt / $Bind_BI)" Silicon
term name = ReactBV add eqn = "$KrateBV * (boronSub * Vac - boronVac / $Bind_BV)" Silicon

# Cluster reaction: boronSub + Inter <-> BIClust (effective 2.5 eV binding)
term name = ReactClust add eqn = "$KrateBI * (boronSub * Inter - BIClust / $Bind_Cl)" Silicon

# Pair diffusivities: (D0 + Dp*Poni) / (Binding * EqDefect * Poni)
# Absorbs charge factor so flux = DiffPair * grad(pair * Poni)
term name = DiffBI add eqn = "($D0_BI + $Dp_BI * Poni) / ($Bind_BI * EqInter * Poni)" Silicon
term name = DiffBV add eqn = "($D0_BV + $Dp_BV * Poni) / ($Bind_BV * EqVac * Poni)" Silicon

# =============================================================================
# DIFFUSION EQUATIONS
# =============================================================================

# --- Point defects: diffusion + IV recomb + pair reaction coupling ---
pdbSetString Silicon Inter Equation \
    "ddt(Inter) - $DiffI*EqInter*grad(ScaleInter) + $kIV*(Inter*Vac - EqInter*EqVac) + ReactBI + ReactClust"
pdbSetString Silicon Vac Equation \
    "ddt(Vac) - $DiffV*EqVac*grad(ScaleVac) + $kIV*(Inter*Vac - EqInter*EqVac) + ReactBV"

# --- Surface recombination at oxide/silicon interface ---
set EqI_surf "($cis * $InumI / $IdenI)"
set EqV_surf "($cvs * $InumV / $IdenV)"
pdbSetString Oxide_Silicon Inter Silicon Equation "-$KsurfI*(Inter(Silicon) - $EqI_surf)"
pdbSetString Oxide_Silicon Vac Silicon Equation "-$KsurfV*(Vac(Silicon) - $EqV_surf)"

# --- Total boron: no direct flux, only pair reaction source/sink ---
pdbSetString Silicon boron Equation \
    "ddt(boron) + ReactBI + ReactBV + ReactClust"

# --- BI pairs: diffusion + reaction ---
pdbSetString Silicon boronInt Equation \
    "ddt(boronInt) - DiffBI * grad(boronInt * Poni) - ReactBI"

# --- BV pairs: diffusion + reaction ---
pdbSetString Silicon boronVac Equation \
    "ddt(boronVac) - DiffBV * grad(boronVac * Poni) - ReactBV"

# --- BIClust: immobile cluster, formed from B+I, dissolved thermally ---
pdbSetString Silicon BIClust Equation \
    "ddt(BIClust) - ReactClust"

# =============================================================================
# OUTPUT: PRE-ANNEAL PROFILE
# =============================================================================

puts "=== PRE-ANNEAL ==="
sel z=boron
foreach v [print.1d] {
    lassign $v d c
    if {$c == "Value"} continue
    if {$d < 0} continue
    puts "$d $c"
}

# =============================================================================
# ANNEAL SEQUENCE (matches TSUPREM-4)
# =============================================================================

set ramprate [expr {(${temp} - 800.0) / 2700.0}]

puts "=== RAMP UP: 45 min, 800C -> ${temp}C ==="
diffuse time=45 temp=800 ramprate=$ramprate init=1e-6 damp.trbdf

puts "=== DWELL: ${time} min at ${temp}C ==="
diffuse time=${time} temp=${temp} init=1e-6 damp.trbdf

puts "=== RAMP DOWN: 45 min, ${temp}C -> 800C ==="
diffuse time=45 temp=${temp} ramprate=-$ramprate init=1e-6 damp.trbdf

# =============================================================================
# OUTPUT: POST-ANNEAL PROFILES
# =============================================================================

puts "=== POST-ANNEAL ==="
sel z=boron
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

puts "=== POST-ANNEAL CLUSTERS ==="
sel z=BIClust
foreach v [print.1d] {
    lassign $v d c
    if {$c == "Value"} continue
    if {$d < 0} continue
    if {$d > 0.5} break
    puts "$d $c"
}

# Junction depth from total boron (unpaired + BI pairs + BV pairs + clusters)
sel z=boron+boronInt+boronVac+BIClust-1.4e15
set xj [interpolate silicon val=0.0]
puts "=== RESULTS ==="
puts "Xj: $xj um"

exit
