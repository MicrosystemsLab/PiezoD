# Ion implantation and diffusion simulation for PiezoD lookup tables
# FLOOXS template with Fermi-level dependent diffusion
#
# Parameters:
#   ${dopant}      - boron, phosphorus, or arsenic
#   ${dose}        - implant dose (cm^-2)
#   ${energy}      - implant energy (keV)
#   ${temp}        - anneal temperature (C)
#   ${time}        - anneal time (minutes)
#
# Physics:
#   - Fermi-level dependent dopant diffusivity (charge neutrality)
#   - Detailed-balance defect diffusion: D * Eq * grad(C/Eq)
#   - Charge-state dependent equilibrium defect concentrations
#   - I-V bulk recombination (diffusion-limited)
#   - Surface recombination at oxide/silicon interface
#   - Dopant diffusion via interstitial and vacancy pair mechanisms
#   - Background doping (phosphorus, 10 ohm-cm n-type)
#   - Temperature ramp matching TSUPREM-4
#
# All parameters from FLOOXS_2026/Params/Silicon/

# scale = row-column equilibration of the Jacobian (A[i,j] *= 1/sqrt(max_row_i * max_col_j))
# Essential for coupled systems with variables spanning many orders of magnitude
# (e.g. dopant ~1e21, defects ~1e15, Poni ~500)
math diffuse dim=1 umf none col scale

# 1D mesh (5um depth, fine near surface for implant resolution)
line x loc=-0.001 spac=0.001 tag=TopOx
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
# Matches TSUPREM-4: initialize <100> impurity=phosphorus i.resistivity=10
sel z=1.4e15 name=Phosphorus

# Ion implantation
implant ${dopant} dose=${dose} energy=${energy} tilt=7

# Initialize interstitials from implant damage (+1 model)
sel z=${dopant} name=Inter

# Set up solutions
options !constdatafields storenodes
solution name=Temp !negative add const val=800

# Solver settings
pdbSetDouble Math iterLimit 500
pdbSetDouble Math updateLimit 1e-8
pdbSetDouble Math rhsLimit 1e-20
pdbSetDouble Math rhsMin 1e-2

# =============================================================================
# THERMAL VOLTAGE AND INTRINSIC CARRIER CONCENTRATION
# =============================================================================

# Vt = kT/q in eV (T in Celsius via FLOOXS Temp variable)
set Vt "(8.617383e-5*(Temp+273.0))"

# ni from Law et al. - accurate at process temperatures (900-1100C)
# ni(T) = 3.87e16 * T_K^1.5 * exp(-7014/T_K) cm^-3
# Gives ni(1000C) ~ 7.1e18, ni(900C) ~ 3.5e18, ni(1100C) ~ 1.4e19
set ni "(3.87e16 * sqrt((Temp+273.0)*(Temp+273.0)*(Temp+273.0)) * exp(-7014.0/(Temp+273.0)))"

# =============================================================================
# POINT DEFECT PARAMETERS (temperature-dependent via Arrhenius)
# Source: FLOOXS_2026/Params/Silicon/Interstitial, Vacancy, Info
# =============================================================================

set lattice 2.714417617e-8

# Equilibrium concentrations: Cstar = prefactor * exp(-Ea/kT)
# Interstitial: 5.0e22*exp(11.2) = 3.6484e27
# Vacancy: 5.0e22*exp(9.0) = 4.0515e26
pdbSetDouble Silicon Inter Cstar {[Arrhenius 3.6484e27 3.7]}
pdbSetDouble Silicon Vac Cstar {[Arrhenius 4.0515e26 3.97]}
set cis [pdbDelayDouble Silicon Inter Cstar]
set cvs [pdbDelayDouble Silicon Vac Cstar]

# Interstitial charge-state contributions (from Params/Silicon/Interstitial)
# neutral=1.0, negative=Arrhenius(5.68, 0.48), positive=Arrhenius(5.68, 0.42)
pdbSetDouble Silicon Inter neutral 1.0
pdbSetDouble Silicon Inter negative {[Arrhenius 5.68 0.48]}
pdbSetDouble Silicon Inter positive {[Arrhenius 5.68 0.42]}
set Ineu [pdbDelayDouble Silicon Inter neutral]
set Ineg [pdbDelayDouble Silicon Inter negative]
set Ipos [pdbDelayDouble Silicon Inter positive]
set IdenI "($Ineu + $Ineg + $Ipos)"

# Vacancy charge-state contributions (from Params/Silicon/Vacancy)
# neutral=1.0, negative=Arrhenius(5.68, 0.145), positive=Arrhenius(5.68, 0.455)
# dnegative=Arrhenius(32.47, 0.62)
pdbSetDouble Silicon Vac neutral 1.0
pdbSetDouble Silicon Vac negative {[Arrhenius 5.68 0.145]}
pdbSetDouble Silicon Vac positive {[Arrhenius 5.68 0.455]}
pdbSetDouble Silicon Vac dnegative {[Arrhenius 32.47 0.62]}
set Vneu [pdbDelayDouble Silicon Vac neutral]
set Vneg [pdbDelayDouble Silicon Vac negative]
set Vpos [pdbDelayDouble Silicon Vac positive]
set Vdng [pdbDelayDouble Silicon Vac dnegative]
set IdenV "($Vneu + $Vneg + $Vdng + $Vpos)"

# Defect diffusivities (D0 in FLOOXS convention)
pdbSetDouble Silicon Inter Diff {[Arrhenius 0.138 1.37]}
pdbSetDouble Silicon Vac Diff {[Arrhenius 1.18e-4 0.1]}
set DiffI [pdbDelayDouble Silicon Inter Diff]
set DiffV [pdbDelayDouble Silicon Vac Diff]

# I-V bulk recombination rate: kIV = 4*pi*(D_I+D_V)*lattice
set kIV "(4.0*3.14159*($DiffI+$DiffV)*$lattice)"

# =============================================================================
# DOPANT DIFFUSION PARAMETERS (Fermi-level dependent)
# Source: FLOOXS_2026/Params/Silicon/{Dopant}/Interstitial, Vacancy
# =============================================================================
#
# Split into D0 (neutral) and Dp/Dn (charge-enhanced) components.
# Acceptors (boron): D = D0 + Dp*Poni, where Poni = p/ni
# Donors (phosphorus, arsenic): D = D0 + Dn*Noni + Dnn*Noni^2, where Noni = n/ni

switch ${dopant} {
    boron {
        # Acceptor: Fermi enhancement via Poni = p/ni
        # B-Interstitial: D0=Arr(0.743,3.56), Dp=Arr(0.617,3.56)
        pdbSetDouble Silicon ${dopant} DiffI_D0 {[Arrhenius 0.743 3.56]}
        pdbSetDouble Silicon ${dopant} DiffI_Dp {[Arrhenius 0.617 3.56]}
        # B-Vacancy: D0=Arr(0.186,3.56), Dp=Arr(0.154,3.56)
        pdbSetDouble Silicon ${dopant} DiffV_D0 {[Arrhenius 0.186 3.56]}
        pdbSetDouble Silicon ${dopant} DiffV_Dp {[Arrhenius 0.154 3.56]}
    }
    phosphorus {
        # Donor: Fermi enhancement via Noni = n/ni
        # P-Interstitial: D0=Arr(5.6,3.71), Dn=Arr(6.38,4.05), Dnn=Arr(0.0245,3.23)
        pdbSetDouble Silicon ${dopant} DiffI_D0 {[Arrhenius 5.6 3.71]}
        pdbSetDouble Silicon ${dopant} DiffI_Dn {[Arrhenius 6.38 4.05]}
        pdbSetDouble Silicon ${dopant} DiffI_Dnn {[Arrhenius 2.45e-2 3.23]}
        # No vacancy path for phosphorus
    }
    arsenic {
        # Donor: Fermi enhancement via Noni = n/ni
        # As-Interstitial: D0=Arr(0.0666,3.45), Dn=0
        pdbSetDouble Silicon ${dopant} DiffI_D0 {[Arrhenius 0.0666 3.45]}
        # As-Vacancy: D0=0, Dn=Arr(12.8,4.05)
        pdbSetDouble Silicon ${dopant} DiffV_Dn {[Arrhenius 12.8 4.05]}
    }
}

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

# Fermi level coupling: Noni and Poni from charge neutrality (no Potential variable)
# n = 0.5 * (Charge + sqrt(Charge^2 + 4*ni^2)), Noni = n/ni, Poni = p/ni
solution name = Charge add const val = 0.0 Silicon
solution name = Noni add const val = "0.5*(Charge + sqrt(Charge*Charge + 4.0*$ni*$ni)) / $ni" Silicon
solution name = Poni add const val = "0.5*(-Charge + sqrt(Charge*Charge + 4.0*$ni*$ni)) / $ni" Silicon

# Initialize vacancies from implant damage (Frenkel pair model)
# Each implanted atom creates one I-V pair, so initial V = initial I
sel z=${dopant} name=Vac

# =============================================================================
# EQUILIBRIUM DEFECT CONCENTRATIONS (Fermi-level dependent)
# Source: FLOOXS_2026/TclLib/Models/Defect.tcl DefectBulk proc
# =============================================================================
#
# EqI = Cstar_I * (neu + Noni*neg + Poni*pos) / (neu + neg + pos)
#
# Noni and Poni are now live solution variables (not hardcoded to 1.0).
# In p-type (boron): Poni >> 1, enhancing Eq_I via positive charge state.
# In n-type (phosphorus/arsenic): Noni >> 1, enhancing Eq_V via negative charge state.

set InumI "($Ineu + Noni*$Ineg + Poni*$Ipos)"
set InumV "($Vneu + Noni*($Vneg + Noni*$Vdng) + Poni*$Vpos)"

term name = EqInter add eqn = "$cis * $InumI / $IdenI" Silicon
term name = EqVac add eqn = "$cvs * $InumV / $IdenV" Silicon

# Scale variables: Sol/Eq (thermodynamic driving force for diffusion)
term name = ScaleInter add eqn = "Inter/EqInter" Silicon
term name = ScaleVac add eqn = "Vac/EqVac" Silicon

# =============================================================================
# CHARGE NEUTRALITY (Fermi level solve)
# Source: FLOOXS_2026/TclLib/Models/Potential.tcl PotentialEqns proc
# =============================================================================
#
# Charge = net ionized donor concentration (positive = n-type, negative = p-type)
# n = 0.5 * (Charge + sqrt(Charge^2 + 4*ni^2))
# Potential / Vt = log(n / ni)
#
# This is the algebraic (non-Poisson) charge neutrality model, equivalent to
# TSUPREM-4's "Fermi" diffusion model.

switch ${dopant} {
    boron {
        # Acceptor: Charge = N_D(background) - N_A(boron)
        term name = Charge add eqn = "1.4e15 - ${dopant}" Silicon
    }
    phosphorus {
        # Donor: Charge = N_D(phosphorus), includes background + implanted
        term name = Charge add eqn = "${dopant}" Silicon
    }
    arsenic {
        # Donor: Charge = N_D(arsenic) + N_D(P_background)
        term name = Charge add eqn = "${dopant} + 1.4e15" Silicon
    }
}

# =============================================================================
# DIFFUSION EQUATIONS
# =============================================================================

# --- Point defect transport (detailed balance with Fermi-corrected equilibrium) ---
# EqInter/EqVac now depend on Noni/Poni, so they vary with local Fermi level
pdbSetString Silicon Inter Equation "ddt(Inter) - $DiffI*EqInter*grad(ScaleInter) + $kIV*(Inter*Vac - EqInter*EqVac)"
pdbSetString Silicon Vac Equation "ddt(Vac) - $DiffV*EqVac*grad(ScaleVac) + $kIV*(Inter*Vac - EqInter*EqVac)"

# --- Surface recombination at oxide/silicon interface ---
# Use Fermi-corrected equilibrium (Noni/Poni are solutions, available at boundary)
set EqI_surf "($cis * $InumI / $IdenI)"
set EqV_surf "($cvs * $InumV / $IdenV)"
pdbSetString Oxide_Silicon Inter Silicon Equation "-$KsurfI*(Inter(Silicon) - $EqI_surf)"
pdbSetString Oxide_Silicon Vac Silicon Equation "-$KsurfV*(Vac(Silicon) - $EqV_surf)"

# --- Dopant transport (Fermi-enhanced with defect pair mechanism) ---
# DopantPair form: flux = (D/chg) * grad(C * ScaleDefect * chg)
# where chg = Poni for acceptors, Noni for donors
# D = D0 + Dp*Poni (acceptors) or D0 + Dn*Noni + Dnn*Noni^2 (donors)

switch ${dopant} {
    boron {
        set D0_I [pdbDelayDouble Silicon ${dopant} DiffI_D0]
        set Dp_I [pdbDelayDouble Silicon ${dopant} DiffI_Dp]
        set D0_V [pdbDelayDouble Silicon ${dopant} DiffV_D0]
        set Dp_V [pdbDelayDouble Silicon ${dopant} DiffV_Dp]
        # Effective diffusivity / Poni for each path
        term name = DiffDopI add eqn = "($D0_I + $Dp_I * Poni) / Poni" Silicon
        term name = DiffDopV add eqn = "($D0_V + $Dp_V * Poni) / Poni" Silicon
        pdbSetString Silicon ${dopant} Equation \
            "ddt(${dopant}) - DiffDopI * grad(${dopant} * ScaleInter * Poni) - DiffDopV * grad(${dopant} * ScaleVac * Poni)"
    }
    phosphorus {
        set D0_I [pdbDelayDouble Silicon ${dopant} DiffI_D0]
        set Dn_I [pdbDelayDouble Silicon ${dopant} DiffI_Dn]
        set Dnn_I [pdbDelayDouble Silicon ${dopant} DiffI_Dnn]
        # Effective diffusivity / Noni
        term name = DiffDopI add eqn = "($D0_I + $Dn_I * Noni + $Dnn_I * Noni * Noni) / Noni" Silicon
        pdbSetString Silicon ${dopant} Equation \
            "ddt(${dopant}) - DiffDopI * grad(${dopant} * ScaleInter * Noni)"
    }
    arsenic {
        set D0_I [pdbDelayDouble Silicon ${dopant} DiffI_D0]
        set Dn_V [pdbDelayDouble Silicon ${dopant} DiffV_Dn]
        # I-path: only D0 (no Fermi enhancement), V-path: only Dn (Noni cancels)
        term name = DiffDopI add eqn = "$D0_I / Noni" Silicon
        term name = DiffDopV add eqn = "$Dn_V" Silicon
        pdbSetString Silicon ${dopant} Equation \
            "ddt(${dopant}) - DiffDopI * grad(${dopant} * ScaleInter * Noni) - DiffDopV * grad(${dopant} * ScaleVac * Noni)"
    }
}

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
diffuse time=45 temp=800 ramprate=$ramprate init=1e-6 damp.trbdf

puts "=== DWELL: ${time} min at ${temp}C ==="
diffuse time=${time} temp=${temp} init=1e-6 damp.trbdf

puts "=== RAMP DOWN: 45 min, ${temp}C -> 800C ==="
diffuse time=45 temp=${temp} ramprate=-$ramprate init=1e-6 damp.trbdf

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

# Junction depth where implanted dopant equals background (1.4e15 cm^-3)
sel z=${dopant}-1.4e15
set xj [interpolate silicon val=0.0]
puts "=== RESULTS ==="
puts "Xj: $xj um"

exit
