# Ion implantation and diffusion simulation for PiezoD lookup tables
# FLOOXS QSS model: quasi-steady-state pair diffusion for B, P, As
#
# Parameters:
#   ${dopant}      - boron, phosphorus, or arsenic
#   ${dose}        - implant dose (cm^-2)
#   ${energy}      - implant energy (keV)
#   ${temp}        - anneal temperature (C)
#   ${time}        - anneal time (minutes)
#
# Physics:
#   - Fermi-level dependent dopant diffusivity (same as fermi model)
#   - Enhanced by ScaleInter = I/I* from explicit interstitial transport
#   - Quasi-steady-state pair approximation: pairs in local equilibrium,
#     folded into ScaleI/ScaleV multipliers on dopant diffusivity
#   - V assumed at equilibrium (ScaleVac = 1): implant creates I, not V
#   - Per-dopant effective diffusivity:
#       Boron:      D = D_I(Poni) * ScaleI + D_V(Poni) (B-V floor stabilizes solver)
#       Phosphorus: D = D_P(Noni) * (fI*ScaleI + 1-fI), fI(C) = 0.17..1.0
#                   (concentration-dependent fI: Jones/LBL kink-and-tail model)
#       Arsenic:    D = D_I + D_V(Noni) (fermi-equivalent, As-I binding=0 eV)
#   - Solved variables: 2 (dopant + Inter) for all dopants
#   - Background doping (~1.4e15 cm^-3, 10 ohm-cm)
#   - 250A protection oxide matching TSUPREM-4
#   - Temperature ramp matching TSUPREM-4
#
# All parameters from FLOOXS_2026/Params/Silicon/

math diffuse dim=1 umf none col scale

# =============================================================================
# 1. MESH
# =============================================================================

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

# =============================================================================
# 2. BACKGROUND DOPING
# =============================================================================

switch ${dopant} {
    boron {
        sel z=1.4e15 name=Phosphorus
    }
    phosphorus {
        sel z=1.4e15 name=Boron
    }
    arsenic {
        sel z=1.4e15 name=Boron
    }
}

# =============================================================================
# 3. ION IMPLANTATION
# =============================================================================

implant ${dopant} dose=${dose} energy=${energy} tilt=7

# =============================================================================
# 4. SOLVER SETTINGS
# =============================================================================

options !constdatafields storenodes
solution name=Temp !negative add const val=800

pdbSetDouble Math iterLimit 500
pdbSetDouble Math updateLimit 1e-8
pdbSetDouble Math rhsLimit 1e-20
pdbSetDouble Math rhsMin 1e-2

# =============================================================================
# 5. INTRINSIC CARRIER CONCENTRATION
# =============================================================================

set ni "(3.87e16 * sqrt((Temp+273.0)*(Temp+273.0)*(Temp+273.0)) * exp(-7014.0/(Temp+273.0)))"

# =============================================================================
# 6. POINT DEFECT PARAMETERS
# Source: FLOOXS_2026/Params/Silicon/Interstitial, Vacancy, Info
# =============================================================================

set lattice 2.714417617e-8

# Equilibrium interstitial concentration and charge states
pdbSetDouble Silicon Inter Cstar {[Arrhenius 3.6484e27 3.7]}
set cis [pdbDelayDouble Silicon Inter Cstar]

pdbSetDouble Silicon Inter neutral 1.0
pdbSetDouble Silicon Inter negative {[Arrhenius 5.68 0.48]}
pdbSetDouble Silicon Inter positive {[Arrhenius 5.68 0.42]}
set Ineu [pdbDelayDouble Silicon Inter neutral]
set Ineg [pdbDelayDouble Silicon Inter negative]
set Ipos [pdbDelayDouble Silicon Inter positive]
set IdenI "($Ineu + $Ineg + $Ipos)"

# Interstitial diffusivity
pdbSetDouble Silicon Inter Diff {[Arrhenius 0.138 1.37]}
set DiffI [pdbDelayDouble Silicon Inter Diff]

# Surface recombination coefficient
set KsurfI "(3.14159*$DiffI*$lattice*1.3e15)"

# =============================================================================
# 7. PER-DOPANT DIFFUSIVITY (split I and V components)
# Source: FLOOXS_2026/Params/Silicon/{Dopant}/Interstitial, Vacancy
# =============================================================================
#
# QSS: D_eff = D_I_components * ScaleInter + D_V_components * 1
# The I-component gets enhanced by ScaleI; V-component stays at equilibrium.

switch ${dopant} {
    boron {
        # B-I path: D0=Arr(0.743,3.56), Dp=Arr(0.617,3.56)
        pdbSetDouble Silicon boron DI_D0 {[Arrhenius 0.743 3.56]}
        pdbSetDouble Silicon boron DI_Dp {[Arrhenius 0.617 3.56]}
        # B-V path: D0=Arr(0.186,3.56), Dp=Arr(0.154,3.56)
        pdbSetDouble Silicon boron DV_D0 {[Arrhenius 0.186 3.56]}
        pdbSetDouble Silicon boron DV_Dp {[Arrhenius 0.154 3.56]}
    }
    phosphorus {
        # P-I path (I-only diffuser, Fi=1.0):
        # D0=Arr(5.6,3.71), Dn=Arr(6.38,4.05), Dnn=Arr(0.0245,3.23)
        pdbSetDouble Silicon phosphorus DI_D0 {[Arrhenius 5.6 3.71]}
        pdbSetDouble Silicon phosphorus DI_Dn {[Arrhenius 6.38 4.05]}
        pdbSetDouble Silicon phosphorus DI_Dnn {[Arrhenius 2.45e-2 3.23]}
    }
    arsenic {
        # As-I path: D0=Arr(0.0666,3.45) (weak, neutral only)
        pdbSetDouble Silicon arsenic DI_D0 {[Arrhenius 0.0666 3.45]}
        # As-V path: Dn=Arr(12.8,4.05) (dominant, electron-enhanced)
        pdbSetDouble Silicon arsenic DV_Dn {[Arrhenius 12.8 4.05]}
    }
}

# =============================================================================
# 8. SOLUTION DEFINITIONS
# =============================================================================

# 2 solved variables for all dopants: dopant + Inter
solution name=${dopant} solve !negative
solution name=Inter solve !negative add

# Fermi level
solution name = Charge add const val = 0.0 Silicon
solution name = Noni add const val = "0.5*(Charge + sqrt(Charge*Charge + 4.0*$ni*$ni)) / $ni" Silicon
solution name = Poni add const val = "0.5*(-Charge + sqrt(Charge*Charge + 4.0*$ni*$ni)) / $ni" Silicon

# =============================================================================
# 9. INITIALIZATION
# =============================================================================

# I from implant (+1 model: each implanted ion creates one interstitial)
sel z=${dopant} name=Inter

# =============================================================================
# 10. EQUILIBRIUM INTERSTITIALS (Fermi-level dependent)
# =============================================================================

set InumI "($Ineu + Noni*$Ineg + Poni*$Ipos)"
term name = EqInter add eqn = "$cis * $InumI / $IdenI" Silicon
term name = ScaleInter add eqn = "Inter/EqInter" Silicon

# =============================================================================
# 11. CHARGE NEUTRALITY
# =============================================================================

switch ${dopant} {
    boron {
        term name = Charge add eqn = "1.4e15 - ${dopant}" Silicon
    }
    phosphorus {
        term name = Charge add eqn = "${dopant} - 1.4e15" Silicon
    }
    arsenic {
        term name = Charge add eqn = "${dopant} - 1.4e15" Silicon
    }
}

# =============================================================================
# 12. DOPANT DIFFUSIVITY AND EQUATION
# =============================================================================
#
# QSS pair approximation: D_eff = D_I * ScaleInter + D_V * 1
# Flux = D_eff / chg * grad(C * chg)
#   chg = Poni (acceptors) or Noni (donors)

switch ${dopant} {
    boron {
        set D0_I [pdbDelayDouble Silicon boron DI_D0]
        set Dp_I [pdbDelayDouble Silicon boron DI_Dp]
        set D0_V [pdbDelayDouble Silicon boron DV_D0]
        set Dp_V [pdbDelayDouble Silicon boron DV_Dp]
        # D_eff = (D0_I + Dp_I*Poni)*ScaleI + (D0_V + Dp_V*Poni)
        term name = DiffDop add eqn = "(($D0_I + $Dp_I * Poni) * ScaleInter + ($D0_V + $Dp_V * Poni)) / Poni" Silicon
        pdbSetString Silicon boron Equation \
            "ddt(boron) - DiffDop * grad(boron * Poni)"
    }
    phosphorus {
        set D0_I [pdbDelayDouble Silicon phosphorus DI_D0]
        set Dn_I [pdbDelayDouble Silicon phosphorus DI_Dn]
        set Dnn_I [pdbDelayDouble Silicon phosphorus DI_Dnn]
        # Concentration-dependent fractional interstitial (Jones/LBL):
        #   fI -> 1.0 at C << ni (intrinsic, full TED)
        #   fI -> 0.17 at C >> ni (extrinsic, PV- E-centers dominate)
        # D_eff = D_P * (fI * ScaleI + (1 - fI))
        term name = FracI add eqn = "0.17 + 0.83 / (1.0 + phosphorus / $ni)" Silicon
        term name = DiffDop add eqn = "($D0_I + $Dn_I * Noni + $Dnn_I * Noni * Noni) * (FracI * ScaleInter + 1.0 - FracI) / Noni" Silicon
        pdbSetString Silicon phosphorus Equation \
            "ddt(phosphorus) - DiffDop * grad(phosphorus * Noni)"
    }
    arsenic {
        set D0_I [pdbDelayDouble Silicon arsenic DI_D0]
        set Dn_V [pdbDelayDouble Silicon arsenic DV_Dn]
        # As-I binding = 0 eV: pairs dissolve in ~6 ns, never reach QSS.
        # ScaleI enhancement is invalid for As. Use equilibrium diffusivity.
        # D_eff = D0_AsI + Dn_AsV * Noni (same as fermi model)
        term name = DiffDop add eqn = "($D0_I + $Dn_V * Noni) / Noni" Silicon
        pdbSetString Silicon arsenic Equation \
            "ddt(arsenic) - DiffDop * grad(arsenic * Noni)"
    }
}

# =============================================================================
# 13. INTERSTITIAL EQUATION
# =============================================================================
#
# Defect transport with surface recombination as the sole I sink.
# No I-V recombination (no V equation): surface recomb dominates for
# near-surface implant damage. No dopant coupling: QSS pairs don't
# create net I source/sink terms.

pdbSetString Silicon Inter Equation \
    "ddt(Inter) - $DiffI*EqInter*grad(ScaleInter)"

# =============================================================================
# 14. SURFACE RECOMBINATION BC
# =============================================================================

set EqI_surf "($cis * $InumI / $IdenI)"
pdbSetString Oxide_Silicon Inter Silicon Equation "-$KsurfI*(Inter(Silicon) - $EqI_surf)"

# =============================================================================
# 15. PRE-ANNEAL OUTPUT
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
# 16. PER-DOPANT SOLVER TUNING
# =============================================================================
#
# B: well-conditioned (D_V floor), tight tolerance OK
# P: stiffer (concentration-dependent fI provides partial floor),
#    slightly relaxed updateLimit for faster timestep growth
# As: fermi-equivalent (no ScaleI coupling), default settings

switch ${dopant} {
    phosphorus {
        pdbSetDouble Math updateLimit 1e-6
    }
}

# =============================================================================
# 17. ANNEAL SEQUENCE (matches TSUPREM-4)
# =============================================================================

set ramprate [expr {(${temp} - 800.0) / 2700.0}]

puts "=== RAMP UP: 45 min, 800C -> ${temp}C ==="
diffuse time=45 temp=800 ramprate=$ramprate init=1e-6 damp.trbdf

puts "=== DWELL: ${time} min at ${temp}C ==="
diffuse time=${time} temp=${temp} init=1e-6 damp.trbdf

puts "=== RAMP DOWN: 45 min, ${temp}C -> 800C ==="
diffuse time=45 temp=${temp} ramprate=-$ramprate init=1e-6 damp.trbdf

# =============================================================================
# 18. POST-ANNEAL OUTPUT
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

# =============================================================================
# 19. JUNCTION DEPTH
# =============================================================================

sel z=${dopant}-1.4e15
set xj [interpolate silicon val=0.0]
puts "=== RESULTS ==="
puts "Xj: $xj um"

exit
