# Ion implantation and diffusion simulation for PiezoD lookup tables
# FLOOXS QSS/Fermi hybrid model
#
# Parameters:
#   ${dopant}      - boron, phosphorus, or arsenic
#   ${dose}        - implant dose (cm^-2)
#   ${energy}      - implant energy (keV)
#   ${temp}        - anneal temperature (C)
#   ${time}        - anneal time (minutes)
#
# Physics:
#   Boron: QSS pair diffusion with explicit interstitial transport
#     - D_eff = D_I(Poni) * ScaleI + D_V(Poni)
#     - ScaleInter = I/I* from solved interstitial PDE
#     - D_V floor stabilizes solver; 2 solved variables (B + Inter)
#   Phosphorus: Fermi-level dependent diffusivity (no TED)
#     - D_eff = D0 + Dn*Noni + Dnn*Noni^2
#     - Inter PDE too stiff at high dose (2e16); fermi-only converges everywhere
#     - 1 solved variable (P only)
#   Arsenic: QSS pair diffusion with explicit interstitial transport (ScaleV=1)
#     - D_eff = D_I * ScaleI + D_V(Noni)
#     - ScaleInter = I/I* from solved interstitial PDE; V at equilibrium
#     - 2 solved variables (As + Inter)
#
# Background doping: ~1.4e15 cm^-3, 10 ohm-cm
# 250A protection oxide matching TSUPREM-4
# Temperature ramp matching TSUPREM-4
# All parameters from FLOOXS_2026/Params/Silicon/
#
# QSS P model (for reference, not used due to Inter PDE stiffness at high dose):
#   Solved 2 variables (P + Inter). ScaleInter soft-capped via Michaelis-Menten
#   (Smax * Inter / (Smax * EqInter + Inter - EqInter)) to bound enhancement.
#   With fI=0.2: D_eff = D_P(Noni) * (0.2*ScaleInter + 0.8).
#   Converges at 2e14 (Smax=5000 gives -7% vs TSUPREM-4) but diverges at 2e16.

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
# 6. POINT DEFECT PARAMETERS (boron and arsenic QSS)
# Source: FLOOXS_2026/Params/Silicon/Interstitial, Vacancy, Info
# =============================================================================

switch ${dopant} {
    boron {
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

        # Surface recombination: temperature-dependent KinkSite via Temp variable
        set KinkI "(0.51 * exp(2.63 / (8.617e-5 * (Temp + 273.15))))"
        set KsurfI "(3.14159*$DiffI*$lattice*$KinkI)"
    }
    arsenic {
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

        # Surface recombination: temperature-dependent KinkSite via Temp variable
        set KinkI "(0.51 * exp(2.63 / (8.617e-5 * (Temp + 273.15))))"
        set KsurfI "(3.14159*$DiffI*$lattice*$KinkI)"
    }
}

# =============================================================================
# 7. PER-DOPANT DIFFUSIVITY PARAMETERS
# Source: FLOOXS_2026/Params/Silicon/{Dopant}/Interstitial, Vacancy
# =============================================================================

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
        # P: I-only diffuser, combined D0+Dn+Dnn (fermi model)
        pdbSetDouble Silicon phosphorus Diff_D0 {[Arrhenius 5.6 3.71]}
        pdbSetDouble Silicon phosphorus Diff_Dn {[Arrhenius 6.38 4.05]}
        pdbSetDouble Silicon phosphorus Diff_Dnn {[Arrhenius 2.45e-2 3.23]}
    }
    arsenic {
        # As-I path: D0=Arr(0.0666,3.45)
        pdbSetDouble Silicon arsenic DI_D0 {[Arrhenius 0.0666 3.45]}
        # As-V path: Dn=Arr(12.8,4.05)
        pdbSetDouble Silicon arsenic DV_Dn {[Arrhenius 12.8 4.05]}
    }
}

# =============================================================================
# 8. SOLUTION DEFINITIONS
# =============================================================================

solution name=${dopant} solve !negative
switch ${dopant} {
    boron {
        solution name=Inter solve !negative add
    }
    arsenic {
        solution name=Inter solve !negative add
    }
}

# Fermi level
solution name = Charge add const val = 0.0 Silicon
solution name = Noni add const val = "0.5*(Charge + sqrt(Charge*Charge + 4.0*$ni*$ni)) / $ni" Silicon
solution name = Poni add const val = "0.5*(-Charge + sqrt(Charge*Charge + 4.0*$ni*$ni)) / $ni" Silicon

# =============================================================================
# 9. INITIALIZATION (QSS: Inter from implant, effective +n model)
# =============================================================================

# Effective +n model (Hobler-Moroz, TSUPREM-4 manual Eq 3-534)
# f_pl = 1 + D.PHDF * m^D.PME * E^(D.PLF * m^D.PLME)
# Default params: D.PHDF=0.0905, D.PME=0.85, D.PLF=-2, D.PLME=-0.5
switch ${dopant} {
    boron       { set mass 10.81 }
    phosphorus  { set mass 30.97 }
    arsenic     { set mass 74.92 }
}
set lambda_inf [expr {-2.0 * pow($mass, -0.5)}]
set f_pl [expr {1.0 + 0.0905 * pow($mass, 0.85) * pow(${energy}, $lambda_inf)}]
puts "f_pl ($mass amu, ${energy} keV) = $f_pl"

switch ${dopant} {
    boron {
        sel z=$f_pl*${dopant} name=Inter
    }
    arsenic {
        sel z=$f_pl*${dopant} name=Inter
    }
}

# =============================================================================
# 10. EQUILIBRIUM INTERSTITIALS (boron and arsenic QSS)
# =============================================================================

switch ${dopant} {
    boron {
        set InumI "($Ineu + Noni*$Ineg + Poni*$Ipos)"
        term name = EqInter add eqn = "$cis * $InumI / $IdenI" Silicon
        term name = ScaleInter add eqn = "Inter/EqInter" Silicon
    }
    arsenic {
        set InumI "($Ineu + Noni*$Ineg + Poni*$Ipos)"
        term name = EqInter add eqn = "$cis * $InumI / $IdenI" Silicon
        term name = ScaleInter add eqn = "Inter/EqInter" Silicon
    }
}

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

switch ${dopant} {
    boron {
        # QSS: D_eff = D_I * ScaleInter + D_V (V at equilibrium)
        set D0_I [pdbDelayDouble Silicon boron DI_D0]
        set Dp_I [pdbDelayDouble Silicon boron DI_Dp]
        set D0_V [pdbDelayDouble Silicon boron DV_D0]
        set Dp_V [pdbDelayDouble Silicon boron DV_Dp]
        term name = DiffDop add eqn = "(($D0_I + $Dp_I * Poni) * ScaleInter + ($D0_V + $Dp_V * Poni)) / Poni" Silicon
        pdbSetString Silicon boron Equation \
            "ddt(boron) - DiffDop * grad(boron * Poni)"
    }
    phosphorus {
        # Fermi: D_eff = D0 + Dn*Noni + Dnn*Noni^2
        set D0 [pdbDelayDouble Silicon phosphorus Diff_D0]
        set Dn [pdbDelayDouble Silicon phosphorus Diff_Dn]
        set Dnn [pdbDelayDouble Silicon phosphorus Diff_Dnn]
        term name = DiffDop add eqn = "($D0 + $Dn * Noni + $Dnn * Noni * Noni) / Noni" Silicon
        pdbSetString Silicon phosphorus Equation \
            "ddt(phosphorus) - DiffDop * grad(phosphorus * Noni)"
    }
    arsenic {
        # QSS: D_eff = D_I * ScaleInter + D_V(Noni) (ScaleV=1)
        set D0_I [pdbDelayDouble Silicon arsenic DI_D0]
        set Dn_V [pdbDelayDouble Silicon arsenic DV_Dn]
        term name = DiffDop add eqn = "($D0_I * ScaleInter + $Dn_V * Noni) / Noni" Silicon
        pdbSetString Silicon arsenic Equation \
            "ddt(arsenic) - DiffDop * grad(arsenic * Noni)"
    }
}

# =============================================================================
# 13. INTERSTITIAL EQUATION (boron and arsenic QSS)
# =============================================================================
#
# Defect transport with surface recombination as the sole I sink.
# No I-V recombination (no V equation): surface recomb dominates for
# near-surface implant damage.

switch ${dopant} {
    boron {
        pdbSetString Silicon Inter Equation \
            "ddt(Inter) - $DiffI*EqInter*grad(ScaleInter)"

        # Surface recombination BC
        set EqI_surf "($cis * $InumI / $IdenI)"
        pdbSetString Oxide_Silicon Inter Silicon Equation "-$KsurfI*(Inter(Silicon) - $EqI_surf)"
    }
    arsenic {
        pdbSetString Silicon Inter Equation \
            "ddt(Inter) - $DiffI*EqInter*grad(ScaleInter)"

        # Surface recombination BC
        set EqI_surf "($cis * $InumI / $IdenI)"
        pdbSetString Oxide_Silicon Inter Silicon Equation "-$KsurfI*(Inter(Silicon) - $EqI_surf)"
    }
}

# =============================================================================
# 14. PRE-ANNEAL OUTPUT
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
# 15. ANNEAL SEQUENCE (matches TSUPREM-4)
# =============================================================================

set ramprate [expr {(${temp} - 800.0) / 2700.0}]

puts "=== RAMP UP: 45 min, 800C -> ${temp}C ==="
diffuse time=45 temp=800 ramprate=$ramprate init=1e-6 damp.trbdf

puts "=== DWELL: ${time} min at ${temp}C ==="
diffuse time=${time} temp=${temp} init=1e-6 damp.trbdf

puts "=== RAMP DOWN: 45 min, ${temp}C -> 800C ==="
diffuse time=45 temp=${temp} ramprate=-$ramprate init=1e-6 damp.trbdf

# =============================================================================
# 16. POST-ANNEAL OUTPUT
# =============================================================================

puts "=== POST-ANNEAL ==="
sel z=${dopant}
foreach v [print.1d] {
    lassign $v d c
    if {$c == "Value"} continue
    if {$d < 0} continue
    puts "$d $c"
}

switch ${dopant} {
    boron {
        puts "=== POST-ANNEAL INTERSTITIALS ==="
        sel z=Inter
        foreach v [print.1d] {
            lassign $v d c
            if {$c == "Value"} continue
            if {$d < 0} continue
            if {$d > 0.5} break
            puts "$d $c"
        }
    }
    arsenic {
        puts "=== POST-ANNEAL INTERSTITIALS ==="
        sel z=Inter
        foreach v [print.1d] {
            lassign $v d c
            if {$c == "Value"} continue
            if {$d < 0} continue
            if {$d > 0.5} break
            puts "$d $c"
        }
    }
}

# =============================================================================
# 17. JUNCTION DEPTH
# =============================================================================

sel z=${dopant}-1.4e15
set xj [interpolate silicon val=0.0]
puts "=== RESULTS ==="
puts "Xj: $xj um"

exit
