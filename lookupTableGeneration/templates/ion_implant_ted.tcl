# Ion implantation and diffusion simulation for PiezoD lookup tables
# FLOOXS TED model: explicit dopant-defect pair diffusion for B, P, As
#
# Parameters:
#   ${dopant}      - boron, phosphorus, or arsenic
#   ${dose}        - implant dose (cm^-2)
#   ${energy}      - implant energy (keV)
#   ${temp}        - anneal temperature (C)
#   ${time}        - anneal time (minutes)
#
# Physics:
#   - Explicit interstitial and vacancy transport (not equilibrium)
#   - Dopant-defect pair formation, dissolution, and diffusion
#   - Fermi-level dependent defect concentrations and pair mobilities
#   - I-V bulk recombination
#   - Surface recombination at oxide/silicon interface
#   - Per-dopant pair kinetics:
#       Boron:      B-I + B-V pairs
#       Phosphorus: P-I pairs only (Fi=1.0, strong binding)
#       Arsenic:    As-I + As-V pairs (V-dominant, electron-enhanced)
#   - Background doping (~1.4e15 cm^-3, 10 ohm-cm)
#   - 250A protection oxide matching TSUPREM-4
#   - Temperature ramp matching TSUPREM-4
#
# All parameters from FLOOXS_2026/Params/Silicon/

math diffuse dim=1 umf none col scale

# =============================================================================
# 1. MESH
# =============================================================================

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

# =============================================================================
# 2. BACKGROUND DOPING
# =============================================================================

# ~1.4e15 cm^-3 (10 ohm-cm)
# Boron implant -> n-type substrate (phosphorus background)
# P/As implant -> p-type substrate (boron background)
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

# ni from Law et al. - accurate at process temperatures (900-1100C)
set ni "(3.87e16 * sqrt((Temp+273.0)*(Temp+273.0)*(Temp+273.0)) * exp(-7014.0/(Temp+273.0)))"

# =============================================================================
# 6. POINT DEFECT PARAMETERS
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

# 7. I-V bulk recombination rate: kIV = 4*pi*(D_I+D_V)*lattice
set kIV "(4.0*3.14159*($DiffI+$DiffV)*$lattice)"

# =============================================================================
# 8. PER-DOPANT PAIR PARAMETERS
# Source: FLOOXS_2026/Params/Silicon/{Dopant}/Interstitial, Vacancy
# =============================================================================

switch ${dopant} {
    boron {
        # BI pair: Binding=Arr(8e-23,-1.0), D0=Arr(0.743,3.56), Dp=Arr(0.617,3.56)
        pdbSetDouble Silicon boron BI_Binding {[Arrhenius 8.0e-23 -1.0]}
        pdbSetDouble Silicon boron BI_D0 {[Arrhenius 0.743 3.56]}
        pdbSetDouble Silicon boron BI_Dp {[Arrhenius 0.617 3.56]}
        set Bind_BI [pdbDelayDouble Silicon boron BI_Binding]
        set D0_BI [pdbDelayDouble Silicon boron BI_D0]
        set Dp_BI [pdbDelayDouble Silicon boron BI_Dp]

        # BV pair: Binding=Arr(8e-23,-0.5), D0=Arr(0.186,3.56), Dp=Arr(0.154,3.56)
        pdbSetDouble Silicon boron BV_Binding {[Arrhenius 8.0e-23 -0.5]}
        pdbSetDouble Silicon boron BV_D0 {[Arrhenius 0.186 3.56]}
        pdbSetDouble Silicon boron BV_Dp {[Arrhenius 0.154 3.56]}
        set Bind_BV [pdbDelayDouble Silicon boron BV_Binding]
        set D0_BV [pdbDelayDouble Silicon boron BV_D0]
        set Dp_BV [pdbDelayDouble Silicon boron BV_Dp]

        # Reaction rates: diffusion-limited
        set KrateBI "(4.0*3.14159*$DiffI*$lattice)"
        set KrateBV "(4.0*3.14159*$DiffV*$lattice)"
    }
    phosphorus {
        # PI pair: Binding=Arr(8e-23,-1.49), D0=Arr(5.6,3.71),
        #          Dn=Arr(6.38,4.05), Dnn=Arr(2.45e-2,3.23)
        pdbSetDouble Silicon phosphorus PI_Binding {[Arrhenius 8.0e-23 -1.49]}
        pdbSetDouble Silicon phosphorus PI_D0 {[Arrhenius 5.6 3.71]}
        pdbSetDouble Silicon phosphorus PI_Dn {[Arrhenius 6.38 4.05]}
        pdbSetDouble Silicon phosphorus PI_Dnn {[Arrhenius 2.45e-2 3.23]}
        set Bind_PI [pdbDelayDouble Silicon phosphorus PI_Binding]
        set D0_PI [pdbDelayDouble Silicon phosphorus PI_D0]
        set Dn_PI [pdbDelayDouble Silicon phosphorus PI_Dn]
        set Dnn_PI [pdbDelayDouble Silicon phosphorus PI_Dnn]

        # Reaction rate: diffusion-limited
        set KratePI "(4.0*3.14159*$DiffI*$lattice)"
    }
    arsenic {
        # AsI pair: Binding=Arr(8e-23,0.0), D0=Arr(0.0666,3.45)
        # Weak binding (0.0 eV), negligible at moderate ScaleI
        pdbSetDouble Silicon arsenic AsI_Binding {[Arrhenius 8.0e-23 0.0]}
        pdbSetDouble Silicon arsenic AsI_D0 {[Arrhenius 0.0666 3.45]}
        set Bind_AsI [pdbDelayDouble Silicon arsenic AsI_Binding]
        set D0_AsI [pdbDelayDouble Silicon arsenic AsI_D0]

        # AsV pair: Binding=Arr(8e-23,-0.5), Dn=Arr(12.8,4.05)
        # Dominant diffusion path (electron-enhanced, D0=0)
        pdbSetDouble Silicon arsenic AsV_Binding {[Arrhenius 8.0e-23 -0.5]}
        pdbSetDouble Silicon arsenic AsV_Dn {[Arrhenius 12.8 4.05]}
        set Bind_AsV [pdbDelayDouble Silicon arsenic AsV_Binding]
        set Dn_AsV [pdbDelayDouble Silicon arsenic AsV_Dn]

        # Reaction rates: diffusion-limited
        set KrateAsI "(4.0*3.14159*$DiffI*$lattice)"
        set KrateAsV "(4.0*3.14159*$DiffV*$lattice)"
    }
}

# =============================================================================
# 9. PER-DOPANT SOLUBILITY
# Source: FLOOXS_2026/Params/Silicon/{Dopant}/Info
# =============================================================================

switch ${dopant} {
    boron {
        pdbSetDouble Silicon boron Solubility {[Arrhenius 7.68e22 0.7086]}
    }
    phosphorus {
        pdbSetDouble Silicon phosphorus Solubility {[Arrhenius 3.89e21 0.265]}
    }
    arsenic {
        pdbSetDouble Silicon arsenic Solubility {[Arrhenius 2.24e22 0.494]}
    }
}
set Css [pdbDelayDouble Silicon ${dopant} Solubility]

# =============================================================================
# 10. SURFACE RECOMBINATION PARAMETERS
# Source: FLOOXS_2026/Params/Oxide_Silicon/Interstitial, Vacancy
# =============================================================================

set KsurfI "(3.14159*$DiffI*$lattice*1.3e15)"
set KsurfV "(3.14159*$DiffV*$lattice*1.0e5)"

# =============================================================================
# 12. SOLUTION DEFINITIONS (per-dopant variable sets)
# =============================================================================

solution name=Inter solve !negative add

switch ${dopant} {
    boron {
        # 5 vars: boron, boronInt, boronVac, Inter, Vac
        solution name=Vac solve !negative add
        solution name=boron solve !negative
        solution name=boronInt solve !negative add
        solution name=boronVac solve !negative add
    }
    phosphorus {
        # 3 vars: phosphorus, phosphorusInt, Inter (no Vac: P is I-only)
        solution name=phosphorus solve !negative
        solution name=phosphorusInt solve !negative add
    }
    arsenic {
        # 5 vars: arsenic, arsenicInt, arsenicVac, Inter, Vac
        solution name=Vac solve !negative add
        solution name=arsenic solve !negative
        solution name=arsenicInt solve !negative add
        solution name=arsenicVac solve !negative add
    }
}

# =============================================================================
# 13. FERMI LEVEL
# =============================================================================

solution name = Charge add const val = 0.0 Silicon
solution name = Noni add const val = "0.5*(Charge + sqrt(Charge*Charge + 4.0*$ni*$ni)) / $ni" Silicon
solution name = Poni add const val = "0.5*(-Charge + sqrt(Charge*Charge + 4.0*$ni*$ni)) / $ni" Silicon

# =============================================================================
# 14. INITIALIZATION
# =============================================================================

# I from implant (+1 model: each implanted ion creates one interstitial)
sel z=${dopant} name=Inter

set kT_init [expr {8.617e-5 * 1073.15}]

# Per-dopant initialization
switch ${dopant} {
    boron {
        # V = 1.0 (near zero): IV recomb negligible initially, allowing
        # excess I to persist and drive TED through BI pair diffusion.
        sel z=1.0 name=Vac
        # Pairs initialized to small values (solver forms them dynamically)
        sel z=1.0 name=boronInt
        sel z=1.0 name=boronVac
    }
    phosphorus {
        # No Vac variable for P (I-only diffuser, no P-V pairs).
        # I sinks: surface recombination + P-I pair formation.
        # Pairs initialized small; solver forms them dynamically.
        sel z=1.0 name=phosphorusInt
    }
    arsenic {
        # V at thermal equilibrium for 800C (ramp start temperature).
        # Physical: implant creates I, not V. V_eq avoids extreme IV stiffness.
        set V_eq_init [expr {4.0515e26 * exp(-3.97 / $kT_init)}]
        sel z=$V_eq_init name=Vac

        # No pre-partitioning needed (weak AsI binding 0 eV, V_eq tiny)
        sel z=1.0 name=arsenicInt
        sel z=1.0 name=arsenicVac
    }
}

# =============================================================================
# 15. EQUILIBRIUM DEFECT CONCENTRATIONS (Fermi-level dependent)
# =============================================================================

set InumI "($Ineu + Noni*$Ineg + Poni*$Ipos)"
set InumV "($Vneu + Noni*($Vneg + Noni*$Vdng) + Poni*$Vpos)"

term name = EqInter add eqn = "$cis * $InumI / $IdenI" Silicon
term name = ScaleInter add eqn = "Inter/EqInter" Silicon

# EqVac and ScaleVac only for dopants that solve Vac (B, As; not P)
switch ${dopant} {
    boron -
    arsenic {
        term name = EqVac add eqn = "$cvs * $InumV / $IdenV" Silicon
        term name = ScaleVac add eqn = "Vac/EqVac" Silicon
    }
}

# =============================================================================
# 16. PER-DOPANT PAIR TERMS
# =============================================================================

# Solubility-limited active dopant
term name = ${dopant}Active add eqn = "($Css) * ${dopant} / (($Css) + ${dopant})" Silicon

switch ${dopant} {
    boron {
        # Substitutional boron: active minus paired
        term name = boronSub add eqn = "boronActive - boronInt - boronVac" Silicon

        # Reaction terms: Krate * (Sub * Defect - Pair / Binding)
        term name = ReactBI add eqn = "$KrateBI * (boronSub * Inter - boronInt / $Bind_BI)" Silicon
        term name = ReactBV add eqn = "$KrateBV * (boronSub * Vac - boronVac / $Bind_BV)" Silicon

        # Pair diffusivities: (D0 + Dp*Poni) / (Binding * EqDefect * Poni)
        term name = DiffBI add eqn = "($D0_BI + $Dp_BI * Poni) / ($Bind_BI * EqInter * Poni)" Silicon
        term name = DiffBV add eqn = "($D0_BV + $Dp_BV * Poni) / ($Bind_BV * EqVac * Poni)" Silicon
    }
    phosphorus {
        # Substitutional phosphorus: active minus paired
        term name = phosphorusSub add eqn = "phosphorusActive - phosphorusInt" Silicon

        # Reaction term: Krate * (Sub * Defect - Pair / Binding)
        term name = ReactPI add eqn = "$KratePI * (phosphorusSub * Inter - phosphorusInt / $Bind_PI)" Silicon

        # Pair diffusivity: (D0 + Dn*Noni + Dnn*Noni^2) / (Binding * EqDefect * Noni)
        term name = DiffPI add eqn = "($D0_PI + $Dn_PI * Noni + $Dnn_PI * Noni * Noni) / ($Bind_PI * EqInter * Noni)" Silicon
    }
    arsenic {
        # Substitutional arsenic: active minus paired
        term name = arsenicSub add eqn = "arsenicActive - arsenicInt - arsenicVac" Silicon

        # Reaction terms
        term name = ReactAsI add eqn = "$KrateAsI * (arsenicSub * Inter - arsenicInt / $Bind_AsI)" Silicon
        term name = ReactAsV add eqn = "$KrateAsV * (arsenicSub * Vac - arsenicVac / $Bind_AsV)" Silicon

        # Pair diffusivities
        # AsI: D0 only (no charge enhancement), acceptor=Noni for donor
        term name = DiffAsI add eqn = "($D0_AsI) / ($Bind_AsI * EqInter * Noni)" Silicon
        # AsV: Dn*Noni only (electron-enhanced, D0=0)
        term name = DiffAsV add eqn = "($Dn_AsV * Noni) / ($Bind_AsV * EqVac * Noni)" Silicon
    }
}

# =============================================================================
# 17. CHARGE NEUTRALITY (using Sub, not total dopant)
# =============================================================================

switch ${dopant} {
    boron {
        # Acceptor on n-type: Charge = N_D - N_A (substitutional)
        term name = Charge add eqn = "1.4e15 - boronSub" Silicon
    }
    phosphorus {
        # Donor on p-type: Charge = N_D - N_A (substitutional)
        term name = Charge add eqn = "phosphorusSub - 1.4e15" Silicon
    }
    arsenic {
        # Donor on p-type: Charge = N_D - N_A (substitutional)
        term name = Charge add eqn = "arsenicSub - 1.4e15" Silicon
    }
}

# =============================================================================
# 18. DIFFUSION EQUATIONS
# =============================================================================

# --- Point defects: diffusion + IV recomb + pair reaction coupling ---

switch ${dopant} {
    boron {
        pdbSetString Silicon Inter Equation \
            "ddt(Inter) - $DiffI*EqInter*grad(ScaleInter) + $kIV*(Inter*Vac - EqInter*EqVac) + ReactBI"
        pdbSetString Silicon Vac Equation \
            "ddt(Vac) - $DiffV*EqVac*grad(ScaleVac) + $kIV*(Inter*Vac - EqInter*EqVac) + ReactBV"

        # Total boron: no direct flux, only pair reaction source/sink
        pdbSetString Silicon boron Equation \
            "ddt(boron) + ReactBI + ReactBV"

        # BI pairs: diffusion + reaction
        pdbSetString Silicon boronInt Equation \
            "ddt(boronInt) - DiffBI * grad(boronInt * Poni) - ReactBI"

        # BV pairs: diffusion + reaction
        pdbSetString Silicon boronVac Equation \
            "ddt(boronVac) - DiffBV * grad(boronVac * Poni) - ReactBV"
    }
    phosphorus {
        # No Vac equation for P (I-only diffuser).
        # I sinks: surface recombination + P-I pair formation.
        pdbSetString Silicon Inter Equation \
            "ddt(Inter) - $DiffI*EqInter*grad(ScaleInter) + ReactPI"

        # Total phosphorus: no direct flux
        pdbSetString Silicon phosphorus Equation \
            "ddt(phosphorus) + ReactPI"

        # PI pairs: diffusion + reaction (donor: use Noni)
        pdbSetString Silicon phosphorusInt Equation \
            "ddt(phosphorusInt) - DiffPI * grad(phosphorusInt * Noni) - ReactPI"
    }
    arsenic {
        pdbSetString Silicon Inter Equation \
            "ddt(Inter) - $DiffI*EqInter*grad(ScaleInter) + $kIV*(Inter*Vac - EqInter*EqVac) + ReactAsI"
        pdbSetString Silicon Vac Equation \
            "ddt(Vac) - $DiffV*EqVac*grad(ScaleVac) + $kIV*(Inter*Vac - EqInter*EqVac) + ReactAsV"

        # Total arsenic: no direct flux
        pdbSetString Silicon arsenic Equation \
            "ddt(arsenic) + ReactAsI + ReactAsV"

        # AsI pairs: diffusion + reaction (donor: use Noni)
        pdbSetString Silicon arsenicInt Equation \
            "ddt(arsenicInt) - DiffAsI * grad(arsenicInt * Noni) - ReactAsI"

        # AsV pairs: diffusion + reaction (donor: use Noni)
        pdbSetString Silicon arsenicVac Equation \
            "ddt(arsenicVac) - DiffAsV * grad(arsenicVac * Noni) - ReactAsV"
    }
}

# =============================================================================
# 19. SURFACE RECOMBINATION BCS
# =============================================================================

set EqI_surf "($cis * $InumI / $IdenI)"
pdbSetString Oxide_Silicon Inter Silicon Equation "-$KsurfI*(Inter(Silicon) - $EqI_surf)"

# Vac surface recombination only for dopants that solve Vac
switch ${dopant} {
    boron -
    arsenic {
        set EqV_surf "($cvs * $InumV / $IdenV)"
        pdbSetString Oxide_Silicon Vac Silicon Equation "-$KsurfV*(Vac(Silicon) - $EqV_surf)"
    }
}

# =============================================================================
# 20. PRE-ANNEAL OUTPUT
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
# 21. ANNEAL SEQUENCE (matches TSUPREM-4)
# =============================================================================

# Ramp rate: (temp - 800) / (45 min * 60 s/min) in C/s
set ramprate [expr {(${temp} - 800.0) / 2700.0}]

puts "=== RAMP UP: 45 min, 800C -> ${temp}C ==="
diffuse time=45 temp=800 ramprate=$ramprate init=1e-6 damp.trbdf

puts "=== DWELL: ${time} min at ${temp}C ==="
diffuse time=${time} temp=${temp} init=1e-6 damp.trbdf

puts "=== RAMP DOWN: 45 min, ${temp}C -> 800C ==="
diffuse time=45 temp=${temp} ramprate=-$ramprate init=1e-6 damp.trbdf

# =============================================================================
# 22. POST-ANNEAL OUTPUT
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

switch ${dopant} {
    boron -
    arsenic {
        puts "=== POST-ANNEAL VACANCIES ==="
        sel z=Vac
        foreach v [print.1d] {
            lassign $v d c
            if {$c == "Value"} continue
            if {$d < 0} continue
            if {$d > 0.5} break
            puts "$d $c"
        }
    }
}

switch ${dopant} {
    boron {
        puts "=== POST-ANNEAL PAIRS-BI ==="
        sel z=boronInt
        foreach v [print.1d] {
            lassign $v d c
            if {$c == "Value"} continue
            if {$d < 0} continue
            if {$d > 0.5} break
            puts "$d $c"
        }
        puts "=== POST-ANNEAL PAIRS-BV ==="
        sel z=boronVac
        foreach v [print.1d] {
            lassign $v d c
            if {$c == "Value"} continue
            if {$d < 0} continue
            if {$d > 0.5} break
            puts "$d $c"
        }
    }
    phosphorus {
        puts "=== POST-ANNEAL PAIRS-PI ==="
        sel z=phosphorusInt
        foreach v [print.1d] {
            lassign $v d c
            if {$c == "Value"} continue
            if {$d < 0} continue
            if {$d > 0.5} break
            puts "$d $c"
        }
    }
    arsenic {
        puts "=== POST-ANNEAL PAIRS-AsI ==="
        sel z=arsenicInt
        foreach v [print.1d] {
            lassign $v d c
            if {$c == "Value"} continue
            if {$d < 0} continue
            if {$d > 0.5} break
            puts "$d $c"
        }
        puts "=== POST-ANNEAL PAIRS-AsV ==="
        sel z=arsenicVac
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
# 23. JUNCTION DEPTH (total dopant = free + pairs)
# =============================================================================

switch ${dopant} {
    boron {
        sel z=boron+boronInt+boronVac-1.4e15
    }
    phosphorus {
        sel z=phosphorus+phosphorusInt-1.4e15
    }
    arsenic {
        sel z=arsenic+arsenicInt+arsenicVac-1.4e15
    }
}
set xj [interpolate silicon val=0.0]
puts "=== RESULTS ==="
puts "Xj: $xj um"

exit
