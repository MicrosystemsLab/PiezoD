# Ion implantation and diffusion simulation for PiezoD lookup tables
# Effective diffusivity model v2: FLOOXS pair D + linear IV relaxation.
#
# Boron: FLOOXS_2026 point defects + v5 {311} clustering + v5 surface recomb
#   - FLOOXS_2026 pair-derived effective D (QSS-equivalent to v5 explicit pairs)
#   - Linear IV relaxation: kIV*EqVac*(Inter-EqInter) replaces stiff Vac PDE
#     (V equilibrates in microseconds due to 235,000x faster D_V vs D_I)
#   - alpha_311=0.1 (matches v5; no recalibration needed with correct D values)
#   - 3 PDEs: boron, Inter, CIc
#
# P/As: TSUPREM-4 Appendix A parameters (unchanged from v1)
#   - Table A-17 point defects, Table A-24 {311}, Tables A-18/A-19 surface recomb
#   - P: 3 PDEs (phosphorus, Inter, CIc)
#   - As: 4 PDEs (arsenic, Inter, Vac, CIc)
#
# Parameters:
#   ${dopant}      - boron, phosphorus, or arsenic
#   ${dose}        - implant dose (cm^-2)
#   ${energy}      - implant energy (keV)
#   ${temp}        - anneal temperature (C)
#   ${time}        - anneal time (minutes)

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
# 5. INTRINSIC CARRIER CONCENTRATION (Law et al.)
# =============================================================================

set ni "(3.87e16 * sqrt((Temp+273.0)*(Temp+273.0)*(Temp+273.0)) * exp(-7014.0/(Temp+273.0)))"

# =============================================================================
# 6. POINT DEFECT PARAMETERS
# =============================================================================

switch ${dopant} {
    boron {
        # FLOOXS_2026 (well-calibrated, matches v5 results)
        set lattice 2.714417617e-8

        pdbSetDouble Silicon Inter Cstar {[Arrhenius 3.6484e27 3.7]}
        pdbSetDouble Silicon Vac Cstar {[Arrhenius 4.0515e26 3.97]}

        pdbSetDouble Silicon Inter neutral 1.0
        pdbSetDouble Silicon Inter negative {[Arrhenius 5.68 0.48]}
        pdbSetDouble Silicon Inter positive {[Arrhenius 5.68 0.42]}

        pdbSetDouble Silicon Vac neutral 1.0
        pdbSetDouble Silicon Vac negative {[Arrhenius 5.68 0.145]}
        pdbSetDouble Silicon Vac positive {[Arrhenius 5.68 0.455]}
        pdbSetDouble Silicon Vac dnegative {[Arrhenius 32.47 0.62]}

        pdbSetDouble Silicon Inter Diff {[Arrhenius 0.138 1.37]}
        pdbSetDouble Silicon Vac Diff {[Arrhenius 1.18e-4 0.1]}
    }
    phosphorus -
    arsenic {
        # Table A-17
        set lattice 2.714e-8

        pdbSetDouble Silicon Inter Cstar {[Arrhenius 1.25e29 3.26]}
        pdbSetDouble Silicon Vac Cstar {[Arrhenius 1.25e29 3.26]}

        pdbSetDouble Silicon Inter neutral 1.0
        pdbSetDouble Silicon Inter negative {[Arrhenius 5.68 0.50]}
        pdbSetDouble Silicon Inter positive {[Arrhenius 5.68 0.26]}

        pdbSetDouble Silicon Vac neutral 1.0
        pdbSetDouble Silicon Vac negative {[Arrhenius 5.68 0.145]}
        pdbSetDouble Silicon Vac positive {[Arrhenius 5.68 0.455]}
        pdbSetDouble Silicon Vac dnegative {[Arrhenius 32.47 0.62]}

        pdbSetDouble Silicon Inter Diff {[Arrhenius 3.65e-4 1.58]}
        pdbSetDouble Silicon Vac Diff {[Arrhenius 3.65e-4 1.58]}
    }
}

set cis [pdbDelayDouble Silicon Inter Cstar]
set cvs [pdbDelayDouble Silicon Vac Cstar]

set Ineu [pdbDelayDouble Silicon Inter neutral]
set Ineg [pdbDelayDouble Silicon Inter negative]
set Ipos [pdbDelayDouble Silicon Inter positive]
set IdenI "($Ineu + $Ineg + $Ipos)"

set Vneu [pdbDelayDouble Silicon Vac neutral]
set Vneg [pdbDelayDouble Silicon Vac negative]
set Vpos [pdbDelayDouble Silicon Vac positive]
set Vdng [pdbDelayDouble Silicon Vac dnegative]
set IdenV "($Vneu + $Vneg + $Vdng + $Vpos)"

set DiffI [pdbDelayDouble Silicon Inter Diff]
set DiffV [pdbDelayDouble Silicon Vac Diff]

# =============================================================================
# 7. {311} CLUSTERING PARAMETERS
# =============================================================================

switch ${dopant} {
    boron {
        # v5 reformulated: R_clust = Kr * (alpha * ScaleI * (CIc+Inter) - CIc)
        # Kr Ea=3.6 eV (Eaglesham et al.), alpha=0.1 (matches v5; FLOOXS pair D
        # gives same D_eff as v5 at QSS, so no recalibration needed)
        set Kr_311 {[Arrhenius 2.0e11 3.6]}
        set alpha_311 0.1
    }
    phosphorus -
    arsenic {
        # Table A-24
        set Kfc_311 {[Arrhenius 5.207e14 3.774]}
        set Kr_311 {[Arrhenius 9.431e13 3.017]}
    }
}

# =============================================================================
# 8. IV BULK RECOMBINATION (Eq 3-286)
# =============================================================================

# kIV = 4*pi*(D_I + D_V)*a0
# B: used in linearized IV relaxation (no Vac PDE)
# As: used in full IV recombination (with Vac PDE)
# P: I-only, no IV term
set kIV "(4.0*3.14159*($DiffI+$DiffV)*$lattice)"

# =============================================================================
# 9. PER-DOPANT EFFECTIVE DIFFUSIVITIES
# =============================================================================

switch ${dopant} {
    boron {
        # FLOOXS_2026 pair-derived (QSS-equivalent to v5 explicit pairs)
        # At QSS, eff model D maps directly to pair diffusivities:
        # DIX=D0_BI, DIP=Dp_BI, DVX=D0_BV, DVP=Dp_BV
        # Binding energy cancels exactly in the derivation.
        # All Ea = 3.56 eV (FLOOXS_2026 pair migration barrier)
        pdbSetDouble Silicon boron DIX {[Arrhenius 0.743 3.56]}
        pdbSetDouble Silicon boron DIP {[Arrhenius 0.617 3.56]}
        pdbSetDouble Silicon boron DVX {[Arrhenius 0.186 3.56]}
        pdbSetDouble Silicon boron DVP {[Arrhenius 0.154 3.56]}
        set DIX_B [pdbDelayDouble Silicon boron DIX]
        set DIP_B [pdbDelayDouble Silicon boron DIP]
        set DVX_B [pdbDelayDouble Silicon boron DVX]
        set DVP_B [pdbDelayDouble Silicon boron DVP]
    }
    phosphorus {
        # I pathway only: neutral (DIX) + minus (DIM * Noni) + double-minus (DIMM * Noni^2)
        pdbSetDouble Silicon phosphorus DIX {[Arrhenius 3.85 3.66]}
        pdbSetDouble Silicon phosphorus DIM {[Arrhenius 4.44 4.00]}
        pdbSetDouble Silicon phosphorus DIMM {[Arrhenius 44.2 4.37]}
        set DIX_P [pdbDelayDouble Silicon phosphorus DIX]
        set DIM_P [pdbDelayDouble Silicon phosphorus DIM]
        set DIMM_P [pdbDelayDouble Silicon phosphorus DIMM]
    }
    arsenic {
        # I pathway: neutral (DIX) + minus (DIM * Noni)
        # V pathway: neutral (DVX) + minus (DVM * Noni)
        pdbSetDouble Silicon arsenic DIX {[Arrhenius 1.37e7 3.44]}
        pdbSetDouble Silicon arsenic DIM {[Arrhenius 6.20 4.15]}
        pdbSetDouble Silicon arsenic DVX {[Arrhenius 5.47e7 3.44]}
        pdbSetDouble Silicon arsenic DVM {[Arrhenius 24.83 4.15]}
        set DIX_As [pdbDelayDouble Silicon arsenic DIX]
        set DIM_As [pdbDelayDouble Silicon arsenic DIM]
        set DVX_As [pdbDelayDouble Silicon arsenic DVX]
        set DVM_As [pdbDelayDouble Silicon arsenic DVM]
    }
}

# =============================================================================
# 10. PER-DOPANT SOLUBILITY (FLOOXS_2026)
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
# 11. SURFACE RECOMBINATION PARAMETERS
# =============================================================================

switch ${dopant} {
    boron {
        # v5: Arrhenius KinkSite, calibrated I sink rate
        set KinkI {[Arrhenius 0.51 -2.63]}
        set KsurfI "(3.14159*$DiffI*$lattice*$KinkI)"
    }
    phosphorus -
    arsenic {
        # Tables A-18/A-19
        set KsurfI {[Arrhenius 1.4e-6 -1.75]}
        set KsurfV {[Arrhenius 4.0e-11 -1.75]}
    }
}

# =============================================================================
# 12. SOLUTION DEFINITIONS
# =============================================================================

solution name=Inter solve !negative add
solution name=CIc solve !negative add

switch ${dopant} {
    boron {
        # 3 PDEs: boron, Inter, CIc (no Vac PDE; IV handled by linear relaxation)
        solution name=boron solve !negative
    }
    phosphorus {
        # 3 PDEs: phosphorus, Inter, CIc
        solution name=phosphorus solve !negative
    }
    arsenic {
        # 4 PDEs: arsenic, Inter, Vac, CIc
        solution name=Vac solve !negative add
        solution name=arsenic solve !negative
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

# Effective +n model (Hobler-Moroz, Eq 3-534)
switch ${dopant} {
    boron       { set mass 10.81 }
    phosphorus  { set mass 30.97 }
    arsenic     { set mass 74.92 }
}
set lambda_inf [expr {-2.0 * pow($mass, -0.5)}]
set f_pl [expr {1.0 + 0.0905 * pow($mass, 0.85) * pow(${energy}, $lambda_inf)}]
puts "f_pl ($mass amu, ${energy} keV) = $f_pl"

set kT_init [expr {8.617e-5 * 1073.15}]

# Pre-cluster: most excess I starts in {311} clusters
set f_cluster 0.9
sel z=(1.0-$f_cluster)*$f_pl*${dopant} name=Inter
sel z=$f_cluster*$f_pl*${dopant} name=CIc
puts "Init: Inter=(1-f)*f_pl*dop, CIc=f*f_pl*dop (f=$f_cluster, f_pl=$f_pl)"

# Vac at thermal equilibrium for 800C (ramp start) - As only
switch ${dopant} {
    arsenic {
        set V_eq_init [expr {1.25e29 * exp(-3.26 / $kT_init)}]
        sel z=$V_eq_init name=Vac
    }
}

# =============================================================================
# 15. EQUILIBRIUM DEFECT CONCENTRATIONS (Fermi-level dependent)
# =============================================================================

set InumI "($Ineu + Noni*$Ineg + Poni*$Ipos)"
set InumV "($Vneu + Noni*($Vneg + Noni*$Vdng) + Poni*$Vpos)"

term name = EqInter add eqn = "$cis * $InumI / $IdenI" Silicon
term name = ScaleInter add eqn = "Inter/EqInter" Silicon

# EqVac: algebraic for B (IV relaxation sink), PDE-coupled for As
switch ${dopant} {
    boron {
        # Algebraic EqVac for linear IV relaxation (no Vac PDE needed)
        term name = EqVac add eqn = "$cvs * $InumV / $IdenV" Silicon
    }
    arsenic {
        term name = EqVac add eqn = "$cvs * $InumV / $IdenV" Silicon
        term name = ScaleVac add eqn = "Vac/EqVac" Silicon
    }
}

# =============================================================================
# 16. SOLUBILITY-LIMITED ACTIVE DOPANT
# =============================================================================

term name = ${dopant}Active add eqn = "($Css) * ${dopant} / (($Css) + ${dopant})" Silicon

# =============================================================================
# 17. CHARGE NEUTRALITY (using active dopant, no pairs tracked)
# =============================================================================

switch ${dopant} {
    boron {
        term name = Charge add eqn = "1.4e15 - boronActive" Silicon
    }
    phosphorus {
        term name = Charge add eqn = "phosphorusActive - 1.4e15" Silicon
    }
    arsenic {
        term name = Charge add eqn = "arsenicActive - 1.4e15" Silicon
    }
}

# =============================================================================
# 18. PER-DOPANT EFFECTIVE D TERMS
# Eq 3-30/31/45/54: DiffDopI = D_I_dop * ScaleI / eta
# where D_I_dop = sum of charge-state components
# and eta = Poni (acceptors) or Noni (donors)
# =============================================================================

switch ${dopant} {
    boron {
        # DiffDopI = (DIX + DIP * Poni) * ScaleInter / Poni
        # DiffDopV = (DVX + DVP * Poni) * 1.0 / Poni  (V at equilibrium, ScaleVac=1)
        term name = DiffDopI add eqn = "($DIX_B + $DIP_B * Poni) * ScaleInter / Poni" Silicon
        term name = DiffDopV add eqn = "($DVX_B + $DVP_B * Poni) * 1.0 / Poni" Silicon
        term name = DiffDop add eqn = "DiffDopI + DiffDopV" Silicon
    }
    phosphorus {
        # DiffDopI = (DIX + DIM * Noni + DIMM * Noni^2) * ScaleInter / Noni
        term name = DiffDopI add eqn = "($DIX_P + $DIM_P * Noni + $DIMM_P * Noni * Noni) * ScaleInter / Noni" Silicon
        term name = DiffDop add eqn = "DiffDopI" Silicon
    }
    arsenic {
        # DiffDopI = (DIX + DIM * Noni) * ScaleInter / Noni
        # DiffDopV = (DVX + DVM * Noni) * ScaleVac / Noni
        term name = DiffDopI add eqn = "($DIX_As + $DIM_As * Noni) * ScaleInter / Noni" Silicon
        term name = DiffDopV add eqn = "($DVX_As + $DVM_As * Noni) * ScaleVac / Noni" Silicon
        term name = DiffDop add eqn = "DiffDopI + DiffDopV" Silicon
    }
}

# =============================================================================
# 19. {311} CLUSTERING
# =============================================================================

switch ${dopant} {
    boron {
        # v5 reformulated form (avoids cancellation stiffness at QSS)
        # R_clust = Kr * (alpha * ScaleInter * (CIc + Inter) - CIc)
        term name = R_clust add eqn = "$Kr_311 * ($alpha_311 * ScaleInter * (CIc + Inter) - CIc)" Silicon
    }
    phosphorus -
    arsenic {
        # Table A-24: R_clust = Kfc * ScaleInter * (CIc + Inter) - Kr * CIc
        term name = R_clust add eqn = "$Kfc_311 * ScaleInter * (CIc + Inter) - $Kr_311 * CIc" Silicon
    }
}

# =============================================================================
# 20. DOPANT EQUATION (Eq 3-30: dC/dt = -div(J_m + J_n))
# =============================================================================

switch ${dopant} {
    boron {
        pdbSetString Silicon boron Equation \
            "ddt(boron) - DiffDop * grad(boron * Poni)"
    }
    phosphorus {
        pdbSetString Silicon phosphorus Equation \
            "ddt(phosphorus) - DiffDop * grad(phosphorus * Noni)"
    }
    arsenic {
        pdbSetString Silicon arsenic Equation \
            "ddt(arsenic) - DiffDop * grad(arsenic * Noni)"
    }
}

# =============================================================================
# 21. INTER EQUATION (Eq 3-280 with effective D coupling)
# =============================================================================

switch ${dopant} {
    boron {
        # Linear IV relaxation: kIV*EqVac*(Inter - EqInter)
        # Since D_V >> D_I (235,000x), V equilibrates in microseconds.
        # Assuming V=EqVac, the IV recomb term kIV*(I*V - I*V*) becomes
        # kIV*EqVac*(I - I*), a linear damping toward equilibrium.
        # tau = 1/(kIV*EqVac): 30s at 900C, 0.75s at 1000C, 0.08s at 1100C.
        pdbSetString Silicon Inter Equation \
            "ddt(Inter) - $DiffI*EqInter*grad(ScaleInter) - DiffDopI*grad(boron*Poni) + $kIV*EqVac*(Inter - EqInter) + R_clust"
    }
    phosphorus {
        # No IV recombination (no Vac PDE)
        pdbSetString Silicon Inter Equation \
            "ddt(Inter) - $DiffI*EqInter*grad(ScaleInter) - DiffDopI*grad(phosphorus*Noni) + R_clust"
    }
    arsenic {
        pdbSetString Silicon Inter Equation \
            "ddt(Inter) - $DiffI*EqInter*grad(ScaleInter) - DiffDopI*grad(arsenic*Noni) + $kIV*(Inter*Vac - EqInter*EqVac) + R_clust"
    }
}

# =============================================================================
# 22. VAC EQUATION (Eq 3-281 with pair coupling, As only)
# =============================================================================

switch ${dopant} {
    arsenic {
        pdbSetString Silicon Vac Equation \
            "ddt(Vac) - $DiffV*EqVac*grad(ScaleVac) - DiffDopV*grad(arsenic*Noni) + $kIV*(Inter*Vac - EqInter*EqVac)"
    }
}

# =============================================================================
# 23. CIc EQUATION (Eq 3-319)
# =============================================================================

pdbSetString Silicon CIc Equation "ddt(CIc) - R_clust"

# =============================================================================
# 24. SURFACE RECOMBINATION BCS
# =============================================================================

set EqI_surf "($cis * $InumI / $IdenI)"
pdbSetString Oxide_Silicon Inter Silicon Equation "-$KsurfI*(Inter(Silicon) - $EqI_surf)"

switch ${dopant} {
    arsenic {
        set EqV_surf "($cvs * $InumV / $IdenV)"
        pdbSetString Oxide_Silicon Vac Silicon Equation "-$KsurfV*(Vac(Silicon) - $EqV_surf)"
    }
}

# =============================================================================
# 25. PRE-ANNEAL OUTPUT
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
# 26. ANNEAL SEQUENCE (matches TSUPREM-4)
# =============================================================================

set ramprate [expr {(${temp} - 800.0) / 2700.0}]

puts "=== RAMP UP: 45 min, 800C -> ${temp}C ==="
diffuse time=45 temp=800 ramprate=$ramprate init=1e-6 damp.trbdf

puts "=== DWELL: ${time} min at ${temp}C ==="
diffuse time=${time} temp=${temp} init=1e-6 damp.trbdf

puts "=== RAMP DOWN: 45 min, ${temp}C -> 800C ==="
diffuse time=45 temp=${temp} ramprate=-$ramprate init=1e-6 damp.trbdf

# =============================================================================
# 27. POST-ANNEAL OUTPUT
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

puts "=== POST-ANNEAL CLUSTERED-I ==="
sel z=CIc
foreach v [print.1d] {
    lassign $v d c
    if {$c == "Value"} continue
    if {$d < 0} continue
    if {$d > 0.5} break
    puts "$d $c"
}

# =============================================================================
# 28. JUNCTION DEPTH (total dopant, no pairs to add)
# =============================================================================

sel z=${dopant}-1.4e15
set xj [interpolate silicon val=0.0]
puts "=== RESULTS ==="
puts "Xj: $xj um"

exit
