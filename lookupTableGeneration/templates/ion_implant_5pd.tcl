# Ion implantation and diffusion simulation for PiezoD lookup tables
# FLOOXS 5-stream + {311} clustering model: explicit dopant-defect pair
# diffusion for B, P, As with {311} interstitial clustering.
#
# CIc tracks clustered interstitials (Eq 3-319 from TSUPREM-4 manual).
#
# Parameter sources (TSUPREM-4 manual Appendix A):
#   Point defects (D, C*, charges): Table A-17 (PDF p.752-753)
#   {311} clustering (Kfc, Kr):     Table A-24 (PDF p.755)
#   Surface recombination (Ksurf):  Tables A-18/A-19 (PDF p.753-754)
#   Pair kinetics (Kf, Kr):         Table A-6 (PDF p.748)
#   Pair diffusivities (D_pair):    Table A-3 (PDF p.746-747)
#   Effective +n (Hobler-Moroz):    Table A-46 (PDF p.761)
#   Solid solubility:               FLOOXS_2026 (not tabulated in TSUPREM-4)
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
#   - {311} interstitial clustering (simplified 1-moment model)
#   - Per-dopant pair kinetics:
#       Boron:      B-I + B-V pairs
#       Phosphorus: P-I pairs only (Fi=1.0, strong binding)
#       Arsenic:    As-I + As-V pairs (V-dominant, electron-enhanced)
#   - Background doping (~1.4e15 cm^-3, 10 ohm-cm)
#   - 250A protection oxide matching TSUPREM-4
#   - Temperature ramp matching TSUPREM-4

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
# Source: TSUPREM-4 manual Appendix A, Table A-17 (PDF p.752-753)
# =============================================================================

# Pre-compute all temperature-dependent Arrhenius parameters at anneal temp.
# This avoids Arrhenius expressions in the PDE (solver stability) and lets us
# use TSUPREM-4 Appendix A values directly as Tcl scalars.
set kT_anneal [expr {8.617e-5 * (${temp} + 273.15)}]

set lattice 2.714417617e-8

# Equilibrium concentrations (Table A-17: CEQUIL.0=1.25e29, CEQUIL.E=3.26)
pdbSetDouble Silicon Inter Cstar {[Arrhenius 1.25e29 3.26]}
pdbSetDouble Silicon Vac Cstar {[Arrhenius 1.25e29 3.26]}
set cis [pdbDelayDouble Silicon Inter Cstar]
set cvs [pdbDelayDouble Silicon Vac Cstar]

# Interstitial charge states (Table A-17)
pdbSetDouble Silicon Inter neutral 1.0
pdbSetDouble Silicon Inter negative {[Arrhenius 5.68 0.50]}
pdbSetDouble Silicon Inter positive {[Arrhenius 5.68 0.26]}
set Ineu [pdbDelayDouble Silicon Inter neutral]
set Ineg [pdbDelayDouble Silicon Inter negative]
set Ipos [pdbDelayDouble Silicon Inter positive]
set IdenI "($Ineu + $Ineg + $Ipos)"

# Vacancy charge states (Table A-17)
pdbSetDouble Silicon Vac neutral 1.0
pdbSetDouble Silicon Vac negative {[Arrhenius 5.68 0.145]}
pdbSetDouble Silicon Vac positive {[Arrhenius 5.68 0.455]}
pdbSetDouble Silicon Vac dnegative {[Arrhenius 32.47 0.62]}
set Vneu [pdbDelayDouble Silicon Vac neutral]
set Vneg [pdbDelayDouble Silicon Vac negative]
set Vpos [pdbDelayDouble Silicon Vac positive]
set Vdng [pdbDelayDouble Silicon Vac dnegative]
set IdenV "($Vneu + $Vneg + $Vdng + $Vpos)"

# Defect diffusivities (Table A-17: D.0=3.65e-4, D.E=1.58 for both I and V)
pdbSetDouble Silicon Inter Diff {[Arrhenius 3.65e-4 1.58]}
pdbSetDouble Silicon Vac Diff {[Arrhenius 3.65e-4 1.58]}
set DiffI [pdbDelayDouble Silicon Inter Diff]
set DiffV [pdbDelayDouble Silicon Vac Diff]

# Pre-computed defect diffusivities at anneal temp (for reaction rates, surface recomb)
set DiffI_val [expr {3.65e-4 * exp(-1.58 / $kT_anneal)}]
set DiffV_val [expr {3.65e-4 * exp(-1.58 / $kT_anneal)}]

# {311} interstitial clustering
# Source: TSUPREM-4 manual Appendix A, Table A-24 (PDF p.755)
#   CL.KFC.0=5.207e14, CL.KFC.E=3.774 (forward clustering)
#   CL.KR.0=9.431e13, CL.KR.E=3.017 (dissolution)
#   CL.CF=0.9398 (approximated as 1.0; FLOOXS lacks pow())
set Kfc_311 [expr {5.207e14 * exp(-3.774 / $kT_anneal)}]
set Kr_311  [expr {9.431e13 * exp(-3.017 / $kT_anneal)}]
puts "311 at ${temp}C: Kfc=$Kfc_311, Kr=$Kr_311, Kfc/Kr=[expr {$Kfc_311/$Kr_311}]"

# 7. I-V bulk recombination rate: kIV = 4*pi*(D_I+D_V)*lattice
set kIV "(4.0*3.14159*($DiffI+$DiffV)*$lattice)"

# =============================================================================
# 8. PER-DOPANT PAIR PARAMETERS
# Source: TSUPREM-4 manual Appendix A
#   Pair kinetics: Table A-6 (PDF p.748)
#   Pair diffusivities: Table A-3 (PDF p.746-747)
#
# Formation rate: Kf = 4*pi*D_defect*a0 (diffusion-limited capture)
# Dissolution rate: Kr = R.I.S = R.V.S = 10 s^-1 (Table A-6, constant)
# Pair diffusivities: pre-computed at anneal temp from Table A-3.
#   Table A-3 units: um^2/min -> cm^2/s via *1.667e-10, except where noted cm^2/s.
# =============================================================================

# Pair formation rates (diffusion-limited)
set Kf_I [expr {4.0 * 3.14159 * $DiffI_val * $lattice}]
set Kf_V [expr {4.0 * 3.14159 * $DiffV_val * $lattice}]
set Kr_pair 10.0

switch ${dopant} {
    boron {
        # BI pair diffusivity (Table A-3, PDF p.746)
        #   DIX.0=2.11e8 um2/min (->0.0352 cm2/s), DIP.0=4.10e9 um2/min (->0.683 cm2/s)
        #   Ea=3.46 eV (same for all charge states)
        set DIX_BI [expr {0.0352 * exp(-3.46 / $kT_anneal)}]
        set DIP_BI [expr {0.683  * exp(-3.46 / $kT_anneal)}]

        # BV pair diffusivity (Table A-3)
        #   DVX.0=1.11e7 um2/min (->1.85e-3 cm2/s), DVP.0=2.16e8 um2/min (->0.0360 cm2/s)
        #   Ea=3.46 eV
        set DVX_BV [expr {1.85e-3 * exp(-3.46 / $kT_anneal)}]
        set DVP_BV [expr {0.0360  * exp(-3.46 / $kT_anneal)}]
    }
    phosphorus {
        # PI pair diffusivity (Table A-3, PDF p.746)
        #   DIX.0=2.31e10 um2/min (->3.85 cm2/s), Ea=3.66
        #   DIM.0=2.664e10 um2/min (->4.44 cm2/s), Ea=4.00
        #   DIMM.0=2.652e11 um2/min (->44.2 cm2/s), Ea=4.37
        set DIX_PI  [expr {3.85  * exp(-3.66 / $kT_anneal)}]
        set DIM_PI  [expr {4.44  * exp(-4.00 / $kT_anneal)}]
        set DIMM_PI [expr {44.2  * exp(-4.37 / $kT_anneal)}]
    }
    arsenic {
        # AsI pair diffusivity (Table A-3, PDF p.746)
        #   DIX.0=1.37e7 cm2/s, Ea=3.44
        #   DIM.0=3.72e10 um2/min (->6.20 cm2/s), Ea=4.15
        set DIX_AsI [expr {1.37e7 * exp(-3.44 / $kT_anneal)}]
        set DIM_AsI [expr {6.20   * exp(-4.15 / $kT_anneal)}]

        # AsV pair diffusivity (Table A-3)
        #   DVX.0=5.47e7 cm2/s, Ea=3.44
        #   DVM.0=1.49e11 um2/min (->24.83 cm2/s), Ea=4.15
        set DVX_AsV [expr {5.47e7 * exp(-3.44 / $kT_anneal)}]
        set DVM_AsV [expr {24.83  * exp(-4.15 / $kT_anneal)}]
    }
}

# =============================================================================
# 9. PER-DOPANT SOLUBILITY
# Source: FLOOXS_2026/Params/Silicon/{Dopant}/Info
# (TSUPREM-4 Table A-14 is tabulated, not Arrhenius; irrelevant at our doses)
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
# Source: TSUPREM-4 manual Appendix A
#   Table A-19 (PDF p.753-754): Inter KSURF.0=1.4e-6, KSURF.E=-1.75
#   Table A-18 (PDF p.753):     Vac  KSURF.0=4.0e-11, KSURF.E=-1.75
# Same for all dopants. Negative Ea -> rate increases with temperature.
# =============================================================================

set KsurfI_val [expr {1.4e-6 * exp(1.75 / $kT_anneal)}]
set KsurfV_val [expr {4.0e-11 * exp(1.75 / $kT_anneal)}]
puts "Ksurf at ${temp}C: I=$KsurfI_val, V=$KsurfV_val"

# =============================================================================
# 12. SOLUTION DEFINITIONS (per-dopant variable sets)
# =============================================================================

solution name=Inter solve !negative add
solution name=CIc solve !negative add

switch ${dopant} {
    boron {
        # 6 vars: boron, boronInt, boronVac, Inter, Vac, CIc
        solution name=Vac solve !negative add
        solution name=boron solve !negative
        solution name=boronInt solve !negative add
        solution name=boronVac solve !negative add
    }
    phosphorus {
        # 4 vars: phosphorus, phosphorusInt, Inter, CIc (no Vac: P is I-only)
        solution name=phosphorus solve !negative
        solution name=phosphorusInt solve !negative add
    }
    arsenic {
        # 6 vars: arsenic, arsenicInt, arsenicVac, Inter, Vac, CIc
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

set kT_init [expr {8.617e-5 * 1073.15}]

# Pre-cluster: most excess I starts in {311} clusters, small fraction free.
# Equivalent to rapid clustering in first microseconds (CL.INI.F analog).
# Avoids initial stiff transient where ScaleInter >> 1.
set f_cluster 0.9
sel z=(1.0-$f_cluster)*$f_pl*${dopant} name=Inter
sel z=$f_cluster*$f_pl*${dopant} name=CIc
puts "Init: Inter=(1-f)*f_pl*dop, CIc=f*f_pl*dop (f=$f_cluster, f_pl=$f_pl)"

# Per-dopant initialization
switch ${dopant} {
    boron {
        # V at thermal equilibrium for 800C (ramp start temperature).
        # With TSUPREM-4 CEQUIL (1.25e29/3.26), V_eq(800C) ~ 5.5e13.
        # Must initialize near equilibrium to avoid extreme kIV stiffness.
        set V_eq_init [expr {1.25e29 * exp(-3.26 / $kT_init)}]
        sel z=$V_eq_init name=Vac
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
        set V_eq_init [expr {1.25e29 * exp(-3.26 / $kT_init)}]
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

        # Reaction: TSUPREM-4 Eq 3-61 with Table A-6 kinetics
        # Formation: Kf * Sub * Defect, Dissolution: Kr_pair * Pair
        term name = ReactBI add eqn = "$Kf_I * boronSub * Inter - $Kr_pair * boronInt" Silicon
        term name = ReactBV add eqn = "$Kf_V * boronSub * Vac - $Kr_pair * boronVac" Silicon

        # Pair diffusivities: Table A-3 (neutral + positive charge state)
        term name = DiffBI add eqn = "$DIX_BI + $DIP_BI * Poni" Silicon
        term name = DiffBV add eqn = "$DVX_BV + $DVP_BV * Poni" Silicon
    }
    phosphorus {
        # Substitutional phosphorus: active minus paired
        term name = phosphorusSub add eqn = "phosphorusActive - phosphorusInt" Silicon

        # Reaction: formation/dissolution
        term name = ReactPI add eqn = "$Kf_I * phosphorusSub * Inter - $Kr_pair * phosphorusInt" Silicon

        # Pair diffusivity: neutral + negative + double-negative
        term name = DiffPI add eqn = "$DIX_PI + $DIM_PI * Noni + $DIMM_PI * Noni * Noni" Silicon
    }
    arsenic {
        # Substitutional arsenic: active minus paired
        term name = arsenicSub add eqn = "arsenicActive - arsenicInt - arsenicVac" Silicon

        # Reaction: formation/dissolution
        term name = ReactAsI add eqn = "$Kf_I * arsenicSub * Inter - $Kr_pair * arsenicInt" Silicon
        term name = ReactAsV add eqn = "$Kf_V * arsenicSub * Vac - $Kr_pair * arsenicVac" Silicon

        # Pair diffusivities: neutral + negative charge state
        term name = DiffAsI add eqn = "$DIX_AsI + $DIM_AsI * Noni" Silicon
        term name = DiffAsV add eqn = "$DVX_AsV + $DVM_AsV * Noni" Silicon
    }
}

# {311} clustering: TSUPREM-4 Eq 3-319 (PDF p.114) with CF=1 approximation
# R_clust > 0: net clustering (I consumed, CIc grows)
# R_clust < 0: net dissolution (I released, CIc shrinks)
term name = R_clust add eqn = "$Kfc_311 * ScaleInter * (CIc + Inter) - $Kr_311 * CIc" Silicon

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
            "ddt(Inter) - $DiffI*EqInter*grad(ScaleInter) + $kIV*(Inter*Vac - EqInter*EqVac) + ReactBI + R_clust"
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
            "ddt(Inter) - $DiffI*EqInter*grad(ScaleInter) + ReactPI + R_clust"

        # Total phosphorus: no direct flux
        pdbSetString Silicon phosphorus Equation \
            "ddt(phosphorus) + ReactPI"

        # PI pairs: diffusion + reaction (donor: use Noni)
        pdbSetString Silicon phosphorusInt Equation \
            "ddt(phosphorusInt) - DiffPI * grad(phosphorusInt * Noni) - ReactPI"
    }
    arsenic {
        pdbSetString Silicon Inter Equation \
            "ddt(Inter) - $DiffI*EqInter*grad(ScaleInter) + $kIV*(Inter*Vac - EqInter*EqVac) + ReactAsI + R_clust"
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

# {311} cluster equation: immobile, grows/shrinks via R_clust (Eq 3-319)
pdbSetString Silicon CIc Equation "ddt(CIc) - R_clust"

# =============================================================================
# 19. SURFACE RECOMBINATION BCS
# =============================================================================

set EqI_surf "($cis * $InumI / $IdenI)"
pdbSetString Oxide_Silicon Inter Silicon Equation "-$KsurfI_val*(Inter(Silicon) - $EqI_surf)"

# Vac surface recombination only for dopants that solve Vac
switch ${dopant} {
    boron -
    arsenic {
        set EqV_surf "($cvs * $InumV / $IdenV)"
        pdbSetString Oxide_Silicon Vac Silicon Equation "-$KsurfV_val*(Vac(Silicon) - $EqV_surf)"
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

puts "=== POST-ANNEAL CLUSTERED-I ==="
sel z=CIc
foreach v [print.1d] {
    lassign $v d c
    if {$c == "Value"} continue
    if {$d < 0} continue
    if {$d > 0.5} break
    puts "$d $c"
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
