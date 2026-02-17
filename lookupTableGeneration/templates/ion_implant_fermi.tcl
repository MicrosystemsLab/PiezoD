# Ion implantation and diffusion simulation for PiezoD lookup tables
# FLOOXS Fermi model: Fermi-level dependent diffusivity, equilibrium defects
#
# Parameters:
#   ${dopant}      - boron, phosphorus, or arsenic
#   ${dose}        - implant dose (cm^-2)
#   ${energy}      - implant energy (keV)
#   ${temp}        - anneal temperature (C)
#   ${time}        - anneal time (minutes)
#
# Physics:
#   - Fermi-level dependent dopant diffusivity via charge neutrality
#   - Equilibrium defects (no TED, no explicit defect transport)
#   - Combined interstitial + vacancy diffusion paths
#   - Background doping (~1.4e15 cm^-3, 10 ohm-cm)
#   - 250A protection oxide matching TSUPREM-4
#   - Temperature ramp matching TSUPREM-4
#
# All diffusivity parameters from FLOOXS_2026/Params/Silicon/

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

# Ion implantation
implant ${dopant} dose=${dose} energy=${energy} tilt=7

# Background doping: ~1.4e15 cm^-3 (10 ohm-cm)
# Boron implant: n-type substrate (phosphorus background)
# P/As implant: p-type substrate (boron background)
switch ${dopant} {
    boron {
        sel z=1.4e15 name=Phosphorus
    }
    phosphorus {
        sel z=${dopant}+1.4e15 name=${dopant}
        sel z=1.4e15 name=Boron
    }
    arsenic {
        sel z=1.4e15 name=Boron
    }
}

options !constdatafields storenodes
solution name=Temp !negative add const val=800

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
# DOPANT DIFFUSION PARAMETERS (Fermi-level dependent)
# Source: FLOOXS_2026/Params/Silicon/{Dopant}/Interstitial, Vacancy
# =============================================================================
#
# Combined I+V paths with equilibrium defects (ScaleI = ScaleV = 1):
#   Acceptors (boron): D = (D0_I+D0_V) + (Dp_I+Dp_V)*Poni
#   Donors (P):        D = D0_I + Dn_I*Noni + Dnn_I*Noni^2
#   Donors (As):       D = D0_I + Dn_V*Noni

switch ${dopant} {
    boron {
        # B-I: D0=Arr(0.743,3.56), Dp=Arr(0.617,3.56)
        # B-V: D0=Arr(0.186,3.56), Dp=Arr(0.154,3.56)
        # Combined: D0=Arr(0.929,3.56), Dp=Arr(0.771,3.56)
        pdbSetDouble Silicon ${dopant} Diff_D0 {[Arrhenius 0.929 3.56]}
        pdbSetDouble Silicon ${dopant} Diff_Dp {[Arrhenius 0.771 3.56]}
    }
    phosphorus {
        # P-I: D0=Arr(5.6,3.71), Dn=Arr(6.38,4.05), Dnn=Arr(0.0245,3.23)
        pdbSetDouble Silicon ${dopant} Diff_D0 {[Arrhenius 5.6 3.71]}
        pdbSetDouble Silicon ${dopant} Diff_Dn {[Arrhenius 6.38 4.05]}
        pdbSetDouble Silicon ${dopant} Diff_Dnn {[Arrhenius 2.45e-2 3.23]}
    }
    arsenic {
        # As-I: D0=Arr(0.0666,3.45); As-V: Dn=Arr(12.8,4.05)
        pdbSetDouble Silicon ${dopant} Diff_D0 {[Arrhenius 0.0666 3.45]}
        pdbSetDouble Silicon ${dopant} Diff_Dn {[Arrhenius 12.8 4.05]}
    }
}

# =============================================================================
# SOLUTION DEFINITIONS
# =============================================================================

solution name=${dopant} solve !negative

# Fermi level: Noni = n/ni, Poni = p/ni from charge neutrality
solution name = Charge add const val = 0.0 Silicon
solution name = Noni add const val = "0.5*(Charge + sqrt(Charge*Charge + 4.0*$ni*$ni)) / $ni" Silicon
solution name = Poni add const val = "0.5*(-Charge + sqrt(Charge*Charge + 4.0*$ni*$ni)) / $ni" Silicon

# =============================================================================
# CHARGE NEUTRALITY
# =============================================================================

switch ${dopant} {
    boron {
        # Acceptor on n-type: Charge = N_D - N_A
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
# DIFFUSION EQUATION
# =============================================================================
#
# Flux = (D_total / chg) * grad(C * chg)
#   chg = Poni (acceptors) or Noni (donors)
# This form gives both concentration gradient and built-in field contributions.

switch ${dopant} {
    boron {
        set D0 [pdbDelayDouble Silicon ${dopant} Diff_D0]
        set Dp [pdbDelayDouble Silicon ${dopant} Diff_Dp]
        term name = DiffDop add eqn = "($D0 + $Dp * Poni) / Poni" Silicon
        pdbSetString Silicon ${dopant} Equation \
            "ddt(${dopant}) - DiffDop * grad(${dopant} * Poni)"
    }
    phosphorus {
        set D0 [pdbDelayDouble Silicon ${dopant} Diff_D0]
        set Dn [pdbDelayDouble Silicon ${dopant} Diff_Dn]
        set Dnn [pdbDelayDouble Silicon ${dopant} Diff_Dnn]
        term name = DiffDop add eqn = "($D0 + $Dn * Noni + $Dnn * Noni * Noni) / Noni" Silicon
        pdbSetString Silicon ${dopant} Equation \
            "ddt(${dopant}) - DiffDop * grad(${dopant} * Noni)"
    }
    arsenic {
        set D0 [pdbDelayDouble Silicon ${dopant} Diff_D0]
        set Dn [pdbDelayDouble Silicon ${dopant} Diff_Dn]
        term name = DiffDop add eqn = "($D0 + $Dn * Noni) / Noni" Silicon
        pdbSetString Silicon ${dopant} Equation \
            "ddt(${dopant}) - DiffDop * grad(${dopant} * Noni)"
    }
}

# =============================================================================
# PRE-ANNEAL PROFILE
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

# Ramp rate: (temp - 800) / (45 min * 60 s/min) in C/s
set ramprate [expr {(${temp} - 800.0) / 2700.0}]

puts "=== RAMP UP: 45 min, 800C -> ${temp}C ==="
diffuse time=45 temp=800 ramprate=$ramprate init=1e-6 damp.trbdf

puts "=== DWELL: ${time} min at ${temp}C ==="
diffuse time=${time} temp=${temp} init=1e-6 damp.trbdf

puts "=== RAMP DOWN: 45 min, ${temp}C -> 800C ==="
diffuse time=45 temp=${temp} ramprate=-$ramprate init=1e-6 damp.trbdf

# =============================================================================
# POST-ANNEAL PROFILE
# =============================================================================

puts "=== POST-ANNEAL ==="
sel z=${dopant}
foreach v [print.1d] {
    lassign $v d c
    if {$c == "Value"} continue
    if {$d < 0} continue
    puts "$d $c"
}

# Junction depth where implanted dopant equals background (1.4e15 cm^-3)
sel z=${dopant}-1.4e15
set xj [interpolate silicon val=0.0]
puts "=== RESULTS ==="
puts "Xj: $xj um"

exit
