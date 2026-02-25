# FLOOXS built-in Fermi diffusion template
# Equation generated entirely by DopantBulk — zero custom equations
# Requires Dopant_patched.tcl (3 bug fixes for FLOOXS 2026 C++ pdb)
#
# Parameters: ${dopant}, ${dose}, ${energy}, ${temp}, ${time}

source /opt/flooxs/TclLib/Name.tcl
source /opt/flooxs/Params/paramFunc
source /work/Dopant_patched.tcl

# =============================================================================
# 1. MESH (silicon only, no oxide)
# =============================================================================

math diffuse dim=1 umf none col scale

line x loc=0.0 spac=0.0002 tag=Top
line x loc=0.05 spac=0.0002
line x loc=0.15 spac=0.001
line x loc=0.5 spac=0.005
line x loc=2.0 spac=0.02
line x loc=5.0 spac=0.05
line x loc=10.0 spac=0.1 tag=Bottom

region Silicon xlo=Top xhi=Bottom
init

# =============================================================================
# 2. BACKGROUND DOPING AND IMPLANT
# =============================================================================

switch ${dopant} {
    boron       { sel z=1.4e15 name=Phosphorus }
    phosphorus  { sel z=1.4e15 name=Boron }
    arsenic     { sel z=1.4e15 name=Boron }
}

implant ${dopant} dose=${dose} energy=${energy} tilt=7

# =============================================================================
# 3. SOLVER SETTINGS
# =============================================================================

options !constdatafields storenodes
solution name = Temp !negative add const val = ${temp}

pdbSetDouble Math iterLimit 500
pdbSetDouble Math updateLimit 1e-8
pdbSetDouble Math rhsLimit 1e-20
pdbSetDouble Math rhsMin 1e-2

solution name = ${dopant} solve !negative

# =============================================================================
# 4. FERMI LEVEL (infrastructure required by DopantFermi)
# =============================================================================

set ni "(3.87e16 * sqrt((Temp+273.0)*(Temp+273.0)*(Temp+273.0)) * exp(-7014.0/(Temp+273.0)))"

solution name = Charge add const val = 0.0 Silicon
solution name = Noni add const val = "0.5*(Charge + sqrt(Charge*Charge + 4.0*$ni*$ni)) / $ni" Silicon
solution name = Poni add const val = "0.5*(-Charge + sqrt(Charge*Charge + 4.0*$ni*$ni)) / $ni" Silicon

switch ${dopant} {
    boron       { term name = Charge add eqn = "1.4e15" Silicon }
    phosphorus  { term name = Charge add eqn = "-1.4e15" Silicon }
    arsenic     { term name = Charge add eqn = "-1.4e15" Silicon }
}

# =============================================================================
# 5. POPULATE C++ PDB (parameters only — DopantBulk generates the equations)
# =============================================================================
# D0/Dp/Dn/Dnn must be single Arrhenius expressions (operator precedence).
# For boron: D0 = BI_D0+BV_D0, Dp = BI_Dp+BV_Dp (same Ea=3.56, combinable)
# For P: I-only diffuser (single component per coeff)
# For As: I component dominates (V negligible), use I-only

switch ${dopant} {
    boron {
        pdbSet Silicon ${dopant} DiffModel 1
        pdbSet Silicon ${dopant} ActiveModel 1
        pdbSet Silicon ${dopant} Charge 1
        pdbSet Silicon ${dopant} D0 {[Arrhenius 0.929 3.56]}
        pdbSet Silicon ${dopant} Dp {[Arrhenius 0.771 3.56]}
        pdbSet Silicon ${dopant} Solubility {[Arrhenius 7.68e22 0.7086]}
    }
    phosphorus {
        pdbSet Silicon ${dopant} DiffModel 1
        pdbSet Silicon ${dopant} ActiveModel 1
        pdbSet Silicon ${dopant} Charge 2
        pdbSet Silicon ${dopant} D0 {[Arrhenius 3.85 3.66]}
        pdbSet Silicon ${dopant} Dn {[Arrhenius 44.2 4.0]}
        pdbSet Silicon ${dopant} Dnn {[Arrhenius 44.2 4.37]}
        pdbSet Silicon ${dopant} Solubility {[Arrhenius 2.45e23 0.518]}
    }
    arsenic {
        pdbSet Silicon ${dopant} DiffModel 1
        pdbSet Silicon ${dopant} ActiveModel 1
        pdbSet Silicon ${dopant} Charge 2
        pdbSet Silicon ${dopant} D0 {[Arrhenius 1.37e7 4.97]}
        pdbSet Silicon ${dopant} Dn {[Arrhenius 1.37e7 4.97]}
        pdbSet Silicon ${dopant} Solubility {[Arrhenius 7.14e22 1.24]}
    }
}

# =============================================================================
# 6. BUILT-IN MODEL (generates diffusion equation — zero custom equations)
# =============================================================================

DopantBulk Silicon ${dopant}

puts "Equation: [pdbGet Silicon ${dopant} Equation]"

# =============================================================================
# 7. DIFFUSE AND EXTRACT Xj
# =============================================================================

diffuse time=${time} temp=${temp} init=1e-6 damp.trbdf

sel z=${dopant}-1.4e15
set xj [interpolate silicon val=0.0]
puts "Xj: $xj um"

exit
