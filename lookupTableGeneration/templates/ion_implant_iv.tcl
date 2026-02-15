# Ion implantation and diffusion simulation for PiezoD lookup tables
# FLOOXS template with I-V recombination - parameters substituted by runner script
#
# Parameters:
#   ${dopant}      - boron, phosphorus, or arsenic
#   ${dose}        - implant dose (cm^-2)
#   ${energy}      - implant energy (keV)
#   ${temp}        - anneal temperature (C)
#   ${time}        - anneal time (minutes)
#
# Physics:
#   - Interstitial and vacancy point defect dynamics
#   - I-V bulk recombination (diffusion-limited)
#   - Surface recombination at oxide/silicon interface
#   - Dopant diffusion via both interstitial and vacancy mechanisms
#   - Temperature ramp matching TSUPREM-4
#
# All parameters from FLOOXS_2026/Params/Silicon/

math diffuse dim=1 sblock scale bcgs row

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

# Ion implantation
implant ${dopant} dose=${dose} energy=${energy} tilt=7

# Initialize interstitials from implant damage (+1 model)
sel z=${dopant} name=Inter

# Set up solutions
options !constdatafields storenodes
solution name=Temp !negative add const val=800

# Temperature-dependent parameters (Arrhenius)
# Vt = kT in eV
set Vt "(8.612e-5*(Temp+273.0))"

# =============================================================================
# POINT DEFECT PARAMETERS
# Source: FLOOXS_2026/Params/Silicon/Interstitial, Vacancy, Info
# =============================================================================

# Lattice spacing (Silicon/Info)
set lattice "2.714417617e-8"

# Equilibrium interstitial concentration (Silicon/Interstitial)
# Cstar = Arrhenius(5.0e22*exp(11.2), 3.7)
set cis "(5.0e22*exp(11.2)*exp(-3.7/$Vt))"

# Equilibrium vacancy concentration (Silicon/Vacancy)
# Cstar = Arrhenius(5.0e22*exp(9.0), 3.97)
set cvs "(5.0e22*exp(9.0)*exp(-3.97/$Vt))"

# Interstitial diffusivity (Silicon/Interstitial)
# D0 = Arrhenius(0.138, 1.37)
set DiffI "0.138*exp(-1.37/$Vt)"

# Vacancy diffusivity (Silicon/Vacancy)
# D0 = Arrhenius(1.18e-4, 0.1)
set DiffV "1.18e-4*exp(-0.1/$Vt)"

# I-V bulk recombination rate - diffusion limited (Silicon/Interstitial, Vacancy)
# Kbulk = 4 * pi * (D_I + D_V) * lattice_spacing
set kIV "(4.0 * 3.14159 * ($DiffI + $DiffV) * $lattice)"

# =============================================================================
# DOPANT DIFFUSION PARAMETERS
# Source: FLOOXS_2026/Params/Silicon/{Dopant}/Interstitial, Vacancy
# =============================================================================

# Dopant-specific diffusivities (interstitial and vacancy mediated)
# Each dopant has different D0, Dp (positive charge state) values
switch ${dopant} {
    boron {
        # Boron/Interstitial: D0=0.743, Dp=0.617, Ea=3.56 eV
        # Boron/Vacancy: D0=0.186, Dp=0.154, Ea=3.56 eV
        set DiffDopantI "(0.743+0.617)*exp(-3.56/$Vt)"
        set DiffDopantV "(0.186+0.154)*exp(-3.56/$Vt)"
    }
    phosphorus {
        # Phosphorus/Interstitial: D0=5.6, Ea=3.71 eV (neutral only)
        # No Phosphorus/Vacancy file - phosphorus is 100% interstitial-mediated
        set DiffDopantI "5.6*exp(-3.71/$Vt)"
        set DiffDopantV "0.0"
    }
    arsenic {
        # Arsenic/Interstitial: D0=0.0666, Ea=3.45 eV
        # Arsenic/Vacancy: Dn=12.8, Ea=4.05 eV (primarily via negative vacancy pairs)
        set DiffDopantI "0.0666*exp(-3.45/$Vt)"
        set DiffDopantV "12.8*exp(-4.05/$Vt)"
    }
}

# =============================================================================
# SURFACE RECOMBINATION PARAMETERS
# Source: FLOOXS_2026/Params/Oxide_Silicon/Interstitial, Vacancy
# =============================================================================

# Temperature-dependent kink site density for interstitials (Oxide_Silicon/Interstitial)
# KinkSite = Arrhenius(0.51, -2.63) * 1e15 / (1e15 + Arrhenius(0.51, -2.63))
# Note: negative activation energy means kink occupation increases with temperature
set KinkArrh "(0.51*exp(2.63/$Vt))"
set KinkSiteI "($KinkArrh * 1.0e15 / (1.0e15 + $KinkArrh))"

# Interstitial surface recombination: Ksurf = pi * D_I * lattice * KinkSite
set KsurfI "(3.14159 * $DiffI * $lattice * $KinkSiteI)"

# Vacancy kink site density is constant (Oxide_Silicon/Vacancy)
# KinkSite = 1.0e5 (much lower than interstitial)
set KinkSiteV "1.0e5"
set KsurfV "(3.14159 * $DiffV * $lattice * $KinkSiteV)"

# =============================================================================
# SOLUTION DEFINITIONS
# =============================================================================

solution name=Inter solve !negative add
solution name=Vac solve !negative add
solution name=${dopant} solve !negative

# Equilibrium concentration solutions (for normalization in equations)
solution name=CIStar add const val="$cis"
solution name=CVStar add const val="$cvs"

# Initialize vacancies from implant damage (Frenkel pair model)
# Each implanted atom creates one I-V pair, so initial V = initial I
# This is more physical than +1 model which assumes V=V* (equilibrium)
# The fast I-V recombination will leave net I that drives TED
sel z=${dopant} name=Vac

# =============================================================================
# DIFFUSION EQUATIONS
# =============================================================================

# Interstitial diffusion with I-V bulk recombination
# dI/dt = D_I * laplacian(I) - k_IV * (I*V - I*V*)
# Simplified form matching reference_iv_recomb.tcl
pdbSetString Silicon Inter Equation "ddt(Inter) - $DiffI*grad(Inter) - $kIV*(Inter*Vac - CIStar*CVStar)"

# Vacancy diffusion with I-V bulk recombination
# dV/dt = D_V * laplacian(V) - k_IV * (I*V - I*V*)
pdbSetString Silicon Vac Equation "ddt(Vac) - $DiffV*grad(Vac) - $kIV*(Inter*Vac - CIStar*CVStar)"

# Surface recombination at oxide/silicon interface
# Interstitials recombine toward equilibrium
pdbSetString Oxide_Silicon Inter Silicon Equation "$KsurfI*(Inter(Silicon)-CIStar)"
# Vacancies recombine toward equilibrium
pdbSetString Oxide_Silicon Vac Silicon Equation "$KsurfV*(Vac(Silicon)-CVStar)"

# Dopant diffusion with interstitial AND vacancy enhancement
# dC/dt = div(D_I * (I/I*) * grad(C)) + div(D_V * (V/V*) * grad(C))
# Combined: dC/dt = div((D_I*(I/I*) + D_V*(V/V*)) * grad(C))
pdbSetString Silicon ${dopant} Equation "ddt(${dopant}) - ($DiffDopantI) * (Inter/CIStar) * grad(${dopant}) - ($DiffDopantV) * (Vac/CVStar) * grad(${dopant})"

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

# Junction depth at 1e17 cm^-3
sel z=${dopant}-1e17
set xj [interpolate silicon val=0.0]
puts "=== RESULTS ==="
puts "Xj: $xj um"

exit
