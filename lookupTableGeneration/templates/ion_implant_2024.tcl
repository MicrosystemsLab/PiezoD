# Ion implantation and diffusion simulation for PiezoD lookup tables
# FLOOXS_2024 template - parameters substituted by runner script
#
# Parameters:
#   ${dopant}      - boron, phosphorus, or arsenic
#   ${dose}        - implant dose (cm^-2)
#   ${energy}      - implant energy (keV)
#   ${temp}        - anneal temperature (C)
#   ${time}        - anneal time (minutes)
#
# Temperature ramp matches TSUPREM-4:
#   - 45 min ramp from 800C to ${temp}
#   - ${time} min dwell at ${temp}
#   - 45 min ramp from ${temp} to 800C

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

# Set up solutions (FLOOXS_2024 syntax)
solution name=Temp !negative add const val=800

# Temperature-dependent parameters (Arrhenius)
# Vt = kT in eV
set Vt "(8.612e-5*(Temp+273.0))"

solution name=Inter solve !negative add

# Equilibrium interstitial concentration
# Source: Params/Silicon/Interstitial
#   Cstar = Arrhenius(5.0e22*exp(11.2), 3.7)
set cis "(5.0e22*exp(11.2)*exp(-3.7/$Vt))"
solution name=CIStar add const val = "$cis"

# TED enhancement factor (required as intermediate solution in 2024)
solution name=TED const silicon val = (Inter/CIStar)

solution name=${dopant} solve !negative

# Interstitial diffusivity
# Source: Params/Silicon/Interstitial
#   D0 = Arrhenius(0.138, 1.37)
set DiffI "0.138*exp(-1.37/$Vt)"

# Interstitial diffusion equation
pdbSetString Silicon Inter Equation "ddt(Inter) - CIStar * $DiffI * grad(Inter/CIStar)"

# Surface recombination at oxide/silicon interface
# Source: Test/Segregation/SurfRecomb.tcl
#   Ksurf = pi * DiffI * lattice_spacing * kink_site_density
#   lattice_spacing = 2.714e-8 cm
#   kink_site_density = 1.3e15 cm^-2
set Ksurf "(3.14159 * $DiffI * 2.714e-8 * 1.3e15)"
pdbSetString Oxide_Silicon Inter Silicon Equation "$Ksurf*(Inter(Silicon)-CIStar)"

# Dopant diffusion with TED enhancement
# Source: Params/Silicon/Boron/Interstitial
#   D0 = Arrhenius(0.743, 3.56)
#   Dp = Arrhenius(0.617, 3.56)
#   Total interstitial-mediated: (D0 + Dp) = 1.36, Ea = 3.56 eV
# Note: Vacancy-mediated diffusion not included (would add 0.34 to prefactor)
set DiffDopant "1.36*exp(-3.56/$Vt)"
pdbSetString Silicon ${dopant} Equation "ddt(${dopant}) - ($DiffDopant) * TED * grad(${dopant})"

# Pre-anneal profile
puts "=== PRE-ANNEAL ==="
sel z=${dopant}
foreach v [print.1d] {
    lassign $v d c
    if {$c == "Value"} continue
    if {$d < 0} continue
    puts "$d $c"
}

# Calculate ramp rate: (temp - 800) / (45 min * 60 s/min) in C/s
set ramprate [expr {(${temp} - 800.0) / 2700.0}]

# TSUPREM-4 anneal sequence
puts "=== RAMP UP: 45 min, 800C -> ${temp}C ==="
diffuse time=45 temp=800 ramprate=$ramprate

puts "=== DWELL: ${time} min at ${temp}C ==="
diffuse time=${time} temp=${temp}

puts "=== RAMP DOWN: 45 min, ${temp}C -> 800C ==="
diffuse time=45 temp=${temp} ramprate=-$ramprate

# Post-anneal profile
puts "=== POST-ANNEAL ==="
sel z=${dopant}
foreach v [print.1d] {
    lassign $v d c
    if {$c == "Value"} continue
    if {$d < 0} continue
    puts "$d $c"
}

# Junction depth at 1e17 cm^-3
sel z=${dopant}-1e17
set xj [interpolate silicon val=0.0]
puts "=== RESULTS ==="
puts "Xj: $xj um"

exit
