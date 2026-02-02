# FLOOXS validation script
# Tests implant and diffusion, outputs concentration vs depth

# Configure 1D diffusion solver
math diffuse dim=1 umf none col scale

# Set up Boron as a solution variable
solution add name=Boron solve !negative

# 1D mesh (5um depth)
line x loc=0.0 spac=0.001 tag=Top
line x loc=0.1 spac=0.001
line x loc=0.5 spac=0.005
line x loc=2.0 spac=0.02
line x loc=5.0 spac=0.05 tag=Bottom

region Silicon xlo=Top xhi=Bottom
init

# Implant boron: 1e15 cm^-2 at 20 keV
implant boron dose=1e15 energy=20 tilt=7

# Output pre-anneal profile
puts "=== PRE-ANNEAL LAYERS ==="
sel z=Boron
puts [layers]

puts "=== PRE-ANNEAL PROFILE (TOTAL) ==="
puts "Depth(um) Concentration(cm^-3)"
sel z=Boron
foreach v [print.1d] {
    lassign $v dist val
    if { $val == "Value" } { continue }
    puts "$dist $val"
}

# Set up Boron diffusion equation
# D = D0 * exp(-Ea/kT), using Boron in Si values
# D0 ~ 0.76 cm^2/s, Ea ~ 3.46 eV for intrinsic diffusion
pdbSetString Si Boron Equation "ddt(Boron) - 0.76*exp(-3.46/(8.62e-5*(Temp+273))) * grad(Boron)"

# Diffuse: 30 min at 1000C
diffuse time=30 temp=1000

# Output post-anneal profiles
puts "=== POST-ANNEAL LAYERS ==="
sel z=Boron
puts [layers]

puts "=== POST-ANNEAL PROFILE (TOTAL) ==="
puts "Depth(um) Concentration(cm^-3)"
sel z=Boron
foreach v [print.1d] {
    lassign $v dist val
    if { $val == "Value" } { continue }
    puts "$dist $val"
}

# Define active concentration using solid solubility model
# Css = C0 * exp(-Ea/kT), typical Boron values at 1000C ~ 2e20 cm^-3
# Active = Css * Total / (Css + Total)
sel z=2e20*Boron/(2e20+Boron) name=BoronActive

puts "=== POST-ANNEAL PROFILE (ACTIVE) ==="
puts "Depth(um) Concentration(cm^-3)"
sel z=BoronActive
foreach v [print.1d] {
    lassign $v dist val
    if { $val == "Value" } { continue }
    puts "$dist $val"
}

puts "=== VALIDATION COMPLETE ==="
exit
