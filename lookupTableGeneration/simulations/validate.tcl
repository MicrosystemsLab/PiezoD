# FLOOXS validation script
# Tests implant and diffusion, outputs concentration vs depth

# Configure 1D diffusion solver
math diffuse dim=1 umf

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

# Set up Boron as a solution variable for diffusion
solution add name=Boron solve !negative

# Set ActiveModel to solid solubility (1 = Solid)
pdbSetSwitch Si Boron ActiveModel 1

# Define BoronActive based on solid solubility
# Solubility at 1000C ~ 2e20 cm^-3
set Css 2.0e20
term name=BoronActive add eqn="($Css) * Boron / (($Css) + Boron)" Silicon

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
