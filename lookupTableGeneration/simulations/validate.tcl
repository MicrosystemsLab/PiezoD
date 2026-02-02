# FLOOXS validation script
# Tests implant and diffusion, outputs concentration vs depth

# 1D mesh - fine near surface, coarser at depth
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

puts "=== PRE-ANNEAL PROFILE ==="
puts "Depth(um) Concentration(cm^-3)"
sel z=Boron
foreach v [print.1d] {
    lassign $v dist val
    if { $val == "Value" } { continue }
    puts "$dist $val"
}

# Set up Boron as PDE solution for diffusion
solution add name=Boron pde solve !negative

# Set diffusion equation: ddt(C) - D*grad(C) = 0
# D ~ 1e-14 cm^2/s for boron at 1000C
pdbSet Si Boron Equation "ddt(Boron) - 1e-14*grad(Boron)"

# Diffuse: 30 min = 1800 seconds
diffuse time=1800

# Output post-anneal profile
puts "=== POST-ANNEAL LAYERS ==="
sel z=Boron
puts [layers]

puts "=== POST-ANNEAL PROFILE ==="
puts "Depth(um) Concentration(cm^-3)"
sel z=Boron
foreach v [print.1d] {
    lassign $v dist val
    if { $val == "Value" } { continue }
    puts "$dist $val"
}

puts "=== VALIDATION COMPLETE ==="
exit
