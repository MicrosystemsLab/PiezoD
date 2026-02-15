# Boron TED with I-V recombination

math diffuse dim=1 umf none col scale

line x loc=-0.001 spac=0.001 tag=TopOx
line x loc=0.00 spac=0.002 tag=TopSi
line x loc=0.15 spac=0.002
line x loc=0.50 spac=0.01
line x loc=2.00 spac=0.02 tag=Bottom

region Oxide xlo=TopOx xhi=TopSi
region Silicon xlo=TopSi xhi=Bottom
init

implant boron dose=2e15 energy=20 tilt=7

# Interstitial damage (+1 model)
sel z=Boron name=Inter

options !constdatafields storenodes
solution name=Temp !negative add const val=1000

set kT [expr {8.612e-5*(1000+273.0)}]

# Equilibrium concentrations
set CIStar [expr {5.0e22*exp(11.2)*exp(-3.7/$kT)}]
set CVStar [expr {1.0e23*exp(7.0)*exp(-3.6/$kT)}]
puts "CIStar: $CIStar cm^-3"
puts "CVStar: $CVStar cm^-3"

solution name=Inter solve !negative add
solution name=Vac solve !negative add

# Initialize vacancies at equilibrium
sel z=$CVStar name=Vac store

# Diffusivities
set DiffI [expr {0.138*exp(-1.37/$kT)}]
set DiffV [expr {0.001*exp(-0.9/$kT)}]
set DiffB [expr {0.76*exp(-3.46/$kT)}]
puts "DiffI: $DiffI"
puts "DiffV: $DiffV"
puts "DiffB: $DiffB"

# I-V recombination rate (diffusion-limited)
set lattice 2.35e-8
set kIV [expr {4 * 3.14159 * ($DiffI + $DiffV) * $lattice}]
set IVprod [expr {$CIStar * $CVStar}]
puts "kIV: $kIV"
puts "I*V*: $IVprod"

# Interstitial equation with bulk recombination
# dI/dt = D_I * grad(I) - k*(I*V - I*V*)
pdbSetString Silicon Inter Equation "ddt(Inter) - $DiffI*grad(Inter) + $kIV*(Inter*Vac - $IVprod)"

# Vacancy equation with bulk recombination
pdbSetString Silicon Vac Equation "ddt(Vac) - $DiffV*grad(Vac) + $kIV*(Inter*Vac - $IVprod)"

# Surface recombination for interstitials
set Ksurf [expr {3.14159 * $DiffI * $lattice * 1.3e15}]
pdbSetString Oxide_Silicon Inter Silicon Equation "-$Ksurf*(Inter(Silicon)-$CIStar)"

# Boron with TED
solution name=Boron solve !negative
pdbSetString Silicon Boron Equation "ddt(Boron) - $DiffB * (Inter/$CIStar) * grad(Boron)"

puts "=== DIFFUSING 30 min ==="
diffuse time=30 temp=1000

puts "=== POST-ANNEAL BORON ==="
sel z=Boron
foreach v [print.1d] {
    lassign $v d c
    if {$c == "Value"} continue
    if {$d < 0} continue
    if {$d > 1.5} break
    puts "$d $c"
}

puts "=== POST-ANNEAL INTERSTITIALS ==="
sel z=Inter
foreach v [print.1d] {
    lassign $v d c
    if {$c == "Value"} continue
    if {$d < 0} continue
    if {$d > 0.3} break
    puts "$d $c"
}

exit
