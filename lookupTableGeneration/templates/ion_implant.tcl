# Ion implantation and diffusion simulation for PiezoD lookup tables
# FLOOXS template - parameters substituted by runner script
#
# Parameters:
#   ${dopant}      - boron, phosphorus, or arsenic
#   ${dose}        - implant dose (cm^-2)
#   ${energy}      - implant energy (keV)
#   ${temp}        - anneal temperature (C)
#   ${time}        - anneal time (minutes)
#   ${outputFile}  - output filename

# 1D mesh setup (5um depth, fine mesh near surface)
line x loc=0.0 spac=0.001 tag=Top
line x loc=0.5 spac=0.005
line x loc=5.0 spac=0.05 tag=Bottom

region Silicon xlo=Top xhi=Bottom
init

# Ion implantation
implant ${dopant} dose=${dose} energy=${energy} tilt=7

# TODO: Verify diffusion setup - may need pdbSetString for dopant equations
# diffuse time=${time} temp=${temp}

# Output dopant profile
sel z=${dopant}
puts [layers]

# TODO: Add file output for post-processing
# print.1d equivalent needed
