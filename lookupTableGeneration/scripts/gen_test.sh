#!/bin/bash
# Generate a test script from a template
# Usage: bash gen_test.sh <dopant> <dose> <energy> <temp> <time> [template] > output.tcl
TEMPLATE="${6:-../templates/ion_implant_fermi.tcl}"
sed -e "s/\${dopant}/$1/g" \
    -e "s/\${dose}/$2/g" \
    -e "s/\${energy}/$3/g" \
    -e "s/\${temp}/$4/g" \
    -e "s/\${time}/$5/g" \
    "$TEMPLATE"
