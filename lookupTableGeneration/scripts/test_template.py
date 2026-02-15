"""Test ion_implant.tcl template with reference case parameters."""

import string
import subprocess
from pathlib import Path

# Reference case: Boron, 2e15 cm^-2, 20 keV, 1000C, 30 min
params = {
    "dopant": "boron",
    "dose": "2e15",
    "energy": "20",
    "temp": "1000",
    "time": "30",
}

# Load template
template_path = Path(__file__).parent.parent / "templates" / "ion_implant.tcl"
template_string = template_path.read_text()

# Substitute parameters using string.Template (same as simulationControl.py)
# Use safe_substitute to leave Tcl variables like $Vt unchanged
s = string.Template(template_string)
file_contents = s.safe_substitute(**params)

# Save to simulations directory
output_path = Path(__file__).parent.parent / "simulations" / "test_template.tcl"
output_path.write_text(file_contents)
print(f"Created {output_path}")

# Run with Docker
cmd = "docker compose run --rm flooxs test_template.tcl"
result = subprocess.run(
    cmd,
    cwd=Path(__file__).parent.parent,
    shell=True,
    capture_output=True,
    text=True,
)

# Show return code
print(f"Return code: {result.returncode}")

# Extract results from output
output = result.stdout + result.stderr
for line in output.split("\n"):
    if "Xj:" in line or "RESULTS" in line or "Target" in line:
        print(line)

if not any("Xj:" in line for line in output.split("\n")):
    print("No Xj found in output. Last 20 lines:")
    for line in output.split("\n")[-20:]:
        print(line)

# Clean up
output_path.unlink()
print("Cleaned up test file")
