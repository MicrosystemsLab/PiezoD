"""Quickstart example for piezod.

Demonstrates basic cantilever modeling based on Harley's 1999 epitaxial cantilever.
Reference: Harley & Kenny, "High-sensitivity piezoresistive cantilevers" (1999)
"""

from piezod import CantileverEpitaxy

# Create an epitaxial cantilever similar to Harley's 89 nm device
c = CantileverEpitaxy()

# Set geometry (matching Harley 1999)
c.l = 300e-6  # length: 300 um
c.w = 44e-6  # width: 44 um
c.t = 89e-9  # thickness: 89 nm
c.l_pr_ratio = 45 / 300  # piezoresistor covers 45 um of 300 um length

# Set electrical parameters
c.v_bridge = 5  # Wheatstone bridge bias: 5V
c.doping_type = "boron"
c.number_of_piezoresistors = 1

# Set operating conditions
c.freq_min = 10  # Hz
c.freq_max = 1000  # Hz
c.fluid = "vacuum"

# Print cantilever properties
print("=== Cantilever Geometry ===")
print(f"Length: {c.l * 1e6:.0f} um")
print(f"Width: {c.w * 1e6:.0f} um")
print(f"Thickness: {c.t * 1e9:.0f} nm")
print(f"Piezoresistor length: {c.l_pr() * 1e6:.0f} um")
print()

print("=== Mechanical Properties ===")
print(f"Stiffness: {c.stiffness() * 1e3:.3f} mN/m")
print(f"Resonant frequency (vacuum): {c.omega_vacuum_hz() / 1e3:.1f} kHz")
print(f"Elastic modulus: {c.modulus() / 1e9:.0f} GPa")
print()

# Compare vacuum vs fluid damping
print("=== Fluid Damping Effects ===")
for fluid in ["vacuum", "air", "water"]:
    c.fluid = fluid
    freq_hz, Q = c.omega_damped_hz_and_Q()
    print(f"{fluid:8s}: f = {freq_hz / 1e3:6.1f} kHz, Q = {Q:6.1f}")
