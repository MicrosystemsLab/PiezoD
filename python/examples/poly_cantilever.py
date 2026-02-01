"""Polycrystalline thin-film cantilever example.

Demonstrates modeling of multi-layer polycrystalline cantilevers with different
piezoresistor materials (poly-Si, Ti, Al). These cantilevers use thin-film
deposition processes rather than bulk silicon micromachining.

The multi-layer beam mechanics account for different elastic moduli and densities
in each layer, calculating the neutral axis position and effective stiffness
using the parallel axis theorem.
"""

import addcopyfighandler  # noqa: F401
import matplotlib.pyplot as plt
import numpy as np

plt.style.use("piezod")

from piezod import CantileverPoly, Material

# =============================================================================
# Example 1: Polysilicon piezoresistor cantilever
# =============================================================================

# Create a 3-layer cantilever: poly-Si / Si / poly-Si
# This symmetric structure places piezoresistors on top and bottom
c_poly = CantileverPoly(
    freq_min=1,
    freq_max=1000,
    l=200e-6,  # 200 um length
    w=20e-6,  # 20 um width
    t_top=200e-9,  # 200 nm poly-Si piezoresistor
    t_mid=2e-6,  # 2 um Si structural layer
    t_bot=200e-9,  # 200 nm poly-Si (symmetric)
    matl_top=Material.POLY,
    matl_mid=Material.SI,
    matl_bot=Material.POLY,
    l_pr_ratio=0.3,  # Piezoresistor covers 30% of length
    v_bridge=5.0,  # 5V bridge bias
    dopant_concentration=1e19,  # 1e19 cm^-3 phosphorus
    number_of_piezoresistors=2,
    number_of_piezoresistors_on_cantilever=2,
)

print("=" * 70)
print("Polycrystalline Thin-Film Cantilever Examples")
print("=" * 70)
print()

print("=== Example 1: Poly-Si Piezoresistor ===")
print(f"Structure: {c_poly.matl_top.value} / {c_poly.matl_mid.value} / {c_poly.matl_bot.value}")
print(f"Thicknesses: {c_poly.t_top * 1e9:.0f} / {c_poly.t_mid * 1e6:.1f}um / {c_poly.t_bot * 1e9:.0f} nm")
print()
print(f"Length: {c_poly.l * 1e6:.0f} um")
print(f"Width: {c_poly.w * 1e6:.0f} um")
print(f"Total thickness: {c_poly.t_total * 1e6:.2f} um")
print()
print(f"Neutral axis: {c_poly.neutral_axis() * 1e6:.3f} um from bottom")
print(f"Stiffness: {c_poly.stiffness():.3f} N/m")
print(f"Resonant frequency: {c_poly.resonant_frequency() / 1e3:.1f} kHz")
print()
print(f"Dopant concentration: {c_poly.dopant_concentration:.0e} cm^-3")
print(f"Sheet resistance: {c_poly.sheet_resistance():.1f} ohm/sq")
print(f"Resistance: {c_poly.resistance():.0f} ohm")
print()
print(f"Force sensitivity: {c_poly.force_sensitivity():.0f} V/N")
print(f"Force resolution: {c_poly.force_resolution() * 1e12:.1f} pN")
print()

# =============================================================================
# Example 2: Metal piezoresistors (Ti and Al)
# =============================================================================

# Titanium piezoresistor on oxide/Si structure
c_ti = CantileverPoly(
    freq_min=1,
    freq_max=1000,
    l=200e-6,
    w=20e-6,
    t_top=100e-9,  # 100 nm Ti piezoresistor
    t_mid=500e-9,  # 500 nm oxide insulator
    t_bot=2e-6,  # 2 um Si structural layer
    matl_top=Material.TI,
    matl_mid=Material.OXIDE,
    matl_bot=Material.SI,
    l_pr_ratio=0.3,
    v_bridge=2.0,  # Lower voltage for metal piezoresistors
    number_of_piezoresistors=2,
    number_of_piezoresistors_on_cantilever=1,  # Single-sided
)

# Aluminum piezoresistor
c_al = CantileverPoly(
    freq_min=1,
    freq_max=1000,
    l=200e-6,
    w=20e-6,
    t_top=100e-9,  # 100 nm Al piezoresistor
    t_mid=500e-9,  # 500 nm oxide insulator
    t_bot=2e-6,  # 2 um Si structural layer
    matl_top=Material.AL,
    matl_mid=Material.OXIDE,
    matl_bot=Material.SI,
    l_pr_ratio=0.3,
    v_bridge=2.0,
    number_of_piezoresistors=2,
    number_of_piezoresistors_on_cantilever=1,
)

print("=== Example 2: Metal Piezoresistors ===")
print(f"{'Parameter':<25} {'Ti':>15} {'Al':>15}")
print("-" * 57)
print(f"{'Sheet resistance (ohm/sq)':<25} {c_ti.sheet_resistance():>15.1f} {c_al.sheet_resistance():>15.1f}")
print(f"{'Resistance (ohm)':<25} {c_ti.resistance():>15.0f} {c_al.resistance():>15.0f}")
print(f"{'Piezo coefficient (1/Pa)':<25} {c_ti.piezo_coefficient():>15.2e} {c_al.piezo_coefficient():>15.2e}")
print(f"{'Force sensitivity (V/N)':<25} {c_ti.force_sensitivity():>15.1f} {c_al.force_sensitivity():>15.1f}")
print(f"{'Force resolution (pN)':<25} {c_ti.force_resolution() * 1e12:>15.1f} {c_al.force_resolution() * 1e12:>15.1f}")
print()

# =============================================================================
# Example 3: Effect of doping concentration on poly-Si performance
# =============================================================================

print("=== Example 3: Doping Concentration Sweep ===")
concentrations = [1e17, 1e18, 1e19, 1e20]

print(f"{'N (cm^-3)':<12} {'Rs (ohm/sq)':<15} {'Sens (V/N)':<15} {'Fmin (pN)':<15}")
print("-" * 57)

for N in concentrations:
    c = CantileverPoly(
        freq_min=1,
        freq_max=1000,
        l=200e-6,
        w=20e-6,
        t_top=200e-9,
        t_mid=2e-6,
        t_bot=200e-9,
        matl_top=Material.POLY,
        matl_mid=Material.SI,
        matl_bot=Material.POLY,
        l_pr_ratio=0.3,
        v_bridge=5.0,
        dopant_concentration=N,
    )
    print(
        f"{N:<12.0e} {c.sheet_resistance():<15.1f} "
        f"{c.force_sensitivity():<15.0f} {c.force_resolution() * 1e12:<15.1f}"
    )
print()

# =============================================================================
# Plot: Noise spectrum comparison
# =============================================================================

fig, ax = plt.subplots(figsize=(10, 6))

freq = np.logspace(0, 3, 500)  # 1 Hz to 1 kHz

# Calculate noise components for poly-Si cantilever
johnson_psd = c_poly.johnson_PSD(freq)
hooge_psd = c_poly.hooge_PSD(freq)
thermo_psd = c_poly.thermo_PSD(freq)
amp_psd = c_poly.amplifier_PSD(freq)
total_psd = johnson_psd + hooge_psd + thermo_psd + amp_psd

ax.loglog(freq, np.sqrt(johnson_psd) * 1e9, "b-", label="Johnson")
ax.loglog(freq, np.sqrt(hooge_psd) * 1e9, "r-", label="1/f (Hooge)")
ax.loglog(freq, np.sqrt(thermo_psd) * 1e9, "g-", label="Thermomechanical")
ax.loglog(freq, np.sqrt(amp_psd) * 1e9, "m-", label="Amplifier")
ax.loglog(freq, np.sqrt(total_psd) * 1e9, "k--", label="Total")

ax.axvline(x=c_poly.knee_frequency(), color="gray", linestyle=":", alpha=0.7)
ax.text(
    c_poly.knee_frequency() * 1.1,
    np.sqrt(johnson_psd[0]) * 1e9 * 0.5,
    f"Knee: {c_poly.knee_frequency():.0f} Hz",
)

ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Voltage Noise (nV/sqrt(Hz))")
ax.set_title("Poly-Si Cantilever Noise Spectrum")
ax.legend(loc="upper right")
ax.set_xlim([1, 1000])

plt.tight_layout()
plt.savefig("poly_cantilever_noise.png")
print("Saved noise spectrum to: poly_cantilever_noise.png")
