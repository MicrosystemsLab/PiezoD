"""Phosphorus-diffused piezoresistive cantilever example.

Demonstrates modeling of POCl3-diffused piezoresistors, which are commonly used
in commercial MEMS fabrication. The phosphorus diffusion model includes the
characteristic "kink" profile due to electric field-enhanced diffusion.

Reference:
    Doll & Bhatti, "The Design and Optimization of Piezoresistive Cantilevers"
    Journal of Applied Physics (2009)
    https://github.com/jcdoll/jcdoll.github.io/blob/master/papers/2009_JAP_CantileverOptimization.pdf

The phosphorus diffusion model is based on:
    Tsai, "POCl3 Diffusion Model for VLSI Process Simulation" (1983)
"""

import addcopyfighandler  # noqa: F401
import matplotlib.pyplot as plt
import numpy as np

from piezod import CantileverDiffusion, CantileverEpitaxy

# Cantilever geometry - typical MEMS force sensor
L = 200e-6  # length: 200 um
W = 20e-6  # width: 20 um
T = 2e-6  # thickness: 2 um

# Piezoresistor covers 30% of cantilever length (U-shaped at base)
L_PR_RATIO = 0.3

# Operating frequency band for force sensing
FREQ_MIN = 1  # Hz
FREQ_MAX = 1000  # Hz

# Electrical parameters
V_BRIDGE = 5.0  # Wheatstone bridge bias (V)

# POCl3 diffusion process parameters
# Typical foundry process: 850C for 30 minutes
DIFFUSION_TEMP = 850 + 273.15  # Temperature in Kelvin
DIFFUSION_TIME = 30 * 60  # Time in seconds (30 minutes)

# Create diffusion cantilever
c_diff = CantileverDiffusion(
    freq_min=FREQ_MIN,
    freq_max=FREQ_MAX,
    l=L,
    w=W,
    t=T,
    l_pr_ratio=L_PR_RATIO,
    v_bridge=V_BRIDGE,
    doping_type="phosphorus",
    diffusion_time=DIFFUSION_TIME,
    diffusion_temp=DIFFUSION_TEMP,
)

# Set operating environment
c_diff.fluid = "air"
c_diff.number_of_piezoresistors = 4  # Full Wheatstone bridge

# Create equivalent epitaxial cantilever for comparison
# Use typical epitaxial doping parameters
c_epi = CantileverEpitaxy()
c_epi.l = L
c_epi.w = W
c_epi.t = T
c_epi.l_pr_ratio = L_PR_RATIO
c_epi.v_bridge = V_BRIDGE
c_epi.doping_type = "phosphorus"
c_epi.dopant_concentration = 1e19  # 1e19 cm^-3
c_epi.t_pr_ratio = 0.1  # 10% of thickness
c_epi.freq_min = FREQ_MIN
c_epi.freq_max = FREQ_MAX
c_epi.fluid = "air"
c_epi.number_of_piezoresistors = 4

print("=" * 70)
print("Phosphorus-Diffused vs Epitaxial Piezoresistive Cantilever Comparison")
print("=" * 70)
print()

# Cantilever geometry
print("=== Cantilever Geometry ===")
print(f"Length: {L * 1e6:.0f} um")
print(f"Width: {W * 1e6:.0f} um")
print(f"Thickness: {T * 1e6:.1f} um")
print(f"Piezoresistor length: {c_diff.l_pr() * 1e6:.0f} um")
print()

# Diffusion process parameters
print("=== Diffusion Process Parameters ===")
print(f"Temperature: {DIFFUSION_TEMP - 273.15:.0f} C")
print(f"Time: {DIFFUSION_TIME / 60:.0f} minutes")
print(f"Junction depth: {c_diff.junction_depth * 1e9:.0f} nm")
print()

# Compare electrical properties
print("=== Electrical Properties Comparison ===")
print(f"{'Parameter':<30} {'Diffusion':>15} {'Epitaxial':>15}")
print("-" * 62)

Rs_diff = c_diff.sheet_resistance()
Rs_epi = c_epi.sheet_resistance()
print(f"{'Sheet resistance (Ohm/sq)':<30} {Rs_diff:>15.1f} {Rs_epi:>15.1f}")

R_diff = c_diff.resistance()
R_epi = c_epi.resistance()
print(f"{'Total resistance (Ohm)':<30} {R_diff:>15.0f} {R_epi:>15.0f}")

Nz_diff = c_diff.Nz()
Nz_epi = c_epi.Nz()
print(f"{'Effective Nz (cm^-2)':<30} {Nz_diff:>15.2e} {Nz_epi:>15.2e}")

alpha_diff = c_diff.alpha()
alpha_epi = c_epi.alpha()
print(f"{'Hooge alpha':<30} {alpha_diff:>15.2e} {alpha_epi:>15.2e}")

print()

# Compare performance metrics
print("=== Performance Metrics ===")
print(f"{'Parameter':<30} {'Diffusion':>15} {'Epitaxial':>15}")
print("-" * 62)

k = c_diff.stiffness()
print(f"{'Stiffness (N/m)':<30} {k:>15.3f} {c_epi.stiffness():>15.3f}")

f0 = c_diff.omega_vacuum_hz()
print(f"{'Resonant freq (kHz)':<30} {f0/1e3:>15.1f} {c_epi.omega_vacuum_hz()/1e3:>15.1f}")

beta_diff = c_diff.beta()
beta_epi = c_epi.beta()
print(f"{'Beta (efficiency)':<30} {beta_diff:>15.3f} {beta_epi:>15.3f}")

sens_diff = c_diff.force_sensitivity()
sens_epi = c_epi.force_sensitivity()
print(f"{'Force sensitivity (V/N)':<30} {sens_diff:>15.0f} {sens_epi:>15.0f}")

res_diff = c_diff.force_resolution() * 1e12
res_epi = c_epi.force_resolution() * 1e12
print(f"{'Force resolution (pN)':<30} {res_diff:>15.1f} {res_epi:>15.1f}")

print()
print("=" * 70)

# Plot doping profiles
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Diffusion profile
z_diff, active_diff, total_diff = c_diff.doping_profile()
ax1 = axes[0]
ax1.semilogy(z_diff * 1e9, active_diff, "b-", linewidth=2, label="Active")
ax1.semilogy(z_diff * 1e9, total_diff, "r--", linewidth=2, label="Total")
ax1.axhline(y=1e15, color="gray", linestyle=":", label="Background (1e15)")
ax1.axvline(x=c_diff.junction_depth * 1e9, color="green", linestyle="--", label="Junction")
ax1.set_xlabel("Depth (nm)")
ax1.set_ylabel("Concentration (cm$^{-3}$)")
ax1.set_title(f"POCl3 Diffusion Profile\n({DIFFUSION_TEMP-273.15:.0f}C, {DIFFUSION_TIME/60:.0f} min)")
ax1.legend()
ax1.set_ylim([1e14, 1e22])
ax1.grid(True, alpha=0.3)

# Epitaxial profile (step function)
z_epi, active_epi, total_epi = c_epi.doping_profile()
ax2 = axes[1]
ax2.semilogy(z_epi * 1e9, active_epi, "b-", linewidth=2, label="Active = Total")
ax2.axhline(y=1e15, color="gray", linestyle=":", label="Background (1e15)")
ax2.axvline(x=c_epi.junction_depth * 1e9, color="green", linestyle="--", label="Junction")
ax2.set_xlabel("Depth (nm)")
ax2.set_ylabel("Concentration (cm$^{-3}$)")
ax2.set_title(f"Epitaxial Profile\n(N = {c_epi.dopant_concentration:.0e} cm$^{{-3}}$)")
ax2.legend()
ax2.set_ylim([1e14, 1e22])
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("doping_profiles_comparison.png", dpi=150)
print("Saved doping profile comparison to: doping_profiles_comparison.png")

# Plot temperature dependence of diffusion
fig2, ax3 = plt.subplots(figsize=(8, 5))

temps = [800, 850, 900, 950, 1000]  # Celsius
colors = plt.cm.hot(np.linspace(0.3, 0.8, len(temps)))

for temp_c, color in zip(temps, colors):
    c_temp = CantileverDiffusion(
        freq_min=FREQ_MIN,
        freq_max=FREQ_MAX,
        l=L,
        w=W,
        t=T,
        l_pr_ratio=L_PR_RATIO,
        v_bridge=V_BRIDGE,
        doping_type="phosphorus",
        diffusion_time=DIFFUSION_TIME,
        diffusion_temp=temp_c + 273.15,
    )
    z, active, total = c_temp.doping_profile()
    ax3.semilogy(z * 1e9, active, color=color, linewidth=1.5, label=f"{temp_c}C")

ax3.axhline(y=1e15, color="gray", linestyle=":", alpha=0.5)
ax3.set_xlabel("Depth (nm)")
ax3.set_ylabel("Active Concentration (cm$^{-3}$)")
ax3.set_title("POCl3 Diffusion Profile vs Temperature\n(30 minute diffusion)")
ax3.legend(title="Temperature")
ax3.set_ylim([1e14, 1e22])
ax3.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("diffusion_temperature_dependence.png", dpi=150)
print("Saved temperature dependence plot to: diffusion_temperature_dependence.png")
