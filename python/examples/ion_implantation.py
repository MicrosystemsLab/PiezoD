"""Ion-implanted piezoresistive cantilever example.

Based on the optimized design from:
  Park, Doll, Rastegar, Pruitt - "Piezoresistive Cantilever Performance"
  Part I: Analytical Model for Sensitivity (JMEMS 2010)
  Part II: Optimization (JMEMS 2010)

This example demonstrates a cantilever optimized for force sensing with
0.05 N/m stiffness and 2.5 kHz resonant frequency, designed for studying
touch sensation in C. elegans.
"""

from pathlib import Path

import addcopyfighandler  # noqa: F401
import scipy.io

from piezod import CantileverImplantation

# Create ion-implanted cantilever based on Park et al. 2010 optimized design
# Parameters from Table III of Part II

# Cantilever geometry
L = 2000e-6  # length: 2000 um
W = 30e-6  # width: 30 um
T = 7e-6  # thickness: 7 um (SOI device layer)

# Piezoresistor geometry: U-shaped with 5 um air gap
# Length: 50-57 um depending on alignment
L_PR_RATIO = 57 / 2000  # piezoresistor length ratio

# Operating frequency band
FREQ_MIN = 1  # Hz
FREQ_MAX = 1000  # Hz

# Electrical parameters
V_BRIDGE = 2.0  # Wheatstone bridge bias (V), constrained by power dissipation

# Ion implantation process parameters (optimal from Fig. 10 of Part II)
IMPLANT_DOSE = 5e15  # 5e15 cm^-2
IMPLANT_ENERGY = 50  # 50 keV

# Annealing parameters (from Table III)
# Paper: Wet oxidation 15 min + N2 anneal 300 min = 315 min total
# Lookup table max is 120 min, so we use that (sqrt(Dt) won't match exactly)
ANNEAL_TEMP = 1000 + 273.15  # 1000C in Kelvin
ANNEAL_TIME = 120 * 60  # 120 min (lookup table max), paper uses 315 min

c = CantileverImplantation(
    freq_min=FREQ_MIN,
    freq_max=FREQ_MAX,
    l=L,
    w=W,
    t=T,
    l_pr_ratio=L_PR_RATIO,
    v_bridge=V_BRIDGE,
    doping_type="boron",
    annealing_time=ANNEAL_TIME,
    annealing_temp=ANNEAL_TEMP,
    annealing_type="oxide",  # wet oxidation for passivation
    implantation_energy=IMPLANT_ENERGY,
    implantation_dose=IMPLANT_DOSE,
)

# Set operating environment
c.fluid = "air"
c.number_of_piezoresistors = 4  # Full bridge for temperature compensation

# Load the lookup table for electrical calculations
lookup_path = Path(__file__).parent.parent.parent / "matlab" / "PiezoD" / "lookupTable.mat"
if lookup_path.exists():
    mat_data = scipy.io.loadmat(str(lookup_path), squeeze_me=True)
    c.load_lookup_table(mat_data)
    lookup_loaded = True
else:
    print(f"Warning: Lookup table not found at {lookup_path}")
    print("Electrical parameters requiring lookup table will not be available.")
    lookup_loaded = False

# Print design summary
print("=== Ion-Implanted Cantilever Design ===")
print("Based on Park et al. JMEMS 2010 (optimized for C. elegans force sensing)")
print()

# Expected values from Table III of Part II paper (Target column)
paper = {
    "stiffness_mN_m": 50,  # mN/m (target constraint)
    "freq_kHz": 2.5,  # kHz (target constraint)
    "resistance_ohm": 500,  # Total piezoresistor resistance (Ohm)
    "beta_star": 0.49,  # efficiency factor
    "power_mW": 2.0,  # mW (constrained by self-heating)
    "force_resolution_pN": 68.1,  # Theoretical F_min (pN)
    "junction_depth_um": 2.17,  # t_p from TSUPREM4 (um)
    "sqrt_Dt_um": 0.0374,  # Diffusion length (um)
    "sensitivity_V_N": 1414,  # S_FV (V/N)
}


def check_match(model_val, paper_val, tolerance=0.15):
    """Check if model value matches paper value within tolerance."""
    if isinstance(paper_val, tuple):
        low, high = paper_val
        return low <= model_val <= high
    return abs(model_val - paper_val) / paper_val < tolerance


def fmt_paper(val):
    """Format paper value for display."""
    if isinstance(val, tuple):
        return f"{val[0]}-{val[1]}"
    return f"{val}"


def print_row(name, model_str, paper_val, tol=0.15):
    """Print a comparison table row."""
    model_num = float(model_str)
    match = check_match(model_num, paper_val, tol)
    status = "OK" if match else "DIFF"
    print(f"  {name:<28} {model_str:>12} {fmt_paper(paper_val):>12} {status:>6}")


print("=" * 62)
print(f"  {'Parameter':<28} {'Model':>12} {'Paper':>12} {'':>6}")
print("=" * 62)

# Mechanical properties
stiffness = c.stiffness() * 1e3  # mN/m
freq = c.omega_vacuum_hz() / 1e3  # kHz
print_row("Stiffness (mN/m)", f"{stiffness:.2f}", paper["stiffness_mN_m"])
print_row("Resonant freq (kHz)", f"{freq:.2f}", paper["freq_kHz"])

# Electrical properties (require lookup table)
if lookup_loaded:
    Rs = c.sheet_resistance()
    beta = c.beta()
    Xj = c.junction_depth * 1e6  # um
    resistance = c.resistance()
    power = c.power_dissipation() * 1e3  # mW
    sqrt_Dt = c.diffusion_length * 1e4  # um

    print_row("Junction depth t_p (um)", f"{Xj:.2f}", paper["junction_depth_um"], tol=0.25)
    print_row("Resistance R_piezo (Ohm)", f"{resistance:.0f}", paper["resistance_ohm"], tol=0.25)
    print_row("Beta* (efficiency)", f"{beta:.2f}", paper["beta_star"])
    print_row("sqrt(Dt) (um)", f"{sqrt_Dt:.4f}", paper["sqrt_Dt_um"])
    print_row("Power dissipation (mW)", f"{power:.2f}", paper["power_mW"])

    # Performance metrics
    print("-" * 62)
    force_sensitivity = c.force_sensitivity()
    force_res = c.force_resolution() * 1e12  # pN

    print_row("Sensitivity S_FV (V/N)", f"{force_sensitivity:.0f}", paper["sensitivity_V_N"])
    print_row("Force resolution (pN)", f"{force_res:.1f}", paper["force_resolution_pN"])

    # Additional derived parameters
    print("-" * 62)
    print("  Additional Parameters:")
    print(f"    Sheet resistance: {Rs:.1f} Ohm/sq")
    print(f"    Hooge parameter alpha: {c.alpha():.2e}")
    print(f"    Number of squares: {c.number_of_squares():.1f}")
    print(f"    Piezoresistor width: {c.w_pr() * 1e6:.1f} um")

print("=" * 62)
