"""Generate process parameter figures from Park et al. JMEMS 2010 Part I.

Plots Beta*, Nz, and Rs vs diffusion length sqrt(Dt).
Note: Lookup table covers anneal times 15-120 min; paper shows 1-900 min.
"""

from pathlib import Path

import addcopyfighandler  # noqa: F401
import matplotlib.pyplot as plt
import numpy as np
import scipy.io

plt.style.use("piezod.default")

# Load lookup table
lookup_path = Path(__file__).parent.parent.parent / "matlab" / "PiezoD" / "lookupTable.mat"
mat_data = scipy.io.loadmat(str(lookup_path), squeeze_me=True)

# Available parameter ranges from lookup table
DOSES = mat_data["ImplantDoses"]  # [2e14, 2e15, 2e16]
TEMPS_C = mat_data["AnnealTemps"]  # [900, 1000, 1100]
TIMES_MIN = mat_data["AnnealTimes"]  # [15, 30, 45, 60, 75, 90, 105, 120]

# Fixed parameters
T_CANTILEVER = 7e-6  # 7 um thick device
ENERGY = 50  # 50 keV implant

# Color map for temperatures
TEMP_COLORS = {900: "blue", 1000: "green", 1100: "red"}

# Boron diffusion parameters for sqrt(Dt) calculation
D0_BORON = 0.76  # cm^2/s
EA_BORON = 3.46  # eV
K_B_EV = 8.617343e-5  # eV/K


def calc_sqrt_Dt(temp_c, time_min):
    """Calculate diffusion length sqrt(Dt) in cm."""
    temp_k = float(temp_c) + 273.15
    time_s = float(time_min) * 60
    D = D0_BORON * np.exp(-EA_BORON / K_B_EV / temp_k)
    return np.sqrt(D * time_s)


def interpolate_value(field, dose, temp_c, time_min):
    """Interpolate a value from the lookup table."""
    from scipy.interpolate import interpn

    return float(
        interpn(
            (
                mat_data["ImplantDopants"].astype(float),
                mat_data["ImplantDoses"],
                mat_data["ImplantEnergies"].astype(float),
                mat_data["AnnealTemps"].astype(float),
                mat_data["AnnealTimes"].astype(float),
                mat_data["AnnealOxidation"].astype(float),
            ),
            mat_data[field],
            np.array([[1, dose, ENERGY, temp_c, time_min, 2]]),  # boron, oxide
            method="linear",
        )[0]
    )


def get_beta(dose, temp_c, time_min):
    """Calculate beta* from lookup table Beta1 and Beta2."""
    Beta1 = interpolate_value("Beta1", dose, temp_c, time_min)
    Beta2 = interpolate_value("Beta2", dose, temp_c, time_min)
    t_um = T_CANTILEVER * 1e6
    return Beta1 - 2 / t_um * Beta2


def sweep_parameters():
    """Sweep all parameter combinations and collect results."""
    results = []
    for dose in DOSES:
        for temp_c in TEMPS_C:
            for time_min in TIMES_MIN:
                sqrt_Dt = calc_sqrt_Dt(temp_c, time_min)
                beta = get_beta(dose, temp_c, time_min)
                Nz = interpolate_value("Nz", dose, temp_c, time_min) / 1e4  # m^-2 to cm^-2
                Rs = interpolate_value("Rs", dose, temp_c, time_min)
                results.append(
                    {
                        "dose": dose,
                        "temp_c": temp_c,
                        "time_min": time_min,
                        "sqrt_Dt": sqrt_Dt,
                        "beta": beta,
                        "Nz": Nz,
                        "Rs": Rs,
                    }
                )
    return results


def plot_beta(results):
    """Plot Beta* vs sqrt(Dt)."""
    fig, ax = plt.subplots(figsize=(6, 5))

    for dose in DOSES:
        for temp_c in TEMPS_C:
            data = [r for r in results if r["dose"] == dose and r["temp_c"] == temp_c]
            sqrt_Dt = [r["sqrt_Dt"] for r in data]
            beta = [r["beta"] for r in data]
            ax.plot(sqrt_Dt, beta, color=TEMP_COLORS[temp_c], linewidth=1)

    ax.set_xscale("log")
    ax.set_xlabel(r"$\sqrt{Dt}$ (cm)")
    ax.set_ylabel(r"$\beta^*$")
    ax.set_ylim(0, 1)
    ax.set_title(r"Efficiency Factor $\beta^*$ (7 $\mu$m device, 50 keV B)")

    for temp_c, color in TEMP_COLORS.items():
        ax.plot([], [], color=color, label=f"{temp_c}C")
    ax.legend(title="Anneal Temp", loc="lower left")

    fig.tight_layout()
    return fig


def plot_Nz(results):
    """Plot Nz vs sqrt(Dt)."""
    fig, ax = plt.subplots(figsize=(6, 5))

    for dose in DOSES:
        for temp_c in TEMPS_C:
            data = [r for r in results if r["dose"] == dose and r["temp_c"] == temp_c]
            sqrt_Dt = [r["sqrt_Dt"] for r in data]
            Nz = [r["Nz"] for r in data]
            ax.plot(sqrt_Dt, Nz, color=TEMP_COLORS[temp_c], linewidth=1)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\sqrt{Dt}$ (cm)")
    ax.set_ylabel(r"$N_z$ (ions/cm$^2$)")
    ax.set_title(r"Effective Doping $N_z$ (50 keV B)")

    for temp_c, color in TEMP_COLORS.items():
        ax.plot([], [], color=color, label=f"{temp_c}C")
    ax.legend(title="Anneal Temp", loc="lower left")

    fig.tight_layout()
    return fig


def plot_Rs(results):
    """Plot Rs vs sqrt(Dt)."""
    fig, ax = plt.subplots(figsize=(6, 5))

    for dose in DOSES:
        for temp_c in TEMPS_C:
            data = [r for r in results if r["dose"] == dose and r["temp_c"] == temp_c]
            sqrt_Dt = [r["sqrt_Dt"] for r in data]
            Rs = [r["Rs"] for r in data]
            ax.plot(sqrt_Dt, Rs, color=TEMP_COLORS[temp_c], linewidth=1)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\sqrt{Dt}$ (cm)")
    ax.set_ylabel(r"$R_s$ ($\Omega$/sq)")
    ax.set_title(r"Sheet Resistance $R_s$ (50 keV B)")

    for temp_c, color in TEMP_COLORS.items():
        ax.plot([], [], color=color, label=f"{temp_c}C")
    ax.legend(title="Anneal Temp", loc="upper right")

    fig.tight_layout()
    return fig


if __name__ == "__main__":
    print("Sweeping parameters...")
    results = sweep_parameters()
    print(f"Collected {len(results)} data points")

    output_dir = Path(__file__).parent

    print("Generating Beta* plot...")
    fig_beta = plot_beta(results)
    fig_beta.savefig(output_dir / "fig_beta.png", dpi=150)

    print("Generating Nz plot...")
    fig_Nz = plot_Nz(results)
    fig_Nz.savefig(output_dir / "fig_Nz.png", dpi=150)

    print("Generating Rs plot...")
    fig_Rs = plot_Rs(results)
    fig_Rs.savefig(output_dir / "fig_Rs.png", dpi=150)

    print("Saved: fig_beta.png, fig_Nz.png, fig_Rs.png")
