"""Compare FLOOXS simulation results to TSUPREM-4 lookup table values.

Parses FLOOXS output, calculates Xj, Rs, Beta1, Beta2, Nz,
and compares to reference values.

Usage:
    python scripts/compare_results.py simulations/output.txt
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from extract_reference import extract_reference_values, load_lookup_table


def parse_flooxs_output(output_path: Path) -> dict[str, NDArray]:
    """Parse FLOOXS output file to extract depth vs concentration profiles.

    Handles both old format (=== PRE-ANNEAL PROFILE ===) and new format
    (=== PRE-ANNEAL ===) with raw depth/concentration pairs.

    Args:
        output_path: Path to FLOOXS output file

    Returns:
        Dictionary with 'pre_anneal' and 'post_anneal' profiles,
        each containing 'depth' and 'concentration' arrays
    """
    with open(output_path) as f:
        lines = f.readlines()

    profiles = {}
    for key, marker in [("pre_anneal", "PRE-ANNEAL"), ("post_anneal", "POST-ANNEAL")]:
        depths = []
        concs = []
        in_section = False
        for line in lines:
            stripped = line.strip()
            if f"=== {marker}" in stripped and "INTERST" not in stripped and "VACAN" not in stripped:
                in_section = True
                continue
            if in_section and stripped.startswith("==="):
                break
            if in_section:
                parts = stripped.split()
                if len(parts) == 2:
                    try:
                        d, c = float(parts[0]), float(parts[1])
                        if d >= 0:
                            depths.append(d)
                            concs.append(c)
                    except ValueError:
                        continue
        if depths:
            profiles[key] = {
                "depth": np.array(depths),
                "concentration": np.array(concs),
            }

    return profiles


def calculate_junction_depth(
    depth: NDArray,
    concentration: NDArray,
    background_doping: float = 1e15,
) -> float:
    """Calculate junction depth where doping equals background level.

    Args:
        depth: Depth array in microns
        concentration: Concentration array in cm^-3
        background_doping: Background doping level in cm^-3

    Returns:
        Junction depth in microns
    """
    for i in range(len(concentration) - 1):
        if concentration[i] > background_doping and concentration[i + 1] <= background_doping:
            frac = (concentration[i] - background_doping) / (concentration[i] - concentration[i + 1])
            return float(depth[i] + frac * (depth[i + 1] - depth[i]))

    idx = np.where(concentration > background_doping)[0]
    if len(idx) > 0:
        return float(depth[idx[-1]])

    return 0.0


def _mobility_and_conductivity(
    concentration: NDArray,
    doping_type: str,
    temperature: float = 300.0,
) -> tuple[NDArray, NDArray]:
    """Compute mobility and conductivity from Reggiani et al. (2002) model.

    Matches PiezoD cantilever.py mobility() implementation.

    Args:
        concentration: Doping concentration in cm^-3
        doping_type: "boron", "phosphorus", or "arsenic"
        temperature: Temperature in Kelvin

    Returns:
        Tuple of (mobility in cm^2/V-s, conductivity in (ohm-cm)^-1)
    """
    q = 1.602e-19
    k_b_eV = 8.617e-5
    T = temperature
    Tnorm = T / 300.0

    Eg = 1.170 - (4.730e-4 * T**2) / (T + 636)
    ni = np.sqrt(2.4e31 * T**3 * np.exp(-Eg / (k_b_eV * T)))

    if doping_type == "boron":
        mumax, c, gamma = 470.5, 0.0, 2.16
        mu0d = 90.0 * Tnorm**-1.3
        mu0a = 44.0 * Tnorm**-0.7
        mu1d = 28.2 * Tnorm**-2.0
        mu1a = 28.2 * Tnorm**-0.8
        Cr1 = 1.3e18 * Tnorm**2.2
        Cr2 = 2.45e17 * Tnorm**3.1
        Cs1 = 1.1e18 * Tnorm**6.2
        Cs2 = 6.1e20
        alpha1, alpha2 = 0.77, 0.719
        p = concentration / 2 + np.sqrt((concentration / 2) ** 2 + ni**2)
        n = ni**2 / p
        ND, NA = n, p
    elif doping_type == "phosphorus":
        mumax, c, gamma = 1441, 0.07, 2.45
        mu0d = 62.2 * Tnorm**-0.7
        mu0a = 132.0 * Tnorm**-1.3
        mu1d = 48.6 * Tnorm**-0.7
        mu1a = 73.5 * Tnorm**-1.25
        Cr1 = 8.5e16 * Tnorm**3.65
        Cr2 = 1.22e17 * Tnorm**2.65
        Cs1, Cs2 = 4e20, 7e20
        alpha1, alpha2 = 0.68, 0.72
        n = concentration / 2 + np.sqrt((concentration / 2) ** 2 + ni**2)
        p = ni**2 / n
        ND, NA = n, p
    else:  # arsenic
        mumax, c, gamma = 1441, 0.07, 2.45
        mu0d = 55.0 * Tnorm**-0.6
        mu0a = 132.0 * Tnorm**-1.3
        mu1d = 42.4 * Tnorm**-0.5
        mu1a = 73.5 * Tnorm**-1.25
        Cr1 = 8.9e16 * Tnorm**3.65
        Cr2 = 1.22e17 * Tnorm**2.65
        Cs1, Cs2 = 2.9e20, 7e20
        alpha1, alpha2 = 0.68, 0.72
        n = concentration / 2 + np.sqrt((concentration / 2) ** 2 + ni**2)
        p = ni**2 / n
        ND, NA = n, p

    mu0 = (mu0d * ND + mu0a * NA) / (ND + NA)
    mu1 = (mu1d * ND + mu1a * NA) / (ND + NA)
    muL = mumax * Tnorm ** (-gamma + c * Tnorm)
    second = (muL - mu0) / (1 + (ND / Cr1) ** alpha1 + (NA / Cr2) ** alpha2)
    third = mu1 / (1 + (ND / Cs1 + NA / Cs2) ** -2)
    mobility = mu0 + second - third

    if doping_type == "boron":
        conductivity = q * mobility * p
    else:
        conductivity = q * mobility * n

    return mobility, conductivity


def _piezoresistance_factor(concentration: NDArray, temperature: float = 300.0) -> NDArray:
    """Piezoresistance factor P(N) from Richter et al.

    Matches PiezoD cantilever.py piezoresistance_factor().

    Args:
        concentration: Doping concentration in cm^-3
        temperature: Temperature in Kelvin

    Returns:
        Piezoresistance factor (dimensionless, 0 to 1)
    """
    Nb = 6e19
    Nc = 7e20
    alpha = 0.43
    gamma = 1.6
    beta = 0.1
    eta = 3
    theta = 0.9

    T0 = 300
    Theta = temperature / T0

    return (
        Theta**-theta
        * (
            1
            + Theta**-beta * (concentration / Nb) ** alpha
            + Theta**-eta * (concentration / Nc) ** gamma
        )
        ** -1
    )


def calculate_sheet_resistance(
    depth: NDArray,
    concentration: NDArray,
    doping_type: str = "boron",
    temperature: float = 300.0,
) -> float:
    """Calculate sheet resistance from concentration profile.

    Args:
        depth: Depth array in microns
        concentration: Concentration array in cm^-3
        doping_type: Dopant type for mobility calculation
        temperature: Temperature in Kelvin

    Returns:
        Sheet resistance in ohm/sq
    """
    _, conductivity = _mobility_and_conductivity(concentration, doping_type, temperature)
    depth_cm = depth * 1e-4
    conductance = np.trapezoid(conductivity, depth_cm)
    if conductance > 0:
        return float(1.0 / conductance)
    return float("inf")


def calculate_beta_nz(
    depth: NDArray,
    concentration: NDArray,
    doping_type: str = "boron",
    temperature: float = 300.0,
) -> tuple[float, float, float]:
    """Calculate Beta1, Beta2, and Nz from concentration profile.

    Beta1 and Beta2 are thickness-independent decompositions of beta*:
        beta*(t) = Beta1 - 2 * Beta2 / t

    Derived from PiezoD cantilever.py beta():
        beta* = 2/(t * S) * integral(sigma * P * z dz)
    where z = t/2 - depth, S = integral(sigma dz).

    Expanding z = t/2 - d:
        Beta1 = integral(sigma * P dz) / integral(sigma dz)
        Beta2 = integral(sigma * P * d dz) / integral(sigma dz)
    where d is depth from surface in microns.

    Nz is mobility-weighted effective carrier density (cm^-2):
        Nz = (integral(N*mu*dz))^2 / integral(N*mu^2*dz)
    Matches PiezoD cantilever_diffusion.py Nz().

    Args:
        depth: Depth array in microns
        concentration: Concentration array in cm^-3
        doping_type: Dopant type
        temperature: Temperature in Kelvin

    Returns:
        Tuple of (Beta1, Beta2 in um, Nz in cm^-2)
    """
    mobility, conductivity = _mobility_and_conductivity(concentration, doping_type, temperature)
    P = _piezoresistance_factor(concentration, temperature)

    depth_cm = depth * 1e-4

    S = np.trapezoid(conductivity, depth_cm)
    sigma_P = np.trapezoid(conductivity * P, depth_cm)
    # depth in cm for integral, but Beta2 is in um so convert back
    sigma_P_d = np.trapezoid(conductivity * P * depth * 1e-4, depth_cm)

    Beta1 = float(sigma_P / S) if S > 0 else 0.0
    # sigma_P_d is in cm^2 units, convert to um for Beta2
    Beta2 = float(sigma_P_d / S * 1e4) if S > 0 else 0.0

    # Mobility-weighted effective carrier density
    num = np.trapezoid(concentration * mobility, depth_cm) ** 2
    den = np.trapezoid(concentration * mobility**2, depth_cm)
    Nz = float(num / den) if den > 0 else 0.0

    return Beta1, Beta2, Nz


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare FLOOXS results to lookup table")
    parser.add_argument("output_file", type=Path, help="FLOOXS output file")
    parser.add_argument("--dopant", type=str, default="boron", choices=["boron", "phosphorus", "arsenic"])
    parser.add_argument("--dose", type=float, default=2e15, help="Implant dose (cm^-2)")
    parser.add_argument("--energy", type=int, default=20, help="Implant energy (keV)")
    parser.add_argument("--temp", type=int, default=1000, help="Anneal temperature (C)")
    parser.add_argument("--time", type=int, default=30, help="Anneal time (min)")
    parser.add_argument("--oxidation", type=str, default="inert", choices=["inert", "oxide"])

    args = parser.parse_args()

    if not args.output_file.exists():
        print(f"Output file not found: {args.output_file}")
        sys.exit(1)

    # Parse FLOOXS output
    print(f"Parsing FLOOXS output: {args.output_file}")
    profiles = parse_flooxs_output(args.output_file)

    if "post_anneal" not in profiles:
        print("ERROR: Could not find post-anneal profile in output")
        sys.exit(1)

    post = profiles["post_anneal"]
    print(f"Found {len(post['depth'])} data points in post-anneal profile")

    # Calculate all metrics from FLOOXS output
    xj = calculate_junction_depth(post["depth"], post["concentration"])
    rs = calculate_sheet_resistance(post["depth"], post["concentration"], args.dopant)
    beta1, beta2, nz = calculate_beta_nz(post["depth"], post["concentration"], args.dopant)

    print(f"\nFLOOXS results:")
    print(f"  Xj:    {xj:.4f} um")
    print(f"  Rs:    {rs:.2f} ohm/sq")
    print(f"  Beta1: {beta1:.4e}")
    print(f"  Beta2: {beta2:.4e} um")
    print(f"  Nz:    {nz:.4e} cm^-2")

    # Load lookup table and get reference values
    script_dir = Path(__file__).parent
    mat_path = script_dir.parent / "lookupTable.mat"

    if not mat_path.exists():
        print(f"\nLookup table not found at: {mat_path}")
        return

    data = load_lookup_table(mat_path)
    dopant_idx = {"boron": 1, "phosphorus": 2, "arsenic": 3}[args.dopant]
    oxidation_idx = {"inert": 1, "oxide": 2}[args.oxidation]

    ref = extract_reference_values(
        data,
        dopant_idx=dopant_idx,
        dose=args.dose,
        energy=args.energy,
        temp=args.temp,
        time=args.time,
        oxidation_idx=oxidation_idx,
    )

    xj_ref = ref["Xj"] * 1e6  # m to um
    rs_ref = ref["Rs"]
    beta1_ref = ref["Beta1"]
    beta2_ref = ref["Beta2"]
    nz_ref = ref["Nz"] * 1e-4  # m^-2 to cm^-2

    print(f"\nTSUPREM-4 reference:")
    print(f"  Xj:    {xj_ref:.4f} um")
    print(f"  Rs:    {rs_ref:.2f} ohm/sq")
    print(f"  Beta1: {beta1_ref:.4e}")
    print(f"  Beta2: {beta2_ref:.4e} um")
    print(f"  Nz:    {nz_ref:.4e} cm^-2")

    def pct_err(val: float, ref_val: float) -> str:
        if ref_val == 0:
            return "N/A"
        return f"{(val - ref_val) / ref_val * 100:+.1f}%"

    print(f"\nComparison:")
    print(f"  Xj:    {pct_err(xj, xj_ref)}")
    print(f"  Rs:    {pct_err(rs, rs_ref)}")
    print(f"  Beta1: {pct_err(beta1, beta1_ref)}")
    print(f"  Beta2: {pct_err(beta2, beta2_ref)}")
    print(f"  Nz:    {pct_err(nz, nz_ref)}")


if __name__ == "__main__":
    main()
