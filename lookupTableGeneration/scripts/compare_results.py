"""Compare FLOOXS simulation results to TSUPREM-4 lookup table values.

Parses FLOOXS output, calculates Xj and Rs, and compares to reference values.

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

    Args:
        output_path: Path to FLOOXS output file

    Returns:
        Dictionary with 'pre_anneal' and 'post_anneal' profiles,
        each containing 'depth' and 'concentration' arrays
    """
    with open(output_path) as f:
        content = f.read()

    profiles = {}

    # Parse pre-anneal profile
    pre_match = re.search(
        r"=== PRE-ANNEAL PROFILE ===\s*\n"
        r"Depth\(um\) Concentration\(cm\^-3\)\s*\n"
        r"([\s\S]*?)(?:===|$)",
        content,
    )

    if pre_match:
        lines = pre_match.group(1).strip().split("\n")
        depths = []
        concs = []
        for line in lines:
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    depths.append(float(parts[0]))
                    concs.append(float(parts[1]))
                except ValueError:
                    continue
        profiles["pre_anneal"] = {
            "depth": np.array(depths),
            "concentration": np.array(concs),
        }

    # Parse post-anneal profile
    post_match = re.search(
        r"=== POST-ANNEAL PROFILE ===\s*\n"
        r"Depth\(um\) Concentration\(cm\^-3\)\s*\n"
        r"([\s\S]*?)(?:===|$)",
        content,
    )

    if post_match:
        lines = post_match.group(1).strip().split("\n")
        depths = []
        concs = []
        for line in lines:
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    depths.append(float(parts[0]))
                    concs.append(float(parts[1]))
                except ValueError:
                    continue
        profiles["post_anneal"] = {
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
    # Find where concentration crosses background doping
    for i in range(len(concentration) - 1):
        if concentration[i] > background_doping and concentration[i + 1] <= background_doping:
            # Linear interpolation to find exact crossing
            frac = (concentration[i] - background_doping) / (concentration[i] - concentration[i + 1])
            return float(depth[i] + frac * (depth[i + 1] - depth[i]))

    # If no crossing found, return last depth where concentration > background
    idx = np.where(concentration > background_doping)[0]
    if len(idx) > 0:
        return float(depth[idx[-1]])

    return 0.0


def calculate_sheet_resistance(
    depth: NDArray,
    concentration: NDArray,
    doping_type: str = "boron",
    temperature: float = 300.0,
) -> float:
    """Calculate sheet resistance from concentration profile.

    Uses Reggiani et al. (2002) mobility model, matching PiezoD MATLAB code.
    Rs = 1 / integral(q * mu(N) * N dx)

    Args:
        depth: Depth array in microns
        concentration: Concentration array in cm^-3
        doping_type: Dopant type for mobility calculation
        temperature: Temperature in Kelvin (default 300K)

    Returns:
        Sheet resistance in ohm/sq
    """
    q = 1.602e-19  # Coulombs
    k_b_eV = 8.617e-5  # eV/K
    T = temperature
    Tnorm = T / 300.0

    # Bandgap and intrinsic carrier concentration
    Eg = 1.170 - (4.730e-4 * T**2) / (T + 636)
    ni = np.sqrt(2.4e31 * T**3 * np.exp(-Eg / (k_b_eV * T)))

    # Reggiani et al. (2002) mobility model parameters
    if doping_type == "boron":
        mumax = 470.5
        c = 0.0
        gamma = 2.16
        mu0d = 90.0 * Tnorm**-1.3
        mu0a = 44.0 * Tnorm**-0.7
        mu1d = 28.2 * Tnorm**-2.0
        mu1a = 28.2 * Tnorm**-0.8
        Cr1 = 1.3e18 * Tnorm**2.2
        Cr2 = 2.45e17 * Tnorm**3.1
        Cs1 = 1.1e18 * Tnorm**6.2
        Cs2 = 6.1e20
        alpha1 = 0.77
        alpha2 = 0.719
        # For p-type (boron), p ≈ dopant concentration
        p = concentration / 2 + np.sqrt((concentration / 2) ** 2 + ni**2)
        n = ni**2 / p
        ND = n
        NA = p
    elif doping_type == "phosphorus":
        mumax = 1441
        c = 0.07
        gamma = 2.45
        mu0d = 62.2 * Tnorm**-0.7
        mu0a = 132.0 * Tnorm**-1.3
        mu1d = 48.6 * Tnorm**-0.7
        mu1a = 73.5 * Tnorm**-1.25
        Cr1 = 8.5e16 * Tnorm**3.65
        Cr2 = 1.22e17 * Tnorm**2.65
        Cs1 = 4e20
        Cs2 = 7e20
        alpha1 = 0.68
        alpha2 = 0.72
        # For n-type (phosphorus), n ≈ dopant concentration
        n = concentration / 2 + np.sqrt((concentration / 2) ** 2 + ni**2)
        p = ni**2 / n
        ND = n
        NA = p
    else:  # arsenic
        mumax = 1441
        c = 0.07
        gamma = 2.45
        mu0d = 55.0 * Tnorm**-0.6
        mu0a = 132.0 * Tnorm**-1.3
        mu1d = 42.4 * Tnorm**-0.5
        mu1a = 73.5 * Tnorm**-1.25
        Cr1 = 8.9e16 * Tnorm**3.65
        Cr2 = 1.22e17 * Tnorm**2.65
        Cs1 = 2.9e20
        Cs2 = 7e20
        alpha1 = 0.68
        alpha2 = 0.72
        # For n-type (arsenic), n ≈ dopant concentration
        n = concentration / 2 + np.sqrt((concentration / 2) ** 2 + ni**2)
        p = ni**2 / n
        ND = n
        NA = p

    # Calculate mobility using Reggiani model
    mu0 = (mu0d * ND + mu0a * NA) / (ND + NA)
    mu1 = (mu1d * ND + mu1a * NA) / (ND + NA)
    muL = mumax * Tnorm ** (-gamma + c * Tnorm)
    second = (muL - mu0) / (1 + (ND / Cr1) ** alpha1 + (NA / Cr2) ** alpha2)
    third = mu1 / (1 + (ND / Cs1 + NA / Cs2) ** -2)
    mobility = mu0 + second - third

    # Calculate conductivity (ohm-cm)^-1
    # For p-type, conductivity dominated by holes
    if doping_type == "boron":
        conductivity = q * mobility * p
    else:
        conductivity = q * mobility * n

    # Integrate using trapezoidal rule (depth in um, convert to cm)
    depth_cm = depth * 1e-4
    conductance = np.trapezoid(conductivity, depth_cm)

    # Sheet resistance
    if conductance > 0:
        return float(1.0 / conductance)
    return float("inf")


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

    # Calculate Xj and Rs from FLOOXS output
    xj_flooxs = calculate_junction_depth(post["depth"], post["concentration"])
    rs_flooxs = calculate_sheet_resistance(post["depth"], post["concentration"], args.dopant)

    print(f"\nFLOOXS results:")
    print(f"  Junction depth (Xj): {xj_flooxs:.4f} um")
    print(f"  Sheet resistance (Rs): {rs_flooxs:.2f} ohm/sq")

    # Load lookup table and get reference values
    script_dir = Path(__file__).parent
    mat_path = script_dir.parent / "lookupTable.mat"

    if not mat_path.exists():
        print(f"\nLookup table not found at: {mat_path}")
        print("Cannot compare to reference values")
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

    xj_ref = ref["Xj"] * 1e6  # Convert m to um
    rs_ref = ref["Rs"]

    print(f"\nTSUPREM-4 reference values:")
    print(f"  Junction depth (Xj): {xj_ref:.4f} um")
    print(f"  Sheet resistance (Rs): {rs_ref:.2f} ohm/sq")

    # Calculate errors
    xj_error = abs(xj_flooxs - xj_ref) / xj_ref * 100
    rs_error = abs(rs_flooxs - rs_ref) / rs_ref * 100

    print(f"\nComparison:")
    print(f"  Xj error: {xj_error:.2f}%")
    print(f"  Rs error: {rs_error:.2f}%")

    if xj_error < 2 and rs_error < 2:
        print("\nVALIDATION PASSED: Both errors < 2%")
    else:
        print("\nVALIDATION FAILED: One or more errors >= 2%")
        sys.exit(1)


if __name__ == "__main__":
    main()
