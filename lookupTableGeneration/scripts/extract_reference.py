"""Extract reference values from TSUPREM-4 lookup table for validation.

Reference case: Boron, 2e15 cm^-2, 20 keV, 1000C, 30 min, inert (no oxide)

Usage:
    python scripts/extract_reference.py
"""

from pathlib import Path

import scipy.io as sio
from scipy.interpolate import interpn

# Lookup table indices
# dopants: 1=B, 2=P, 3=As
# doses: [2e14, 2e15, 2e16]
# energies: [20, 50, 80] keV
# temps: [900, 1000, 1100] C
# times: [15, 30, 45, 60, 75, 90, 105, 120] min
# oxidation: 1=inert, 2=oxide


def load_lookup_table(mat_path: Path) -> dict:
    """Load the MATLAB lookup table."""
    data = sio.loadmat(mat_path)
    return data


def extract_reference_values(
    data: dict,
    dopant_idx: int = 1,  # 1=Boron
    dose: float = 2e15,
    energy: int = 20,
    temp: int = 1000,
    time: int = 30,
    oxidation_idx: int = 1,  # 1=inert
) -> dict[str, float]:
    """Extract Xj and Rs for a specific condition.

    Args:
        data: Lookup table data from scipy.io.loadmat
        dopant_idx: Dopant index (1=B, 2=P, 3=As)
        dose: Implant dose (cm^-2)
        energy: Implant energy (keV)
        temp: Anneal temperature (C)
        time: Anneal time (minutes)
        oxidation_idx: Oxidation type (1=inert, 2=oxide)

    Returns:
        Dictionary with Xj (m) and Rs (ohm/sq)
    """
    # Grid points for each dimension
    dopants = data["ImplantDopants"].flatten()
    doses = data["ImplantDoses"].flatten()
    energies = data["ImplantEnergies"].flatten()
    temps = data["AnnealTemps"].flatten()
    times = data["AnnealTimes"].flatten()
    oxidations = data["AnnealOxidation"].flatten()

    # Create interpolation grids
    points = (dopants, doses, energies, temps, times, oxidations)
    query = [[dopant_idx, dose, energy, temp, time, oxidation_idx]]

    # Extract Xj (junction depth in meters)
    xj = interpn(points, data["Xj"], query, method="linear")[0]

    # Extract Rs (sheet resistance in ohm/sq)
    rs = interpn(points, data["Rs"], query, method="linear")[0]

    # Extract Beta1 and Beta2 for piezoresistive coefficient
    beta1 = interpn(points, data["Beta1"], query, method="linear")[0]
    beta2 = interpn(points, data["Beta2"], query, method="linear")[0]

    # Extract Nz (effective doping concentration in cm^-2)
    nz = interpn(points, data["Nz"], query, method="linear")[0]

    return {
        "Xj": float(xj),
        "Rs": float(rs),
        "Beta1": float(beta1),
        "Beta2": float(beta2),
        "Nz": float(nz),
    }


def main() -> None:
    # Find the lookup table
    script_dir = Path(__file__).parent
    mat_path = script_dir.parent / "lookupTable.mat"

    if not mat_path.exists():
        print(f"Lookup table not found at: {mat_path}")
        return

    print(f"Loading lookup table from: {mat_path}")
    data = load_lookup_table(mat_path)

    # Reference case: Boron, 2e15 cm^-2, 20 keV, 1000C, 30 min, inert
    print("\nReference case: Boron, 2e15 cm^-2, 20 keV, 1000C, 30 min, inert")
    print("-" * 60)

    ref = extract_reference_values(
        data,
        dopant_idx=1,  # Boron
        dose=2e15,
        energy=20,
        temp=1000,
        time=30,
        oxidation_idx=1,  # inert
    )

    print(f"Junction depth (Xj): {ref['Xj'] * 1e6:.4f} um")
    print(f"Sheet resistance (Rs): {ref['Rs']:.2f} ohm/sq")
    print(f"Beta1: {ref['Beta1']:.4e} Pa^-1")
    print(f"Beta2: {ref['Beta2']:.4e} Pa^-1 um")
    print(f"Effective doping (Nz): {ref['Nz']:.4e} cm^-2")

    # Also print values for comparison with other conditions
    print("\n\nAdditional reference cases:")
    print("=" * 60)

    test_cases = [
        {"dose": 2e14, "energy": 20, "temp": 900, "time": 15},
        {"dose": 2e15, "energy": 50, "temp": 1000, "time": 60},
        {"dose": 2e16, "energy": 80, "temp": 1100, "time": 120},
    ]

    for case in test_cases:
        print(f"\nBoron, {case['dose']:.0e} cm^-2, {case['energy']} keV, {case['temp']}C, {case['time']} min, inert")
        ref = extract_reference_values(data, dopant_idx=1, oxidation_idx=1, **case)
        print(f"  Xj: {ref['Xj'] * 1e6:.4f} um, Rs: {ref['Rs']:.2f} ohm/sq")


if __name__ == "__main__":
    main()
