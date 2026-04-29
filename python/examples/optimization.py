"""Cantilever optimization example.

Mirrors the optimization section of the MATLAB tutorial: an epitaxial
piezoresistive cantilever is optimized for force resolution subject to
power-dissipation and resonant-frequency constraints, with the bridge
voltage capped at 10 V and the thickness floor at 1 um.

Run:
    uv run python examples/optimization.py
"""

from __future__ import annotations

import warnings

from piezod import (
    CantileverEpitaxy,
    CantileverMetric,
    CantileverMetricConstraint,
    force_resolution_goal,
    optimize_performance,
)


def main() -> None:
    # Starting design: 200 um x 20 um x 2 um epitaxial boron piezoresistor.
    c = CantileverEpitaxy(
        freq_min=1,
        freq_max=1000,
        l=200e-6,
        w=20e-6,
        t=2e-6,
        l_pr_ratio=0.3,
        v_bridge=2.0,
        doping_type="boron",
        dopant_concentration=1e19,
        t_pr_ratio=0.3,
    )
    c.fluid = "vacuum"
    c.number_of_piezoresistors = 4

    print("=== Starting design ===")
    print(f"  geometry: l={c.l * 1e6:.1f} um, w={c.w * 1e6:.1f} um, t={c.t * 1e6:.2f} um")
    print(f"  v_bridge: {c.v_bridge:.2f} V")
    print(f"  N_dopant: {c.dopant_concentration:.2e} cm^-3")
    print(f"  t_pr_ratio: {c.t_pr_ratio:.3f}")
    print(f"  force_resolution: {c.force_resolution() * 1e12:.2f} pN")
    print(f"  power: {c.power_dissipation() * 1e3:.3f} mW")
    print(f"  vacuum freq: {c.omega_vacuum_hz() / 1e3:.1f} kHz")
    print(f"  stiffness: {c.stiffness():.3g} N/m")
    print()

    constraints = [
        CantileverMetricConstraint(CantileverMetric.POWER_DISSIPATION, maximum=2e-3),
        CantileverMetricConstraint(CantileverMetric.OMEGA_VACUUM_HZ, minimum=5 * 1000),
        CantileverMetricConstraint(CantileverMetric.STIFFNESS, minimum=1e-3, maximum=1e1),
    ]

    print("Optimizing for force resolution with multi-start (this can take ~10 s)...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = optimize_performance(
            c,
            force_resolution_goal(),
            parameter_constraints={"max_v_bridge": 10.0},
            metric_constraints=constraints,
            n_starts=5,
            max_iterations=10,
            random_seed=0,
        )

    opt = result.optimized
    print(f"  ran {len(result.all_results)} starts; best objective {result.objective_value:.3f} pN")
    print()

    print("=== Optimized design ===")
    print(f"  geometry: l={opt.l * 1e6:.1f} um, w={opt.w * 1e6:.1f} um, t={opt.t * 1e6:.3f} um")
    print(f"  v_bridge: {opt.v_bridge:.2f} V")
    print(f"  N_dopant: {opt.dopant_concentration:.2e} cm^-3")
    print(f"  t_pr_ratio: {opt.t_pr_ratio:.3f}")
    print(f"  force_resolution: {opt.force_resolution() * 1e12:.3f} pN")
    print(f"  power: {opt.power_dissipation() * 1e3:.3f} mW (max {2.0:.1f})")
    print(f"  vacuum freq: {opt.omega_vacuum_hz() / 1e3:.2f} kHz (min {5.0:.1f})")
    print(f"  stiffness: {opt.stiffness():.3g} N/m")
    print(f"  l/w ratio: {opt.l / opt.w:.2f} (default >= 2)")
    print(f"  w/t ratio: {opt.w / opt.t:.2f} (default >= 2)")
    print()
    improvement = c.force_resolution() / opt.force_resolution()
    print(f"Force resolution improved by {improvement:.1f}x")


if __name__ == "__main__":
    main()
