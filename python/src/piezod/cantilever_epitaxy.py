"""Cantilever with epitaxial doping profile."""

from piezod.cantilever import Cantilever


class CantileverEpitaxy(Cantilever):
    """Cantilever with epitaxial piezoresistor doping."""

    def __init__(self) -> None:
        """Initialize epitaxial cantilever with default parameters."""
        super().__init__()
        self.dopant_concentration = 1e19
        self.t_pr_ratio = 0.3

    def doping_profile(self) -> None:
        """Calculate the doping profile."""
        return

    def doping_optimization_scaling(self) -> None:
        """Get optimization scaling for doping parameters."""
        return

    def doping_cantilever_from_state(self) -> None:
        """Reconstruct cantilever from optimization state."""
        return

    def doping_current_state(self) -> None:
        """Get current doping state."""
        return

    def doping_initial_conditions_random(self) -> None:
        """Generate random initial conditions for optimization."""
        return

    def doping_optimization_bounds(self, parameter_constraints: dict) -> None:
        """Get optimization bounds for doping parameters."""
        return

    def Nz(self) -> None:
        """Calculate carrier concentration profile."""
        return

    def alpha(self) -> None:
        """Calculate piezoresistive coefficient."""
        return

    def sheet_resistance(self) -> None:
        """Calculate sheet resistance."""
        return
