from Cantilever import *


class CantileverEpitaxy(Cantilever):

    def __init__(self):
        Cantilever.__init__(self)
        self.dopant_concentration = 1e19
        self.t_pr_ratio = 0.3

    def doping_profile(self):
        return

    def doping_optimization_scaling(self):
        return

    def doping_cantilever_from_state(self):
        return

    def doping_current_state(self):
        return

    def doping_initial_conditions_random(self):
        return

    def doping_optimization_bounds(self, parameter_constraints):
        return

    def Nz(self):
        return

    def alpha(self):
        return

    def sheet_resistance(self):
        return