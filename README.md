# Introduction
PiezoD is a set of Matlab classes that implement analytical models for piezoresistive and piezoelectric sensor cantilever performance. Several types of piezoresistive cantilevers are supported: 1) single crystal silicon doped via diffusion, epitaxy or ion implantation, 2) polycrystalline silicon and metals, and 3) piezoelectric cantilevers.

# Overview of Analytical Models
There are four basic parts to the analytical models: 1) mechanics, 2) electrical, 3) thermal. 

For each dopant type, the elastic modulus is calculated assuming that the piezoresistors are oriented along the optimal longitudinal axis (<110> for p-type, <100> for n-type). The piezoresistance factor is calculated as a function of dopant concentration using the model presented by Richter et al in 2008.

There are two models for the cantilever thermal conductivity. The first (simpler) is to calculate the temperature-independent thermal conductivity. If cantilever self-heating is small, then this is a good approximation. If the self-heating is substantial (e.g. temperature rises of >200C), then a more complicated k(T) model can be used.

# Circuit Assumptions
We assume that you're using a 1/4-active Wheatstone bridge and measuring the signal using an instrumentation amplifier. The noise characteristics of the INA103 and the AD8221 are currently included in the code and can be chosen by setting a flag.

# Optimization
Optimization is performed using fmincon, a Matlab function that is part of the Optimization Toolbox. It implements an l-bfgs-b optimizer that handles nonlinear constraints and uses finite differences to approximate the Hessian (matrix of second derivatives used in the optimization). The code supplies the model functions to the optimizer, allowing the optimization target (e.g. force resolution) to be calculated and its derivatives with respect to each of the active design parameters to be calculated. The design paramters vary by about 29 orders of magnitude in size (thickness ~ 100nm to dopant concentration ~ 1e20/cc), so the parameters are scaled to O(1) each iteration. There are several possible optimization goals that can be set by a flag - integrated force resolution, integrated displacement resolution, and force noise at a particular frequency (i.e. resonant sensing).

# Usage
Examples are provided in sampleCode.m.

# Code Structure and Extending the Code
Cantilever is the abstract base class that implements most of the details. In order to extend the code (i.e. alternative doping methods), you need to extend cantilever.m and implement these methods:
  doping_profile
  optimization_scaling
  cantilever_from_state
  optimize_performance
  optimization_constraints

You probably also would want to override
  print_performance
  print_performance_for_excel
  optimization_constraints

For an example of how to do this, take a look at cantileverDiffusion.m and cantileverEpitaxy.m.
