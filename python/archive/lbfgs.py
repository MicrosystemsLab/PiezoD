from cantilever_divingboard import *

# We need to scale the parameters before applying the optimization algorithm
# Normally there are about 20 orders of magnitude between the dimensions and 
# the doping concentration, so this is a critical step

        

# Run the script
freq_min = 1e3
freq_max = 1e5
omega_min = 100e3
initial_guess = (50e-6, 1e-6, 1e-6, 
                30e-6, 1e-6, 1e-6, 500e-9, 5., 1e15)
constraints = ((30e-6, 100e-6), (500e-9, 20e-6), (1e-6, 10e-6),
               (2e-6, 100e-6), (500e-9, 5e-6), (500e-9, 20e-6), (30e-9, 10e-6),
               (1., 10.), (1e15, 4e19))
x = optimize_cantilever(initial_guess, constraints, freq_min, freq_max, omega_min)
c = cantilever_divingboard(freq_min, freq_max, x)
c.print_performance()