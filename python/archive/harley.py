from cantilever_divingboard import *

t = 89e-9
w = 44e-6
l = 300e-6

l_pr = 45e-6/l # estimated from SEM in paper
w_pr = 44e-6/w
t_pr = 30e-9/t
w_gap = 3e-6/w

freq_min = 1e1
freq_max = 1e3

V_bias = 5
N = 4e19
beta = 0.7

c = cantilever_divingboard(freq_min, freq_max, 
                           l, w, t, 
                           l_pr, w_pr, w_gap, t_pr, 
                           V_bias, N)

c.print_performance()

print 'Force resolution (pN): %g calculated, %g actual' % (c.force_resolution(), 500e-15)
print 'Integrated noise (V): %g calculated, %g actual' % (c.integrated_noise(), 1.14e-6)
print 'Stiffness (N/m): %g calculated, %g actual' % (c.stiffness(), 3e-5)
print 'Knee frequency (Hz): %g calculated, %g actual' % (c.knee_frequency(), 1e3)
print 'Piezo factor(m^2/N): %g calculated, %g actual' % (c.piezoresistance_factor(), 4e-10)
print 'Johnson noise level: %g calculated, %g actual' % (c.johnson_noise_voltage(), 1.5e-8)
print ""

print c.omega_damped()

# # Harley noise spectrum (optional)
# c.plot_noise_spectrum()
# slope_of_harley_noise = (log10(1.75e-7)-log10(5e-8))/(log10(10)-log10(100))



# Design time
omega_min = c.omega_vacuum_hz()
lower_bounds = (10e-6, 1e-6, 89e-9, 0., 0., 0., 0., 0., 1e15)
upper_bounds = (300e-6, 100e-6, 30e-6, 1., 1., 1., 1., 15., 8e19)
max_power = 12.4755e-3
fluid_type = 'vacuum'
c.freq_min = 1e3
c.freq_max = 1e5
c.t = 1e-6
c.t_pr_ratio = 0.5

c.optimize_openopt(omega_min, lower_bounds, upper_bounds, max_power, fluid_type)
c.print_performance()




# Optimizer junkyard (thinks that haven't been working)

# # cobyla
# # Pack all of our constraints into a single tuple
# omega_min = c.omega_vacuum_hz()
# max_voltage = 15.
# max_power = 12.4755e-3
# min_N = 1e15
# max_N = 4e19
# fluid_type = 'vacuum'
# 
# c.l = 200e-6
# c.l_pr = 80e-6
# c.w = 10e-6
# c.w_pr = 10e-6
# c.t = 1e-6
# c.t_pr = 100e-9
# x = c.optimize_cantilever_new(omega_min, min_dimensions, max_voltage, max_power, min_N, max_N, fluid_type)
# c.print_performance()




# # lbfgs
# constraints = ((30e-6, 400e-6), (5e-6, 100e-6), (89e-9, 10e-6),
#                (0., 1.), (0., 1.), (0., 1.), (0., 1.),
#                (1., 15.), (1e15, 4e19))
# omega_min = c.omega_vacuum_hz()
# 
# c.freq_min = 1e1
# c.freq_max = 1e3
# c.V_bias = 10.
# c.optimize_cantilever(constraints, omega_min, fluid_type = 'water')
# c.print_performance()

# Harley:
# L = 350um, t = 89nm, w = 44um
# t_pr = 30nm, w_pr = 44um, gap is small
# 
# deltaR/R = beta*6*l*pi_l/w/t^2*deltaF = beta*3*E*t*pi_l/2l^2*deltaX
# deflection sensitivity = 50 ppm/micron
# beta = 0.7 (pi_l = 4e-10 for 4e19 doping
# k = 3e-5 N/m
# knee @ 1kHz
# total noise 10Hz - 1kHz = 1.14uV => force resolution = 500 fN
# at 1kHz the force resolution is 8.6 fN/rtHz
# print c.deflection_sensitivity() # expected = 10 V/m
# issues: twisting and too soft to slide on the surface