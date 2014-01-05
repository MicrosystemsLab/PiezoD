from cantilever_divingboard import *

def optimize_cantilever(x0):
    c = cantilever_divingboard(*x0)
    return c.force_resolution()
    
# def cantilever_constraints(foo):
#     c = cantilever_divingboard(*foo)
#         
#     min_thickness = 1e-6
#     min_length = 20e-6
#     min_width = 5e-6
#     
#     enforce_positive_dimensions = min(c.l, c.w, c.t, c.l_pr, c.w_pr, c.w_gap, c.t_pr)
#     enforce_minimum_dimensions = min(c.l - min_length, c.w - min_width, c.t - min_thickness)
#     enforce_maximum_piezo_dimensions = min(c.l - c.l_pr, c.w - c.w_pr, c.t - c.t_pr)
#     enforce_positive_resolution = c.force_resolution()
#     
#     min_natural_frequency_hz = 50e3
#     enforce_natural_frequency = c.omega_vacuum_hz() - min_natural_frequency_hz
#     
#     min_constraint = min(enforce_positive_dimensions,
#                 enforce_minimum_dimensions,
#                enforce_maximum_piezo_dimensions,
#                enforce_natural_frequency,
#                enforce_positive_resolution)
#                
#     # print min_constraint
#     return min_constraint
    
def enforce_positive_dimensions(x0):
    c = cantilever_divingboard(*x0)
    return min(c.l, c.w, c.t, c.l_pr, c.w_pr, c.w_gap, c.t_pr)
    
def enforce_minimum_dimensions(x0):
    c = cantilever_divingboard(*x0)
    
    min_thickness = 1e-6
    min_length = 20e-6
    min_width = 500e-9    
    
    return min(c.l - min_length, c.w - min_width, c.t - min_thickness)
    
def enforce_maximum_piezo_dimensions(x0):
    c = cantilever_divingboard(*x0)
    return min(c.l - c.l_pr, c.w - c.w_pr, c.t - c.t_pr)
    
def enforce_minimum_piezo_dimensions(x0):
    c = cantilever_divingboard(*x0)
    
    min_pr_thickness = 30e-9
    min_pr_length = 1e-6
    min_pr_width = 500e-9
    
    return min(c.l_pr - min_pr_length, c.w_pr - min_pr_width, c.t_pr - min_pr_thickness)

def enforce_positive_resolution(x0):
    c = cantilever_divingboard(*x0)
    return c.force_resolution()
    
def enforce_minimum_natural_frequency(x0):
    c = cantilever_divingboard(*x0)
    min_natural_frequency_hz = 50e3
    print 'Force resolution: %f' % c.force_resolution()
    return (c.omega_vacuum_hz() - min_natural_frequency_hz)


c1 = cantilever_divingboard(*x)
print 'Cantilever L/W/T: %f %f %f' % (c1.l*1e6, c1.w*1e6, c1.t*1e6)
print 'Piezo L/W/T: %f %f %f' % (c1.l_pr*1e6, c1.w_pr*1e6, c1.t_pr*1e6)
print 'Resistance: %f' % log10(c1.resistance())
print 'Doping Concentration: %f' % c1.N
print 'Force resolution: %f' % c1.force_resolution()
print 'Natural freq: %f' % c1.omega_vacuum_hz()

# fmin_tnc and blfs won't work, because they don't provide enough constraint flexibility