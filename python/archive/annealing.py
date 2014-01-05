from cantilever_divingboard import *

def optimize_cantilever(x0):
    c = cantilever_divingboard(*x0)
    print x0
    return c.force_resolution()

x = optimize.anneal(func = optimize_cantilever,
                    x0 = ( 150e-6, 1e-6, 5e-6, 
                            10e-6, 1e-6, 1e-6, 1e-6, 3., 1e15))
                            
c1 = cantilever_divingboard(*x)
print 'Cantilever L/W/T: %f %f %f' % (c1.l*1e6, c1.w*1e6, c1.t*1e6)
print 'Piezo L/W/T: %f %f %f' % (c1.l_pr*1e6, c1.w_pr*1e6, c1.t_pr*1e6)
print 'Resistance: %f' % log10(c1.resistance())
print 'Doping Concentration: %f' % c1.N
print 'Force resolution: %f' % c1.force_resolution()
print 'Natural freq: %f' % c1.omega_vacuum_hz()