from cantilever_divingboard import *



# Agressive
# Minimal thickness and width

# Moderative

# Conservative


v_all = linspace(1, 100, 101)
for v_test in v_all:
    c = cantilever_divingboard( l = 80e-6, w = 500e-9, t = 1e-6,
                            l_pr = 10e-6, w_pr = 500e-9, w_gap = 1e-6, t_pr = 100e-9,
                            V_bias = v_test, N = 4e17, freq_min = 1e0, freq_max = 1e2)

    print "Resolution (pN): %f Freq (kHz): %f" % (c.force_resolution(), c.omega_damped_hz()/1000)