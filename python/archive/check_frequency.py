from cantilever_divingboard import *

t = 1e-6
w = 1e-6
l = 400e-6

l_pr = 1 # estimated from SEM in paper
w_pr = 1
t_pr = 1
w_gap = 1

freq_min = 1e1
freq_max = 1e5

V_bias = 5
N = 4e19
beta = 0.7

c = cantilever_divingboard(freq_min, freq_max, 
                           l, w, t, 
                           l_pr, w_pr, w_gap, t_pr, 
                           V_bias, N)

# For L = 400e-6
c.l = 400e-6

c.w = 0.1e-6
print c.omega_damped_hz()

c.w = 0.5e-6
print c.omega_damped_hz()

c.w = 1e-6
print c.omega_damped_hz()

c.w = 5e-6
print c.omega_damped_hz()

c.w = 10e-6
print c.omega_damped_hz()


# For L = 100e-6
c.l = 100e-6

c.w = 0.1e-6
print c.omega_damped_hz()

c.w = 0.5e-6
print c.omega_damped_hz()

c.w = 1e-6
print c.omega_damped_hz()

c.w = 5e-6
print c.omega_damped_hz()

c.w = 10e-6
print c.omega_damped_hz()

