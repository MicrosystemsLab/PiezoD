from cvxmod import *
from cvxmod.atoms import norm1
from cvxmod.sets import probsimp
from cantilever_divingboard import *


w = optvar('w', 1)
l = optvar('l', 1)
t = optvar('t', 1)

def frequency(w, l, t):
    return l/(square(t))
    
p = problem(maximize(l+square(t)), [l <= 5, t <= 5, w <= 5])
p.solve()

print "Optimal problem value is %.4f." % p.value
printval(w)


