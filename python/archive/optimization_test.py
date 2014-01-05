from cvxmod import *
from cvxmod.atoms import norm1
from cvxmod.sets import probsimp

A = randn(10, 5)
b = randn(10, 1)
x = optvar('x', 1)
y = optvar('y', 5)

def foo(x):
    return x

p = problem(minimize(foo(x)), [x >= 1])
p.solve()

print "Optimal problem value is %.4f." % p.value
printval(x)
