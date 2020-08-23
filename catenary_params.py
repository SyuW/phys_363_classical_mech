from scipy.optimize import fsolve
import math

'''
For determining parameters for catenary curve
'''

def equations(p):
    a, b = p
    return (a*math.cosh(b/a) - 1, a*math.cosh((1-b)/a) - 1)

x, y = fsolve(equations, (0.1, 1))

print(x)
print(y)

print(equations((x, y)))