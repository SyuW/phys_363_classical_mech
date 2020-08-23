import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from mpl_toolkits.mplot3d import Axes3D

'''
Plot the locus of points traced by the spin axis of a spinning top under gravity
'''

parameters_list = [{"a" : 2.0, "b" : 1.0, "alpha" : 1.0, "beta" : 2.0},
                   {"a" : 2.0, "b" : 2.0, "alpha" : 2.0, "beta" : 1.0},
                   {"a" : 5.0, "b" : 2.0, "alpha" : 2.0, "beta" : 3.0},
                   {"a" : 2.0, "b" : 1.0, "alpha" : 1.0, "beta" : 0.0},
                   {"a" : 1.0, "b" : 1.0, "alpha" : 1.0, "beta" : 0.0}]

def u_dot_squared_func(u, a, b, alpha, beta):
    return (1 - u**2)*(alpha - beta*u) - (b - a*u)**2

def u_double_dot_func(u, a, b, alpha, beta):
    return (3/2)*beta*(u**2) - (alpha + a**2)*u + a*b - (beta/2)

def phi_dot_func(u, a, b, alpha, beta):
    return (b - a*u)/(1 - u**2)

def find_roots(params):
    return fsolve(u_dot_squared_func, x0=[-1, 1], args=tuple(params.values()))

def main():
    ## Initial conditions
    delta_t = 0.001
    ## Misc
    current_plot_num = 1

    for parameters in parameters_list:
        fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')

        u = sum(find_roots(parameters)) / 3
        u_dot = np.sqrt(u_dot_squared_func(u, **parameters))
        phi = 0

        U = [u]
        PHI = [phi]

        for n in range(1, 10000):
            u_dot += u_double_dot_func(u, **parameters)*delta_t
            phi += phi_dot_func(u, **parameters)*delta_t

            u += u_dot*delta_t
            U += [u]; PHI += [phi]
        
        print(max(U))
        THETA = np.arccos(U)
        
        fig.suptitle(f"Phys 363 PS5 4({chr(96 + current_plot_num)})"); current_plot_num += 1
        plt.xlabel(r"$\theta(t)$"); plt.ylabel(r"$\phi(t)$", rotation=0)
        plt.plot(THETA, PHI)
        plt.show()

main()
