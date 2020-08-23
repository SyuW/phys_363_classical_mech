from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt
import numpy as np

initial_position = (4, 0)
initial_velocity=(0, 0.5001)
g = 1
m = 0.01
delta_t = 0.01

def r_double_dot_func(r, theta_dot):
    return (r*theta_dot**2 - g)/2

def theta_dot_func(r):
    r_i = initial_position[0]
    theta_dot_i = initial_velocity[1]
    return ((theta_dot_i)*(r_i**2))/(r**2)

def main():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    r, theta = initial_position
    r_dot, theta_dot = initial_velocity
    R, THETA = [r], [theta]
    for n in range(1, 15000):
        r_dot += delta_t*r_double_dot_func(r, theta_dot)
        theta += delta_t*theta_dot

        r += delta_t*r_dot
        theta_dot = theta_dot_func(r)
        R += [r]; THETA += [theta]
    
    print(np.array(R)); print(np.array(THETA))
    
    X = np.array(R)*np.cos(THETA); Y = np.array(R)*np.sin(THETA)
    Z = np.array(R)

    ax.plot(X, Y, Z, label=r"($r_i$, $\dot{r}_i$, $\theta_i$, $\dot{\theta}_i$) =" + f"({initial_position[0]}, {initial_velocity[0]}, {initial_position[1]}, {initial_velocity[1]})")
    ax.legend()

    plt.show()


main()