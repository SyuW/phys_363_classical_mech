import numpy as np
import matplotlib.pyplot as plt


# Constants
a = 1
G = 1
m = 1
M_2 = 1
delta_t = 0.001
mass_ratios = {"i": 1, "ii": 0.041, "iii": 0.039, "iv": 19/20000}

def determine_mass_of_M_1(mass_ratio, M_2):
    return mass_ratio * M_2

def compute_acceleration(velocity, position, M_1, M_2, M):
    x, y = position
    x_dot, y_dot = velocity

    omega = np.sqrt(G * M / a**3)

    helper_x_1 = a*M_2/M - x
    helper_x_2 = a*M_1/M + x
    helper_y_1 = a*M_2/M - x
    helper_y_2 = a*M_1/M + x

    #Terms:         Coriolis          Centrifugal    Gravitational Potential
    x_double_dot = -2*omega*(y_dot) + (omega**2)*x + (G*M_1*helper_x_1)/((helper_x_1**2 + y**2)**(3/2)) - (G*M_2*helper_x_2)/((helper_x_2**2 + y**2)**(3/2))
    y_double_dot =  2*omega*(x_dot) + (omega**2)*y - (G*M_1*y)         /((helper_y_1**2 + y**2)**(3/2)) - (G*M_2*y)         /((helper_y_2**2 + y**2)**(3/2))

    return np.array([x_double_dot, y_double_dot])

def main():
    # Set masses
    for part in mass_ratios:
        M_1 = determine_mass_of_M_1(mass_ratios[part], M_2)
        M = M_1 + M_2
        delta = (M_2 - M_1) / M

        # Lagrange Points
        L_4 = (a / 2) * np.array([delta, np.sqrt(3)])
        L_5 = (a / 2) * np.array([delta, -np.sqrt(3)])

        # Initial conditions
        perturbation = np.array([0, 0.001*a])

        positions  = [L_4 + perturbation]
        velocities = [np.array([0, 0])]
        
        # For comparing nature of eigenfrequencies to observed trajectory
        eigenfreqs_disc = 27*delta**2 - 23
        print(f"Eigenfrequencies discriminant: {eigenfreqs_disc}")
	
	# Solve DE using Euler's Method
        for n in range(1, 100000):
            velocities += [velocities[-1] + delta_t*compute_acceleration(velocities[-1], positions[-1], M_1, M_2, M)]
            positions  += [positions[-1] + delta_t*velocities[-1]]

        x = [coords[0] for coords in positions]; y = [coords[1] for coords in positions]
        plt.plot(x, y)
        plt.xlabel(r"$y_{x}^{\prime}$")
        plt.ylabel(r"$y_{y}^{\prime}$", rotation=0)
        plt.show()

main()
