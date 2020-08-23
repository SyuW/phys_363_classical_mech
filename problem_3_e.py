import numpy as np

import matplotlib.pyplot as plt

M_2 = 1
M_1 = 2*M_2
a = 1
G = 1
m = 1

angular_freq = np.sqrt(G*(M_1+M_2)*m/(a**3))
fig, ax = plt.subplots(1,1)

limit = 5
x_list = np.linspace(-limit, limit, 1000)
y_list = np.linspace(-limit, limit, 1000)

X, Y = np.meshgrid(x_list, y_list)
Z = -(1/2)*m*(angular_freq**2)*(X**2 + Y**2) - (G*M_1*m)/(np.sqrt((a*M_2/(M_1 + M_2) - X)**2 + Y**2)) - (G*M_2*m)/(np.sqrt((a*M_1/(M_1 + M_2) + X)**2 + Y**2))

cp = ax.contour(X, Y, Z, levels=np.linspace(-30, 30, 200))
plt.xlabel(r"${y_x}'$")
plt.ylabel(r"${y_y}'$")
plt.show()


