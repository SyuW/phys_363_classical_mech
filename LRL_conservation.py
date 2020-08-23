def main():
    '''
    Show that the Laplace-Runge-Lenz vector is conserved
    '''
    # Masses of m_1 and m_2, m_1 set as the focus
    m_1 = 1
    m_2 = 1
    mu = 0.5     # Reduced mass: (m_1*m_2)/(m_1+m_2)

    # Orbital params
    k = 1        # (G*m_1*m_2)
    e = 0.5       # Eccentricity
    a = 1        # Semi-major axis

    # Angular momentum
    l_z = np.sqrt(a*(1 - e**2)*mu)

    # Sample xi
    xi = np.arange(0, 2*np.pi, 0.01)

    # Trajectory
    t = np.sqrt((mu*a**3)/k)*(xi - e*np.sin(xi))
    r = a*(1 - e*np.cos(xi))
    cosine_phi = (np.cos(xi) - e)/(1 - e*np.cos(xi))
    sine_phi = (np.sqrt(1 - e**2)*np.sin(xi))/(1 - e*np.cos(xi))
    r_hat = np.array(list(zip(cosine_phi, sine_phi, np.zeros(len(sine_phi)))))
    phi_hat = np.array(list(zip(-sine_phi, cosine_phi)))

    phi_dot = l_z/(mu*r**2)

    # Calculate momentum using finite differences
    r_dot = []
    for i, t_i in enumerate(t[:-1]):
        dt = t[i+1] - t[i]
        dr = r[i+1] - r[i]
        r_dot += [dr/dt]
    
    velocity_radial = np.array(list(zip(r_dot*(cosine_phi[:-1]), r_dot*(sine_phi[:-1]), np.zeros(len(r_dot)))))
    velocity_azimuthal = np.array(list(zip(-r*phi_dot*(sine_phi), r*phi_dot*(cosine_phi), np.zeros(len(r_dot)))))

    momentum = m_2 * (velocity_radial + velocity_azimuthal)
    position = np.array(list(zip(r*cosine_phi, r*sine_phi, np.zeros(len(r)))))

    for i, x in enumerate(position):
        print(np.cross(momentum[i], np.cross(position[i], momentum[i])) - m_2*k*r_hat)

main()