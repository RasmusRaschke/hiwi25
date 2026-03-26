### Written by Rasmus Raschke, Universität Hamburg, 2026 ###
## Contains functions returning analytical solutions for different magnetic balls ##

import numpy as np

def magnetic_top(x0: float, y0: float, z0: float, mu: list[float], spin: float, t: list[float], B: float, M: float, R: float) -> float:
    """Calculates analytic trajectory of a magnetic sphere on a flat surface with initial normal spin.

    Parameters
    ----------
    x0 : float
        Initial x position
    y0 : float
        Initial y position
    z0 : float
        Initial z position
    mu : list[float]
        Initial magnetic dipole moment of the sphere
    spin : float
        Initial spin about the surface normal
    t : list[float]
        Time discretization
    B : float
        Magnetic field strength in z direction
    M : float
        Mass of the sphere
    R : float
        Radius of the sphere
    
    Returns
    -------
    x : float
        List containing x projection of trajectory
    y : float
        List containing y projection of trajectory
    z : float
        List containing z projection of trajectory
    """

    mu_norm = np.linalg.norm(mu)
    omega0 = np.sqrt((5*B*mu_norm)/(7*R**2*M))
    delta = np.sqrt(spin**2 + 4*omega0**2)
    etap = 0.5 * (spin + delta)
    etam = 0.5 * (spin - delta)
    w0 = (mu[0] + 1j * mu[1]) / mu_norm
    q = (x0 + 1j * y0) + (w0 * R / delta) * (- etam * (np.exp(1j*etap*t)-1) + etap * (np.exp(1j*etam*t) - 1))
    x = np.real(q)
    y = np.imag(q)
    z = np.full_like(t, z0, dtype=float)
    return x,y,z


def get_beta(system: int, nu: float, mu_norm: float, beta0: float, k: float, t: float, B: float, M: float, R: float):
    """Calculates generalized angle for a magnetic sphere on a flat surface with small magnetic misalignment; normally called by get_pos()

    Parameters
    ----------
    system : int
        1 for harmonic oscillator, 2, 3, and 4 for underdamped, critial, and overdamped harmonic oscillator, and 5 for dry friction
    nu : float
        Rolling resistance coefficient
    mu_norm : float
        Norm of initial magnetic dipole moment
    beta0 : float
        Initial angle between B and mu
    k : float
        Smoothing coefficient in C++ code for tanh
    t : list[float]
        Time discretization
    B : float
        Magnetic field strength in z direction
    M : float
        Mass of the sphere
    R : float
        Radius of the sphere
    
    Returns
    -------
    beta : list[float]
        List containing the angle beta at each t
    """
    if system == 1:
        print("System is harmonic oscillator")
        omega = np.sqrt((5 * mu_norm * B) / (7 * R**2 * M))
        return beta0 * np.cos(omega * t)
    elif system == 2:
        print("System is underdamped harmonic oscillator.")
        omega0 = np.sqrt((5 * mu_norm * B) / (7 * R**2 * M))
        zeta = (5 * nu * g * k) / (14 * omega0) 
        print(omega0)
        print(zeta)
        if zeta > 1:
            print("Warning: omega_0 will be imaginary.")
        omegad = omega0 * np.sqrt(1.0 - zeta**2)
        print(np.exp(-zeta * omega0 * t))
        return beta0 * np.exp(-zeta * omega0 * t) * (np.cos(omegad * t) + (zeta * omega0 / omegad) * np.sin(omegad * t))
    elif system == 3:
        print("System is critically damped hamonic oscillator.")
        omega0 = np.sqrt((5 * mu_norm * B) / (7 * R**2 * M))
        return beta0 * (1 + omega0 * t) * np.exp(- omega0 * t)
    elif system == 4:
        print("System is overdamped harmonic oscillator.")
        omega0 = np.sqrt((5 * mu_norm * B) / (7 * R**2 * M))
        zeta = (5 * nu * g * k) / (14 * omega0) 
        if zeta < 1:
            print("Warning: lambda will be imaginary.")
        lambda1 = -zeta * omega0 + omega0 * np.sqrt(zeta**2 - 1)
        lambda2 = -zeta * omega0 - omega0 * np.sqrt(zeta**2 - 1)
        return (beta0 / (lambda1 - lambda2)) * (-lambda2 * np.exp(lambda1 * t) + lambda1 * np.exp(lambda2 * t))
    elif system == 5:
        print("System is dry harmonic oscillator.")
        omega0 = np.sqrt((5 * mu_norm * B) / (7 * R**2 * M))
        zeta = (5 * nu * g) / (14 * R * omega0)
        T_half = np.pi / omega0
        stick = 2.0 * zeta / omega0
        t_arr = np.atleast_1d(t)
        beta = np.zeros_like(t_arr, dtype=float)
        n = np.floor(t_arr / T_half).astype(int)
        tau = t_arr - n * T_half
        A_n = beta0 - 2.0 * stick * n
        mask_stick = A_n <= stick
        beta[mask_stick] = 0.0
        mask = ~mask_stick
        sign = (-1)**n[mask]
        beta[mask] = sign * (stick + (A_n[mask] - stick) * np.cos(omega0 * tau[mask]))
        return beta if np.ndim(t) else beta[0]
    else:
        raise ValueError("Invalid System!")
            

def get_pos(system: int, x0: float, y0: float, z0: float, vx0: float, vy0: float, vz0: float, mu: list[float], nu: float, k: float, t: list[float], B: float, M: float, R: float) -> float:
    """Calculates analytic trajectory of a magnetic sphere on a flat surface initially slightly misaligned with the outer magnetic field.

    Parameters
    ----------
    system : int
        1 for harmonic oscillator, 2, 3, and 4 for underdamped, critial, and overdamped harmonic oscillator, and 5 for dry friction
    x0 : float
        Initial x position
    y0 : float
        Initial y position
    z0 : float
        Initial z position
    vx0 : float
        Initial x velocity
    vy0 : float
        Initial y velocity
    vz0 : float
        Initial z velocity
    mu : list[float]
        Initial magnetic dipole moment of the sphere
    nu : float
        Rolling resistance coefficient
    k : float
        Smoothing coefficient in C++ code for tanh
    t : list[float]
        Time discretization
    B : float
        Magnetic field strength in z direction
    M : float
        Mass of the sphere
    R : float
        Radius of the sphere
    
    Returns
    -------
    x : float
        List containing x projection of trajectory
    y : float
        List containing y projection of trajectory
    z : float
        List containing z projection of trajectory
    """

    mu_norm = np.linalg.norm(mu)
    beta = get_beta(system, nu, mu_norm, beta0, k, t, B, M, R)
    x = vx0 * t + x0
    y = y0 + R * (beta - beta0)
    z = vz0 * t + z0
    return x, y, z

def exp_envelope(system: int, y0: float, mu_norm: float, beta0: float, nu: float, k: float, t: list[float], B: float, M: float, R: float, g:float) -> float:
    """Calculates the envelope of an underdamped or dry harmonic oscillator

    Parameters
    ----------
    system : int
        1 for underdamped harmonic oscillator, 2 for dry friction
    y0 : float
        Initial y position
    mu_norm : float
        Norm of initial magnetic dipole moment
    beta0 : float
        Initial angle between B and mu
    nu : float
        Rolling resistance coefficient
    k : float
        Smoothing coefficient in C++ code for tanh
    t : list[float]
        Time discretization
    B : float
        Magnetic field strength in z direction
    M : float
        Mass of the sphere
    R : float
        Radius of the sphere
    g : float
        Gravitational acceleration

    Returns
    -------
    upper : list[float]
        List containing maximal envelope for all times t
    lower : list[float]
        List containing minimal envelope for all times t
    """

    omega0 = np.sqrt((5 * mu_norm * B) / (7 * R**2 * M))
    if system == 1:
        zeta = (5 * nu * g * k) / (14 * omega0) 
        A = beta0 / (np.sqrt(1-zeta**2)) * np.exp(- zeta * omega0 * t)
        upper = y0 + R * (A - beta0)
        lower = y0 + R * (-A - beta0)
        return upper, lower
    elif system == 2:
        zeta = (5 * nu * g * k) / (14 * omega0) 
        A = np.maximum(beta0 - (4.0 * zeta / np.pi) * t, 0.0)
        upper = y0 + R * (A - beta0)
        lower = y0 + R * (-A - beta0)
        return upper, lower
    else:
        raise ValueError("Invalid System!")

def incline(x0: float, y0: float, z0: float, vx0: float, vy0: float, vz0: float, phi: float, t: list[float], g: float) -> float:
    """Calculates analytic trajectory of a sphere rolling on an incline.

    Parameters
    ----------
    x0 : float
        Initial x position
    y0 : float
        Initial y position
    z0 : float
        Initial z position
    vx0 : float
        Initial x velocity
    vy0 : float
        Initial y velocity
    vz0 : float
        Initial z velocity
    phi : float
        Incline angle in degrees
    t : list[float]
        Time discretization
    g : float
        Gravitational acceleration

    Returns
    -------
    x : float
        List containing x projection of trajectory
    y : float
        List containing y projection of trajectory
    z : float
        List containing z projection of trajectory
    """

    phi = np.deg2rad(phi)
    x = vx0 * t + x0
    y = (-5*g / 14) * np.sin(phi) * np.cos(phi) * t**2 + vy0 * t + y0
    z = (5*g / 14) * (1-np.cos(phi)) * t**2 + vz0 * t + z0
    return x, y, z