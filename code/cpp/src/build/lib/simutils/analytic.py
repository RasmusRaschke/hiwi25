### Written by Rasmus Raschke, Universität Hamburg, 2026 ###
## Contains functions returning analytical solutions for different magnetic balls ##

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import cumulative_trapezoid, solve_ivp
from scipy.special import ellipj, ellipk
from scipy.spatial.transform import Rotation
from scipy.spatial import cKDTree
import mpmath as mp

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
    z = -(5*g / 14) * (1-np.cos(phi)**2) * t**2 + vz0 * t + z0
    return x, y, z

### Quartic Oscillator ###

def reduced_constants(M: float, R: float, I: float, mu_norm: float, B: float, spin: float, C: float, E: float, L: float):
    """Returns the constants of the exact quartic reduction.

    The exact normal-field rolling problem uses
        J = M*R**2 + I,
        kappa = mu_norm*B/J,
        a = spin**2/4 - C/2,
    and physical solutions satisfy
        2*E + spin*L + C**2/4 = kappa**2.

    Returns
    -------
    constants : dict
        Dictionary containing J, kappa, a, E, L, C, spin,
        physical_constraint_residual, and omega_xi.
    """

    J = M * R**2 + I
    kappa = mu_norm * B / J
    a = spin**2 / 4.0 - C / 2.0
    residual = 2.0 * E + spin * L + C**2 / 4.0 - kappa**2
    omega_xi = -(L + 0.5 * spin * C) / kappa
    return {
        "J": J,
        "kappa": kappa,
        "a": a,
        "E": E,
        "L": L,
        "C": C,
        "spin": spin,
        "physical_constraint_residual": residual,
        "omega_xi": omega_xi,
    }


def quartic_roots(a: float, E: float, L: float):
    """Returns the three roots of the radial cubic P3(q).

    The squared radius q=|y|**2 obeys
        q_dot**2 = -q**3 - 4*a*q**2 + 8*E*q - 4*L**2.
    """

    roots = np.roots([-1.0, -4.0 * a, 8.0 * E, -4.0 * L**2])
    if np.max(np.abs(np.imag(roots))) > 1e-9:
        raise ValueError("The radial cubic does not have three real roots.")
    return np.sort(np.real(roots))


def _elliptic_third_kind_accumulated(u, characteristic: float, parameter: float):
    """Calculates integral_0^u dv/(1-n*sn(v|m)^2) with period continuation."""

    u = np.asarray(u, dtype=float)
    K = float(ellipk(parameter))
    period = 2.0 * K
    complete = float(mp.ellippi(characteristic, parameter))
    periods = np.floor(u / period).astype(int)
    remainder = u - periods * period
    amplitude = ellipj(remainder, parameter)[3]
    incomplete = np.array([
        float(mp.ellippi(characteristic, float(phi), parameter))
        for phi in np.atleast_1d(amplitude)
    ])
    return periods * (2.0 * complete) + incomplete


def _complex_integral_from_first_time(t, values):
    """Cumulative trapezoidal integral of complex samples starting at t[0]."""



    t = np.asarray(t, dtype=float)
    values = np.asarray(values, dtype=complex)
    if t.ndim != 1 or values.shape != t.shape:
        raise ValueError("t and values must be one-dimensional arrays of equal length.")
    if np.any(np.diff(t) <= 0.0):
        raise ValueError("t must be strictly increasing.")
    return cumulative_trapezoid(values, t, initial=0.0)


def _quat_multiply(q1, q2):
    """Quaternion Product"""

    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return np.array([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2,
    ], dtype=float)


def _quat_from_axis_direction(xi, body_phase: float = 0.0):
    """Quaternion mapping body b_z onto xi, followed by a body-z phase."""

    xi = np.asarray(xi, dtype=float)
    xi = xi / np.linalg.norm(xi)
    if xi[2] > -1.0 + 1e-12:
        w = np.sqrt(0.5 * (1.0 + xi[2]))
        q_align = np.array([
            w,
            -xi[1] / (2.0 * w),
            xi[0] / (2.0 * w),
            0.0,
        ])
    else:
        q_align = np.array([0.0, 1.0, 0.0, 0.0])

    q_phase = np.array([
        np.cos(0.5 * body_phase),
        0.0,
        0.0,
        np.sin(0.5 * body_phase),
    ])
    q = _quat_multiply(q_align, q_phase)
    return q / np.linalg.norm(q)


def _reconstruct_quaternion(t, omega_world, xi0, body_phase0: float = 0.0):
    """Reconstructs attitude by quadrature."""

    t = np.asarray(t, dtype=float)
    omega_world = np.asarray(omega_world, dtype=float)
    q0 = _quat_from_axis_direction(xi0, body_phase0)
    splines = [CubicSpline(t, omega_world[:, j]) for j in range(3)]

    def rhs(time, q):
        omega = np.array([splines[j](time) for j in range(3)])
        pure = np.array([0.0, omega[0], omega[1], omega[2]])
        return 0.5 * _quat_multiply(pure, q)

    solution = solve_ivp(
        rhs,
        (float(t[0]), float(t[-1])),
        q0,
        t_eval=t,
        rtol=2e-11,
        atol=2e-13,
    )
    quaternion = solution.y.T
    quaternion /= np.linalg.norm(quaternion, axis=1)[:, None]
    return unwrap_quaternion_sign(quaternion)


def magnetic_quartic_motion(
    spin: float,
    C: float,
    E: float,
    L: float,
    t: list[float],
    B: float,
    M: float,
    R: float,
    mu_norm: float,
    I: float,
    x0: float = 0.0,
    y0: float = 0.0,
    z0: float | None = None,
    theta0: float = 0.0,
    radial_phase_time: float = 0.0,
    body_phase0: float = 0.0,
):
    """Calculates the full bounded exact motion of the magnetic rolling sphere.

    Returns
    -------
    state : dict
        Arrays for x,y,z,vx,vy,vz,qw,qx,qy,qz, omega_world,
        xi, w, y_reduced, y_dot_reduced, radial_q, radial_q_dot,
        theta, normal_spin and the reduction constants.
    """

    t = np.asarray(t, dtype=float)
    if t.ndim != 1 or len(t) < 2 or np.any(np.diff(t) <= 0.0):
        raise ValueError("t must be a strictly increasing one-dimensional array.")

    constants = reduced_constants(M, R, I, mu_norm, B, spin, C, E, L)
    kappa = constants["kappa"]
    a = constants["a"]
    if abs(constants["physical_constraint_residual"]) > 1e-7 * max(1.0, kappa**2):
        raise ValueError("E,L,C,kappa,spin do not satisfy the physical unit-vector constraint.")

    q1, q2, q3 = quartic_roots(a, E, L)
    if not (q1 < q2 < q3 and q2 > 0.0):
        raise ValueError("This function implements the bounded branch q2 <= q <= q3 with q2>0.")

    omega_e = 0.5 * np.sqrt(q3 - q1)
    parameter = (q3 - q2) / (q3 - q1)
    characteristic = (q3 - q2) / q3
    u = omega_e * (t - radial_phase_time)
    sn, cn, dn, _ = ellipj(u, parameter)
    radial_q = q3 - (q3 - q2) * sn**2
    radial_q_dot = -2.0 * omega_e * (q3 - q2) * sn * cn * dn
    rho = np.sqrt(radial_q)
    rho_dot = radial_q_dot / (2.0 * rho)

    Pi_acc = _elliptic_third_kind_accumulated(u, characteristic, parameter)
    u_ref = omega_e * (t[0] - radial_phase_time)
    Pi_ref = float(_elliptic_third_kind_accumulated(np.array([u_ref]), characteristic, parameter)[0])
    theta = theta0 + (L / (q3 * omega_e)) * (Pi_acc - Pi_ref)
    y_reduced = rho * np.exp(1j * theta)
    y_dot_reduced = (rho_dot + 1j * L / rho) * np.exp(1j * theta)

    rotating_phase = 0.5 * spin * (t - t[0])
    phase_factor = np.exp(1j * rotating_phase)
    w = phase_factor * y_reduced
    w_dot = phase_factor * (y_dot_reduced + 0.5j * spin * y_reduced)
    xi_complex = 1j * w_dot / kappa
    xi3 = (radial_q - C) / (2.0 * kappa)
    xi = np.column_stack((np.real(xi_complex), np.imag(xi_complex), xi3))

    xi_norm_error = np.max(np.abs(np.linalg.norm(xi, axis=1) - 1.0))
    if xi_norm_error > 5e-7:
        raise RuntimeError(f"Reconstructed magnetic axis is not unit length; max error={xi_norm_error:.3e}.")

    omega_world = np.column_stack((np.real(w), np.imag(w), np.full_like(t, spin)))
    velocity_complex = -1j * R * w
    displacement = _complex_integral_from_first_time(t, velocity_complex)
    position_complex = (x0 + 1j * y0) + displacement
    z_value = R if z0 is None else z0

    quaternion = _reconstruct_quaternion(t, omega_world, xi[0], body_phase0)

    return {
        "t": t,
        "x": np.real(position_complex),
        "y": np.imag(position_complex),
        "z": np.full_like(t, z_value),
        "vx": np.real(velocity_complex),
        "vy": np.imag(velocity_complex),
        "vz": np.zeros_like(t),
        "qw": quaternion[:, 0],
        "qx": quaternion[:, 1],
        "qy": quaternion[:, 2],
        "qz": quaternion[:, 3],
        "quaternion": quaternion,
        "omega_world": omega_world,
        "normal_spin": omega_world[:, 2],
        "xi": xi,
        "w": w,
        "y_reduced": y_reduced,
        "y_dot_reduced": y_dot_reduced,
        "radial_q": radial_q,
        "radial_q_dot": radial_q_dot,
        "theta": theta,
        "roots": np.array([q1, q2, q3]),
        "constants": constants,
    }


def magnetic_homoclinic_motion(
    spin: float,
    t: list[float],
    B: float,
    M: float,
    R: float,
    mu_norm: float,
    I: float,
    centre_time: float,
    x0: float = 0.0,
    y0: float = 0.0,
    z0: float | None = None,
    phase0: float = 0.0,
    body_phase0: float = 0.0,
):
    """Calculates the complete homoclinic motion including path and quaternion.

    The exact reduced orbit is
        w = 2*lambda*sech(lambda*(t-centre_time))
            * exp(i*(phase0 + spin*(t-centre_time)/2)),
    where lambda**2 = kappa - spin**2/4.
    """

    t = np.asarray(t, dtype=float)
    if t.ndim != 1 or len(t) < 2 or np.any(np.diff(t) <= 0.0):
        raise ValueError("t must be a strictly increasing one-dimensional array.")

    J = M * R**2 + I
    kappa = mu_norm * B / J
    lambda2 = kappa - spin**2 / 4.0
    if lambda2 <= 0.0:
        raise ValueError("Homoclinic motion requires kappa > spin**2/4.")
    lambda0 = np.sqrt(lambda2)
    tau = t - centre_time
    sech = 1.0 / np.cosh(lambda0 * tau)
    rho = 2.0 * lambda0 * sech
    rho_dot = -2.0 * lambda0**2 * sech * np.tanh(lambda0 * tau)
    phase = phase0 + 0.5 * spin * tau
    w = rho * np.exp(1j * phase)
    w_dot = (rho_dot + 0.5j * spin * rho) * np.exp(1j * phase)
    xi_complex = 1j * w_dot / kappa
    xi3 = -1.0 + 2.0 * lambda0**2 * sech**2 / kappa
    xi = np.column_stack((np.real(xi_complex), np.imag(xi_complex), xi3))
    omega_world = np.column_stack((np.real(w), np.imag(w), np.full_like(t, spin)))

    velocity_complex = -1j * R * w
    displacement = _complex_integral_from_first_time(t, velocity_complex)
    position_complex = (x0 + 1j * y0) + displacement
    z_value = R if z0 is None else z0
    quaternion = _reconstruct_quaternion(t, omega_world, xi[0], body_phase0)

    return {
        "t": t,
        "x": np.real(position_complex),
        "y": np.imag(position_complex),
        "z": np.full_like(t, z_value),
        "vx": np.real(velocity_complex),
        "vy": np.imag(velocity_complex),
        "vz": np.zeros_like(t),
        "qw": quaternion[:, 0],
        "qx": quaternion[:, 1],
        "qy": quaternion[:, 2],
        "qz": quaternion[:, 3],
        "quaternion": quaternion,
        "omega_world": omega_world,
        "normal_spin": omega_world[:, 2],
        "xi": xi,
        "w": w,
        "lambda": lambda0,
        "kappa": kappa,
        "C": 2.0 * kappa,
        "E": 0.0,
        "L": 0.0,
        "a": -lambda0**2,
    }


def magnetic_weierstrass(a: float, E: float, L: float, t: list[float], theta0: float = 0.0, t0: float = 0.0):
    """Reduced bounded elliptic solution y(t), y_dot(t), q(t)."""

    from scipy.special import ellipj

    t = np.asarray(t, dtype=float)
    q1, q2, q3 = quartic_roots(a, E, L)
    if not (q1 < q2 < q3 and q2 > 0.0):
        raise ValueError("Expected q1 < q2 < q3 with q2>0 for the bounded branch.")
    omega_e = 0.5 * np.sqrt(q3 - q1)
    parameter = (q3 - q2) / (q3 - q1)
    characteristic = (q3 - q2) / q3
    u = omega_e * (t - t0)
    sn, cn, dn, _ = ellipj(u, parameter)
    q = q3 - (q3 - q2) * sn**2
    q_dot = -2.0 * omega_e * (q3 - q2) * sn * cn * dn
    rho = np.sqrt(q)
    rho_dot = q_dot / (2.0 * rho)
    Pi_acc = _elliptic_third_kind_accumulated(u, characteristic, parameter)
    theta = theta0 + (L / (q3 * omega_e)) * Pi_acc
    y = rho * np.exp(1j * theta)
    y_dot = (rho_dot + 1j * L / rho) * np.exp(1j * theta)
    return y, y_dot, q


def unwrap_quaternion_sign(quaternion):
    """Removes the q <-> -q sign ambiguity from a quaternion time series."""

    q = np.asarray(quaternion, dtype=float).copy()
    if q.ndim == 1:
        q = q[None, :]
    q /= np.linalg.norm(q, axis=1)[:, None]
    for i in range(1, len(q)):
        if np.dot(q[i], q[i-1]) < 0.0:
            q[i] *= -1.0
    return q


def world_omega_from_output(qw, qx, qy, qz, Ox, Oy, Oz):
    """Converts the C++ body angular velocity Omega to world angular velocity."""

    q = unwrap_quaternion_sign(np.column_stack((qw, qx, qy, qz)))
    rotation = Rotation.from_quat(np.column_stack((q[:, 1], q[:, 2], q[:, 3], q[:, 0])))
    Omega = np.column_stack((Ox, Oy, Oz))
    return rotation.apply(Omega)


def magnetic_axis_from_quaternion(qw, qx, qy, qz):
    """Returns xi=P_t*b_z from a C++ quaternion time series."""

    q = unwrap_quaternion_sign(np.column_stack((qw, qx, qy, qz)))
    rotation = Rotation.from_quat(np.column_stack((q[:, 1], q[:, 2], q[:, 3], q[:, 0])))
    return rotation.apply(np.tile(np.array([0.0, 0.0, 1.0]), (len(q), 1)))


def reduced_state_from_output(data, drive_frequency: float | None = None):
    """Builds a bounded reduced state vector from a solver output structured array."""

    names = data.dtype.names
    omega = world_omega_from_output(
        data["qw"], data["qx"], data["qy"], data["qz"],
        data["Ox"], data["Oy"], data["Oz"],
    )
    xi = magnetic_axis_from_quaternion(data["qw"], data["qx"], data["qy"], data["qz"])
    state = np.column_stack((omega[:, 0], omega[:, 1], xi))
    if drive_frequency is not None:
        phase = drive_frequency * data["t"]
        state = np.column_stack((state, np.cos(phase), np.sin(phase)))
    return state


def largest_lyapunov_rosenstein(
    t,
    state,
    theiler_time: float | None = None,
    max_divergence_time: float | None = None,
    fit_window: tuple[float, float] | None = None,
    max_points: int = 12000,
    neighbor_candidates: int = 64,
):
    """Estimates the largest Lyapunov exponent from one sampled trajectory.
    
    Returns
    -------
    result : dict
        lag_time, mean_log_divergence, fit_mask, exponent, intercept,
        pair_count and sampling_stride.
    """

    t = np.asarray(t, dtype=float)
    X = np.asarray(state, dtype=float)
    if X.ndim == 1:
        X = X[:, None]
    if len(t) != len(X) or len(t) < 100:
        raise ValueError("Need at least 100 state samples with the same length as t.")

    stride = max(1, int(np.ceil(len(t) / max_points)))
    t = t[::stride]
    X = X[::stride]
    dt = float(np.median(np.diff(t)))

    scale = np.std(X, axis=0)
    scale[scale < 1e-14] = 1.0
    Xn = (X - np.mean(X, axis=0)) / scale

    if theiler_time is None:
        theiler_time = max(10.0 * dt, 0.05 * (t[-1] - t[0]))
    theiler = max(1, int(round(theiler_time / dt)))

    if max_divergence_time is None:
        max_divergence_time = min(0.15 * (t[-1] - t[0]), 200.0 * dt)
    max_lag = max(2, int(round(max_divergence_time / dt)))
    usable = len(t) - max_lag
    if usable <= 10:
        raise ValueError("Trajectory is too short for the requested divergence horizon.")

    tree = cKDTree(Xn[:usable])
    k_query = min(neighbor_candidates + 1, usable)
    _, candidates = tree.query(Xn[:usable], k=k_query)
    if candidates.ndim == 1:
        candidates = candidates[:, None]

    pairs = []
    for i in range(usable):
        for j in np.atleast_1d(candidates[i]):
            j = int(j)
            if j == i:
                continue
            if abs(j - i) <= theiler:
                continue
            if j + max_lag >= len(Xn):
                continue
            pairs.append((i, j))
            break

    if len(pairs) < 20:
        raise RuntimeError("Could not find enough temporally separated neighbour pairs.")

    mean_log = np.empty(max_lag + 1)
    for lag in range(max_lag + 1):
        distances = np.array([
            np.linalg.norm(Xn[i + lag] - Xn[j + lag])
            for i, j in pairs
        ])
        distances = distances[distances > 1e-14]
        mean_log[lag] = np.mean(np.log(distances)) if len(distances) else np.nan

    lag_time = np.arange(max_lag + 1) * dt
    finite = np.isfinite(mean_log)
    if fit_window is None:
        fit_start = 0.1 * max_divergence_time
        fit_end = 0.45 * max_divergence_time
    else:
        fit_start, fit_end = fit_window
    fit_mask = finite & (lag_time >= fit_start) & (lag_time <= fit_end)
    if np.sum(fit_mask) < 3:
        raise RuntimeError("Lyapunov fit window contains too few points.")
    exponent, intercept = np.polyfit(lag_time[fit_mask], mean_log[fit_mask], 1)

    return {
        "lag_time": lag_time,
        "mean_log_divergence": mean_log,
        "fit_mask": fit_mask,
        "exponent": float(exponent),
        "intercept": float(intercept),
        "pair_count": len(pairs),
        "sampling_stride": stride,
    }