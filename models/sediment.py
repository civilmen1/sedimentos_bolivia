import numpy as np
from scipy.optimize import fsolve

def calculate_water_density(temp_c):
    """
    Calculates water density (rho_w) in kg/m3 based on temperature in Celsius.
    Using Eq 1.1.2 from Ramirez Quispe (2021).
    """
    rho_w = 1000 * (1 - ((temp_c + 288.941) * (temp_c - 3.986)**2) / (508929.2 * (temp_c + 68.13)))
    return rho_w

def calculate_kinematic_viscosity(temp_c):
    """
    Calculates kinematic viscosity (nu) in m2/s based on temperature in Celsius.
    Using Eq 1.1.3 from Ramirez Quispe (2021).
    """
    nu = (1.14 - 0.031 * (temp_c - 15) + 0.00068 * (temp_c - 15)**2) * 1e-6
    return nu

def calculate_specific_gravity(rho_s, rho_w):
    """
    Calculates specific gravity (s) = rho_s / rho_w.
    """
    return rho_s / rho_w

def calculate_fall_velocity(d50_mm, s, nu, g=9.807):
    """
    Calculates fall velocity (ws) in m/s using Van Rijn (1993) - Cap 1.2.11.
    """
    d50_m = d50_mm / 1000.0
    if d50_mm <= 0.1:
        ws = ((s - 1) * g * d50_m**2) / (18 * nu)
    elif d50_mm <= 1.0:
        ws = (10 * nu / d50_m) * (np.sqrt(1 + (0.01 * (s - 1) * g * d50_m**3) / nu**2) - 1)
    else:
        ws = 1.1 * np.sqrt((s - 1) * g * d50_m)
    return ws

def calculate_dimensionless_particle_parameter(d50_mm, s, nu, g=9.807):
    """
    Calculates dimensionless particle parameter (D*) - Cap 2.4.
    """
    d50_m = d50_mm / 1000.0
    d_star = (((s - 1) * g) / nu**2)**(1/3) * d50_m
    return d_star

def calculate_critical_shear_stress_shields(dstar, s, rho_w, d50_mm, g=9.807):
    """
    Calculates critical shear stress (tau_c) in Pa using Shields curve approximation.
    """
    d50_m = d50_mm / 1000.0
    # Approximation of Shields parameter (theta_c)
    if dstar <= 4:
        theta_c = 0.24 / dstar
    elif dstar <= 10:
        theta_c = 0.14 * dstar**(-0.64)
    elif dstar <= 20:
        theta_c = 0.04 * dstar**(-0.1)
    elif dstar <= 150:
        theta_c = 0.013 * dstar**0.29
    else:
        theta_c = 0.055

    tau_c = theta_c * (s - 1) * rho_w * g * d50_m
    return tau_c, theta_c

def meyer_peter_muller(s, d50_mm, slope, depth, rho_w, g=9.807):
    """
    Meyer-Peter & Müller (1948) for bed load.
    Returns qb in kg/(m*s).
    """
    d50_m = d50_mm / 1000.0
    tau = rho_w * g * depth * slope
    theta = tau / ((s - 1) * rho_w * g * d50_m)

    theta_c = 0.047
    if theta > theta_c:
        phi = 8 * (theta - theta_c)**1.5
        qb_vol = phi * np.sqrt((s - 1) * g * d50_m**3) # m3/(m*s)
        qb_mass = qb_vol * s * rho_w # kg/(m*s)
        return qb_mass
    return 0.0

def engelund_hansen(v, depth, slope, d50_mm, s, rho_w, g=9.807):
    """
    Engelund-Hansen (1967) for total load.
    Returns qt in kg/(m*s).
    """
    d50_m = d50_mm / 1000.0
    tau = rho_w * g * depth * slope
    theta = tau / ((s - 1) * rho_w * g * d50_m)
    f = 2 * g * depth * slope / v**2 # friction factor

    phi = 0.1 * theta**2.5 / f
    qt_vol = phi * np.sqrt((s - 1) * g * d50_m**3)
    qt_mass = qt_vol * s * rho_w
    return qt_mass

def van_rijn_bedload(v, depth, d50_mm, dstar, s, rho_w, nu, g=9.807):
    """
    Van Rijn (1984) for bed load.
    """
    d50_m = d50_mm / 1000.0
    # Critical velocity
    if 1 < dstar <= 10:
        v_c = 0.19 * d50_m**0.1 * np.log10(12 * depth / (3 * d50_m))
    elif dstar > 10:
        v_c = 0.19 * d50_m**0.1 * np.log10(12 * depth / (3 * d50_m)) # simplified for this task
    else:
        v_c = 0.0

    if v > v_c:
        t = (v**2 - v_c**2) / v_c**2
        if t < 0: t = 0
        qb_vol = 0.053 * ((s - 1) * g)**0.5 * d50_m**1.5 * t**2.1 / dstar**0.3
        return qb_vol * s * rho_w
    return 0.0

def solve_manning_y(q, b, s, n):
    """
    Solves for normal depth (y) in a rectangular channel using Manning's equation.
    Q = (1/n) * A * R^(2/3) * S^(1/2)
    """
    if s <= 0:
        return 0.1 # Fallback

    def func(y):
        if y <= 0: return -q
        a = b * y
        p = b + 2 * y
        r = a / p
        return (1/n) * a * (r**(2/3)) * (s**0.5) - q

    # Initial guess: assume wide channel (R approx y)
    # q = (1/n) * b * y * y^(2/3) * s^0.5 => y^(5/3) = (q * n) / (b * s^0.5)
    y_guess = ((q * n) / (b * s**0.5))**(3/5)

    y_sol = fsolve(func, y_guess)
    return float(y_sol[0])

def calculate_hydraulic_downscaling(q, b, s, n):
    """
    Downscaling analysis to calculate y and v from Q, B, S, n.
    """
    y = solve_manning_y(q, b, s, n)
    v = q / (b * y) if y > 0 else 0
    return y, v
