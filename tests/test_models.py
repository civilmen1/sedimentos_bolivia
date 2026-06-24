import pytest
import numpy as np
from models.sediment import (
    calculate_water_density,
    calculate_kinematic_viscosity,
    calculate_specific_gravity,
    calculate_fall_velocity,
    calculate_dimensionless_particle_parameter,
    meyer_peter_muller,
    engelund_hansen
)

def test_water_properties():
    # At 20C, rho_w should be approx 998.2
    rho_w = calculate_water_density(20)
    assert pytest.approx(rho_w, rel=1e-3) == 998.2

    # At 15C, nu should be approx 1.14e-6
    nu = calculate_kinematic_viscosity(15)
    assert pytest.approx(nu, rel=1e-3) == 1.14e-6

def test_sediment_parameters():
    rho_s = 2650
    rho_w = 1000
    s = calculate_specific_gravity(rho_s, rho_w)
    assert s == 2.65

    # Test fall velocity for d50=0.45mm
    nu = 1.0e-6
    ws = calculate_fall_velocity(0.45, 2.65, nu)
    assert ws > 0

    dstar = calculate_dimensionless_particle_parameter(0.45, 2.65, nu)
    assert dstar > 0

def test_transport_models():
    s = 2.65
    d50_mm = 0.45
    slope = 0.01
    depth = 1.0
    rho_w = 1000
    v = 1.5

    qb_mpm = meyer_peter_muller(s, d50_mm, slope, depth, rho_w)
    assert qb_mpm >= 0

    qt_eh = engelund_hansen(v, depth, slope, d50_mm, s, rho_w)
    assert qt_eh >= 0
