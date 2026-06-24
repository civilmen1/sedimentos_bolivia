import pytest
from models.sediment import calculate_hydraulic_downscaling

def test_hec_ras_logic():
    # Test Normal Depth iteration
    q = 10.0
    b = 5.0
    s = 0.01
    n = 0.03

    results = calculate_hydraulic_downscaling(q, b, s, n)

    assert 'depth' in results
    assert 'critical_depth' in results
    assert 'froude' in results
    assert 'regime' in results

    # Check Q = V * A
    yn = results['depth']
    v = results['velocity']
    assert abs(v * b * yn - q) < 1e-4

    # Critical depth for Q=10, B=5
    # unit_q = 2.0
    # yc = (2^2 / 9.807)^(1/3) = (0.4078)^(1/3) = 0.7416
    assert abs(results['critical_depth'] - 0.7416) < 0.01

if __name__ == "__main__":
    test_hec_ras_logic()
    print("HEC-RAS logic test passed!")
