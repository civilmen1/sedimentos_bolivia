import pytest
from models.sediment import calculate_hydraulic_downscaling

def test_hydraulic_downscaling():
    # Q = 10, B = 5, S = 0.01, n = 0.03
    q = 10.0
    b = 5.0
    s = 0.01
    n = 0.03

    results = calculate_hydraulic_downscaling(q, b, s, n)

    y = results['depth']
    v = results['velocity']

    assert y > 0
    assert v > 0
    # Q = v * A = v * b * y
    assert abs(v * b * y - q) < 1e-3

if __name__ == "__main__":
    test_hydraulic_downscaling()
    print("Test passed!")
