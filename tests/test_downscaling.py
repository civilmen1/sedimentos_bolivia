import pytest
from models.sediment import calculate_hydraulic_downscaling

def test_hydraulic_downscaling():
    # Q = 10, B = 5, S = 0.01, n = 0.03
    # Manning's: Q = (1/n) * A * R^(2/3) * S^(1/2)
    # If y = 1.0: A = 5, P = 7, R = 5/7 = 0.714
    # Q = (1/0.03) * 5 * (0.714)**(2/3) * (0.01)**0.5
    # Q = 33.33 * 5 * 0.798 * 0.1 = 13.3
    # So for Q=10, y should be less than 1.0

    q = 10.0
    b = 5.0
    s = 0.01
    n = 0.03

    y, v = calculate_hydraulic_downscaling(q, b, s, n)

    assert y > 0
    assert v > 0
    # Q = v * A = v * b * y
    assert abs(v * b * y - q) < 1e-3

if __name__ == "__main__":
    test_hydraulic_downscaling()
    print("Test passed!")
