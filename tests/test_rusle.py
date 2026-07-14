import math
import pytest
from models import rusle as rs


def test_r_factor_formulas():
    Pa = 500.0
    assert rs.r_factor(Pa, "hurni") == pytest.approx(-8.12 + 0.562 * Pa, rel=1e-6)
    assert rs.r_factor(Pa, "lineal") == pytest.approx(79 + 0.363 * Pa, rel=1e-6)
    assert rs.r_factor(Pa, "potencia") == pytest.approx(0.0483 * Pa ** 1.610, rel=1e-6)
    # nunca negativo
    assert rs.r_factor(1.0, "hurni") >= 0.0


def test_k_factor_textura():
    assert rs.k_factor_textura(12) == pytest.approx(0.0053, rel=1e-9)
    assert rs.k_factor_textura(9) == pytest.approx(0.0500, rel=1e-9)
    assert rs.k_factor_textura(99) == 0.0   # clase fuera de rango


def test_ls_factor_mccool_creciente():
    # LS crece con la pendiente y da valores físicamente sanos
    ls5 = rs.ls_factor(5)
    ls20 = rs.ls_factor(20)
    assert 0 < ls5 < 2
    assert ls20 > ls5
    # con longitud de ladera > 0, L amplifica
    assert rs.ls_factor(10, slope_length_m=100) > rs.ls_factor(10)


def test_c_factor_decrece_con_vegetacion():
    assert rs.c_factor(0.2) > rs.c_factor(0.8)
    assert 0.0 <= rs.c_factor(0.9) <= 1.0
    # NDVI cercano a 1 no rompe (singularidad acotada)
    assert rs.c_factor(0.999) >= 0.0


def test_p_factor():
    assert rs.p_factor(10, es_cultivo=False) == 1.0
    assert rs.p_factor(1, es_cultivo=True) == 0.6
    assert rs.p_factor(25, es_cultivo=True) == 0.9


def test_compute_rusle_desglose():
    r = rs.compute_rusle(precip_anual_mm=500, clase_textura=6, slope_pct=15,
                         ndvi=0.3, es_cultivo=True, r_formula="hurni")
    for k in ("R", "K", "LS", "C", "P", "A", "severidad", "color"):
        assert k in r
    assert r["A"] == pytest.approx(r["R"] * r["K"] * r["LS"] * r["C"] * r["P"], rel=1e-2)
    assert r["A"] >= 0


def test_clase_severidad():
    assert rs.clase_severidad(3)[0] == "Ligera"
    assert rs.clase_severidad(50)[0] == "Severa"
    assert rs.clase_severidad(100)[0] == "Crítica"
