"""
Modelo RUSLE (Revised Universal Soil Loss Equation) — pérdida de suelo.

    A = R · K · LS · C · P     [t·ha⁻¹·año⁻¹]

Funciones puras (números → números). La orquestación con Google Earth Engine
(obtención de Pₐ, textura, pendiente, NDVI, cobertura en el punto) vive en
utils/gee_handler.compute_rusle_point, que llama a estas funciones.

Referencias: Wischmeier & Smith (1978); Renard et al. RUSLE (1997);
van der Knijff (C); McCool (factor S); Sharpley & Williams EPIC (K);
Hurni (1985) y fórmulas empíricas regionales para R.
"""

import math

# K por clase de textura USDA (OpenLandMap USDA texture class 1..12).
# Valores del script RUSLE-GEE de referencia (t·ha·h·ha⁻¹·MJ⁻¹·mm⁻¹).
_K_TEXTURE = {
    1: 0.0288, 2: 0.0341, 3: 0.0360, 4: 0.0394, 5: 0.0423, 6: 0.0264,
    7: 0.0394, 8: 0.0499, 9: 0.0500, 10: 0.0450, 11: 0.0170, 12: 0.0053,
}

R_FORMULAS = ("hurni", "lineal", "potencia")


def r_factor(precip_anual_mm, formula="hurni"):
    """
    Erosividad de la lluvia R [MJ·mm·ha⁻¹·h⁻¹·año⁻¹] a partir de la
    precipitación media anual Pₐ (mm).
        hurni    : R = −8.12 + 0.562·Pₐ   (montaña árida/semiárida — Bolivia)
        lineal   : R = 79 + 0.363·Pₐ      (templado/subtropical)
        potencia : R = 0.0483·Pₐ^1.610    (intensidades extremas)
    """
    Pa = float(precip_anual_mm)
    if formula == "lineal":
        R = 79.0 + 0.363 * Pa
    elif formula == "potencia":
        R = 0.0483 * Pa ** 1.610
    else:  # hurni (default)
        R = -8.12 + 0.562 * Pa
    return max(R, 0.0)


def k_factor_textura(clase_textura):
    """Erodibilidad K según la clase de textura USDA (1..12)."""
    try:
        tc = int(round(float(clase_textura)))
    except (TypeError, ValueError):
        return 0.0
    return _K_TEXTURE.get(tc, 0.0)


def ls_factor(slope_pct, slope_length_m=None):
    """
    Factor topográfico LS. Se usa el factor de pendiente S de McCool (RUSLE),
    dimensionalmente sano, con L≈1 para una estimación puntual (sin longitud de
    ladera se asume una unidad de referencia de 22.13 m):
        θ = atan(slope_pct/100)
        S = 10.8·senθ + 0.03      si slope < 9 %
        S = 16.8·senθ − 0.50      si slope ≥ 9 %
        L = (λ/22.13)^m           con m según la pendiente (si se da λ)
    """
    s = max(float(slope_pct), 0.0)
    theta = math.atan(s / 100.0)
    sin_t = math.sin(theta)
    if s < 9.0:
        S = 10.8 * sin_t + 0.03
    else:
        S = 16.8 * sin_t - 0.50
    S = max(S, 0.0)
    if slope_length_m and slope_length_m > 0:
        beta = (sin_t / 0.0896) / (3.0 * sin_t ** 0.8 + 0.56)
        m = beta / (1.0 + beta)
        L = (slope_length_m / 22.13) ** m
    else:
        L = 1.0
    return L * S


def c_factor(ndvi, alpha=2.0, beta=1.0):
    """
    Factor de cobertura C (van der Knijff) a partir del NDVI:
        C = exp(−α·NDVI/(β−NDVI))     acotado a [0, 1].
    NDVI se recorta a < β para evitar la singularidad.
    """
    nd = min(float(ndvi), beta - 1e-6)
    C = math.exp(-alpha * nd / (beta - nd))
    return min(max(C, 0.0), 1.0)


def p_factor(slope_pct, es_cultivo=False):
    """
    Factor de prácticas de soporte P.
    Sin intervención (bosque, agua, urbano) → P = 1.0. En cultivo, P crece con
    la pendiente (siembras en contorno menos efectivas en laderas empinadas).
    """
    if not es_cultivo:
        return 1.0
    s = float(slope_pct)
    if s < 2:
        return 0.6
    if s < 5:
        return 0.5
    if s < 8:
        return 0.5
    if s < 12:
        return 0.6
    if s < 16:
        return 0.7
    if s < 20:
        return 0.8
    return 0.9


def soil_loss(R, K, LS, C, P):
    """Pérdida de suelo A = R·K·LS·C·P [t·ha⁻¹·año⁻¹]."""
    return max(R * K * LS * C * P, 0.0)


def clase_severidad(A):
    """Clasifica la tasa de pérdida de suelo (t·ha⁻¹·año⁻¹)."""
    if A < 5:
        return "Ligera", "#1a9850"
    if A < 10:
        return "Moderada", "#a6d96a"
    if A < 20:
        return "Alta", "#fee08b"
    if A < 40:
        return "Muy alta", "#fdae61"
    if A < 80:
        return "Severa", "#f46d43"
    return "Crítica", "#d73027"


def compute_rusle(precip_anual_mm, clase_textura, slope_pct, ndvi,
                  es_cultivo=False, r_formula="hurni", slope_length_m=None):
    """
    Combina los cinco factores y devuelve el desglose + la pérdida de suelo.
    Todas las entradas son valores puntuales (ya extraídos de GEE).
    """
    R = r_factor(precip_anual_mm, r_formula)
    K = k_factor_textura(clase_textura)
    LS = ls_factor(slope_pct, slope_length_m)
    C = c_factor(ndvi)
    P = p_factor(slope_pct, es_cultivo)
    A = soil_loss(R, K, LS, C, P)
    sev, color = clase_severidad(A)
    return {
        "R": round(R, 2), "K": round(K, 4), "LS": round(LS, 3),
        "C": round(C, 4), "P": round(P, 3),
        "A": round(A, 2), "severidad": sev, "color": color,
        "r_formula": r_formula,
    }
