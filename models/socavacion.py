"""
Módulo de socavación fluvial.

Fórmulas y coeficientes verificados contra las FUENTES PRIMARIAS del cuaderno
HIDRÁULICA FLUVIAL (NotebookLM), aislando cada PDF para evitar recuperación
cruzada:
    - Maza Álvarez, J.A. (1968). Socavación en Cauces Naturales. SOP, México.
      → Lischtvan-Lebediev completo, tablas β (p.9), x (p.10) y μ (p.13), ψ.
    - Manual de Hidrología y Drenaje, Administradora Boliviana de Carreteras
      (ABC) → §6.5: socavación en pilas (Tablas 6.5-1 a 6.5-4) y estribos
      (§6.5.3, Froehlich/HEC-18, tablas Ks y Kθ). REFERENCIA NORMATIVA BOLIVIA.
    - Hydraulic Design of Safe Bridges (FHWA, HEC-18, 2012).

Cuatro tipos:
    1. GENERAL       (descenso generalizado del fondo)
    2. CONTRACCIÓN   (estrechamiento de la sección)
    3. LOCAL EN PILAS
    4. LOCAL EN ESTRIBOS

OJO con las unidades del diámetro: cambian de método a método
(mm en Lacey, Blench y HEC-18 forma "1.48"; m en Lischtvan-Lebediev, Kellerhals
y V_c).

Todas las funciones son puras (números → números, SI). La función
`compute_socavacion(params)` arma el resultado estructurado para la app.
"""

import numpy as np

G = 9.807


# ---------------------------------------------------------------------------
# TABLAS DE COEFICIENTES  (verificadas contra Maza 1968 y Manual ABC)
# ---------------------------------------------------------------------------

# β — coeficiente de paso según periodo de retorno (Maza 1968, p.9).
# Prob. excedencia 100%→Tr1 ... 0.1%→Tr1000.
_BETA_TR = {
    1: 0.77, 2: 0.82, 5: 0.86, 10: 0.90, 20: 0.94,
    50: 0.97, 100: 1.00, 300: 1.03, 500: 1.05, 1000: 1.07,
}

# x — exponente, suelos NO cohesivos vs diámetro medio d_m (mm). Maza p.10.
_X_NO_COHESIVO = {
    0.05: 0.43, 0.15: 0.42, 0.50: 0.41, 1.0: 0.40, 1.5: 0.39,
    2.5: 0.38, 4.0: 0.37, 6.0: 0.36, 8.0: 0.35, 10: 0.34,
    15: 0.33, 20: 0.32, 25: 0.31, 40: 0.30, 60: 0.29,
    90: 0.28, 140: 0.27, 190: 0.26, 250: 0.25, 310: 0.24,
    370: 0.23, 450: 0.22, 570: 0.21, 750: 0.20, 1000: 0.19,
}

# x — exponente, suelos COHESIVOS vs peso específico seco γs (t/m³). Maza p.10.
_X_COHESIVO = {
    0.80: 0.52, 0.83: 0.51, 0.86: 0.50, 0.88: 0.49, 0.90: 0.48,
    0.93: 0.47, 0.96: 0.46, 0.98: 0.45, 1.00: 0.44, 1.04: 0.43,
    1.08: 0.42, 1.12: 0.41, 1.16: 0.40, 1.20: 0.39, 1.24: 0.38,
    1.28: 0.37, 1.34: 0.36, 1.40: 0.35, 1.46: 0.34, 1.52: 0.33,
    1.58: 0.32, 1.64: 0.31, 1.71: 0.30, 1.80: 0.29, 1.89: 0.28,
    2.00: 0.27,
}

# ψ — corrección por sedimento en suspensión (Maza p.11, Tabla 6.5-14 ABC).
# Entrada: peso específico de la mezcla agua-sedimento γ_mezcla (t/m³).
_PSI_SUSPENSION = {
    1.00: 1.00, 1.05: 1.08, 1.10: 1.13, 1.15: 1.20, 1.20: 1.27,
    1.25: 1.34, 1.30: 1.42, 1.35: 1.50, 1.40: 1.60,
}

# μ — factor de contracción por pilas (Maza 1968, p.13).
# Doble entrada: velocidad media V (filas) × claro libre entre pilas (columnas).
# La celda (2.0 m/s, 16 m) del PDF sale 0.97 por error de OCR; se corrige a 0.95
# para respetar la monotonía de la tabla.
_MU_CLAROS = [10, 13, 16, 18, 21, 25, 30, 42, 52, 63, 106, 124, 200]
_MU_VELOCIDADES = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
_MU_MATRIZ = [
    [0.96, 0.97, 0.98, 0.98, 0.99, 0.99, 0.99, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00],
    [0.94, 0.96, 0.97, 0.97, 0.97, 0.98, 0.99, 0.99, 0.99, 0.99, 1.00, 1.00, 1.00],
    [0.93, 0.94, 0.95, 0.96, 0.97, 0.97, 0.98, 0.98, 0.99, 0.99, 0.99, 0.99, 1.00],
    [0.90, 0.93, 0.94, 0.95, 0.96, 0.96, 0.97, 0.98, 0.98, 0.99, 0.99, 0.99, 1.00],
    [0.89, 0.91, 0.93, 0.94, 0.95, 0.96, 0.96, 0.97, 0.98, 0.98, 0.99, 0.99, 0.99],
    [0.87, 0.90, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.98, 0.99, 0.99, 0.99],
    [0.85, 0.89, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.99, 0.99],
]

# k1 — exponente de Laursen para contracción, según u*/w (HEC-18).
# (se resuelve en k1_contraccion)

# Ks — forma de la pila (Manual ABC, Tabla 6.5-1). Valor representativo.
PIER_SHAPE_K1 = {
    "circular": 1.00,
    "lenticular": 0.75,
    "eliptica": 0.70,
    "rectangular": 1.10,
    "rect_semicircular": 0.90,       # rectangular con extremo semicircular
    "rect_redondeado": 1.01,
    "semicircular_triangular": 0.86,  # nariz semicircular y cola triangular
    "rect_triangular": 0.70,          # rectangular con nariz triangular
    # compatibilidad con nombres previos del formulario:
    "cuadrada": 1.10,
    "triangular": 0.70,
    "grupo_cilindros": 1.00,
}

# K3 — condición del lecho (CSU / HEC-18).
PIER_BED_K3 = {
    "plano": 1.1, "dunas_medianas": 1.1, "dunas_grandes": 1.3, "antidunas": 1.1,
}

# Kd (Ettema 1980) — tamaño relativo del sedimento en pilas (ABC Tabla 6.5-4).
_KD_ETTEMA = {8: 0.72, 10: 0.77, 15: 0.87, 20: 0.94, 25: 1.00, 30: 1.00, 50: 1.00}

# Ks — forma del estribo (Manual ABC, §6.5.3.2).
ABUT_SHAPE_KS = {
    "vertical": 1.00,
    "vertical_semicircular": 0.75,
    "vertical_aletas": 0.75,
    "talud_1_2": 0.60,   # pared inclinada H:V = 1:2
    "talud_1_1": 0.50,
    "talud_15_1": 0.45,
    "talud_2_1": 0.30,
    # compatibilidad: "spill-through" con talud ≈ 0.55, vertical = 1.00
    "spill": 0.55,
}

# Kθ — ángulo de esviaje del estribo (Manual ABC). φ medido entre el eje del
# estribo y la ribera aguas arriba.
_KTHETA_ESTRIBO = {30: 1.10, 60: 1.05, 90: 1.00, 120: 0.98, 150: 0.90}

# ── ARTAMONOV — estribos (Maza 1968, ec. 36, p.92) ──────────────────────────
# S_T = Pα · Pq · Pk · Ho  (medido desde la superficie libre del agua).
# Pα — ángulo de esviaje del estribo. (El PDF empieza en 30°; la tabla clásica
# de Artamonov arranca en 20° para 0.84 — posible error de lectura, ver aviso.)
_ARTAMONOV_PALPHA = {30: 0.84, 60: 0.94, 90: 1.00, 120: 1.07, 150: 1.18}
# Pq — relación de gastos Q1/Q (Q1 = gasto teórico por la zona del estribo).
_ARTAMONOV_PQ = {0.10: 2.00, 0.20: 2.65, 0.30: 3.22, 0.40: 3.45,
                 0.50: 3.67, 0.60: 3.87, 0.70: 4.06, 0.80: 4.20}
# Pk — talud del estribo k (horizontal:vertical).
_ARTAMONOV_PK = {0: 1.00, 0.5: 0.91, 1.0: 0.85, 1.5: 0.83, 2.0: 0.61, 3.0: 0.50}

# ── LAURSEN-TOCH — pilas (Maza 1968, ec. 33-34, pp.59-64) ───────────────────
# K2 — forma de la nariz (Tabla XII, p.59). TABLA REAL citable.
# K1 (Fig.21) y K3 (Fig.22) son GRÁFICAS: el usuario debe ingresarlas leídas
# del PDF (el cuaderno advirtió que los valores leídos por IA no cuadran).
LAURSEN_TOCH_K2 = {
    "rectangular": 1.00,          # a/b = 4
    "semicircular": 0.90,
    "circular": 0.90,             # nariz circular ≈ semicircular
    "eliptica_2_1": 0.81,         # P/r = 2/1
    "eliptica_3_1": 0.75,         # P/r = 3/1
    "lenticular_2_1": 0.81,       # P/r = 2/1
    "lenticular_3_1": 0.69,       # P/r = 3/1
    "biselada": 0.78,             # a/b = 4 (Tison)
    "hidrodinamico": 0.75,        # perfil hidrodinámico a/b = 4 (Tison)
}
# Mapeo desde las formas del formulario a la K2 de Laursen-Toch.
_LT_K2_FROM_FORMA = {
    "circular": 0.90, "rectangular": 1.00, "cuadrada": 1.00,
    "rect_semicircular": 0.90, "rect_redondeado": 0.90,
    "semicircular_triangular": 0.90, "rect_triangular": 0.78,
    "lenticular": 0.81, "eliptica": 0.81, "grupo_cilindros": 0.90,
}

# ── YAROSLAVTZIEV — pilas (Maza 1968, ec. 35, pp.65-74) ─────────────────────
# Kf — factor de forma de pila (valores tabulados del libro).
_YAROS_KF = {
    "rectangular": 12.4,       # Tipo I
    "cuadrada": 12.4,
    "circular": 10.0,          # Tipo II
    "semicircular_triangular": 8.7,   # Tipo III frentes semicirc. (φ≈10°)
    "rect_semicircular": 8.5,
    "eliptica": 8.5,
    "lenticular": 8.5,
    "rect_triangular": 10.0,   # nariz apuntada, tajamar 90°
}

# Coeficientes de Yarnell — obstrucción por pilas (K y C_D).
YARNELL_K = {
    "semicircular": {"K": 0.90, "Cd": 0.90},
    "cilindros_diafragma": {"K": 0.95, "Cd": 0.95},
    "cilindros_sin_diafragma": {"K": 1.05, "Cd": 1.05},
    "triangular_90": {"K": 1.05, "Cd": 1.05},
    "rectangular_plana": {"K": 1.25, "Cd": 1.25},
    "caballete_10_pilotes": {"K": 2.50, "Cd": 2.50},
}


def _interp_table(table, key):
    """Interpolación lineal sobre una tabla {x: y} (extremos planos)."""
    xs = sorted(table.keys())
    if key <= xs[0]:
        return table[xs[0]]
    if key >= xs[-1]:
        return table[xs[-1]]
    return float(np.interp(key, xs, [table[x] for x in xs]))


def beta_periodo_retorno(tr_anios):
    """Coeficiente β de Lischtvan-Lebediev según periodo de retorno (años)."""
    return _interp_table(_BETA_TR, float(tr_anios))


def exponente_x(dm_mm=None, gamma_s=None, cohesivo=False):
    """Exponente x. No cohesivo: f(d_m mm). Cohesivo: f(γs t/m³)."""
    if cohesivo:
        if gamma_s is None:
            raise ValueError("Suelo cohesivo requiere gamma_s (t/m³)")
        return _interp_table(_X_COHESIVO, float(gamma_s))
    if dm_mm is None:
        raise ValueError("Suelo no cohesivo requiere dm (mm)")
    return _interp_table(_X_NO_COHESIVO, float(dm_mm))


def psi_suspension(gamma_mezcla):
    """Corrección ψ por sedimento en suspensión, según γ de la mezcla (t/m³)."""
    return _interp_table(_PSI_SUSPENSION, float(gamma_mezcla))


def mu_contraccion(velocidad, claro_libre):
    """
    Factor de contracción μ por pilas (Maza 1968, p.13), interpolación bilineal.
        velocidad    : velocidad media de la corriente (m/s)
        claro_libre  : claro libre entre pilas / luz del vano (m)
    Para V < 1 m/s, μ = 1.00.
    """
    v = float(velocidad)
    c = float(claro_libre)
    if v < 1.0:
        return 1.00
    # columna interpolada para cada fila de velocidad, luego interpolar en V
    col_vals = [float(np.interp(c, _MU_CLAROS, fila)) for fila in _MU_MATRIZ]
    return float(np.interp(v, _MU_VELOCIDADES, col_vals))


def velocidad_critica(y1, d50_mm):
    """
    Velocidad crítica de arrastre (Laursen, S.I.):
        V_c = 6.2 · y1^(1/6) · d50^(1/3)   [d50 en METROS]
    Sirve para decidir el régimen de contracción:
        V1 < V_c -> agua clara ; V1 >= V_c -> lecho vivo.
    """
    d50_m = d50_mm / 1000.0
    return 6.2 * y1 ** (1.0 / 6.0) * d50_m ** (1.0 / 3.0)


# ===========================================================================
# 1. SOCAVACIÓN GENERAL
# ===========================================================================

def lacey(Q, dm_mm):
    """Lacey (1930). Q (m³/s), d_m (mm). h_ms = 0.389·Q^(1/3)/d_m^(1/6)."""
    return 0.389 * (Q ** (1.0 / 3.0) / dm_mm ** (1.0 / 6.0))


def blench(q, d50_mm):
    """
    Blench (1969). q caudal unitario (m³/s·m), d50 (mm).
        arenas (d50<2mm): 1.200·q^(2/3)/d50^(1/6)
        grava  (d50>2mm): 1.230·q^(2/3)/d50^(1/12)
    """
    if d50_mm > 2.0:
        return 1.230 * q ** (2.0 / 3.0) / d50_mm ** (1.0 / 12.0)
    return 1.200 * q ** (2.0 / 3.0) / d50_mm ** (1.0 / 6.0)


def maza_echavarria(Q, B, d50_m):
    """Maza & Echavarría (1973). Q (m³/s), B (m), d50 en METROS."""
    return 0.365 * (Q ** 0.784) / (B ** 0.784 * d50_m ** 0.157)


def lischtvan_lebediev_reducida(q, d50_m):
    """L-L forma reducida (arenas). q (m³/s·m), d50 en METROS."""
    return 0.333 * q ** 0.710 / d50_m ** 0.199


def kellerhals(q, d90_m):
    """Kellerhals (grava gruesa). q (m³/s·m), d90 en METROS."""
    return 0.470 * q ** 0.800 / d90_m ** 0.120


def lischtvan_lebediev_maza(Ho, Hm, Q, We, mu, beta, cohesivo=False,
                            dm_mm=None, gamma_s=None, psi=1.0):
    """
    Lischtvan-Lebediev / Maza (1968) — método completo por franjas.

    α = Q / (Hm^(5/3) · We · μ)
        No cohesivo: Hs = [ α·Ho^(5/3) / (0.68·β·ψ·d_m^0.28) ]^(1/(1+x))
        Cohesivo:    Hs = [ α·Ho^(5/3) / (0.60·β·ψ·γs^1.18) ]^(1/(1+x))

    ψ (sedimento en suspensión) va en el denominador junto a β (Maza).
    Devuelve (Hs desde la superficie, socavación = Hs - Ho).
    """
    alpha = Q / (Hm ** (5.0 / 3.0) * We * mu)
    x = exponente_x(dm_mm=dm_mm, gamma_s=gamma_s, cohesivo=cohesivo)
    exp = 1.0 / (1.0 + x)
    if cohesivo:
        denom = 0.60 * beta * psi * gamma_s ** 1.18
    else:
        denom = 0.68 * beta * psi * dm_mm ** 0.28
    Hs = (alpha * Ho ** (5.0 / 3.0) / denom) ** exp
    return Hs, max(Hs - Ho, 0.0)


# ===========================================================================
# 2. SOCAVACIÓN POR CONTRACCIÓN
# ===========================================================================

def k1_contraccion(u_star, w):
    """Exponente k1 de Laursen según u*/w (modo de transporte)."""
    if w <= 0:
        return 0.64
    r = u_star / w
    if r < 0.5:
        return 0.59
    if r <= 2.0:
        return 0.64
    return 0.69


def laursen_contraccion_lecho_vivo(y1, W1, W2, k1):
    """Laursen — lecho vivo. y2 = y1·(W1/W2)^k1. Devuelve (y2, y2-y1)."""
    y2 = y1 * (W1 / W2) ** k1
    return y2, max(y2 - y1, 0.0)


def laursen_contraccion_agua_clara(y1, W1, W2, V1, d50_mm):
    """
    Laursen — agua clara (forma completa, S.I.):
        y2/y1 = (W1/W2)^(6/7) · [ V1² / (36·y1^(1/3)·d50^(2/3)) ]^(3/7)
    d50 en METROS. Devuelve (y2, y2-y1).
    """
    d50_m = d50_mm / 1000.0
    ratio = (W1 / W2) ** (6.0 / 7.0) * (
        V1 ** 2 / (36.0 * y1 ** (1.0 / 3.0) * d50_m ** (2.0 / 3.0))) ** (3.0 / 7.0)
    y2 = y1 * ratio
    return y2, max(y2 - y1, 0.0)


def parker_contraccion(y1, W1, W2):
    """Parker (1985) — grava sin acorazamiento. y2 = y1·(W1/W2)^0.825."""
    y2 = y1 * (W1 / W2) ** 0.825
    return y2, max(y2 - y1, 0.0)


def hec18_contraccion_agua_clara(Q2, dm_m, W2, Ku=0.025):
    """
    Richardson & Davis / HEC-18 — agua clara (forma canónica, S.I.):
        y2 = [ Ku · Q2² / (D_m^(2/3) · W2²) ]^(3/7),  Ku = 0.025, D_m en METROS.
    Equivale a y2 = 1.48·[Q2/(d_m^(1/3)·W2)]^(6/7) con d_m en MILÍMETROS.
    """
    return (Ku * Q2 ** 2 / (dm_m ** (2.0 / 3.0) * W2 ** 2)) ** (3.0 / 7.0)


# ===========================================================================
# 3. SOCAVACIÓN LOCAL EN PILAS
# ===========================================================================

def kd_pila(b, d50_mm):
    """Kd (Ettema 1980) — tamaño relativo del sedimento, según b/d50.
    b en metros, d50 en mm → relación adimensional b/d50 = (b·1000)/d50_mm."""
    if d50_mm <= 0:
        return 1.0
    return _interp_table(_KD_ETTEMA, (b * 1000.0) / d50_mm)


def k_esviaje_pila(theta_deg, largo_ancho):
    """
    Factor de esviaje Kω (Froehlich, Manual ABC):
        Kω = (cos θ + (L/b)·sen θ)^0.62
    theta en grados ; largo_ancho = L/b.
    """
    th = np.radians(theta_deg)
    return (np.cos(th) + largo_ancho * np.sin(th)) ** 0.62


def breusers_nicollet_shen(b, y, Ks=1.0, Ktheta=1.0):
    """Breusers-Nicollet-Shen (1977): d_s = b·2.00·tanh(y/b)·Ks·Kθ."""
    return b * 2.00 * np.tanh(y / b) * Ks * Ktheta


def csu_hec18_pila(a, y1, Fr1, K1=1.0, K2=1.0, K3=1.1, K4=1.0, aplicar_tope=True):
    """
    CSU / HEC-18 (Richardson & Davis, 1995):
        y_s = 2.0·K1·K2·K3·K4·a^0.65·y1^0.35·Fr1^0.43
    Topes (Manual ABC): d_s,máx = 2.4a si Fr≤0.8 ; 3.0a si Fr>0.8.
    """
    ys = 2.0 * K1 * K2 * K3 * K4 * a ** 0.65 * y1 ** 0.35 * Fr1 ** 0.43
    if aplicar_tope:
        tope = 2.4 * a if Fr1 <= 0.8 else 3.0 * a
        ys = min(ys, tope)
    return ys


def laursen_toch(b, K1, K2=1.0, K3=None, esviaje=False):
    """
    Laursen-Toch (Maza 1968, ec.33-34) — socavación local en pila (medida
    desde el fondo ya erosionado por socavación general):
        pila alineada:  S_o = K1·K2·b
        pila esviajada: S   = K1·K3·b
    b = ancho de pila (m). K2 = forma de la nariz (tabla real). K1 (Fig.21) y
    K3 (Fig.22) son GRÁFICAS → deben ingresarse leídas del PDF de Maza.
    """
    if esviaje and K3 is not None:
        return K1 * K3 * b
    return K1 * K2 * b


def yaroslavtziev(V, b1, H, Kf, D85_m=0.0, C=0.6):
    """
    Yaroslavtziev (Maza 1968, ec.35) — socavación local en pila (desde el
    fondo ya erosionado por socavación general):
        S_o = Kf·Kv·(C+KH)·V²/g − 30·D85
    con las funciones ANALÍTICAS (el libro advierte usar la ecuación, no la
    figura):
        log Kv = −0.28·√(V²/(g·b1))
        log KH = 0.17 − 0.35·(H/b1)

    Unidades: V (m/s), b1 (m) proyección de la pila ⊥ al flujo, H (m),
    D85 en METROS. C = 0.6 (cauce principal) ó 1.0 (cauce de inundación).
    Si D85 < 0.005 m (arenas < 0.5 cm) el término −30·D85 se omite (sin
    acorazamiento).

    NOTA de unidades: el cuaderno transcribió D85 "en cm" con el término −30·D85,
    lo que es dimensionalmente inconsistente (restaría decenas de metros). Se usa
    la forma publicada estándar con D85 en METROS; conviene verificar contra la
    figura original del PDF de Maza.
    """
    Kv = 10.0 ** (-0.28 * np.sqrt(V ** 2 / (G * b1)))
    KH = 10.0 ** (0.17 - 0.35 * (H / b1))
    termino_d = 30.0 * D85_m if D85_m >= 0.005 else 0.0
    So = Kf * Kv * (C + KH) * V ** 2 / G - termino_d
    return max(So, 0.0)


# ===========================================================================
# 4. SOCAVACIÓN LOCAL EN ESTRIBOS
# ===========================================================================

def froehlich_hec18_estribo(L, ya, Fr, Ks=1.0, Ktheta=1.0):
    """
    Froehlich (1989a) / HEC-18 — el método recomendado (lecho móvil y agua
    clara, no sobreestima):
        d_s = L · 2.27·Ks·Kθ·(ya/L)^0.57·Fr^0.61
    L longitud del estribo (m) ; ya tirante en la llanura (m) ; Fr = V/√(g·ya).
    """
    return L * 2.27 * Ks * Ktheta * (ya / L) ** 0.57 * Fr ** 0.61


def liu_chang_skinner(L, y, Fr, spill_through=True):
    """
    Liu, Chang & Skinner (1961) — lecho móvil.
        spill-through: d_s = L·1.10·(y/L)^0.6·Fr^0.33
        vertical:      d_s = L·2.15·(y/L)^0.6·Fr^0.33
    """
    c = 1.10 if spill_through else 2.15
    return L * c * (y / L) ** 0.6 * Fr ** 0.33


def laursen_estribo(L, y, agua_clara=False):
    """Laursen (1962/63): d_s = L·(1.57 lecho móvil | 1.89 agua clara)·(y/L)^0.5."""
    c = 1.89 if agua_clara else 1.57
    return L * c * (y / L) ** 0.5


def ktheta_estribo(theta_deg):
    """Factor Kθ por esviaje del estribo (Manual ABC, interpolado)."""
    return _interp_table(_KTHETA_ESTRIBO, float(theta_deg))


def artamonov(Ho, Q1_Q, alpha=90, k=0.0, espigones=False):
    """
    Artamonov (Maza 1968, ec.36) — socavación al pie del estribo:
        S_T = Pα · Pq · Pk · Ho     (medido desde la SUPERFICIE del agua)
    Pα = ángulo de esviaje α ; Pq = relación de gastos Q1/Q ; Pk = talud k.
    `espigones=True` (espigones enfrentados en ambas orillas) aplica −25 %.
    Devuelve (S_T desde superficie, socavación local = S_T − Ho).
    """
    Pa = _interp_table(_ARTAMONOV_PALPHA, float(alpha))
    Pq = _interp_table(_ARTAMONOV_PQ, float(Q1_Q))
    Pk = _interp_table(_ARTAMONOV_PK, float(k))
    St = Pa * Pq * Pk * Ho
    if espigones:
        St *= 0.75
    return St, max(St - Ho, 0.0)


# ===========================================================================
# ORQUESTADOR
# ===========================================================================

def _row(metodo, valor, unidad, nota="", referencia=""):
    return {
        "metodo": metodo,
        "valor": round(float(valor), 4) if valor is not None else None,
        "unidad": unidad,
        "nota": nota,
        "referencia": referencia,
    }


def compute_socavacion(params):
    """
    Calcula todos los métodos de socavación cuyos datos estén disponibles.
    Ver README del módulo para la lista de claves de `params`.
    """
    d50 = params.get("d50_mm")
    d90 = params.get("d90_mm")
    dm = params.get("dm_mm", d50)
    y = params.get("depth")
    v = params.get("velocity")
    Q = params.get("Q")
    B = params.get("B") or params.get("W1")
    W2 = params.get("W2")
    tr = params.get("tr", 100)
    cohesivo = bool(params.get("cohesivo", False))
    gamma_s = params.get("gamma_s")
    u_star = params.get("u_star")
    w = params.get("w")
    claro = params.get("claro_pilas")
    gamma_mezcla = params.get("gamma_mezcla")

    out = {"general": [], "contraccion": [], "pila": [], "estribo": [],
           "yarnell": YARNELL_K, "avisos": []}

    # μ: de la tabla de Maza si hay V y claro entre pilas; si no, el dado o 1.0
    if v and claro:
        mu = mu_contraccion(v, claro)
    else:
        mu = params.get("mu", 1.0) or 1.0
    # ψ: de la tabla si se da el peso de la mezcla; si no, 1.0
    psi = psi_suspension(gamma_mezcla) if gamma_mezcla else 1.0

    # Caudal unitario q (m³/s·m)
    q = None
    if Q and B and B > 0:
        q = Q / B
    elif v and y:
        q = v * y

    # --- 1. GENERAL ---
    if Q and dm:
        out["general"].append(_row(
            "Lacey (1930)", lacey(Q, dm), "m",
            "Prof. media desde superficie · d_m en mm", "Lacey 1930"))
    if q and d50:
        out["general"].append(_row(
            "Blench (1969)", blench(q, d50), "m",
            "Sin acorazamiento · d50 en mm", "Blench 1969"))
    if Q and B and d50:
        out["general"].append(_row(
            "Maza & Echavarría (1973)", maza_echavarria(Q, B, d50 / 1000.0),
            "m", "Calibrada en ríos sudamericanos", "Maza-Echavarría 1973"))
    if q and d50:
        out["general"].append(_row(
            "Lischtvan-Lebediev (reducida)",
            lischtvan_lebediev_reducida(q, d50 / 1000.0), "m",
            "Rango de arenas", "Lischtvan-Lebediev"))
    if q and d90:
        out["general"].append(_row(
            "Kellerhals", kellerhals(q, d90 / 1000.0), "m",
            "Grava gruesa · d90", "Kellerhals"))

    # Lischtvan-Lebediev / Maza completo (por franjas)
    if y and Q and B:
        try:
            beta = beta_periodo_retorno(tr)
            Hs, sc = lischtvan_lebediev_maza(
                Ho=y, Hm=y, Q=Q, We=B, mu=mu, beta=beta,
                cohesivo=cohesivo, dm_mm=dm, gamma_s=gamma_s, psi=psi)
            tipo = "cohesivo" if cohesivo else "no cohesivo"
            extra = f" · μ={mu:.2f}" + (f" · ψ={psi:.2f}" if psi != 1.0 else "")
            out["general"].append(_row(
                "Lischtvan-Lebediev / Maza (completo)", sc, "m",
                f"Hs={Hs:.2f} m desde superficie · suelo {tipo} · β={beta:.2f}"
                f" (Tr={tr} a){extra}", "Maza 1968"))
        except Exception as e:
            out["avisos"].append(f"L-L/Maza no calculable: {e}")

    # --- 2. CONTRACCIÓN ---
    if y and B and W2 and W2 > 0:
        # Régimen: comparar V1 con V_c
        regimen = None
        if v and d50:
            vc = velocidad_critica(y, d50)
            regimen = "lecho vivo" if v >= vc else "agua clara"
            out["contraccion"].append(_row(
                "Régimen (V_c Laursen)", vc, "m/s",
                f"V1={v:.2f} m/s vs V_c={vc:.2f} → {regimen}", "HEC-18"))
        # Lecho vivo (Laursen)
        k1 = k1_contraccion(u_star or 0.0, w or 1.0)
        y2, sc = laursen_contraccion_lecho_vivo(y, B, W2, k1)
        out["contraccion"].append(_row(
            "Laursen (lecho vivo)", sc, "m",
            f"y2={y2:.2f} m · k1={k1}"
            + (" ✓ aplica" if regimen == "lecho vivo" else ""), "Laursen"))
        # Agua clara (Laursen, forma completa)
        if v and d50:
            y2c, scc = laursen_contraccion_agua_clara(y, B, W2, v, d50)
            out["contraccion"].append(_row(
                "Laursen (agua clara)", scc, "m",
                f"y2={y2c:.2f} m"
                + (" ✓ aplica" if regimen == "agua clara" else ""), "Laursen"))
        out["contraccion"].append(_row(
            "Parker (1985)", parker_contraccion(y, B, W2)[1], "m",
            "Grava sin acorazar", "Parker 1985"))
        if Q and dm:
            y2h = hec18_contraccion_agua_clara(Q, dm / 1000.0, W2)
            out["contraccion"].append(_row(
                "HEC-18 (agua clara)", max(y2h - y, 0.0), "m",
                f"y2={y2h:.2f} m · Ku=0.025", "Richardson & Davis / HEC-18"))

    # --- 3. PILAS ---
    b = params.get("pila_ancho")
    if b and y:
        forma = params.get("pila_forma", "circular")
        Ks = PIER_SHAPE_K1.get(forma, 1.0)
        theta = params.get("pila_theta", 0.0)
        largo = params.get("pila_largo")
        Kw = 1.0
        if theta and largo and b > 0:
            Kw = k_esviaje_pila(theta, largo / b)
        K3 = PIER_BED_K3.get(params.get("pila_lecho", "plano"), 1.1)
        out["pila"].append(_row(
            "Breusers-Nicollet-Shen (1977)",
            breusers_nicollet_shen(b, y, Ks=Ks, Ktheta=Kw), "m",
            f"Ks={Ks} · Kω={Kw:.2f} · forma {forma}", "Breusers et al. 1977"))
        if params.get("froude"):
            fr = params["froude"]
            ys = csu_hec18_pila(b, y, fr, K1=Ks, K2=Kw, K3=K3)
            tope = 2.4 * b if fr <= 0.8 else 3.0 * b
            capped = " (limitado al tope)" if ys >= tope - 1e-9 else ""
            out["pila"].append(_row(
                "CSU / HEC-18 (normativo)", ys, "m",
                f"Ks={Ks} Kω={Kw:.2f} K3={K3} · Fr={fr} · tope {tope:.2f} m{capped}",
                "Richardson & Davis 1995 / HEC-18"))
        # Yaroslavtziev (Maza) — requiere velocidad
        if v:
            Kf = _YAROS_KF.get(forma, 10.0)
            # b1 = proyección ⊥ al flujo; con largo y ángulo, para rectangular
            b1 = b
            if theta and largo:
                b1 = largo * np.sin(np.radians(theta)) + b * np.cos(np.radians(theta))
            # D85 ≈ d90 (proxy), de mm a m
            d85_m = (d90 / 1000.0) if d90 else 0.0
            C = 1.0 if params.get("pila_inundacion") else 0.6
            sy = yaroslavtziev(v, b1, y, Kf, D85_m=d85_m, C=C)
            zona = "cauce de inundación" if C == 1.0 else "cauce principal"
            out["pila"].append(_row(
                "Yaroslavtziev (Maza)", sy, "m",
                f"Kf={Kf} · b1={b1:.2f} m · C={C} ({zona})"
                + (" · D85≈d90" if d85_m >= 0.005 else " · arena (sin −30·D85)"),
                "Yaroslavtziev / Maza 1968"))
        # Laursen-Toch (Maza) — requiere K1 leído de la Fig.21 (gráfica)
        k1_lt = params.get("pila_k1")
        if k1_lt:
            K2 = _LT_K2_FROM_FORMA.get(forma, 1.0)
            k3_lt = params.get("pila_k3")
            esv = bool(k3_lt) and bool(theta)
            slt = laursen_toch(b, k1_lt, K2=K2, K3=k3_lt, esviaje=esv)
            modo = "esviajada (K1·K3·b)" if esv else "alineada (K1·K2·b)"
            out["pila"].append(_row(
                "Laursen-Toch (Maza)", slt, "m",
                f"K1={k1_lt} (Fig.21) · K2={K2} · forma {forma} · {modo}",
                "Laursen-Toch / Maza 1968"))

    # --- 4. ESTRIBOS ---
    L = params.get("estribo_L")
    if L and y and L > 0:
        spill = bool(params.get("estribo_spill", True))
        agua_clara = bool(params.get("estribo_agua_clara", False))
        forma_e = params.get("estribo_forma", "spill" if spill else "vertical")
        Kse = ABUT_SHAPE_KS.get(forma_e, 0.55 if spill else 1.0)
        theta_e = params.get("estribo_theta", 90)
        Kte = ktheta_estribo(theta_e)
        # Froehlich / HEC-18 (recomendado)
        if params.get("froude"):
            fr = params["froude"]
            out["estribo"].append(_row(
                "Froehlich / HEC-18 (recomendado)",
                froehlich_hec18_estribo(L, y, fr, Ks=Kse, Ktheta=Kte), "m",
                f"Ks={Kse} · Kθ={Kte:.2f} · forma {forma_e}", "Froehlich 1989 / HEC-18"))
            out["estribo"].append(_row(
                "Liu-Chang-Skinner (1961)",
                liu_chang_skinner(L, y, fr, spill_through=spill), "m",
                f"{'spill-through' if spill else 'vertical'} · lecho móvil",
                "Liu, Chang & Skinner 1961"))
        out["estribo"].append(_row(
            "Laursen (estribo)", laursen_estribo(L, y, agua_clara=agua_clara),
            "m", f"{'agua clara' if agua_clara else 'lecho móvil'}", "Laursen 1962/63"))
        # Artamonov (Maza) — requiere la relación de gastos Q1/Q
        q1q = params.get("estribo_q1q")
        if q1q:
            # talud k desde la forma del estribo (H:V)
            k_map = {"talud_1_2": 0.5, "talud_1_1": 1.0, "talud_15_1": 1.5,
                     "talud_2_1": 2.0, "spill": 1.5}
            k = k_map.get(forma_e, 0.0)
            alpha_a = params.get("estribo_theta", 90)
            espig = bool(params.get("estribo_espigones", False))
            St, sc_a = artamonov(y, q1q, alpha=alpha_a, k=k, espigones=espig)
            out["estribo"].append(_row(
                "Artamonov (Maza)", sc_a, "m",
                f"S_T={St:.2f} m desde superficie · Q1/Q={q1q} · k={k} · α={alpha_a}°"
                + (" · espigones −25%" if espig else ""), "Artamonov / Maza 1968"))

    return out
