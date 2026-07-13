"""
Módulo de socavación fluvial.

Fórmulas tomadas del cuaderno HIDRAULICA FLUVIAL (NotebookLM), con las fuentes:
    - AnalisisdelaSocavacionenCaucesNaturales.pdf
    - Socavacion cauces naturales_Maza.pdf
    - Hidrología e Hidráulica Aplicada a puentes.pdf
    - Hydraulic Design of Safe Bridges (FHWA / HEC-18)
    - manual_de_hidrologia_y_drenaje.pdf

Se distinguen cuatro tipos:
    1. Socavación GENERAL      (descenso generalizado del fondo)
    2. Socavación por CONTRACCIÓN (estrechamiento de la sección)
    3. Socavación LOCAL EN PILAS  (pie de pila de puente)
    4. Socavación LOCAL EN ESTRIBOS

OJO con las unidades del diámetro: cambian de método a método
(mm en Lacey y Blench; m en Lischtvan-Lebediev y Kellerhals).

Todas las funciones son puras: reciben números y devuelven números (SI).
La función orquestadora `compute_socavacion(params)` arma el resultado
estructurado para la app, calculando solo los métodos cuyos datos existen.
"""

import numpy as np

G = 9.807


# ---------------------------------------------------------------------------
# TABLAS DE COEFICIENTES
# ---------------------------------------------------------------------------

# Coeficiente β de Lischtvan-Lebediev según periodo de retorno de la avenida.
# (Tablas estándar de Maza Álvarez — confirmar valores exactos contra la fuente.)
_BETA_TR = {
    1: 0.77, 2: 0.82, 5: 0.86, 10: 0.90, 20: 0.94,
    50: 0.97, 100: 1.00, 300: 1.03, 500: 1.05, 1000: 1.07,
}

# Exponente x — suelos NO cohesivos, en función del diámetro medio d_m (mm).
_X_NO_COHESIVO = {
    0.05: 0.43, 0.15: 0.42, 0.50: 0.41, 1.0: 0.40, 1.5: 0.39,
    2.5: 0.38, 4.0: 0.37, 6.0: 0.36, 8.0: 0.35, 10: 0.34,
    15: 0.33, 20: 0.32, 25: 0.31, 40: 0.30, 60: 0.29,
    90: 0.28, 140: 0.27, 190: 0.26, 250: 0.25, 310: 0.24,
    370: 0.23, 450: 0.22, 570: 0.21, 750: 0.20, 1000: 0.19,
}

# Exponente x — suelos COHESIVOS, en función del peso específico seco γs (t/m³).
_X_COHESIVO = {
    0.80: 0.52, 0.86: 0.50, 0.90: 0.48, 0.96: 0.46, 1.00: 0.44,
    1.08: 0.42, 1.16: 0.40, 1.24: 0.38, 1.34: 0.36, 1.40: 0.35,
    1.46: 0.34, 1.52: 0.33, 1.58: 0.32, 1.64: 0.31, 1.71: 0.30,
    1.80: 0.29, 1.89: 0.28, 2.00: 0.27,
}

# Factor K1 por forma de la nariz de la pila (CSU / HEC-18).
PIER_SHAPE_K1 = {
    "cuadrada": 1.1,
    "redondeada": 1.0,
    "circular": 1.0,
    "grupo_cilindros": 1.0,
    "triangular": 0.9,   # nariz triangular afilada
}

# Factor K3 por condición del lecho (CSU / HEC-18).
PIER_BED_K3 = {
    "plano": 1.1,        # lecho plano o dunas pequeñas
    "dunas_medianas": 1.1,
    "dunas_grandes": 1.3,
    "antidunas": 1.1,
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
    """Interpolación lineal sobre una tabla {x: y} ordenada por clave."""
    xs = sorted(table.keys())
    if key <= xs[0]:
        return table[xs[0]]
    if key >= xs[-1]:
        return table[xs[-1]]
    return float(np.interp(key, xs, [table[x] for x in xs]))


def beta_periodo_retorno(tr_anios):
    """Coeficiente β de Lischtvan-Lebediev según el periodo de retorno (años)."""
    return _interp_table(_BETA_TR, float(tr_anios))


def exponente_x(dm_mm=None, gamma_s=None, cohesivo=False):
    """
    Exponente x de Lischtvan-Lebediev.
    - No cohesivo: función de d_m (mm).
    - Cohesivo:    función de γs (t/m³).
    """
    if cohesivo:
        if gamma_s is None:
            raise ValueError("Suelo cohesivo requiere gamma_s (t/m³)")
        return _interp_table(_X_COHESIVO, float(gamma_s))
    if dm_mm is None:
        raise ValueError("Suelo no cohesivo requiere dm (mm)")
    return _interp_table(_X_NO_COHESIVO, float(dm_mm))


# ===========================================================================
# 1. SOCAVACIÓN GENERAL
# ===========================================================================

def lacey(Q, dm_mm):
    """
    Lacey (1930) — cauces aluviales con lecho de arena.
        h_ms = 0.389 * (Q^(1/3) / d_m^(1/6))
    Q en m³/s ; d_m en mm ; devuelve profundidad media de socavación (m).
    """
    return 0.389 * (Q ** (1.0 / 3.0) / dm_mm ** (1.0 / 6.0))


def blench(q, d50_mm):
    """
    Blench (1969). q = caudal unitario (m³/s·m) ; d50 en mm.
        arenas (0.06 < d50 < 2 mm): h_ms = 1.200 * q^(2/3) / d50^(1/6)
        grava gruesa (d50 > 2 mm):  h_ms = 1.230 * q^(2/3) / d50^(1/12)
    """
    if d50_mm > 2.0:
        return 1.230 * q ** (2.0 / 3.0) / d50_mm ** (1.0 / 12.0)
    return 1.200 * q ** (2.0 / 3.0) / d50_mm ** (1.0 / 6.0)


def maza_echavarria(Q, B, d50_m):
    """
    Maza & Echavarría (1973) — calibrada con ríos sudamericanos.
        h_ms = 0.365 * Q^0.784 / (B^0.784 * d50^0.157)
    Q en m³/s ; B (ancho del espejo) en m ; d50 en METROS.
    """
    return 0.365 * (Q ** 0.784) / (B ** 0.784 * d50_m ** 0.157)


def lischtvan_lebediev_reducida(q, d50_m):
    """
    Lischtvan-Lebediev (forma reducida, rango de arenas).
        h_ms = 0.333 * q^0.710 / d50^0.199
    q = caudal unitario (m³/s·m) ; d50 en METROS.
    """
    return 0.333 * q ** 0.710 / d50_m ** 0.199


def kellerhals(q, d90_m):
    """
    Kellerhals — lechos de grava gruesa.
        h_ms = 0.470 * q^0.800 / d90^0.120
    q = caudal unitario local (m³/s·m) ; d90 en METROS.
    """
    return 0.470 * q ** 0.800 / d90_m ** 0.120


def lischtvan_lebediev_maza(Ho, Hm, Q, We, mu, beta, cohesivo=False,
                            dm_mm=None, gamma_s=None):
    """
    Lischtvan-Lebediev / Maza (1968) — método completo por franjas.

    Coeficiente de distribución de gasto:
        α = Q / (Hm^(5/3) · We · μ)

    Profundidad total socavada Hs (medida desde la superficie):
        No cohesivo: Hs = [ α·Ho^(5/3) / (0.68·β·d_m^0.28) ]^(1/(1+x))
        Cohesivo:    Hs = [ α·Ho^(5/3) / (0.60·β·γs^1.18) ]^(1/(1+x))

    Parámetros:
        Ho     tirante inicial de la franja (m)
        Hm     tirante medio de aproximación (m)
        Q      caudal de diseño (m³/s)
        We     ancho efectivo (m)
        mu     factor de contracción (μ)
        beta   coeficiente por periodo de retorno
        dm_mm  d_m en mm (suelo no cohesivo)
        gamma_s peso específico seco (t/m³) (suelo cohesivo)

    Devuelve (Hs, socavacion) con socavacion = Hs - Ho (m).
    """
    alpha = Q / (Hm ** (5.0 / 3.0) * We * mu)
    x = exponente_x(dm_mm=dm_mm, gamma_s=gamma_s, cohesivo=cohesivo)
    exp = 1.0 / (1.0 + x)
    if cohesivo:
        denom = 0.60 * beta * gamma_s ** 1.18
    else:
        denom = 0.68 * beta * dm_mm ** 0.28
    Hs = (alpha * Ho ** (5.0 / 3.0) / denom) ** exp
    return Hs, max(Hs - Ho, 0.0)


# ===========================================================================
# 2. SOCAVACIÓN POR CONTRACCIÓN
# ===========================================================================

def k1_contraccion(u_star, w):
    """
    Exponente k1 de Laursen según la relación u*/w (modo de transporte).
        < 0.5   -> 0.59 (arrastre de fondo)
        0.5-2.0 -> 0.64 (inicio de suspensión)
        > 2.0   -> 0.69 (suspensión)
    """
    if w <= 0:
        return 0.64
    r = u_star / w
    if r < 0.5:
        return 0.59
    if r <= 2.0:
        return 0.64
    return 0.69


def laursen_contraccion_lecho_vivo(y1, W1, W2, k1):
    """
    Laursen — contracción en lecho móvil (live-bed).
        y2 = y1 · (W1/W2)^k1
    Devuelve (y2, socavacion = y2 - y1).
    """
    y2 = y1 * (W1 / W2) ** k1
    return y2, max(y2 - y1, 0.0)


def parker_contraccion(y1, W1, W2):
    """
    Parker (1985) — cauces de grava, sin acorazamiento (conservadora).
        y2 = y1 · (W1/W2)^0.825
    """
    y2 = y1 * (W1 / W2) ** 0.825
    return y2, max(y2 - y1, 0.0)


def hec18_contraccion_agua_clara(Q2, dm_m, W2, Ku=0.025):
    """
    Richardson & Davis / HEC-18 — contracción en agua clara (S.I.).

    Forma completa de HEC-18 (la que devolvió el cuaderno vino truncada):
        y2 = [ Ku · Q2² / (D_m^(2/3) · W2²) ]^(3/7)
    con Ku = 0.025 (S.I.).
    Q2 caudal en la contracción (m³/s) ; d_m en METROS ; W2 ancho (m).
    Devuelve y2 (tirante tras la socavación, m).

    NOTA: confirmar coeficiente y forma exacta contra el PDF de HEC-18.
    """
    return (Ku * Q2 ** 2 / (dm_m ** (2.0 / 3.0) * W2 ** 2)) ** (3.0 / 7.0)


# ===========================================================================
# 3. SOCAVACIÓN LOCAL EN PILAS
# ===========================================================================

def breusers_nicollet_shen(b, y, Ks=1.0, Ktheta=1.0):
    """
    Breusers, Nicollet & Shen (1977) — inicio del movimiento.
        d_s = b · 2.00 · tanh(y/b) · Ks · Kθ
    b = ancho de la pila (m) ; y = tirante (m).
    """
    return b * 2.00 * np.tanh(y / b) * Ks * Ktheta


def csu_hec18_pila(a, y1, Fr1, K1=1.0, K2=1.0, K3=1.1, K4=1.0):
    """
    CSU / HEC-18 (Richardson & Davis, 1995) — método normativo.
        y_s = 2.0 · K1·K2·K3·K4 · a^0.65 · y1^0.35 · Fr1^0.43
    a = ancho de pila (m) ; y1 = calado aguas arriba (m) ; Fr1 = Froude.
    """
    return 2.0 * K1 * K2 * K3 * K4 * a ** 0.65 * y1 ** 0.35 * Fr1 ** 0.43


def k2_angulo_ataque(theta_deg, largo_ancho):
    """
    Factor K2 por ángulo de ataque (HEC-18):
        K2 = (cos θ + (L/a)·sen θ)^0.65
    theta en grados ; largo_ancho = L/a (largo de pila / ancho).
    """
    th = np.radians(theta_deg)
    return (np.cos(th) + largo_ancho * np.sin(th)) ** 0.65


# ===========================================================================
# 4. SOCAVACIÓN LOCAL EN ESTRIBOS
# ===========================================================================

def liu_dodge_skinner(L, y, Fr, spill_through=True):
    """
    Liu, Dodge & Skinner (1961) — lecho móvil.
        spill-through (con talud): d_s = L · 1.10 · (y/L)^0.6 · Fr^0.33
        vertical / con aletas:     d_s = L · 2.15 · (y/L)^0.6 · Fr^0.33
    L = longitud del estribo proyectada al flujo (m).
    """
    c = 1.10 if spill_through else 2.15
    return L * c * (y / L) ** 0.6 * Fr ** 0.33


def laursen_estribo(L, y, agua_clara=False):
    """
    Laursen (1962/1963) — estribo sobre el cauce principal.
        lecho móvil (1962): d_s = L · 1.57 · (y/L)^0.5
        agua clara  (1963): d_s = L · 1.89 · (y/L)^0.5
    """
    c = 1.89 if agua_clara else 1.57
    return L * c * (y / L) ** 0.5


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

    `params` (dict) puede contener:
        d50_mm, d90_mm, dm_mm      granulometría
        depth (y = Ho), velocity   hidráulica de campo
        u_star, w (ws)             para el exponente de contracción
        froude                     Froude de aproximación
        Q                          caudal de diseño (m³/s)
        B / W1                     ancho del espejo (m)
        W2                         ancho contraído (m)
        tr                         periodo de retorno (años)
        cohesivo (bool), gamma_s   tipo de suelo
        mu                         factor de contracción μ (default 1.0)
        pila_ancho, pila_forma, pila_theta, pila_largo
        estribo_L, estribo_spill (bool), estribo_agua_clara (bool)

    Devuelve un dict con listas por tipo, listo para el frontend.
    """
    d50 = params.get("d50_mm")
    d90 = params.get("d90_mm")
    dm = params.get("dm_mm", d50)          # d_m ≈ d50 si no se da explícito
    y = params.get("depth")
    v = params.get("velocity")
    Q = params.get("Q")
    B = params.get("B") or params.get("W1")
    W2 = params.get("W2")
    tr = params.get("tr", 100)
    cohesivo = bool(params.get("cohesivo", False))
    gamma_s = params.get("gamma_s")
    mu = params.get("mu", 1.0)
    u_star = params.get("u_star")
    w = params.get("w")

    out = {"general": [], "contraccion": [], "pila": [], "estribo": [],
           "yarnell": YARNELL_K, "avisos": []}

    # Caudal unitario q (m³/s·m) — de Q/B si hay ambos, si no de v·y
    q = None
    if Q and B and B > 0:
        q = Q / B
    elif v and y:
        q = v * y

    # --- 1. GENERAL ---
    if Q and dm:
        out["general"].append(_row(
            "Lacey (1930)", lacey(Q, dm), "m",
            "Prof. media desde la superficie · d_m en mm", "Lacey 1930"))
    if q and d50:
        out["general"].append(_row(
            "Blench (1969)", blench(q, d50), "m",
            "Sin acorazamiento · d50 en mm", "Blench 1969"))
    if Q and B and d50:
        out["general"].append(_row(
            "Maza & Echavarría (1973)", maza_echavarria(Q, B, d50 / 1000.0),
            "m", "Calibrada en ríos sudamericanos · d50 en m", "Maza-Echavarría 1973"))
    if q and d50:
        out["general"].append(_row(
            "Lischtvan-Lebediev (reducida)",
            lischtvan_lebediev_reducida(q, d50 / 1000.0), "m",
            "Rango de arenas · d50 en m", "Lischtvan-Lebediev"))
    if q and d90:
        out["general"].append(_row(
            "Kellerhals", kellerhals(q, d90 / 1000.0), "m",
            "Grava gruesa · d90 en m", "Kellerhals"))

    # Lischtvan-Lebediev / Maza completo (por franjas)
    if y and Q and B:
        try:
            beta = beta_periodo_retorno(tr)
            Hs, sc = lischtvan_lebediev_maza(
                Ho=y, Hm=y, Q=Q, We=B, mu=mu, beta=beta,
                cohesivo=cohesivo, dm_mm=dm, gamma_s=gamma_s)
            tipo = "cohesivo" if cohesivo else "no cohesivo"
            out["general"].append(_row(
                "Lischtvan-Lebediev / Maza (completo)", sc, "m",
                f"Hs={Hs:.2f} m desde superficie · suelo {tipo} · β={beta:.2f}"
                f" (Tr={tr} a) · μ={mu}", "Maza 1968"))
        except Exception as e:
            out["avisos"].append(f"L-L/Maza no calculable: {e}")

    # --- 2. CONTRACCIÓN ---
    if y and B and W2 and W2 > 0:
        k1 = k1_contraccion(u_star or 0.0, w or 1.0)
        y2, sc = laursen_contraccion_lecho_vivo(y, B, W2, k1)
        out["contraccion"].append(_row(
            "Laursen (lecho vivo)", sc, "m",
            f"y2={y2:.2f} m · k1={k1}", "Laursen"))
        y2p, scp = parker_contraccion(y, B, W2)
        out["contraccion"].append(_row(
            "Parker (1985)", scp, "m",
            f"y2={y2p:.2f} m · grava sin acorazar", "Parker 1985"))
        if Q and dm:
            y2c = hec18_contraccion_agua_clara(Q, dm / 1000.0, W2)
            out["contraccion"].append(_row(
                "HEC-18 (agua clara)", max(y2c - y, 0.0), "m",
                f"y2={y2c:.2f} m · d_m en m", "Richardson & Davis / HEC-18"))

    # --- 3. PILAS ---
    b = params.get("pila_ancho")
    if b and y:
        forma = params.get("pila_forma", "circular")
        K1 = PIER_SHAPE_K1.get(forma, 1.0)
        theta = params.get("pila_theta", 0.0)
        largo = params.get("pila_largo")
        K2 = 1.0
        if theta and largo and b > 0:
            K2 = k2_angulo_ataque(theta, largo / b)
        K3 = PIER_BED_K3.get(params.get("pila_lecho", "plano"), 1.1)
        out["pila"].append(_row(
            "Breusers-Nicollet-Shen (1977)",
            breusers_nicollet_shen(b, y, Ks=K1, Ktheta=1.0), "m",
            f"Ks={K1} · forma {forma}", "Breusers et al. 1977"))
        if params.get("froude"):
            out["pila"].append(_row(
                "CSU / HEC-18 (normativo)",
                csu_hec18_pila(b, y, params["froude"], K1=K1, K2=K2, K3=K3),
                "m", f"K1={K1} K2={K2:.2f} K3={K3} · Fr={params['froude']}",
                "Richardson & Davis 1995 / HEC-18"))

    # --- 4. ESTRIBOS ---
    L = params.get("estribo_L")
    if L and y and L > 0:
        spill = bool(params.get("estribo_spill", True))
        agua_clara = bool(params.get("estribo_agua_clara", False))
        if params.get("froude"):
            out["estribo"].append(_row(
                "Liu-Dodge-Skinner (1961)",
                liu_dodge_skinner(L, y, params["froude"], spill_through=spill),
                "m", f"{'spill-through' if spill else 'vertical'} · lecho móvil",
                "Liu et al. 1961"))
        out["estribo"].append(_row(
            "Laursen (estribo)", laursen_estribo(L, y, agua_clara=agua_clara),
            "m", f"{'agua clara' if agua_clara else 'lecho móvil'}", "Laursen 1962/63"))

    return out
