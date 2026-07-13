import pytest
from models import socavacion as sc


# ── SOCAVACIÓN GENERAL ──────────────────────────────────────────────────────

def test_lacey():
    # 0.389 * (120^(1/3) / 0.5^(1/6))
    assert sc.lacey(120.0, 0.5) == pytest.approx(2.154, rel=1e-3)


def test_blench_arena_vs_grava():
    arena = sc.blench(3.0, 0.5)      # d50 < 2 mm -> exponente 1/6
    grava = sc.blench(3.0, 5.0)      # d50 > 2 mm -> exponente 1/12
    assert arena == pytest.approx(2.802, rel=1e-3)
    assert grava > 0
    # La grava (mayor d50) socava menos que la arena con igual q
    assert grava < arena


def test_maza_echavarria_unidades_metros():
    # d50 en metros
    h = sc.maza_echavarria(120.0, 40.0, 0.0005)
    assert h > 0


def test_kellerhals():
    assert sc.kellerhals(3.0, 0.002) > 0


def test_lischtvan_lebediev_maza_no_cohesivo():
    Hs, socav = sc.lischtvan_lebediev_maza(
        Ho=2.0, Hm=2.0, Q=120.0, We=40.0, mu=1.0,
        beta=1.0, cohesivo=False, dm_mm=0.5)
    assert Hs > 2.0            # el fondo baja respecto al tirante inicial
    assert socav == pytest.approx(Hs - 2.0, rel=1e-6)


def test_lischtvan_lebediev_maza_cohesivo_requiere_gamma():
    with pytest.raises(ValueError):
        sc.lischtvan_lebediev_maza(
            Ho=2.0, Hm=2.0, Q=120.0, We=40.0, mu=1.0,
            beta=1.0, cohesivo=True, gamma_s=None)


# ── COEFICIENTES ────────────────────────────────────────────────────────────

def test_beta_interpolacion_periodo_retorno():
    assert sc.beta_periodo_retorno(100) == pytest.approx(1.00, rel=1e-6)
    assert sc.beta_periodo_retorno(2) == pytest.approx(0.82, rel=1e-6)
    # interpolación entre 100 y 300
    assert 1.00 < sc.beta_periodo_retorno(200) < 1.03


def test_exponente_x_no_cohesivo_decrece_con_diametro():
    x_fino = sc.exponente_x(dm_mm=0.5)
    x_grueso = sc.exponente_x(dm_mm=100)
    assert x_fino > x_grueso


def test_k1_contraccion_por_modo_transporte():
    assert sc.k1_contraccion(0.02, 0.10) == 0.59   # u*/w < 0.5
    assert sc.k1_contraccion(0.06, 0.06) == 0.64   # 0.5-2.0
    assert sc.k1_contraccion(0.30, 0.10) == 0.69   # > 2.0


# ── CONTRACCIÓN ─────────────────────────────────────────────────────────────

def test_laursen_contraccion_estrecha_aumenta_socavacion():
    y2, socav = sc.laursen_contraccion_lecho_vivo(2.0, 40.0, 30.0, 0.64)
    assert y2 > 2.0
    assert socav == pytest.approx(y2 - 2.0, rel=1e-6)


def test_hec18_agua_clara_valor_razonable():
    # Debe dar metros de tirante, no decenas de metros
    y2 = sc.hec18_contraccion_agua_clara(120.0, 0.0005, 30.0)
    assert 3.0 < y2 < 8.0


# ── PILAS ───────────────────────────────────────────────────────────────────

def test_csu_hec18_pila():
    ys = sc.csu_hec18_pila(a=1.5, y1=2.0, Fr1=0.34, K1=1.0, K2=1.0, K3=1.1)
    assert ys > 0


def test_forma_pila_cuadrada_socava_mas_que_triangular():
    cuad = sc.csu_hec18_pila(1.5, 2.0, 0.34, K1=sc.PIER_SHAPE_K1["cuadrada"])
    tri = sc.csu_hec18_pila(1.5, 2.0, 0.34, K1=sc.PIER_SHAPE_K1["triangular"])
    assert cuad > tri


def test_breusers_nicollet_shen():
    ds = sc.breusers_nicollet_shen(b=1.5, y=2.0)
    assert ds > 0


# ── ESTRIBOS ────────────────────────────────────────────────────────────────

def test_liu_vertical_socava_mas_que_spill():
    vert = sc.liu_dodge_skinner(8.0, 2.0, 0.34, spill_through=False)
    spill = sc.liu_dodge_skinner(8.0, 2.0, 0.34, spill_through=True)
    assert vert > spill


def test_laursen_estribo_agua_clara_mayor():
    clara = sc.laursen_estribo(8.0, 2.0, agua_clara=True)
    movil = sc.laursen_estribo(8.0, 2.0, agua_clara=False)
    assert clara > movil


# ── ORQUESTADOR ─────────────────────────────────────────────────────────────

def test_compute_socavacion_completo():
    params = dict(
        d50_mm=0.5, d90_mm=2.0, depth=2.0, velocity=1.5,
        Q=120.0, B=40.0, W2=30.0, tr=100, u_star=0.12, w=0.06,
        froude=0.34, pila_ancho=1.5, pila_forma="circular",
        pila_theta=10, pila_largo=6.0, estribo_L=8.0, estribo_spill=True)
    r = sc.compute_socavacion(params)
    assert len(r["general"]) >= 5
    assert len(r["contraccion"]) >= 2
    assert len(r["pila"]) >= 2
    assert len(r["estribo"]) >= 2
    assert r["avisos"] == []
    # todos los valores positivos y finitos
    for grupo in ("general", "contraccion", "pila", "estribo"):
        for row in r[grupo]:
            assert row["valor"] is not None and row["valor"] > 0


def test_compute_socavacion_sin_datos_puente():
    # Solo granulometría e hidráulica -> métodos generales por q, nada de puente
    params = dict(d50_mm=0.5, d90_mm=2.0, depth=2.0, velocity=1.5)
    r = sc.compute_socavacion(params)
    assert r["pila"] == []
    assert r["estribo"] == []
    # q = v*y disponible -> Blench, L-L reducida, Kellerhals sí salen
    assert len(r["general"]) >= 1
