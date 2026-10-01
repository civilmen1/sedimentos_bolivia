from datetime import date

from utils.gee_handler import (HIST_YEARS, hist_epochs, hist_last_year,
                               hist_sensors, mann_kendall_sen)


def test_last_year_uses_complete_dry_season():
    assert hist_last_year(date(2026, 10, 1)) == 2026
    assert hist_last_year(date(2026, 7, 15)) == 2025


def test_epochs_span_twenty_years():
    ep = hist_epochs(2026)
    assert ep == [2006, 2011, 2016, 2021, 2026]
    assert ep[-1] - ep[0] == HIST_YEARS == 20


def test_sensors_by_year():
    assert hist_sensors(2006)[0] == "TM 5"
    assert hist_sensors(2012) == ["ETM+ 7 (SLC-off)"]
    assert hist_sensors(2016) == ["OLI 8"]
    assert hist_sensors(2024) == ["OLI 8", "OLI-2 9"]


def test_mann_kendall_detects_decline():
    yrs = list(range(2006, 2027))
    t = mann_kendall_sen(yrs, [0.5 - 0.01 * (y - 2006) for y in yrs])
    assert t["tendencia"] == "decreciente"
    assert abs(t["sen"] + 0.01) < 1e-9
    assert t["p"] < 0.05


def test_mann_kendall_no_trend_and_short_series():
    yrs = list(range(2006, 2026))
    vals = [0.3 + (0.02 if i % 2 else -0.02) for i in range(len(yrs))]
    assert mann_kendall_sen(yrs, vals)["tendencia"] == "sin tendencia significativa"
    assert mann_kendall_sen([2020, 2021, 2022], [1, 2, 3]) is None
    # años faltantes (None) se ignoran
    assert mann_kendall_sen(yrs, [None] * 15 + vals[:5]) is None


def test_summarize_trend_cells_aggregates_by_band_trend_priority():
    from utils.gee_handler import summarize_trend_cells
    # banda*100 + tendencia*10 + prioridad -> ha
    cells = {0 * 100 + 0 * 10 + 0: 30.0,   # baja, disminución fuerte, crítica
             2 * 100 + 1 * 10 + 1: 20.0,   # alta, disminución, alta
             1 * 100 + 2 * 10 + 4: 40.0,   # media, estable, baja
             2 * 100 + 4 * 10 + 3: 10.0}   # alta, incremento fuerte, conservación
    r = summarize_trend_cells(cells)
    assert r["total_ha"] == 100.0
    assert r["decl_pct"] == 50.0 and r["stab_pct"] == 40.0 and r["incr_pct"] == 10.0
    assert r["crit_pct"] == 50.0 and r["cons_pct"] == 10.0
    assert r["band_trend"][0][0] == 30.0 and r["band_prio"][2][3] == 10.0
    assert r["band_max_decl"] == "Cuenca baja"
    assert summarize_trend_cells({}) is None


def test_ndvi_report_section_renders_without_data():
    import app
    with app.app.test_request_context():
        html = app.app.jinja_env.get_template("_ndvi_informe.html").render(
            results={"lat": -17.4, "lon": -66.1, "area_name": ""})
    assert "8.8.13 Recomendaciones" in html
    assert "no se obtuvieron datos de NDVI interanual" in html
