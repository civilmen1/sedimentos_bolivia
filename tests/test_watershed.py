"""
Pruebas para la generación de la cuenca hidrográfica, la red de drenaje
y los mapas temáticos con overlay de cuenca.

GEE no está autenticado en el entorno de pruebas, por lo que se ejercita
la ruta de fallback sintético (la que se ejecuta cuando gee_ready() == False).
"""
import numpy as np
import pytest

import app as app_module


def test_synthetic_watershed_structure():
    """El fallback sintético devuelve la estructura esperada."""
    wd = app_module._synthetic_watershed(-16.5, -68.15, 15.0)

    assert wd["is_real"] is False
    # Límite de cuenca cerrado (primer punto == último)
    assert len(wd["boundary"]) >= 4
    assert wd["boundary"][0] == wd["boundary"][-1]
    # Cada vértice es [lon, lat]
    for lon, lat in wd["boundary"]:
        assert isinstance(lon, float)
        assert isinstance(lat, float)
    # Cauce principal y afluentes como polilíneas lon/lat
    assert len(wd["channel"]) >= 2
    for lon, lat in wd["channel"]:
        assert isinstance(lon, float) and isinstance(lat, float)
    assert isinstance(wd["tributaries"], list)
    assert len(wd["tributaries"]) >= 1
    for trib in wd["tributaries"]:
        assert len(trib) >= 2


def test_synthetic_watershed_is_deterministic():
    """La misma coordenada produce la misma cuenca (semilla estable)."""
    a = app_module._synthetic_watershed(-17.0, -66.0, 15.0)
    b = app_module._synthetic_watershed(-17.0, -66.0, 15.0)
    assert a["boundary"] == b["boundary"]
    assert a["channel"] == b["channel"]
    assert a["tributaries"] == b["tributaries"]


def test_get_watershed_overlay_caches(monkeypatch):
    """get_watershed_overlay usa caché y fallback sintético sin GEE."""
    monkeypatch.setattr(app_module, "gee_ready", lambda: False)
    app_module._WATERSHED_CACHE.clear()

    wd1 = app_module.get_watershed_overlay(-15.0, -67.0, 15.0)
    assert wd1["is_real"] is False
    # Segunda llamada devuelve exactamente el mismo objeto (caché)
    wd2 = app_module.get_watershed_overlay(-15.0, -67.0, 15.0)
    assert wd1 is wd2


def test_generate_watershed_map_returns_png(monkeypatch):
    """El mapa de cuenca se genera como PNG base64."""
    monkeypatch.setattr(app_module, "gee_ready", lambda: False)
    app_module._MAP_CACHE.clear()
    app_module._WATERSHED_CACHE.clear()

    png = app_module.generate_watershed_map(-16.5, -68.15, 15.0)
    assert png.startswith("data:image/png;base64,")
    assert len(png) > 5000   # imagen no trivial


def test_generate_all_thematic_maps_has_watershed_first(monkeypatch):
    """Los 9 mapas se generan y 'watershed' es el primero."""
    monkeypatch.setattr(app_module, "gee_ready", lambda: False)
    app_module._MAP_CACHE.clear()
    app_module._WATERSHED_CACHE.clear()

    maps = app_module.generate_all_thematic_maps(-16.5, -68.15, 15.0)
    keys = list(maps.keys())

    assert keys[0] == "watershed"
    expected = {"watershed", "dem", "slope", "ndvi", "ndwi",
                "ndti", "manning", "risk", "jrc"}
    assert set(keys) == expected
    for mt, png in maps.items():
        assert png.startswith("data:image/png;base64,"), mt


def test_map_metadata_includes_watershed():
    """Los diccionarios de metadatos incluyen la entrada watershed."""
    assert "watershed" in app_module.MAP_TITLES
    assert "watershed" in app_module.MAP_SOURCES
    assert "watershed" in app_module.MAP_LEGEND_LABELS
    # watershed debe ser la primera clave en MAP_TITLES
    assert list(app_module.MAP_TITLES.keys())[0] == "watershed"


def test_cartographic_map_accepts_watershed_overlay(monkeypatch):
    """generate_cartographic_map dibuja el overlay sin errores."""
    monkeypatch.setattr(app_module, "gee_ready", lambda: False)
    wd = app_module._synthetic_watershed(-16.5, -68.15, 15.0)
    png = app_module.generate_cartographic_map(
        -16.5, -68.15, "dem", radius_km=15.0,
        data_array=app_module._synthetic_data("dem", lat=-16.5, lon=-68.15),
        watershed_data=wd,
    )
    assert png.startswith("data:image/png;base64,")


def test_utm_epsg_southern_and_northern():
    """utm_epsg devuelve el huso correcto en ambos hemisferios."""
    # La Paz, Bolivia → UTM 19S → EPSG:32719
    assert app_module.utm_epsg(-16.5, -68.15) == 32719
    # Hemisferio norte mismo huso → EPSG:32619
    assert app_module.utm_epsg(16.5, -68.15) == 32619


def test_delineate_watershed_from_dem_synthetic():
    """
    El geoproceso pysheds delimita una cuenca a partir de un DEM en UTM y
    devuelve geometrías en lon/lat alrededor del punto de muestreo.
    """
    import numpy as _np
    import pyproj
    from affine import Affine

    lat, lon = -16.5, -68.15
    epsg = app_module.utm_epsg(lat, lon)
    to_utm = pyproj.Transformer.from_crs(4326, epsg, always_xy=True)
    xc, yc = to_utm.transform(lon, lat)

    scale, n = 18.75, 600
    half = n * scale / 2.0
    transform = Affine(scale, 0, xc - half, 0, -scale, yc + half)

    cols, rows = _np.meshgrid(_np.arange(n), _np.arange(n))
    chan_col = n / 2 + 40 * _np.sin(rows / 90.0)
    dem = (_np.abs(cols - chan_col) * 1.5 + (n - rows) * 0.8 + 3000.0)

    res = app_module.delineate_watershed_from_dem(dem, transform, epsg, lat, lon)
    assert res is not None
    assert res["is_real"] is True
    assert len(res["boundary"]) >= 4
    assert len(res["channel"]) >= 2
    # El parteaguas debe rodear aproximadamente el punto de muestreo
    lons = [p[0] for p in res["boundary"]]
    lats = [p[1] for p in res["boundary"]]
    assert min(lons) - 0.05 <= lon <= max(lons) + 0.05
    assert min(lats) - 0.05 <= lat <= max(lats) + 0.05

    # Morfometría calculada con el hydro-tool
    m = res["morphometry"]
    assert m is not None
    for k in ("area_km2", "perimeter_km", "channel_length_km",
              "channel_slope_mm", "drainage_density", "gravelius_kc",
              "elev_min_m", "elev_max_m", "tc_kirpich_min"):
        assert k in m and m[k] is not None, k
    assert m["area_km2"] > 0
    assert m["elev_max_m"] >= m["elev_min_m"]
    assert m["channel_slope_mm"] > 0


def test_watershed_model_inputs_fallback(monkeypatch):
    """Sin GEE, watershed_model_inputs conserva la pendiente del usuario."""
    monkeypatch.setattr(app_module, "gee_ready", lambda: False)
    app_module._WATERSHED_CACHE.clear()
    wmi = app_module.watershed_model_inputs(-16.5, -68.15, 0.005)
    assert wmi["slope"] == 0.005
    assert wmi["watershed_is_real"] is False
    assert wmi["morphometry"] is None
    assert "usuario" in wmi["slope_source"]


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
