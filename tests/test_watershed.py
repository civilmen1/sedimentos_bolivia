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
    # Máscara de drenaje booleana con al menos un cauce
    assert wd["stream_mask"].dtype == bool
    assert wd["stream_mask"].any()
    # Extent (lon_min, lon_max, lat_min, lat_max)
    lon_min, lon_max, lat_min, lat_max = wd["extent"]
    assert lon_min < lon_max
    assert lat_min < lat_max


def test_synthetic_watershed_is_deterministic():
    """La misma coordenada produce la misma cuenca (semilla estable)."""
    a = app_module._synthetic_watershed(-17.0, -66.0, 15.0)
    b = app_module._synthetic_watershed(-17.0, -66.0, 15.0)
    assert a["boundary"] == b["boundary"]
    assert np.array_equal(a["stream_mask"], b["stream_mask"])


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


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
