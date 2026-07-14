"""
Google Earth Engine handler — autenticación por cuenta de servicio y
obtención de mapas temáticos renderizados (getThumbURL) como arrays RGB.

Autenticación (en orden de prioridad):
  1. EE_SERVICE_ACCOUNT_JSON  → contenido JSON completo de la clave de cuenta
                                de servicio (recomendado, se configura como
                                "secret" del Space).
  2. EE_SERVICE_ACCOUNT_FILE + EE_SERVICE_ACCOUNT_EMAIL → ruta a archivo .json.
  3. ee.Initialize() por defecto (credenciales personales locales).
"""
import ee
import os
import io
import json
import math

_GEE_READY = False
_GEE_ERROR = None          # último error de inicialización (para diagnóstico)
_GEE_EMAIL = None          # cuenta de servicio usada
_GEE_PROJECT = None        # proyecto de Google Cloud usado
_GEE_METHOD = None         # método de autenticación que tuvo éxito


def initialize_gee():
    """Inicializa GEE con cuenta de servicio o credenciales por defecto."""
    global _GEE_READY, _GEE_ERROR, _GEE_EMAIL, _GEE_PROJECT, _GEE_METHOD
    _GEE_ERROR = None
    try:
        sa_json = (os.environ.get("EE_SERVICE_ACCOUNT_JSON")
                   or os.environ.get("GEE_SERVICE_ACCOUNT_JSON"))
        if sa_json:
            try:
                info = json.loads(sa_json)
            except json.JSONDecodeError as je:
                raise ValueError(
                    "EE_SERVICE_ACCOUNT_JSON no es un JSON válido: "
                    f"{je}. Pega el contenido completo del archivo .json de la "
                    "cuenta de servicio (incluyendo las llaves { }).")
            email = info.get("client_email")
            project = info.get("project_id")
            if not email:
                raise ValueError(
                    "El JSON de la cuenta de servicio no contiene 'client_email'.")
            creds = ee.ServiceAccountCredentials(email, key_data=sa_json)
            if project:
                ee.Initialize(creds, project=project)
            else:
                ee.Initialize(creds)
            print(f"GEE inicializado con cuenta de servicio: {email}")
            _GEE_READY = True
            _GEE_EMAIL = email
            _GEE_PROJECT = project
            _GEE_METHOD = "service_account_json"
            return True

        key_file = os.environ.get("EE_SERVICE_ACCOUNT_FILE")
        sa_email = os.environ.get("EE_SERVICE_ACCOUNT_EMAIL")
        if key_file and sa_email:
            creds = ee.ServiceAccountCredentials(sa_email, key_file=key_file)
            ee.Initialize(creds)
            print(f"GEE inicializado con archivo de cuenta de servicio: {sa_email}")
            _GEE_READY = True
            _GEE_EMAIL = sa_email
            _GEE_METHOD = "service_account_file"
            return True

        project = os.environ.get("GEE_PROJECT")
        if project:
            ee.Initialize(project=project)
        else:
            ee.Initialize()
        print("GEE inicializado con credenciales por defecto.")
        _GEE_READY = True
        _GEE_PROJECT = project
        _GEE_METHOD = "default_credentials"
        return True
    except Exception as e:
        print(f"GEE Initialization failed: {e}")
        _GEE_READY = False
        _GEE_ERROR = str(e)
        return False


def gee_ready():
    return _GEE_READY


def gee_status(probe=False):
    """
    Diagnóstico del estado de Google Earth Engine.

    Si probe=True realiza una consulta mínima en vivo para verificar que la
    cuenta tiene acceso real a datos (no solo que la inicialización no falló).
    """
    has_json = bool(os.environ.get("EE_SERVICE_ACCOUNT_JSON")
                    or os.environ.get("GEE_SERVICE_ACCOUNT_JSON"))
    has_file = bool(os.environ.get("EE_SERVICE_ACCOUNT_FILE")
                    and os.environ.get("EE_SERVICE_ACCOUNT_EMAIL"))
    status = {
        "ready": _GEE_READY,
        "method": _GEE_METHOD,
        "service_account": _GEE_EMAIL,
        "project": _GEE_PROJECT,
        "init_error": _GEE_ERROR,
        "env_EE_SERVICE_ACCOUNT_JSON_present": has_json,
        "env_service_account_file_present": has_file,
        "probe": None,
        "probe_error": None,
    }
    if probe and _GEE_READY:
        try:
            val = (ee.Image("USGS/SRTMGL1_003")
                   .reduceRegion(
                       reducer=ee.Reducer.first(),
                       geometry=ee.Geometry.Point([-68.15, -16.5]),
                       scale=30)
                   .getInfo())
            status["probe"] = "ok"
            status["probe_sample"] = val
        except Exception as e:
            status["probe"] = "failed"
            status["probe_error"] = str(e)
    return status


# ════════════════════════════════════════════════════════════════════════════
# DEFINICIÓN DE CAPAS TEMÁTICAS
# ════════════════════════════════════════════════════════════════════════════

def _build_region(lat, lon, radius_km):
    deg_lat = radius_km / 111.0
    deg_lon = radius_km / (111.0 * math.cos(math.radians(lat)))
    return ee.Geometry.Rectangle(
        [lon - deg_lon, lat - deg_lat, lon + deg_lon, lat + deg_lat]
    )


def _mask_s2_clouds(img):
    """Enmascara nubes/sombras/cirros por píxel con la banda SCL (Scene
    Classification): clases 3 = sombra de nube, 8/9 = nubes, 10 = cirro.
    Sin esto, la mediana queda contaminada y los índices (NDVI/NDWI/NDTI)
    no corresponden a valores reales, sobre todo en la Amazonía."""
    scl = img.select("SCL")
    bad = scl.eq(3).Or(scl.eq(8)).Or(scl.eq(9)).Or(scl.eq(10))
    return img.updateMask(bad.Not())


def _s2_median(region):
    """Mediana Sentinel-2 SR 2020–2023 con máscara de nubes por píxel (SCL)."""
    return (
        ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
        .filterBounds(region)
        .filterDate("2020-01-01", "2023-12-31")
        .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 40))
        .map(_mask_s2_clouds)
        .median()
    )


def _l89_col(region):
    """Landsat 8/9 Collection 2 L2 (SR) 2020–2024, nubes/sombras/cirros
    enmascarados por píxel con QA_PIXEL (bits 2,3,4)."""
    def _mask(img):
        qa = img.select("QA_PIXEL")
        bad = (qa.bitwiseAnd(1 << 3).neq(0)      # nube
               .Or(qa.bitwiseAnd(1 << 4).neq(0))  # sombra de nube
               .Or(qa.bitwiseAnd(1 << 2).neq(0)))  # cirro
        return img.updateMask(bad.Not())
    l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
    l9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
    return (l8.merge(l9)
            .filterBounds(region)
            .filterDate("2020-01-01", "2024-12-31")
            .filter(ee.Filter.lt("CLOUD_COVER", 60))
            .map(_mask))


def _index_multisource(region, s2_pair, l89_pair, name):
    """
    Índice normalizado (A−B)/(A+B) multi-fuente: Sentinel-2 (máscara SCL) +
    Landsat 8/9 C2 L2 (reflectancia = DN×0.0000275−0.2; máscara QA_PIXEL),
    mediana 2020–2024. Dos constelaciones → mediana más poblada y estable
    (menos huecos por nubes), valores físicamente reales.

    Mapeo de bandas (fórmulas estándar):
      NDVI (Rouse 1974):    S2 (B8,B4)  | L8/9 (SR_B5,SR_B4)  = (NIR−Rojo)
      NDWI (McFeeters 1996):S2 (B3,B8)  | L8/9 (SR_B3,SR_B5)  = (Verde−NIR)
      NDTI (Lacaux 2007):   S2 (B4,B3)  | L8/9 (SR_B4,SR_B3)  = (Rojo−Verde)
    """
    s2 = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
          .filterBounds(region)
          .filterDate("2020-01-01", "2024-12-31")
          .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 40))
          .map(_mask_s2_clouds)
          .map(lambda i: i.normalizedDifference(list(s2_pair)).rename(name)))

    def _l_index(img):
        a = img.select(l89_pair[0]).multiply(0.0000275).add(-0.2)
        b = img.select(l89_pair[1]).multiply(0.0000275).add(-0.2)
        return a.subtract(b).divide(a.add(b)).rename(name)

    l89 = _l89_col(region).map(_l_index)
    return s2.merge(l89).median().clip(region)


def _ndvi_multisource(region):
    return _index_multisource(region, ("B8", "B4"), ("SR_B5", "SR_B4"), "NDVI")


def _hls_composite(region, start="2020-01-01", end="2025-01-01"):
    """
    Composite HLS (Harmonized Landsat Sentinel-2, NASA) — L30 (Landsat 8/9) +
    S30 (Sentinel-2), armonizado a una respuesta espectral común, NBAR (nadir
    BRDF) y grilla 30 m. Al estar armonizado NO deja costuras entre sensores
    (que era el problema de mezclar S2 crudo + Landsat C2 en una mediana).

    Bandas renombradas a green/red/nir; NIR angosto (B05 en L30, B8A en S30)
    para consistencia espectral; escala de reflectancia ×0.0001; enmascarado
    por Fmask (bit 1 nube, 2 adyacente, 3 sombra). Mediana 2020–2024.
    """
    def _mask_hls(img):
        f = img.select("Fmask")
        clear = (f.bitwiseAnd(1 << 1).eq(0)        # no nube
                 .And(f.bitwiseAnd(1 << 2).eq(0))  # no adyacente a nube
                 .And(f.bitwiseAnd(1 << 3).eq(0)))  # no sombra de nube
        return img.updateMask(clear).multiply(0.0001)

    def _prep(col_id, nir_band):
        return (ee.ImageCollection(col_id)
                .filterBounds(region).filterDate(start, end)
                .map(_mask_hls)
                .select(["B03", "B04", nir_band], ["green", "red", "nir"]))

    l30 = _prep("NASA/HLS/HLSL30/v002", "B05")
    s30 = _prep("NASA/HLS/HLSS30/v002", "B8A")
    return l30.merge(s30).median().clip(region)


def _s2_index(region, kind, start="2021-01-01", end="2025-01-01"):
    """
    Índice espectral desde Sentinel-2 SR SOLO (10 m) — un único sensor, sin
    costuras. Máscara de nubes con Cloud Score+ (banda cs_cdf, umbral 0.60;
    recomendado por el catálogo GEE, superior a SCL en trópicos). Mediana.

    Fórmulas verificadas (citas originales):
      NDVI  (Rouse 1974):  (B8−B4)/(B8+B4)        = (NIR−Rojo)
      MNDWI (Xu 2006):     (B3−B11)/(B3+B11)       = (Verde−SWIR1)
      NDTI  (Lacaux 2007): (B4−B3)/(B4+B3)         = (Rojo−Verde), turbidez

    Para el AGUA se usa MNDWI (Xu 2006) en vez de NDWI-McFeeters: Satgé et al.
    (2017, Lago Poopó) mostró que NDWI subestima el agua somera y turbia del
    Altiplano boliviano; MNDWI/AWEI/WRI son más exactos. El NDTI (turbidez) se
    enmascara a agua (JRC occurrence>30% ∪ MNDWI>0), pues no tiene sentido físico
    sobre tierra. NOTA: NDTI no está calibrado en Bolivia — es un proxy relativo.
    """
    csp = ee.ImageCollection("GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED")
    s2 = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
          .filterBounds(region).filterDate(start, end)
          .linkCollection(csp, ["cs_cdf"])
          .map(lambda i: i.updateMask(i.select("cs_cdf").gte(0.60)))
          .median().clip(region))
    if kind == "ndvi":
        return s2.normalizedDifference(["B8", "B4"]).rename("NDVI")
    # Índice de agua = MNDWI (verde−SWIR1), óptimo para Bolivia
    mndwi = s2.normalizedDifference(["B3", "B11"]).rename("MNDWI")
    if kind == "ndwi":
        return mndwi
    ndti = s2.normalizedDifference(["B4", "B3"]).rename("NDTI")
    jrc = (ee.Image("JRC/GSW1_4/GlobalSurfaceWater")
           .select("occurrence").unmask(0))
    return ndti.updateMask(jrc.gt(30).Or(mndwi.gt(0.0)))


def _hls_index(region, kind):
    """
    Índice espectral desde el composite HLS armonizado (sin costuras):
      NDVI = (NIR−Rojo)/(NIR+Rojo)
      NDWI = (Verde−NIR)/(Verde+NIR)   (McFeeters)
      NDTI = (Rojo−Verde)/(Rojo+Verde) (Lacaux, turbidez)
    El NDTI se enmascara a AGUA (JRC occurrence > 30% ∪ NDWI > 0), porque la
    turbidez solo tiene sentido físico sobre agua (sobre tierra el índice está
    dominado por suelo/vegetación y se ve "mal").
    """
    hls = _hls_composite(region)
    if kind == "ndvi":
        return hls.normalizedDifference(["nir", "red"]).rename("NDVI")
    ndwi = hls.normalizedDifference(["green", "nir"]).rename("NDWI")
    if kind == "ndwi":
        return ndwi
    ndti = hls.normalizedDifference(["red", "green"]).rename("NDTI")
    jrc = (ee.Image("JRC/GSW1_4/GlobalSurfaceWater")
           .select("occurrence").unmask(0))
    water = jrc.gt(30).Or(ndwi.gt(0.0))
    return ndti.updateMask(water)


def _modis_index(region, kind):
    """
    Índices a ESCALA DE CUENCA GRANDE con MODIS (diseñado para escala
    continental; una mediana S2/Landsat de 10–30 m sobre >100 000 km² excede
    los límites de cómputo de GEE y el mapa caía al sintético).
      - NDVI: producto oficial MOD13Q1 (250 m, ×0.0001)
      - NDWI (McFeeters): MOD09A1 500 m → (b04 verde − b02 NIR)/(suma)
      - NDTI (Lacaux):    MOD09A1 500 m → (b01 rojo − b04 verde)/(suma)
    """
    if kind == "ndvi":
        return (ee.ImageCollection("MODIS/061/MOD13Q1")
                .filterBounds(region)
                .filterDate("2020-01-01", "2024-12-31")
                .select("NDVI")
                .median().multiply(0.0001)
                .rename("NDVI").clip(region))

    col = (ee.ImageCollection("MODIS/061/MOD09A1")
           .filterBounds(region)
           .filterDate("2020-01-01", "2024-12-31"))

    def _idx(img):
        r = img.select("sur_refl_b01").multiply(0.0001)   # rojo
        n = img.select("sur_refl_b02").multiply(0.0001)   # NIR
        g = img.select("sur_refl_b04").multiply(0.0001)   # verde
        if kind == "ndwi":
            v = g.subtract(n).divide(g.add(n))
        else:                                             # ndti
            v = r.subtract(g).divide(r.add(g))
        return v.rename(kind.upper())

    return col.map(_idx).median().clip(region)


def fetch_rivers_hydrosheds(boundary_lonlat, max_order=6, max_feats=400):
    """
    Red de drenaje VECTORIAL (WWF HydroSHEDS Free Flowing Rivers) dentro del
    polígono de cuenca — para mega-cuencas donde la red ráster local (ventana
    MERIT de ~74 km) solo cubre el entorno de la salida. RIV_ORD: 1 = río más
    grande. Devuelve lista de polilíneas [[lon,lat],...] ordenadas por longitud.
    """
    if not _GEE_READY or not boundary_lonlat or len(boundary_lonlat) < 4:
        return None
    try:
        ring = [[float(p[0]), float(p[1])] for p in boundary_lonlat]
        poly = ee.Geometry.Polygon([ring])
        fc = (ee.FeatureCollection("WWF/HydroSHEDS/v1/FreeFlowingRivers")
              .filterBounds(poly)
              .filter(ee.Filter.lte("RIV_ORD", max_order))
              .limit(max_feats))
        info = fc.getInfo()
        lines = []
        for f in info.get("features", []):
            g = f.get("geometry", {})
            if g.get("type") == "LineString":
                lines.append(g["coordinates"])
            elif g.get("type") == "MultiLineString":
                lines.extend(g["coordinates"])
        lines = [ln for ln in lines if len(ln) >= 2]
        lines.sort(key=len, reverse=True)
        return lines or None
    except Exception as e:
        print(f"fetch_rivers_hydrosheds failed: {e}")
        return None


# Paletas y rangos por capa (también usadas para construir la leyenda)
LAYER_META = {
    "dem": {
        "vmin": 1500, "vmax": 5500,
        "palette": ["#0b6623", "#a4d65e", "#f2e394", "#b87333", "#7f5539", "#ffffff"],
        "source": "SRTM v3 / NASA (2000) — 30 m",
        "title": "Modelo Digital de Elevación (DEM) — SRTM 30 m",
        "legend": "Elevación (m s.n.m.)",
    },
    "slope": {
        "vmin": 0, "vmax": 40,
        "palette": ["#ffffb2", "#fecc5c", "#fd8d3c", "#f03b20", "#bd0026"],
        "source": "SRTM v3 + ee.Terrain.slope() / GEE — 30 m",
        "title": "Pendiente del Terreno (grados)",
        "legend": "Pendiente (°)",
    },
    "ndvi": {
        "vmin": -0.2, "vmax": 0.85,
        "palette": ["#a50026", "#d73027", "#fdae61", "#a6d96a", "#1a9850", "#006837"],
        "source": "Sentinel-2 L2A (10 m) — mediana 2021–2024, nubes Cloud Score+",
        "title": "Índice de Vegetación Normalizado (NDVI)",
        "legend": "NDVI (−1 a +1)",
    },
    "ndwi": {
        "vmin": -0.5, "vmax": 0.5,
        "palette": ["#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e"],
        "source": "Sentinel-2 L2A (10 m) — MNDWI (Verde−SWIR1, Xu 2006), mediana 2021–2024, Cloud Score+",
        "title": "Índice de Agua Modificado (MNDWI)",
        "legend": "MNDWI (−1 a +1)",
    },
    "ndti": {
        "vmin": -0.3, "vmax": 0.45,
        "palette": ["#ffffe5", "#fff7bc", "#fee391", "#fec44f", "#fe9929", "#cc4c02"],
        "source": "Sentinel-2 L2A (10 m) — mediana 2021–2024, Cloud Score+ · enmascarado a agua (JRC/NDWI)",
        "title": "Índice de Turbidez Normalizado (NDTI)",
        "legend": "NDTI (−1 a +1)",
    },
    "manning": {
        "vmin": 0.02, "vmax": 0.12,
        "palette": ["#fff7fb", "#d0d1e6", "#67a9cf", "#02818a", "#016450"],
        "source": "ESA WorldCover 2021 — 10 m | Reclasificación Manning's n",
        "title": "Coeficiente de Manning (n) — ESA WorldCover",
        "legend": "n de Manning",
    },
    "risk": {
        "vmin": 0, "vmax": 1,
        "palette": ["#1a9850", "#91cf60", "#d9ef8b", "#fee08b", "#fc8d59", "#d73027"],
        "source": "Multi-fuente GEE (SRTM + JRC) — 30 m",
        "title": "Índice Compuesto de Riesgo Hidrosedimentológico",
        "legend": "Índice de Riesgo (0–1)",
    },
    "jrc": {
        "vmin": 0, "vmax": 100,
        "palette": ["#ffffff", "#deebf7", "#9ecae1", "#3182bd", "#08519c"],
        "source": "JRC Global Surface Water 1984–2021 / Landsat — 30 m",
        "title": "Frecuencia de Inundación — JRC Global Surface Water",
        "legend": "Frecuencia de inundación (%)",
    },
    "watershed": {
        "vmin": 0, "vmax": 1,
        "palette": ["#cce5ff", "#4a90d9", "#003399"],
        "source": "Copernicus DEM GLO-30 (remuestreo bicúbico → 12.5 m) + geoproceso pyflwdir (D8) + Sentinel-2 L2A",
        "title": "Cuenca Hidrográfica y Red de Drenaje",
        "legend": "Elevación (m s.n.m.)",
    },
}


# Por encima de este radio de ventana, los índices se calculan con MODIS:
# una mediana S2/Landsat de 10–30 m sobre cientos de miles de km² excede los
# límites de cómputo de GEE (y el mapa caía en silencio al fondo sintético).
_LARGE_WINDOW_KM = 60.0


def _layer_image(map_type, region, radius_km=15.0):
    """Devuelve la ee.Image de banda única para la capa indicada, eligiendo
    la fuente según la escala de la ventana (S2+Landsat vs MODIS)."""
    large = radius_km > _LARGE_WINDOW_KM

    if map_type == "dem":
        return ee.Image("USGS/SRTMGL1_003").clip(region)

    if map_type == "slope":
        dem = ee.Image("USGS/SRTMGL1_003")
        return ee.Terrain.slope(dem).clip(region)

    if map_type == "ndvi":
        if large:
            return _modis_index(region, "ndvi")
        return _s2_index(region, "ndvi")

    if map_type == "ndwi":
        if large:
            return _modis_index(region, "ndwi")
        return _s2_index(region, "ndwi")

    if map_type == "ndti":
        if large:
            return _modis_index(region, "ndti")
        return _s2_index(region, "ndti")

    if map_type == "manning":
        lc = ee.Image("ESA/WorldCover/v100/2020").select("Map").clip(region)
        from_classes = [10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100]
        to_n_milli = [120, 80, 50, 40, 25, 22, 10, 30, 80, 120, 50]
        return lc.remap(from_classes, to_n_milli).divide(1000)

    if map_type == "risk":
        dem = ee.Image("USGS/SRTMGL1_003")
        slope = ee.Terrain.slope(dem)
        jrc = (ee.Image("JRC/GSW1_4/GlobalSurfaceWater")
               .select("occurrence").unmask(0))
        slope_factor = slope.multiply(-1).add(40).divide(40).clamp(0, 1)
        risk = (jrc.divide(100).multiply(0.6)
                .add(slope_factor.multiply(0.4)))
        return risk.clip(region)

    if map_type == "jrc":
        return (ee.Image("JRC/GSW1_4/GlobalSurfaceWater")
                .select("occurrence").clip(region))

    raise ValueError(f"Capa desconocida: {map_type}")


# ════════════════════════════════════════════════════════════════════════════
# DEM COPERNICUS DE ALTA RESOLUCIÓN (para geoproceso hidrológico)
# ════════════════════════════════════════════════════════════════════════════

def utm_epsg(lat, lon):
    """EPSG del huso UTM WGS84 que contiene (lat, lon)."""
    zone = int((lon + 180) / 6) + 1
    return (32600 + zone) if lat >= 0 else (32700 + zone)


_LAST_DEM_ERROR = None


def last_dem_error():
    return _LAST_DEM_ERROR


def fetch_copernicus_dem(lat, lon, radius_km=15.0, scale_m=12.5):
    """
    Descarga el DEM Copernicus GLO-30 con remuestreo bicúbico (downscaling) a
    `scale_m` metros, reproyectado al huso UTM local, como array float.

    Retorna (dem_array, transform_affine, epsg, extent_lonlat) o None si GEE
    no está disponible o la descarga falla.
      - dem_array  : numpy 2-D float (metros), nodata → np.nan
      - transform  : affine.Affine del raster en UTM
      - epsg       : int del huso UTM (p.ej. 32719)
      - extent     : (lon_min, lon_max, lat_min, lat_max)
    """
    global _LAST_DEM_ERROR
    _LAST_DEM_ERROR = None
    if not _GEE_READY:
        _LAST_DEM_ERROR = "GEE no inicializado"
        return None
    try:
        import requests
        import numpy as np
        from rasterio.io import MemoryFile

        region = _build_region(lat, lon, radius_km)
        epsg = utm_epsg(lat, lon)

        dem = (ee.ImageCollection("COPERNICUS/DEM/GLO30")
               .select("DEM").mosaic()
               .resample("bicubic"))

        url = dem.getDownloadURL({
            "region": region,
            "scale": scale_m,
            "crs": f"EPSG:{epsg}",
            "format": "GEO_TIFF",
        })
        resp = requests.get(url, timeout=180)
        resp.raise_for_status()

        # GEE puede responder 200 con un cuerpo de error (no es un GeoTIFF)
        ctype = resp.headers.get("Content-Type", "")
        if "tif" not in ctype and "octet-stream" not in ctype and "zip" not in ctype:
            snippet = resp.content[:300].decode("utf-8", "replace")
            raise RuntimeError(f"respuesta no-GeoTIFF (Content-Type={ctype}): {snippet}")

        with MemoryFile(resp.content) as mf:
            with mf.open() as ds:
                arr = ds.read(1).astype("float64")
                transform = ds.transform
                nod = ds.nodata
        if nod is not None:
            arr[arr == nod] = np.nan
        arr[arr < -1000] = np.nan

        deg_lat = radius_km / 111.0
        deg_lon = radius_km / (111.0 * math.cos(math.radians(lat)))
        extent = (lon - deg_lon, lon + deg_lon, lat - deg_lat, lat + deg_lat)
        return arr, transform, epsg, extent
    except Exception as e:
        _LAST_DEM_ERROR = f"{type(e).__name__}: {e}"
        print(f"fetch_copernicus_dem failed: {_LAST_DEM_ERROR}")
        return None


def fetch_merit_hydro(lat, lon, radius_km=40.0, max_dim=1600):
    """
    Descarga MERIT Hydro (hidrografía global pre-acondicionada, ~90 m): bandas
    'dir' (dirección de flujo D8 ESRI), 'upa' (área acumulada km²) y 'elv'
    (elevación m). Ideal para zonas planas (Amazonía) y cuencas grandes, donde
    el geoproceso sobre DEM crudo falla.

    Retorna (dir_arr, upa_arr, elv_arr, transform, extent) en EPSG:4326,
    o None si falla.
    """
    global _LAST_DEM_ERROR
    _LAST_DEM_ERROR = None
    if not _GEE_READY:
        _LAST_DEM_ERROR = "GEE no inicializado"
        return None
    try:
        import requests
        import numpy as np
        from rasterio.io import MemoryFile

        # La dirección de flujo D8 NO admite remuestreo (romper la malla nativa
        # invalida la topología). Se trabaja SIEMPRE a la escala nativa de MERIT
        # (~92.77 m) y se recorta el radio para respetar el límite de descarga.
        scale_m = 92.77
        max_radius = (max_dim * scale_m) / 2000.0   # km por lado/2
        radius_km = min(radius_km, max_radius)
        region = _build_region(lat, lon, radius_km)
        # .toFloat(): 'dir' es entero y 'upa'/'elv' float — GEE rechaza GeoTIFF
        # multibanda con tipos mezclados ("bands must have compatible types").
        img = (ee.Image("MERIT/Hydro/v1_0_1")
               .select(["dir", "upa", "elv"]).toFloat())
        url = img.getDownloadURL({
            "region": region,
            "scale": scale_m,
            "crs": "EPSG:4326",
            "format": "GEO_TIFF",
        })
        resp = requests.get(url, timeout=180)
        resp.raise_for_status()
        ctype = resp.headers.get("Content-Type", "")
        if "tif" not in ctype and "octet-stream" not in ctype:
            snippet = resp.content[:300].decode("utf-8", "replace")
            raise RuntimeError(f"respuesta no-GeoTIFF ({ctype}): {snippet}")

        with MemoryFile(resp.content) as mf:
            with mf.open() as ds:
                dir_arr = np.rint(ds.read(1)).astype("int16")
                upa_arr = ds.read(2).astype("float64")
                elv_arr = ds.read(3).astype("float64")
                transform = ds.transform

        deg_lat = radius_km / 111.0
        deg_lon = radius_km / (111.0 * math.cos(math.radians(lat)))
        extent = (lon - deg_lon, lon + deg_lon, lat - deg_lat, lat + deg_lat)
        return dir_arr, upa_arr, elv_arr, transform, extent
    except Exception as e:
        _LAST_DEM_ERROR = f"MERIT: {type(e).__name__}: {e}"
        print(f"fetch_merit_hydro failed: {_LAST_DEM_ERROR}")
        return None


def fetch_s2_rgb(lat, lon, radius_km=15.0, dimensions=1100):
    """Fondo satelital Sentinel-2 color verdadero (B4-B3-B2) como array RGB."""
    if not _GEE_READY:
        return None
    try:
        import requests
        from matplotlib import image as mpimg

        region = _build_region(lat, lon, radius_km)
        s2 = _s2_median(region).select(["B4", "B3", "B2"])
        url = s2.visualize(min=0, max=3000, gamma=1.4).getThumbURL({
            "region": region, "dimensions": dimensions, "format": "png",
        })
        resp = requests.get(url, timeout=120)
        resp.raise_for_status()
        return mpimg.imread(io.BytesIO(resp.content))
    except Exception as e:
        print(f"fetch_s2_rgb failed: {e}")
        return None


def fetch_hydrobasins_upstream(lat, lon, level=7):
    """
    Cuenca de aporte con HydroSHEDS/HydroBASINS (vectorial, pre-delimitada
    global): la herramienta correcta para cuencas GRANDES/internacionales que
    ningún ráster local puede contener (p.ej. Iténez/Mamoré, >300 000 km²).

    Toma la sub-cuenca `level` que contiene el punto y agrega aguas arriba por
    'NEXT_DOWN' (traversal iterativo en GEE). Nivel 7 (~sub-cuencas de cientos
    de km²) mantiene pocas features → la unión de geometrías no expira.

    Retorna dict {'boundary': [[lon,lat],...], 'area_km2': UP_AREA del punto}
    o None. El área proviene del campo UP_AREA de HydroBASINS nivel 12 (área
    de drenaje total aguas arriba, exacta e independiente de la ventana).
    """
    if not _GEE_READY:
        return None
    try:
        point = ee.Geometry.Point([lon, lat])

        # Área de drenaje real en el punto (campo UP_AREA, nivel 12 = preciso)
        up_area = None
        try:
            b12 = ee.FeatureCollection("WWF/HydroSHEDS/v1/Basins/hybas_12")
            seed12 = b12.filterBounds(point).first()
            up_area = ee.Number(seed12.get("UP_AREA")).getInfo()
        except Exception as ae:
            print(f"UP_AREA lookup failed: {ae}")

        basins = ee.FeatureCollection(
            f"WWF/HydroSHEDS/v1/Basins/hybas_{level}")
        seed = basins.filterBounds(point).first()
        if seed is None:
            return None
        seed_id = ee.Number(seed.get("HYBAS_ID"))

        def _step(_, state):
            state = ee.Dictionary(state)
            ids = ee.List(state.get("ids"))
            ups = basins.filter(ee.Filter.inList("NEXT_DOWN", ids))
            new_ids = ups.aggregate_array("HYBAS_ID")
            merged = ids.cat(new_ids).distinct()
            return ee.Dictionary({"ids": merged})

        init = ee.Dictionary({"ids": ee.List([seed_id])})
        result = ee.Dictionary(
            ee.List.sequence(1, 25).iterate(_step, init))
        all_ids = ee.List(result.get("ids"))
        upstream = basins.filter(ee.Filter.inList("HYBAS_ID", all_ids))
        geom = upstream.union(500).geometry().simplify(maxError=500)
        info = geom.getInfo()
        boundary = _extract_polygon_coords(info)
        if not boundary or len(boundary) < 4:
            return None
        return {"boundary": boundary, "area_km2": up_area}
    except Exception as e:
        print(f"fetch_hydrobasins_upstream failed: {e}")
        return None


# n de Manning por clase ESA WorldCover (mismo criterio que la capa 'manning')
_WORLDCOVER_MANNING = {
    10: 0.120, 20: 0.080, 30: 0.050, 40: 0.040, 50: 0.025,
    60: 0.022, 70: 0.010, 80: 0.030, 90: 0.080, 95: 0.120, 100: 0.050,
}
_WORLDCOVER_NAMES = {
    10: "Bosque", 20: "Arbustos", 30: "Pastizal", 40: "Cultivos",
    50: "Urbano", 60: "Suelo desnudo", 70: "Nieve/Hielo", 80: "Agua",
    90: "Humedal", 95: "Mangle", 100: "Musgo/Liquen",
}


def compute_basin_weighted_manning(boundary_lonlat):
    """
    Coeficiente de Manning (n) ponderado por área dentro de la cuenca, a partir
    de la cobertura ESA WorldCover 2021 (10 m). `boundary_lonlat` es la lista
    [[lon,lat], ...] del parteaguas.

    Retorna dict {'n_weighted', 'classes':[{code,name,n,area_pct}], 'source'}
    o None si GEE no está disponible o falla.
    """
    if not _GEE_READY or not boundary_lonlat or len(boundary_lonlat) < 4:
        return None
    try:
        ring = [[float(p[0]), float(p[1])] for p in boundary_lonlat]
        basin = ee.Geometry.Polygon([ring])
        lc = ee.Image("ESA/WorldCover/v100/2020").select("Map").clip(basin)
        hist = lc.reduceRegion(
            reducer=ee.Reducer.frequencyHistogram(),
            geometry=basin, scale=10, maxPixels=1e10, bestEffort=True,
        ).get("Map").getInfo()
        if not hist:
            return None
        total = sum(hist.values())
        if total <= 0:
            return None
        n_weighted = 0.0
        classes = []
        for code_str, count in sorted(hist.items(), key=lambda kv: -kv[1]):
            code = int(float(code_str))
            n_cls = _WORLDCOVER_MANNING.get(code, 0.05)
            frac = count / total
            n_weighted += n_cls * frac
            classes.append({
                "code": code,
                "name": _WORLDCOVER_NAMES.get(code, str(code)),
                "n": round(n_cls, 3),
                "area_pct": round(frac * 100, 1),
            })
        return {
            "n_weighted": round(n_weighted, 4),
            "classes": classes,
            "source": "ESA WorldCover 2021 (10 m) — n ponderado por área de cuenca",
        }
    except Exception as e:
        print(f"compute_basin_weighted_manning failed: {e}")
        return None


# ════════════════════════════════════════════════════════════════════════════
# CUENCA HIDROGRÁFICA Y RED DE DRENAJE (HydroSHEDS — respaldo)
# ════════════════════════════════════════════════════════════════════════════

def _extract_polygon_coords(geojson):
    """Extrae el anillo exterior de un GeoJSON Polygon o MultiPolygon."""
    if not geojson:
        return []
    t = geojson.get("type", "")
    if t == "Polygon":
        return geojson["coordinates"][0]
    if t == "MultiPolygon":
        rings = [poly[0] for poly in geojson["coordinates"]]
        return max(rings, key=len)
    return []


def fetch_watershed_data(lat, lon, radius_km=15.0):
    """
    Delimita la cuenca hidrográfica y extrae la red de drenaje usando GEE.

    Retorna dict con:
      'boundary'    : lista de [lon, lat] del polígono de cuenca
      'stream_mask' : numpy bool array (H×W) — True donde hay cauce
      'rgb'         : numpy RGBA — fondo satelital Sentinel-2 color verdadero
      'extent'      : (lon_min, lon_max, lat_min, lat_max)
      'is_real'     : True
    Retorna None si GEE no está disponible o falla.
    """
    if not _GEE_READY:
        return None
    try:
        import requests
        from matplotlib import image as mpimg
        import numpy as np

        point = ee.Geometry.Point([lon, lat])
        region = _build_region(lat, lon, radius_km)

        deg_lat = radius_km / 111.0
        deg_lon = radius_km / (111.0 * math.cos(math.radians(lat)))
        extent = (lon - deg_lon, lon + deg_lon, lat - deg_lat, lat + deg_lat)

        # 1. Límite de cuenca — HydroSHEDS nivel 12 (sub-cuencas)
        basins = ee.FeatureCollection("WWF/HydroSHEDS/v1/Basins/hybas_12")
        basin_feat = basins.filterBounds(point).first()
        basin_geom = basin_feat.geometry()
        clipped = basin_geom.intersection(region).simplify(maxError=100)
        basin_info = clipped.getInfo()
        boundary = _extract_polygon_coords(basin_info)

        # 2. Red de drenaje — umbral de acumulación de flujo HydroSHEDS 15 seg
        flow_acc = ee.Image("WWF/HydroSHEDS/15ACC").clip(region)
        streams = flow_acc.gt(300).selfMask()
        stream_url = streams.visualize(
            min=0, max=1, palette=["000000", "1155bb"]
        ).getThumbURL({"region": region, "dimensions": 512, "format": "png"})
        resp_s = requests.get(stream_url, timeout=60)
        resp_s.raise_for_status()
        sarr = mpimg.imread(io.BytesIO(resp_s.content))
        if sarr.ndim == 3 and sarr.shape[2] >= 3:
            stream_mask = (sarr[:, :, 2] > 0.25) & (sarr[:, :, 0] < 0.15)
        else:
            stream_mask = sarr[:, :, 0] > 0.5

        # 3. Fondo satelital — Sentinel-2 color verdadero (B4-B3-B2)
        s2 = _s2_median(region).select(["B4", "B3", "B2"])
        rgb_url = s2.visualize(min=0, max=3000, gamma=1.4).getThumbURL({
            "region": region, "dimensions": 700, "format": "png"
        })
        resp_rgb = requests.get(rgb_url, timeout=90)
        resp_rgb.raise_for_status()
        rgb_arr = mpimg.imread(io.BytesIO(resp_rgb.content))

        return {
            "boundary": boundary,
            "stream_mask": stream_mask,
            "rgb": rgb_arr,
            "extent": extent,
            "is_real": True,
        }
    except Exception as e:
        print(f"fetch_watershed_data failed: {e}")
        return None


# Capas con rango dependiente del terreno/región: la paleta se estira a los
# percentiles 2–98 de la REGIÓN (como los ejemplos de cuenca del Code Editor
# de GEE, p.ej. "NDVI from 0.056 to 0.776"). Con rangos fijos, el DEM amazónico
# quedaba en un solo color y los índices desperdiciaban la paleta en valores
# inexistentes en la zona.
_DYNAMIC_STRETCH = {"dem", "slope", "ndvi", "ndwi", "ndti"}


def fetch_gee_thumbnail(map_type, lat, lon, radius_km=15.0, dimensions=1024):
    """
    Obtiene la imagen renderizada real desde GEE (getThumbURL).

    Retorna (arr, vmin, vmax): array numpy RGB(A) float [0..1] y el rango de
    valores realmente usado para la paleta (dinámico para dem/slope, fijo para
    el resto). Devuelve None si GEE no está disponible o falla la petición.
    """
    if not _GEE_READY:
        return None
    try:
        import requests
        from matplotlib import image as mpimg

        region = _build_region(lat, lon, radius_km)
        meta = LAYER_META[map_type]
        img = _layer_image(map_type, region, radius_km=radius_km)

        vmin, vmax = meta["vmin"], meta["vmax"]
        if map_type in _DYNAMIC_STRETCH:
            try:
                stats = img.reduceRegion(
                    reducer=ee.Reducer.percentile([2, 98]),
                    geometry=region, scale=150,
                    maxPixels=1e9, bestEffort=True,
                ).getInfo() or {}
                vals = [v for v in stats.values() if v is not None]
                if len(vals) >= 2:
                    lo, hi = min(vals), max(vals)
                    if hi - lo > 1e-6:
                        vmin, vmax = float(lo), float(hi)
            except Exception as se:
                print(f"dynamic stretch failed for '{map_type}': {se}")

        vis = img.visualize(min=vmin, max=vmax, palette=meta["palette"])
        url = vis.getThumbURL({
            "region": region,
            "dimensions": dimensions,
            "format": "png",
        })
        resp = requests.get(url, timeout=120)
        resp.raise_for_status()
        arr = mpimg.imread(io.BytesIO(resp.content))
        return arr, vmin, vmax
    except Exception as e:
        print(f"GEE thumbnail fetch failed for '{map_type}': {e}")
        return None


# ════════════════════════════════════════════════════════════════════════════
# FUNCIONES PUNTUALES PARA LA CALCULADORA (sin cambios funcionales)
# ════════════════════════════════════════════════════════════════════════════

def get_slope_from_dem(lat, lon, buffer_m=100):
    point = ee.Geometry.Point([lon, lat])
    dem = ee.Image("USGS/SRTMGL1_003")
    slope_img = ee.Terrain.slope(dem)
    stats = slope_img.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=point.buffer(buffer_m),
        scale=30,
    ).getInfo()
    slope_deg = stats.get("slope", 0)
    return math.tan(slope_deg * math.pi / 180.0)


def get_map_url(lat, lon, layer_type="slope", zoom=15):
    point = ee.Geometry.Point([lon, lat])
    region = point.buffer(2000).bounds()
    if layer_type == "slope":
        dem = ee.Image("USGS/SRTMGL1_003")
        img = ee.Terrain.slope(dem)
        vis_params = {"min": 0, "max": 30, "palette": ["blue", "green", "red"]}
    elif layer_type == "landcover":
        img = ee.Image("ESA/WorldCover/v100/2020").select("Map")
        vis_params = {}
    else:
        img = ee.Image("USGS/SRTMGL1_003")
        vis_params = {}
    params = {
        "region": region,
        "dimensions": 512,
        "format": "png",
        "min": vis_params.get("min", 0),
        "max": vis_params.get("max", 30),
        "palette": ",".join(vis_params.get("palette", ["blue", "green", "red"])),
    }
    return img.getThumbURL(params)


def get_landcover_at_point(lat, lon):
    point = ee.Geometry.Point([lon, lat])
    lc = ee.Image("ESA/WorldCover/v100/2020").select("Map")
    val = lc.reduceRegion(
        reducer=ee.Reducer.first(), geometry=point, scale=10
    ).getInfo()
    return val.get("Map", "Unknown")


def get_ndti_turbidity(lat, lon):
    point = ee.Geometry.Point([lon, lat])
    s2 = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
          .filterBounds(point)
          .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 20))
          .sort("system:time_start", False)
          .first())
    if not s2:
        return 0.0
    ndti = s2.normalizedDifference(["B4", "B3"])
    val = ndti.reduceRegion(
        reducer=ee.Reducer.mean(), geometry=point.buffer(50), scale=10
    ).getInfo()
    return val.get("nd", 0.0)


def _measure_halfwidth(water, point, radius_m, scale):
    """Semi-anchura (m) = máx. distancia a la orilla dentro del radio, o None
    si no hay agua en la ventana. `water` es máscara 0/1 (1 = agua)."""
    land = water.Not()
    dist_m = (land.fastDistanceTransform(256).sqrt()
              .multiply(scale)
              .updateMask(water)
              .rename("halfw"))
    stats = dist_m.reduceRegion(
        reducer=ee.Reducer.max(),
        geometry=point.buffer(radius_m),
        scale=scale,
        maxPixels=int(1e8),
        bestEffort=True,
    ).getInfo()
    return stats.get("halfw")


def get_channel_width(lat, lon, search_m=300, occurrence_pct=50, scale=30):
    """
    Estima el ancho del cauce B (m) desde una máscara de agua.

    Método: se calcula la transformada de distancia hacia la orilla más cercana
    (banco) y se toma el MÁXIMO dentro de un radio alrededor del punto (línea
    central = SEMI-anchura); B = 2× esa distancia.

    Robustez: los ríos bolivianos (andinos, trenzados, estacionales) a menudo no
    alcanzan una permanencia alta en el JRC Global Surface Water, así que se
    prueba una CASCADA de umbrales de `occurrence` y radios crecientes, y si nada
    aparece, se usa un respaldo con NDWI de Sentinel-2 (mediana reciente). La
    respuesta indica qué método/umbral/radio funcionó.

    Limita: un píxel es de ~30 m (cauces < ~30 m no se resuelven); cerca de lagos
    puede sobrestimar; no ve bajo el agua (no reemplaza batimetría).

    Devuelve {width_m, half_width_m, source, note} o {error}.
    """
    try:
        point = ee.Geometry.Point([lon, lat])
        gsw = ee.Image("JRC/GSW1_4/GlobalSurfaceWater")
        occ = gsw.select("occurrence").unmask(0)

        # Cascada JRC: (umbral %, radio m). De permanente/estrecho a
        # ocasional/amplio, para captar ríos estacionales o mal geolocalizados.
        for thr, rad in [(occurrence_pct, search_m), (25, search_m),
                         (25, search_m * 2), (10, search_m * 2),
                         (10, search_m * 4)]:
            half = _measure_halfwidth(occ.gte(thr), point, rad, scale)
            if half is not None:
                w = round(2.0 * float(half), 1)
                return {
                    "width_m": w, "half_width_m": round(float(half), 1),
                    "source": f"JRC Global Surface Water 1984–2021 "
                              f"(occurrence ≥ {thr}%), radio {rad} m, {scale} m",
                    "note": "Aproximación multitemporal; verificar con "
                            "levantamiento. No resuelve cauces < ~30 m.",
                }

        # Respaldo NDWI (McFeeters) con Sentinel-2 (mediana 2020–2024).
        rad = search_m * 4
        s2 = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
              .filterBounds(point.buffer(rad))
              .filterDate("2020-01-01", "2024-12-31")
              .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 30))
              .median())
        ndwi = s2.normalizedDifference(["B3", "B8"])   # verde, NIR
        half = _measure_halfwidth(ndwi.gt(0.0), point, rad, 10)
        if half is not None:
            w = round(2.0 * float(half), 1)
            return {
                "width_m": w, "half_width_m": round(float(half), 1),
                "source": f"NDWI Sentinel-2 (mediana 2020–2024), radio {rad} m, 10 m",
                "note": "Respaldo NDWI (no había agua en JRC). Verificar con "
                        "levantamiento; cerca de lagos puede sobrestimar.",
            }

        return {"error": "No se encontró agua cerca del punto ni en JRC Global "
                         "Surface Water ni en NDWI Sentinel-2. Verifique que el "
                         "punto caiga sobre el cauce (haga clic sobre el río)."}
    except Exception as e:
        return {"error": f"Error GEE al estimar el ancho: {e}"}


def compute_rusle_point(lat, lon, r_formula="hurni", buffer_m=150,
                        years=("2014-01-01", "2024-01-01")):
    """
    Estima la pérdida de suelo RUSLE (A = R·K·LS·C·P) en el punto, obteniendo
    los cinco factores desde GEE y combinándolos con models.rusle.

    Fuentes:
        R  → CHIRPS (precip. media anual) + fórmula empírica (Hurni por defecto)
        K  → OpenLandMap clase de textura USDA
        LS → Copernicus/SRTM (pendiente) → factor S de McCool
        C  → NDVI Sentinel-2 (van der Knijff)
        P  → ESA WorldCover (cultivo) × pendiente

    Devuelve el desglose de factores + A, o {error}. Verificable en el Space.
    """
    from models.rusle import compute_rusle
    try:
        point = ee.Geometry.Point([lon, lat])
        region = point.buffer(buffer_m)
        start, end = years
        n_years = max(1, (int(end[:4]) - int(start[:4])))

        # R — precipitación media anual (mm) desde CHIRPS
        chirps = (ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY")
                  .filterDate(start, end).select("precipitation").sum())
        pa_total = chirps.reduceRegion(
            reducer=ee.Reducer.mean(), geometry=region, scale=5000,
            maxPixels=int(1e8), bestEffort=True).getInfo().get("precipitation")
        precip_anual = (pa_total / n_years) if pa_total else 0.0

        # K — clase de textura USDA (OpenLandMap, 0 cm)
        textura = (ee.Image("OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02")
                   .select("b0"))
        clase_tex = textura.reduceRegion(
            reducer=ee.Reducer.mode(), geometry=region, scale=250,
            maxPixels=int(1e8), bestEffort=True).getInfo().get("b0") or 6

        # LS — pendiente (%) desde SRTM
        dem = ee.Image("USGS/SRTMGL1_003")
        slope_deg = ee.Terrain.slope(dem).reduceRegion(
            reducer=ee.Reducer.mean(), geometry=region, scale=30,
            maxPixels=int(1e8), bestEffort=True).getInfo().get("slope") or 0.0
        slope_pct = math.tan(math.radians(slope_deg)) * 100.0

        # C — NDVI Sentinel-2 (mediana anual, nubes < 20 %)
        s2 = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
              .filterBounds(region).filterDate(start, end)
              .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 20)).median())
        ndvi = s2.normalizedDifference(["B8", "B4"]).reduceRegion(
            reducer=ee.Reducer.mean(), geometry=region, scale=10,
            maxPixels=int(1e8), bestEffort=True).getInfo().get("nd") or 0.0

        # P — ¿cultivo? (ESA WorldCover clase 40 = cropland)
        lc = (ee.Image("ESA/WorldCover/v100/2020").select("Map")
              .reduceRegion(reducer=ee.Reducer.mode(), geometry=region,
                            scale=10, maxPixels=int(1e8), bestEffort=True)
              .getInfo().get("Map"))
        es_cultivo = (lc == 40)

        res = compute_rusle(precip_anual, clase_tex, slope_pct, ndvi,
                            es_cultivo=es_cultivo, r_formula=r_formula)
        res.update({
            "precip_anual_mm": round(precip_anual, 1),
            "clase_textura": int(clase_tex),
            "slope_pct": round(slope_pct, 2),
            "ndvi": round(ndvi, 3),
            "es_cultivo": bool(es_cultivo),
            "source": "CHIRPS + OpenLandMap + SRTM + Sentinel-2 + ESA WorldCover",
            "note": "Estimación puntual (buffer local). RUSLE da erosión "
                    "potencial en ladera, no el aporte real al cauce (usar SDR).",
        })
        return res
    except Exception as e:
        return {"error": f"Error GEE al estimar RUSLE: {e}"}
