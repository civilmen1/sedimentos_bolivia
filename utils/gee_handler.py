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


def initialize_gee():
    """Inicializa GEE con cuenta de servicio o credenciales por defecto."""
    global _GEE_READY
    try:
        sa_json = (os.environ.get("EE_SERVICE_ACCOUNT_JSON")
                   or os.environ.get("GEE_SERVICE_ACCOUNT_JSON"))
        if sa_json:
            info = json.loads(sa_json)
            email = info.get("client_email")
            project = info.get("project_id")
            creds = ee.ServiceAccountCredentials(email, key_data=sa_json)
            if project:
                ee.Initialize(creds, project=project)
            else:
                ee.Initialize(creds)
            print(f"GEE inicializado con cuenta de servicio: {email}")
            _GEE_READY = True
            return True

        key_file = os.environ.get("EE_SERVICE_ACCOUNT_FILE")
        sa_email = os.environ.get("EE_SERVICE_ACCOUNT_EMAIL")
        if key_file and sa_email:
            creds = ee.ServiceAccountCredentials(sa_email, key_file=key_file)
            ee.Initialize(creds)
            print(f"GEE inicializado con archivo de cuenta de servicio: {sa_email}")
            _GEE_READY = True
            return True

        project = os.environ.get("GEE_PROJECT")
        if project:
            ee.Initialize(project=project)
        else:
            ee.Initialize()
        print("GEE inicializado con credenciales por defecto.")
        _GEE_READY = True
        return True
    except Exception as e:
        print(f"GEE Initialization failed: {e}")
        _GEE_READY = False
        return False


def gee_ready():
    return _GEE_READY


# ════════════════════════════════════════════════════════════════════════════
# DEFINICIÓN DE CAPAS TEMÁTICAS
# ════════════════════════════════════════════════════════════════════════════

def _build_region(lat, lon, radius_km):
    deg_lat = radius_km / 111.0
    deg_lon = radius_km / (111.0 * math.cos(math.radians(lat)))
    return ee.Geometry.Rectangle(
        [lon - deg_lon, lat - deg_lat, lon + deg_lon, lat + deg_lat]
    )


def _s2_median(region):
    """Mediana Sentinel-2 SR libre de nubes 2020–2023."""
    return (
        ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
        .filterBounds(region)
        .filterDate("2020-01-01", "2023-12-31")
        .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 25))
        .median()
    )


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
        "source": "Sentinel-2 L2A / ESA — 10 m | Mediana 2020–2023",
        "title": "Índice de Vegetación Normalizado (NDVI)",
        "legend": "NDVI (−1 a +1)",
    },
    "ndwi": {
        "vmin": -0.5, "vmax": 0.5,
        "palette": ["#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e"],
        "source": "Sentinel-2 L2A / ESA — 10 m | Mediana 2020–2023",
        "title": "Índice de Agua Normalizado (NDWI)",
        "legend": "NDWI (−1 a +1)",
    },
    "ndti": {
        "vmin": -0.3, "vmax": 0.45,
        "palette": ["#ffffe5", "#fff7bc", "#fee391", "#fec44f", "#fe9929", "#cc4c02"],
        "source": "Sentinel-2 L2A / ESA — 10 m | Mediana 2020–2023",
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
        "source": "Copernicus DEM GLO-30 (remuestreo bicúbico → 12.5 m) + geoproceso pysheds + Sentinel-2 L2A",
        "title": "Cuenca Hidrográfica y Red de Drenaje",
        "legend": "Elevación (m s.n.m.)",
    },
}


def _layer_image(map_type, region):
    """Devuelve la ee.Image de banda única para la capa indicada."""
    if map_type == "dem":
        return ee.Image("USGS/SRTMGL1_003").clip(region)

    if map_type == "slope":
        dem = ee.Image("USGS/SRTMGL1_003")
        return ee.Terrain.slope(dem).clip(region)

    if map_type == "ndvi":
        return _s2_median(region).normalizedDifference(["B8", "B4"]).clip(region)

    if map_type == "ndwi":
        return _s2_median(region).normalizedDifference(["B3", "B8"]).clip(region)

    if map_type == "ndti":
        return _s2_median(region).normalizedDifference(["B4", "B3"]).clip(region)

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
    if not _GEE_READY:
        return None
    try:
        import requests
        import numpy as np
        import rasterio
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
        print(f"fetch_copernicus_dem failed: {e}")
        return None


def fetch_s2_rgb(lat, lon, radius_km=15.0, dimensions=700):
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
        resp = requests.get(url, timeout=90)
        resp.raise_for_status()
        return mpimg.imread(io.BytesIO(resp.content))
    except Exception as e:
        print(f"fetch_s2_rgb failed: {e}")
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


def fetch_gee_thumbnail(map_type, lat, lon, radius_km=15.0, dimensions=700):
    """
    Obtiene la imagen renderizada real desde GEE (getThumbURL) y la devuelve
    como array numpy RGB(A) float [0..1]. Devuelve None si GEE no está
    disponible o falla la petición.
    """
    if not _GEE_READY:
        return None
    try:
        import requests
        from matplotlib import image as mpimg

        region = _build_region(lat, lon, radius_km)
        meta = LAYER_META[map_type]
        img = _layer_image(map_type, region)
        vis = img.visualize(min=meta["vmin"], max=meta["vmax"],
                            palette=meta["palette"])
        url = vis.getThumbURL({
            "region": region,
            "dimensions": dimensions,
            "format": "png",
        })
        resp = requests.get(url, timeout=90)
        resp.raise_for_status()
        arr = mpimg.imread(io.BytesIO(resp.content))
        return arr
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
