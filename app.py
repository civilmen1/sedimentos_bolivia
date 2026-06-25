from flask import Flask, render_template, request, jsonify
import numpy as np
import math
from datetime import datetime
from models.sediment import (
    calculate_water_density,
    calculate_kinematic_viscosity,
    calculate_specific_gravity,
    calculate_fall_velocity,
    calculate_dimensionless_particle_parameter,
    calculate_critical_shear_stress_shields,
    meyer_peter_muller,
    engelund_hansen,
    van_rijn_bedload
)
from utils.gee_handler import (
    initialize_gee,
    get_slope_from_dem,
    get_map_url,
    get_landcover_at_point,
    get_ndti_turbidity
)

app = Flask(__name__)
GEE_AVAILABLE = initialize_gee()
G = 9.807

WORLDCOVER = {
    10: ("Bosque", 0.12),
    20: ("Arbustos", 0.08),
    30: ("Pastizal", 0.05),
    40: ("Cultivos", 0.04),
    50: ("Urbano", 0.025),
    60: ("Suelo desnudo", 0.022),
    70: ("Nieve / Hielo", 0.010),
    80: ("Agua abierta", 0.030),
    90: ("Humedal", 0.08),
    95: ("Mangle", 0.12),
    100: ("Musgo / Líquen", 0.05),
}


def classify_particle(d50_mm):
    if d50_mm < 0.004:
        return ("Arcilla", "#c0392b")
    elif d50_mm < 0.0625:
        return ("Limo", "#e67e22")
    elif d50_mm < 0.25:
        return ("Arena fina", "#f1c40f")
    elif d50_mm < 0.5:
        return ("Arena media", "#f39c12")
    elif d50_mm < 2.0:
        return ("Arena gruesa", "#d35400")
    elif d50_mm < 16:
        return ("Grava fina", "#7f8c8d")
    elif d50_mm < 64:
        return ("Grava gruesa", "#566573")
    else:
        return ("Canto rodado", "#2c3e50")


def rouse_mode(z):
    if z > 7.5:
        return "Arrastre de fondo exclusivo"
    elif z > 2.5:
        return "Fondo dominante"
    elif z > 1.2:
        return "Transporte mixto"
    elif z > 0.8:
        return "Suspensión dominante"
    else:
        return "Suspensión / Washload"


def rosgen_classify(slope, d50_mm, froude):
    if slope > 0.10:
        return ("A+", "Muy empinado (S > 10%), confinado, materiales gruesos, cascadas dominantes")
    elif slope > 0.04:
        return ("A", "Empinado (S = 4–10%), ligeramente confinado, rápidos y pozas")
    elif slope > 0.02:
        return ("B", "Moderadamente empinado (S = 2–4%), pocas barras laterales")
    elif slope > 0.005:
        if d50_mm >= 2.0:
            return ("B", "Moderado (S = 0.5–2%), grava/canto, sinuosidad baja–media")
        else:
            return ("C", "Bajo gradiente (S = 0.5–2%), sinuoso, arena/grava, planicie de inundación")
    elif slope > 0.001:
        if froude < 0.3:
            return ("E", "Muy bajo gradiente (S < 0.5%), alta sinuosidad, canal estable")
        else:
            return ("C", "Bajo gradiente, meandriforme, amplia planicie de inundación")
    else:
        return ("D", "Muy bajo gradiente (S < 0.1%), potencial trenzamiento (braided)")


def generate_gee_code(lat, lon, d50, d90):
    template = '''"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   ANÁLISIS MORFOLÓGICO FLUVIAL — GOOGLE EARTH ENGINE PYTHON API             ║
║   Sedimentos Bolivia | Punto: LAT=<<LAT>>, LON=<<LON>>                      ║
║   d₅₀=<<D50>> mm | d₉₀=<<D90>> mm                                          ║
╠══════════════════════════════════════════════════════════════════════════════╣
║   Ecosistema de herramientas integradas:                                     ║
║   ├─ GEE Python API  → procesamiento y análisis satelital                   ║
║   ├─ QGIS / ArcGIS   → cartografía y edición vectorial                      ║
║   ├─ HEC-RAS / HMS   → modelación hidráulica e hidrológica                  ║
║   ├─ iRIC            → morfodinámica fluvial 2D                              ║
║   ├─ SWAT+           → modelación de cuencas a largo plazo                  ║
║   ├─ xarray / pandas → análisis de series temporales                        ║
║   └─ TensorFlow/RF   → clasificación con machine learning                   ║
╚══════════════════════════════════════════════════════════════════════════════╝

ENTORNO RECOMENDADO: Google Colab
  !pip install earthengine-api geemap pandas xarray matplotlib numpy
  import ee; ee.Authenticate(auth_mode='notebook'); ee.Initialize()
"""

import ee
import math

# ════════════════════════════════════════════════════════════════════════════
# PARTE I — CONFIGURACIÓN GENERAL
# ════════════════════════════════════════════════════════════════════════════
LAT    = <<LAT>>     # Latitud del punto de muestreo
LON    = <<LON>>     # Longitud del punto de muestreo
D50_MM = <<D50>>     # Diámetro mediano d₅₀ (mm)
D90_MM = <<D90>>     # Diámetro percentil 90 d₉₀ (mm)

BUFFER_CAUCE  = 2000    # Radio análisis cauce (m)
BUFFER_CUENCA = 15000   # Radio análisis cuenca (m)
EXPORT_FOLDER = "GEE_Sedimentos_Bolivia"
CRS_WGS84     = "EPSG:4326"
SCALE_10      = 10
SCALE_30      = 30
SCALE_90      = 90

ee.Authenticate()
ee.Initialize(project="YOUR_GEE_PROJECT_ID")   # ← Reemplazar con tu proyecto GEE

punto         = ee.Geometry.Point([LON, LAT])
buffer_cauce  = punto.buffer(BUFFER_CAUCE).bounds()
buffer_cuenca = punto.buffer(BUFFER_CUENCA).bounds()

print("=" * 70)
print(f"  ANÁLISIS MORFOLÓGICO FLUVIAL — GEE PYTHON API")
print(f"  LAT={LAT}°  LON={LON}°  d₅₀={D50_MM} mm  d₉₀={D90_MM} mm")
print(f"  Exportación: Google Drive / {EXPORT_FOLDER}")
print("=" * 70)

# ════════════════════════════════════════════════════════════════════════════
# PARTE II — DEM Y MORFOLOGÍA DEL TERRENO
# → QGIS/ArcGIS, HEC-RAS (RAS Mapper), iRIC
# ════════════════════════════════════════════════════════════════════════════
dem       = ee.Image("USGS/SRTMGL1_003").clip(buffer_cuenca)
hillshade = ee.Terrain.hillshade(dem)
pendiente = ee.Terrain.slope(dem)       # grados
aspecto   = ee.Terrain.aspect(dem)
terrain   = ee.Terrain.products(dem)

slope_rad = pendiente.multiply(math.pi / 180)
slope_mm  = slope_rad.tan().rename("slope_m_m")
slope_pct = slope_mm.multiply(100).rename("slope_pct")

# TPI: Topographic Position Index (identifica lechos en depresiones)
dem_focal = dem.focal_mean(radius=300, kernelType="circle", units="meters")
tpi = dem.subtract(dem_focal).rename("TPI")

vis_dem  = dict(min=0, max=5000, palette=["006633","E5FFCC","662A00","D8D8D8","F5F5F5"])
dem_vis  = dem.visualize(**vis_dem)
hill_vis = hillshade.visualize(min=0, max=255, gamma=1.3)
dem_blend = dem_vis.blend(hill_vis.updateMask(ee.Image(0.4)))

# Exportar DEM crudo para HEC-RAS (RAS Mapper) e iRIC
ee.batch.Export.image.toDrive(
    image=dem.rename("elevation"),
    description="Morf01_DEM_HEC-RAS_iRIC",
    folder=EXPORT_FOLDER, region=buffer_cauce,
    scale=SCALE_30, maxPixels=1e13, crs=CRS_WGS84
).start()

# Exportar hillshade + DEM blend para QGIS/ArcGIS
ee.batch.Export.image.toDrive(
    image=dem_blend,
    description="Morf02_DEM_Hillshade_QGIS",
    folder=EXPORT_FOLDER, region=buffer_cuenca,
    scale=SCALE_30, maxPixels=1e13, crs=CRS_WGS84
).start()

# Exportar stack topográfico completo (slope/aspect/TPI) para iRIC
topo_stack = (dem.rename("dem")
              .addBands(slope_pct)
              .addBands(aspecto.rename("aspecto"))
              .addBands(tpi))
ee.batch.Export.image.toDrive(
    image=topo_stack,
    description="Morf03_Topografia_Stack_iRIC",
    folder=EXPORT_FOLDER, region=buffer_cauce,
    scale=SCALE_30, maxPixels=1e13, crs=CRS_WGS84
).start()
print("✓ Parte II: DEM + derivados morfológicos → HEC-RAS, iRIC, QGIS")

# ════════════════════════════════════════════════════════════════════════════
# PARTE III — RED HIDROGRÁFICA Y CUENCA (HydroSHEDS)
# → SWAT+ (delimitación), HEC-HMS (parámetros de cuenca)
# ════════════════════════════════════════════════════════════════════════════
flow_acc = ee.Image("WWF/HydroSHEDS/15ACC").clip(buffer_cuenca)

cauce_menor = flow_acc.gte(500).selfMask().rename("cauce_menor")
cauce_mayor = flow_acc.gte(2000).selfMask().rename("cauce_mayor")
cauce_ppal  = flow_acc.gte(5000).selfMask().rename("cauce_principal")

basins = ee.FeatureCollection("WWF/HydroSHEDS/v1/Basins/hybas_sa_15s_v1c")
cuenca_punto = basins.filterBounds(punto).first()
area_cuenca_km2 = cuenca_punto.geometry().area().divide(1e6)

# Exportar cuenca como parámetros CSV para HEC-HMS
params_cuenca = ee.Feature(None, {
    "area_km2":    area_cuenca_km2,
    "lat_outlet":  LAT,
    "lon_outlet":  LON,
    "d50_mm":      D50_MM,
    "d90_mm":      D90_MM,
})
ee.batch.Export.table.toDrive(
    collection=ee.FeatureCollection([params_cuenca]),
    description="Tabla1_Cuenca_HEC-HMS",
    folder=EXPORT_FOLDER, fileFormat="CSV"
).start()

ee.batch.Export.image.toDrive(
    image=cauce_ppal.addBands(cauce_mayor).addBands(cauce_menor),
    description="Morf04_Red_Hidrografica_SWAT",
    folder=EXPORT_FOLDER, region=buffer_cuenca,
    scale=SCALE_30, maxPixels=1e13, crs=CRS_WGS84
).start()
print("✓ Parte III: Red hidrográfica y cuenca → SWAT+, HEC-HMS")

# ════════════════════════════════════════════════════════════════════════════
# PARTE IV — ANCHO DEL CAUCE (MNDWI Sentinel-2)
# → HEC-RAS (geometría del cauce), iRIC (condición de borde)
# ════════════════════════════════════════════════════════════════════════════
s2_ancho = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
            .filterBounds(punto)
            .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 5))
            .filterDate("2022-06-01", "2024-09-30")
            .median()
            .clip(buffer_cauce))

# MNDWI = (Green - SWIR) / (Green + SWIR) — mejor que NDWI para ancho del cauce
mndwi   = s2_ancho.normalizedDifference(["B3", "B11"]).rename("MNDWI")
canal   = mndwi.gt(0.10).rename("canal_agua")   # máscara binaria del cauce

ee.batch.Export.image.toDrive(
    image=canal,
    description="Morf05_Canal_Mascara_HEC-RAS",
    folder=EXPORT_FOLDER, region=buffer_cauce,
    scale=SCALE_10, maxPixels=1e13, crs=CRS_WGS84
).start()
print("✓ Parte IV: Máscara MNDWI del cauce → HEC-RAS")

# ════════════════════════════════════════════════════════════════════════════
# PARTE V — COBERTURA Y USO DEL SUELO + MANNING'S n + CN
# → HEC-RAS (rugosidad), HEC-HMS (número de curva), SWAT+ (LULC)
# ════════════════════════════════════════════════════════════════════════════
worldcover = ee.Image("ESA/WorldCover/v200/2021").select("Map").clip(buffer_cuenca)

# Tablas de reclasificación WorldCover → Manning's n y CN (AMC-II Tipo B)
clases   = [10,    20,    30,    40,    50,    60,    70,    80,    90,    95,    100 ]
manning  = [0.120, 0.080, 0.050, 0.040, 0.025, 0.022, 0.010, 0.030, 0.080, 0.120, 0.050]
cn_vals  = [60,    65,    74,    81,    91,    91,    30,    100,   85,    78,    68  ]

manning_raster = worldcover.remap(clases, manning).rename("manning_n")
cn_raster      = worldcover.remap(clases, cn_vals ).rename("CN_AMCII")

vis_wc = dict(
    min=10, max=100,
    palette=["006400","FFBB22","FFFF4C","F096FF","FA0000",
             "B4B4B4","F0F0F0","0064C8","0096A0","00CF75","FAE6A0"]
)

ee.batch.Export.image.toDrive(
    image=manning_raster,
    description="Morf06_Manning_n_HEC-RAS",
    folder=EXPORT_FOLDER, region=buffer_cauce,
    scale=SCALE_10, maxPixels=1e13, crs=CRS_WGS84
).start()

ee.batch.Export.image.toDrive(
    image=cn_raster,
    description="Morf07_CN_HEC-HMS_SWAT",
    folder=EXPORT_FOLDER, region=buffer_cuenca,
    scale=SCALE_10, maxPixels=1e13, crs=CRS_WGS84
).start()

ee.batch.Export.image.toDrive(
    image=worldcover.visualize(**vis_wc),
    description="Morf08_WorldCover_QGIS_SWAT",
    folder=EXPORT_FOLDER, region=buffer_cuenca,
    scale=SCALE_10, maxPixels=1e13, crs=CRS_WGS84
).start()
print("✓ Parte V: Manning's n, CN, WorldCover → HEC-RAS, HEC-HMS, SWAT+")

# ════════════════════════════════════════════════════════════════════════════
# PARTE VI — ÍNDICES ESPECTRALES SENTINEL-2
# (NDVI, NDWI, MNDWI, NDTI, EVI, BSI, TSS estimado)
# → QGIS/ArcGIS, análisis multivariado, xarray
# ════════════════════════════════════════════════════════════════════════════
s2 = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
      .filterBounds(punto)
      .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 10))
      .filterDate("2022-01-01", "2024-12-31")
      .median()
      .clip(buffer_cuenca))

ndvi  = s2.normalizedDifference(["B8",  "B4" ]).rename("NDVI")
ndwi  = s2.normalizedDifference(["B3",  "B8" ]).rename("NDWI")
ndti  = s2.normalizedDifference(["B4",  "B3" ]).rename("NDTI")

# EVI = 2.5 * (NIR - Red) / (NIR + 6*Red - 7.5*Blue + 1)
evi = s2.expression(
    "2.5 * ((NIR - Red) / (NIR + 6.0 * Red - 7.5 * Blue + 1.0))",
    {"NIR": s2.select("B8"), "Red": s2.select("B4"), "Blue": s2.select("B2")}
).rename("EVI")

# BSI = Bare Soil Index
bsi = s2.expression(
    "((Red + SWIR) - (NIR + Blue)) / ((Red + SWIR) + (NIR + Blue))",
    {"Red": s2.select("B4"), "SWIR": s2.select("B11"),
     "NIR": s2.select("B8"), "Blue": s2.select("B2")}
).rename("BSI")

# TSS (mg/L) estimado: modelo empírico Chen et al. (2015)
tss = ndti.multiply(3.554).exp().multiply(3.958).rename("TSS_mg_L")

indices_stack = (ndvi.addBands(ndwi).addBands(mndwi)
                 .addBands(ndti).addBands(evi).addBands(bsi).addBands(tss))

ee.batch.Export.image.toDrive(
    image=indices_stack,
    description="Morf09_Indices_S2_xarray",
    folder=EXPORT_FOLDER, region=buffer_cuenca,
    scale=SCALE_10, maxPixels=1e13, crs=CRS_WGS84
).start()
print("✓ Parte VI: NDVI/NDWI/MNDWI/NDTI/EVI/BSI/TSS → QGIS, xarray")

# ════════════════════════════════════════════════════════════════════════════
# PARTE VII — DINÁMICA HIDROLÓGICA: JRC GLOBAL SURFACE WATER
# → HEC-RAS (validación planimetría), SWAT+
# ════════════════════════════════════════════════════════════════════════════
jrc         = ee.Image("JRC/GSW1_4/GlobalSurfaceWater").clip(buffer_cuenca)
ocurrencia  = jrc.select("occurrence")    # % tiempo con agua 1984-2021
recurrencia = jrc.select("recurrence")    # % años húmedos con agua

zona_perm    = ocurrencia.gte(90).rename("zona_permanente")
zona_frec    = ocurrencia.gte(50).And(ocurrencia.lt(90)).rename("zona_frecuente")
zona_estac   = ocurrencia.gte(10).And(ocurrencia.lt(50)).rename("zona_estacional")
zona_excep   = ocurrencia.gt(0).And(ocurrencia.lt(10)).rename("zona_excepcional")

zona_inund = (zona_perm.multiply(4).add(zona_frec.multiply(3))
              .add(zona_estac.multiply(2)).add(zona_excep)
              .rename("zona_inundacion"))

ee.batch.Export.image.toDrive(
    image=ocurrencia.addBands(recurrencia).addBands(zona_inund),
    description="Morf10_JRC_Inundacion_HEC-RAS",
    folder=EXPORT_FOLDER, region=buffer_cuenca,
    scale=SCALE_30, maxPixels=1e13, crs=CRS_WGS84
).start()
print("✓ Parte VII: Zonas JRC GSW (ocurrencia/recurrencia) → HEC-RAS")

# ════════════════════════════════════════════════════════════════════════════
# PARTE VIII — MIGRACIÓN LATERAL MULTITEMPORAL (LANDSAT 8/9)
# Composición RGB: R=2024, G=2015, B=2010
# → QGIS/ArcGIS, xarray para análisis de tendencia
# ════════════════════════════════════════════════════════════════════════════
def get_landsat_sr(year, aoi):
    return (ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
            .merge(ee.ImageCollection("LANDSAT/LC08/C02/T1_L2"))
            .filterBounds(aoi)
            .filterDate(str(year) + "-06-01", str(year) + "-09-30")
            .filter(ee.Filter.lt("CLOUD_COVER", 20))
            .map(lambda img: img.multiply(0.0000275).add(-0.2).clip(aoi))
            .median()
            .select(["SR_B4","SR_B3","SR_B2","SR_B6"], ["R","G","B","SWIR"]))

def get_mndwi_landsat(year, aoi):
    col = get_landsat_sr(year, aoi)
    return (col.normalizedDifference(["G","SWIR"]).gt(0.10)
            .rename("agua_" + str(year)))

anos_ls = [2010, 2015, 2020, 2024]
sc_2010, sc_2015, sc_2020, sc_2024 = [get_landsat_sr(a, buffer_cauce) for a in anos_ls]
agua_ls = [get_mndwi_landsat(a, buffer_cauce) for a in anos_ls]

agua_stack  = agua_ls[0].addBands(agua_ls[1]).addBands(agua_ls[2]).addBands(agua_ls[3])
frec_agua   = agua_stack.reduce(ee.Reducer.sum()).rename("frec_agua")
rgb_cambio  = ee.Image.cat([sc_2024.select("R"), sc_2015.select("G"), sc_2010.select("B")])

ee.batch.Export.image.toDrive(
    image=rgb_cambio.rename(["R_2024","G_2015","B_2010"]),
    description="Morf11_Migracion_Lateral_RGB",
    folder=EXPORT_FOLDER, region=buffer_cauce,
    scale=SCALE_30, maxPixels=1e13, crs=CRS_WGS84
).start()

ee.batch.Export.image.toDrive(
    image=agua_stack.addBands(frec_agua),
    description="Morf12_Agua_Binario_xarray",
    folder=EXPORT_FOLDER, region=buffer_cauce,
    scale=SCALE_30, maxPixels=1e13, crs=CRS_WGS84
).start()
print("✓ Parte VIII: Migración lateral Landsat 2010-2024 → QGIS, xarray")

# ════════════════════════════════════════════════════════════════════════════
# PARTE IX — PRECIPITACIÓN CHIRPS: SERIE MENSUAL 2010-2024
# → CSV para calibración HEC-HMS y análisis con xarray/pandas
# ════════════════════════════════════════════════════════════════════════════
chirps = (ee.ImageCollection("UCSB-CHG/CHIRPS/PENTAD")
          .filterBounds(punto)
          .filterDate("2010-01-01", "2024-12-31")
          .select("precipitation"))

def get_precip_mes(year, month):
    start = ee.Date.fromYMD(year, month, 1)
    pcp   = chirps.filterDate(start, start.advance(1, "month")).sum()
    val   = pcp.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=punto.buffer(5000),
        scale=5566, maxPixels=1e6
    )
    return ee.Feature(None, {"year": year, "month": month,
                              "precip_mm": val.get("precipitation")})

serie_precip = ee.FeatureCollection([
    get_precip_mes(y, m) for y in range(2010, 2025) for m in range(1, 13)
])

ee.batch.Export.table.toDrive(
    collection=serie_precip,
    description="Tabla2_Precipitacion_CHIRPS_HEC-HMS",
    folder=EXPORT_FOLDER, fileFormat="CSV"
).start()

chirps_anual = chirps.filterDate("2015-01-01","2024-12-31").sum().divide(10)
ee.batch.Export.image.toDrive(
    image=chirps_anual.rename("precip_mm_anual"),
    description="Morf13_Precipitacion_CHIRPS",
    folder=EXPORT_FOLDER, region=buffer_cuenca,
    scale=5566, maxPixels=1e13, crs=CRS_WGS84
).start()
print("✓ Parte IX: Precipitación CHIRPS mensual 2010-2024 → HEC-HMS, xarray")

# ════════════════════════════════════════════════════════════════════════════
# PARTE X — SERIE TEMPORAL NDTI: MONITOREO DE TURBIDEZ
# → CSV para análisis con xarray, pandas, matplotlib (Google Colab)
# ════════════════════════════════════════════════════════════════════════════
s2_col = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
          .filterBounds(punto)
          .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 20))
          .filterDate("2019-01-01", "2024-12-31"))

def extraer_ndti(img):
    val = img.normalizedDifference(["B4","B3"]).reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=punto.buffer(100),
        scale=SCALE_10, maxPixels=1000
    )
    return ee.Feature(None, {
        "fecha":   img.date().format("YYYY-MM-dd"),
        "ndti":    val.get("nd"),
        "lat":     LAT,
        "lon":     LON,
        "d50_mm":  D50_MM,
    })

serie_ndti = s2_col.map(extraer_ndti)
ee.batch.Export.table.toDrive(
    collection=serie_ndti,
    description="Tabla3_Serie_NDTI_2019_2024",
    folder=EXPORT_FOLDER, fileFormat="CSV"
).start()

# En Google Colab, analizar con:
#   import pandas as pd, matplotlib.pyplot as plt
#   df = pd.read_csv("Tabla3_Serie_NDTI_2019_2024.csv")
#   df["fecha"] = pd.to_datetime(df["fecha"])
#   df.set_index("fecha").ndti.resample("ME").mean().plot()
print("✓ Parte X: Serie temporal NDTI 2019-2024 → pandas/xarray/Colab")

# ════════════════════════════════════════════════════════════════════════════
# PARTE XI — FEATURES PARA MACHINE LEARNING (RF / TensorFlow)
# → Identificación de zonas fuente de sedimento
# ════════════════════════════════════════════════════════════════════════════
s2_rf = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
         .filterBounds(buffer_cuenca)
         .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 5))
         .filterDate("2023-01-01", "2024-12-31")
         .median().clip(buffer_cuenca))

features_rf = (s2_rf.select(["B2","B3","B4","B5","B6","B7","B8","B8A","B11","B12"])
               .addBands(ndvi).addBands(ndwi).addBands(mndwi)
               .addBands(bsi).addBands(slope_pct).addBands(dem.rename("elevacion")))

# Muestrar píxeles para entrenamiento supervisado
muestra_rf = features_rf.sample(
    region=buffer_cuenca, scale=20, numPixels=2000, seed=42
).map(lambda f: f.set("clase_worldcover",
      worldcover.reduceRegion(ee.Reducer.first(), f.geometry(), 10).get("Map")))

# Random Forest GEE (100 árboles) — ajustar 'inputProperties' según bandNames()
# rf = ee.Classifier.smileRandomForest(numberOfTrees=100, seed=42)
# rf_trained = rf.train(muestra_rf, "clase_worldcover", features_rf.bandNames())
# clasificacion = features_rf.classify(rf_trained).rename("clase_rf")

# Exportar features para TensorFlow externo
ee.batch.Export.table.toDrive(
    collection=muestra_rf.limit(2000),
    description="Tabla4_Features_TensorFlow_RF",
    folder=EXPORT_FOLDER, fileFormat="CSV"
).start()

ee.batch.Export.image.toDrive(
    image=features_rf,
    description="Morf14_Features_ML_Stack",
    folder=EXPORT_FOLDER, region=buffer_cuenca,
    scale=20, maxPixels=1e13, crs=CRS_WGS84
).start()

# En Google Colab con TensorFlow/scikit-learn:
#   import pandas as pd; from sklearn.ensemble import RandomForestClassifier
#   df = pd.read_csv("Tabla4_Features_TensorFlow_RF.csv").dropna()
#   X = df.drop(columns=["clase_worldcover",".geo","system:index"])
#   y = df["clase_worldcover"].astype(int)
#   rf = RandomForestClassifier(n_estimators=100, random_state=42).fit(X, y)
print("✓ Parte XI: Features ML para TensorFlow/RF → Google Colab")

# ════════════════════════════════════════════════════════════════════════════
# PARTE XII — ÍNDICE DE RIESGO INTEGRADO DE MOVILIDAD DE SEDIMENTOS
# Variables: pendiente (30%) + JRC ocurrencia (25%) + NDWI (20%)
#            + BSI (15%) + TPI negativo (10%)
# ════════════════════════════════════════════════════════════════════════════
r_slope  = slope_pct.subtract(2).divide(18).clamp(0, 1).multiply(30)
r_jrc    = ocurrencia.divide(100).clamp(0, 1).multiply(25)
r_ndwi_r = ndwi.add(0.3).divide(0.8).clamp(0, 1).multiply(20)
r_bsi    = bsi.add(0.5).divide(1.0).clamp(0, 1).multiply(15)
r_tpi    = tpi.multiply(-1).add(20).divide(40).clamp(0, 1).multiply(10)

indice_riesgo = (r_slope.add(r_jrc).add(r_ndwi_r).add(r_bsi).add(r_tpi)
                 .clamp(0, 100).rename("riesgo_movilidad"))

ee.batch.Export.image.toDrive(
    image=indice_riesgo,
    description="Morf15_Riesgo_Movilidad_Integrado",
    folder=EXPORT_FOLDER, region=buffer_cuenca,
    scale=SCALE_30, maxPixels=1e13, crs=CRS_WGS84
).start()
print("✓ Parte XII: Índice de riesgo integrado (5 variables, 0-100)")

# ════════════════════════════════════════════════════════════════════════════
# RESUMEN FINAL
# ════════════════════════════════════════════════════════════════════════════
print()
print("=" * 70)
print("  EXPORTACIÓN COMPLETA — 15 CAPAS RASTER + 4 TABLAS CSV")
print(f"  Carpeta: Google Drive / {EXPORT_FOLDER}")
print()
print("  CAPAS GeoTIFF (QGIS / ArcGIS):")
print("  01. DEM 30m crudo              → HEC-RAS RAS Mapper, iRIC")
print("  02. DEM + Hillshade visual     → QGIS cartografía")
print("  03. Stack topo (slope/asp/TPI) → iRIC morfodinámica")
print("  04. Red hidrográfica           → SWAT+ delimitación")
print("  05. Canal MNDWI (máscara)      → HEC-RAS geometría")
print("  06. Manning's n por clase      → HEC-RAS rugosidad")
print("  07. CN AMC-II                  → HEC-HMS infiltración")
print("  08. WorldCover 2021 visual     → QGIS, SWAT+ LULC")
print("  09. Índices S2 multibanda      → QGIS, xarray/rioxarray")
print("  10. JRC zonas de inundación    → HEC-RAS validación planimetría")
print("  11. RGB Landsat 2024/2015/2010 → QGIS migración lateral")
print("  12. Agua binaria 2010-2024     → xarray análisis temporal")
print("  13. Precipitación CHIRPS       → HEC-HMS área de estudio")
print("  14. Features ML multibanda     → TensorFlow / scikit-learn")
print("  15. Índice de riesgo integrado → zonificación, informes")
print()
print("  TABLAS CSV (xarray / pandas / modelos):")
print("  T1. Parámetros morfométricos   → HEC-HMS input")
print("  T2. Precip. CHIRPS mensual     → HEC-HMS calibración")
print("  T3. Serie NDTI 2019-2024       → pandas/matplotlib/Colab")
print("  T4. Features RF/TensorFlow     → clasificación supervisada")
print()
print("  Estado en GEE: https://code.earthengine.google.com/tasks")
print("=" * 70)

# ════════════════════════════════════════════════════════════════════════════
# NOTAS DE INTEGRACIÓN CON HERRAMIENTAS EXTERNAS
# ════════════════════════════════════════════════════════════════════════════
#
# QGIS / ArcGIS:
#   Importar GeoTIFFs como capas raster. Usar GDAL/ogr2ogr para conversión.
#   Proyección: reproyectar a UTM local si es necesario para mediciones.
#
# HEC-RAS (versión 6.x con RAS Mapper):
#   1. New Terrain: importar Morf01_DEM_HEC-RAS_iRIC.tif
#   2. Land Cover: importar Morf06_Manning_n_HEC-RAS.tif como manning raster
#   3. Validar extensión del cauce con Morf05_Canal_Mascara_HEC-RAS.tif
#   4. Comparar calado simulado con Morf10_JRC_Inundacion_HEC-RAS.tif (banda 1)
#
# HEC-HMS:
#   1. Basin file: área de cuenca desde Tabla1_Cuenca_HEC-HMS.csv
#   2. Loss: CN desde Morf07_CN_HEC-HMS_SWAT.tif (promedio de cuenca)
#   3. Precipitation: cargar Tabla2_Precipitacion_CHIRPS_HEC-HMS.csv
#
# SWAT+ (ArcSWAT / QSWAT+):
#   1. DEM: Morf01 como terreno base
#   2. Land Use: Morf08_WorldCover_QGIS_SWAT.tif (reclasificar a códigos SWAT)
#   3. Stream Network: Morf04_Red_Hidrografica_SWAT.tif
#
# iRIC (Nays2DH / MFlow2D):
#   1. Terrain: Morf03_Topografia_Stack_iRIC.tif (banda "dem")
#   2. Roughness: Morf06_Manning_n_HEC-RAS.tif
#   3. Channel mask: Morf05_Canal_Mascara_HEC-RAS.tif para contorno
#
# xarray / rioxarray (Google Colab):
#   import rioxarray as rxr, xarray as xr
#   ds = rxr.open_rasterio("Morf09_Indices_S2_xarray.tif")
#   agua_ts = rxr.open_rasterio("Morf12_Agua_Binario_xarray.tif")
#
# TensorFlow / scikit-learn (Google Colab):
#   df = pd.read_csv("Tabla4_Features_TensorFlow_RF.csv").dropna()
#   from sklearn.ensemble import RandomForestClassifier
#   rf = RandomForestClassifier(n_estimators=200).fit(X_train, y_train)
'''
    return (template
            .replace("<<LAT>>", str(lat))
            .replace("<<LON>>", str(lon))
            .replace("<<D50>>", str(d50))
            .replace("<<D90>>", str(d90)))


@app.route("/")
def index():
    return render_template("index.html", gee_available=GEE_AVAILABLE)


@app.route("/calculate", methods=["POST"])
def calculate():
    data = request.json
    try:
        lat = float(data["lat"])
        lon = float(data["lon"])
        d50_mm = float(data["d50"])
        d90_mm = float(data.get("d90", d50_mm * 2))
        rho_s = float(data.get("rho_s", 2650))
        temp = float(data.get("temp", 20))
        depth = float(data.get("depth", 1.0))
        velocity = float(data.get("velocity", 1.0))

        if depth <= 0 or velocity <= 0 or d50_mm <= 0:
            return jsonify({"error": "Tirante, velocidad y d₅₀ deben ser valores positivos"}), 400

        # GEE or manual slope
        slope = float(data.get("slope", 0.005))
        lc_code = None
        lc_label = "Sin datos GEE"
        ndti = 0.0
        maps = {}
        gee_active = False

        if data.get("use_gee") and GEE_AVAILABLE:
            try:
                slope = get_slope_from_dem(lat, lon)
                lc_code = get_landcover_at_point(lat, lon)
                ndti = get_ndti_turbidity(lat, lon)
                lc_label = WORLDCOVER.get(lc_code, ("Sin clasificar", 0.035))[0]
                maps = {
                    "slope": get_map_url(lat, lon, "slope"),
                    "landcover": get_map_url(lat, lon, "landcover"),
                }
                gee_active = True
            except Exception as e:
                print(f"GEE error: {e}")

        manning_n = WORLDCOVER.get(lc_code, ("", 0.035))[1]

        # Physical properties
        rho_w = calculate_water_density(temp)
        nu = calculate_kinematic_viscosity(temp)
        s = calculate_specific_gravity(rho_s, rho_w)

        # Particle parameters
        ws = calculate_fall_velocity(d50_mm, s, nu)
        dstar = calculate_dimensionless_particle_parameter(d50_mm, s, nu)
        tau_c, theta_c = calculate_critical_shear_stress_shields(dstar, s, rho_w, d50_mm)
        classification, class_color = classify_particle(d50_mm)

        # Hydraulic state
        tau_0 = rho_w * G * depth * slope
        theta_0 = tau_0 / ((s - 1) * rho_w * G * (d50_mm / 1000))
        u_star = math.sqrt(max(tau_0 / rho_w, 0))
        froude = velocity / math.sqrt(G * depth)
        reynolds = velocity * depth / nu
        z_rouse = ws / (0.41 * u_star) if u_star > 1e-9 else 999

        # Sediment transport
        qb_mpm = meyer_peter_muller(s, d50_mm, slope, depth, rho_w)
        qt_eh = engelund_hansen(velocity, depth, slope, d50_mm, s, rho_w)
        qb_vr = van_rijn_bedload(velocity, depth, d50_mm, dstar, s, rho_w, nu)

        # Sensitivity: transport vs depth
        depths_range = [0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0]
        sensitivity = [
            {
                "depth": d,
                "mpm": round(meyer_peter_muller(s, d50_mm, slope, d, rho_w), 8),
                "eh": round(engelund_hansen(velocity, d, slope, d50_mm, s, rho_w), 8),
                "vr": round(van_rijn_bedload(velocity, d, d50_mm, dstar, s, rho_w, nu), 8),
            }
            for d in depths_range
        ]

        return jsonify(
            {
                "physical": {
                    "rho_w": round(rho_w, 2),
                    "nu": f"{nu:.4e}",
                    "s": round(s, 4),
                },
                "particle": {
                    "ws": round(ws, 6),
                    "dstar": round(dstar, 2),
                    "tau_c": round(tau_c, 4),
                    "theta_c": round(theta_c, 4),
                    "classification": classification,
                    "class_color": class_color,
                },
                "hydraulic": {
                    "tau_0": round(tau_0, 4),
                    "theta_0": round(theta_0, 4),
                    "mobile": bool(theta_0 > theta_c),
                    "u_star": round(u_star, 5),
                    "froude": round(froude, 4),
                    "reynolds": f"{reynolds:.2e}",
                    "z_rouse": round(z_rouse, 2),
                    "transport_mode": rouse_mode(z_rouse),
                    "manning_n": round(manning_n, 4),
                },
                "gee": {
                    "slope": round(slope, 8),
                    "slope_pct": round(slope * 100, 5),
                    "landcover": lc_label,
                    "ndti": round(ndti, 4),
                    "gee_active": gee_active,
                },
                "transport": {
                    "meyer_peter_muller": round(qb_mpm, 8),
                    "engelund_hansen": round(qt_eh, 8),
                    "van_rijn": round(qb_vr, 8),
                },
                "sensitivity": sensitivity,
                "maps": maps,
            }
        )

    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({"error": str(e)}), 400


@app.route("/report")
def report():
    try:
        lat = float(request.args.get("lat", -16.5))
        lon = float(request.args.get("lon", -68.15))
        d50 = float(request.args.get("d50", 0.45))
        d90 = float(request.args.get("d90", 0.9))
        rho_s = float(request.args.get("rho_s", 2650))
        depth = float(request.args.get("depth", 1.0))
        velocity = float(request.args.get("velocity", 1.0))
        temp = float(request.args.get("temp", 20))
        slope = float(request.args.get("slope", 0.005))

        rho_w = calculate_water_density(temp)
        nu = calculate_kinematic_viscosity(temp)
        s = calculate_specific_gravity(rho_s, rho_w)
        ws = calculate_fall_velocity(d50, s, nu)
        dstar = calculate_dimensionless_particle_parameter(d50, s, nu)
        tau_c, theta_c = calculate_critical_shear_stress_shields(dstar, s, rho_w, d50)
        tau_0 = rho_w * G * depth * slope
        theta_0 = tau_0 / ((s - 1) * rho_w * G * (d50 / 1000))
        u_star = math.sqrt(max(tau_0 / rho_w, 0))
        froude = velocity / math.sqrt(G * depth)
        classification, _ = classify_particle(d50)
        qb_mpm = meyer_peter_muller(s, d50, slope, depth, rho_w)
        qt_eh = engelund_hansen(velocity, depth, slope, d50, s, rho_w)
        qb_vr = van_rijn_bedload(velocity, depth, d50, dstar, s, rho_w, nu)
        z_rouse = ws / (0.41 * u_star) if u_star > 1e-9 else 999
        ratio = round(theta_0 / theta_c, 3) if theta_c > 0 else "∞"
        rosgen_type, rosgen_desc = rosgen_classify(slope, d50, froude)

        depths_range = [0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0]
        sensitivity = [
            {
                "depth": d,
                "mpm": round(meyer_peter_muller(s, d50, slope, d, rho_w), 8),
                "eh": round(engelund_hansen(velocity, d, slope, d50, s, rho_w), 8),
                "vr": round(van_rijn_bedload(velocity, d, d50, dstar, s, rho_w, nu), 8),
            }
            for d in depths_range
        ]

        return render_template(
            "report.html",
            results={
                "date": datetime.now().strftime("%d/%m/%Y %H:%M"),
                "lat": lat, "lon": lon,
                "d50": d50, "d90": d90, "rho_s": rho_s,
                "temp": temp, "depth": depth, "velocity": velocity,
                "slope": round(slope, 8),
                "slope_pct": round(slope * 100, 5),
                "rho_w": round(rho_w, 2),
                "nu": f"{nu:.4e}",
                "s": round(s, 4),
                "ws": round(ws, 6),
                "dstar": round(dstar, 2),
                "tau_c": round(tau_c, 4),
                "theta_c": round(theta_c, 4),
                "tau_0": round(tau_0, 4),
                "theta_0": round(theta_0, 4),
                "u_star": round(u_star, 5),
                "froude": round(froude, 4),
                "z_rouse": round(z_rouse, 2),
                "transport_mode": rouse_mode(z_rouse),
                "classification": classification,
                "mobile": theta_0 > theta_c,
                "ratio": ratio,
                "rosgen_type": rosgen_type,
                "rosgen_desc": rosgen_desc,
                "sensitivity": sensitivity,
                "transport": {
                    "meyer_peter_muller": round(qb_mpm, 8),
                    "engelund_hansen": round(qt_eh, 8),
                    "van_rijn": round(qb_vr, 8),
                },
            },
        )
    except Exception as e:
        return str(e), 400


@app.route("/gee_code")
def gee_code():
    lat = float(request.args.get("lat", -16.5))
    lon = float(request.args.get("lon", -68.15))
    d50 = float(request.args.get("d50", 0.45))
    d90 = float(request.args.get("d90", 0.9))
    code = generate_gee_code(lat, lon, d50, d90)
    return render_template("gee_maps.html", lat=lat, lon=lon, d50=d50, d90=d90, gee_code=code)


if __name__ == "__main__":
    app.run(debug=True, port=5000)
