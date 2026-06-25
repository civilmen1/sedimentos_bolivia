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
    template = '''import ee
import math

# ─── Configuración del área de estudio ──────────────────────────────────────
LAT    = <<LAT>>    # Latitud del punto de muestreo
LON    = <<LON>>    # Longitud del punto de muestreo
D50_MM = <<D50>>    # Diámetro mediano d₅₀ (mm)
D90_MM = <<D90>>    # Diámetro percentil 90 d₉₀ (mm)
BUFFER_M = 5000     # Radio de análisis (m)
EXPORT_FOLDER = "GEE_Sedimentos_Bolivia"
SCALE = 30          # Resolución de exportación (m)

# Inicializar Google Earth Engine
ee.Authenticate()
ee.Initialize()

punto = ee.Geometry.Point([LON, LAT])
area  = punto.buffer(BUFFER_M).bounds()

print(f"Área de estudio: LAT={LAT}, LON={LON}")
print(f"Buffer: {BUFFER_M} m | Resolución: {SCALE} m")
print("=" * 60)

# ────────────────────────────────────────────────────────────────────────────
# MAPA 1: Modelo de Elevación Digital (SRTM 30 m) + Hillshade
# ────────────────────────────────────────────────────────────────────────────
dem       = ee.Image("USGS/SRTMGL1_003")
hillshade = ee.Terrain.hillshade(dem)
pendiente = ee.Terrain.slope(dem)

vis_dem  = dict(min=0, max=5000, palette=["006633","E5FFCC","662A00","D8D8D8","F5F5F5"])
dem_vis  = dem.visualize(**vis_dem)
hill_vis = hillshade.visualize(min=0, max=255, gamma=1.3)
dem_blend = dem_vis.blend(hill_vis.updateMask(ee.Image(0.4)))

ee.batch.Export.image.toDrive(
    image=dem_blend,
    description="Mapa1_DEM_Hillshade",
    folder=EXPORT_FOLDER, region=area, scale=SCALE, maxPixels=1e13
).start()
print("✓ Mapa 1: DEM + Hillshade exportado")

# ────────────────────────────────────────────────────────────────────────────
# MAPA 2: Red Hidrográfica (HydroSHEDS Flow Accumulation 15s)
# ────────────────────────────────────────────────────────────────────────────
flow_acc = ee.Image("WWF/HydroSHEDS/15ACC")
rios     = flow_acc.gte(1000).selfMask()  # umbral: 1 000 celdas
rios_sobre_dem = dem_blend.blend(rios.visualize(min=0, max=1, palette=["0000FF"]))

ee.batch.Export.image.toDrive(
    image=rios_sobre_dem,
    description="Mapa2_Red_Hidrografica",
    folder=EXPORT_FOLDER, region=area, scale=SCALE, maxPixels=1e13
).start()
print("✓ Mapa 2: Red Hidrográfica (HydroSHEDS) exportado")

# ────────────────────────────────────────────────────────────────────────────
# MAPA 3: Pendiente del Lecho (%)
# ────────────────────────────────────────────────────────────────────────────
pct_slope = pendiente.multiply(math.pi / 180).tan().multiply(100).rename("pendiente_pct")
vis_slope = dict(
    min=0, max=20,
    palette=["FFFDE7","FFF59D","FFEE58","FDD835","F9A825","F57F17","E65100","BF360C","7F0000"]
)

ee.batch.Export.image.toDrive(
    image=pct_slope,
    description="Mapa3_Pendiente_Pct",
    folder=EXPORT_FOLDER, region=area, scale=SCALE, maxPixels=1e13
).start()
print("✓ Mapa 3: Pendiente del Lecho (%) exportado")

# ────────────────────────────────────────────────────────────────────────────
# MAPA 4: Cobertura y Uso del Suelo — ESA WorldCover 2021
# ────────────────────────────────────────────────────────────────────────────
worldcover = ee.Image("ESA/WorldCover/v200/2021").select("Map").clip(area)
vis_wc = dict(
    min=10, max=100,
    palette=["006400","FFBB22","FFFF4C","F096FF","FA0000",
             "B4B4B4","F0F0F0","0064C8","0096A0","00CF75","FAE6A0"]
)

ee.batch.Export.image.toDrive(
    image=worldcover.visualize(**vis_wc),
    description="Mapa4_Cobertura_USO_Suelo",
    folder=EXPORT_FOLDER, region=area, scale=10, maxPixels=1e13
).start()
print("✓ Mapa 4: Cobertura y Uso del Suelo (ESA WorldCover 2021) exportado")

# ────────────────────────────────────────────────────────────────────────────
# Sentinel-2 SR — mediana libre de nubes (2022–2024)
# ────────────────────────────────────────────────────────────────────────────
s2 = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
      .filterBounds(punto)
      .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 10))
      .filterDate("2022-01-01", "2024-12-31")
      .median()
      .clip(area))

# ────────────────────────────────────────────────────────────────────────────
# MAPA 5: NDVI — Vegetación Ribereña
# NDVI = (NIR − Red) / (NIR + Red)  ←  B8 = NIR, B4 = Red
# ────────────────────────────────────────────────────────────────────────────
ndvi = s2.normalizedDifference(["B8", "B4"]).rename("NDVI")
vis_ndvi = dict(
    min=-0.2, max=0.8,
    palette=["d73027","f46d43","fdae61","fee08b","d9ef8b","a6d96a","66bd63","1a9850"]
)

ee.batch.Export.image.toDrive(
    image=ndvi,
    description="Mapa5_NDVI_Vegetacion",
    folder=EXPORT_FOLDER, region=area, scale=10, maxPixels=1e13
).start()
print("✓ Mapa 5: NDVI — Vegetación Ribereña exportado")

# ────────────────────────────────────────────────────────────────────────────
# MAPA 6: NDWI — Cuerpos de Agua
# NDWI = (Green − NIR) / (Green + NIR)  ←  B3 = Green, B8 = NIR
# ────────────────────────────────────────────────────────────────────────────
ndwi = s2.normalizedDifference(["B3", "B8"]).rename("NDWI")
vis_ndwi = dict(
    min=-0.3, max=0.5,
    palette=["8B4513","DEB887","FFFFFF","87CEEB","1E90FF","000080"]
)

ee.batch.Export.image.toDrive(
    image=ndwi,
    description="Mapa6_NDWI_Agua",
    folder=EXPORT_FOLDER, region=area, scale=10, maxPixels=1e13
).start()
print("✓ Mapa 6: NDWI — Cuerpos de Agua exportado")

# ────────────────────────────────────────────────────────────────────────────
# MAPA 7: NDTI — Índice de Turbidez / Sedimentos en Suspensión
# NDTI = (Red − Green) / (Red + Green)  ←  B4 = Red, B3 = Green
# ────────────────────────────────────────────────────────────────────────────
ndti = s2.normalizedDifference(["B4", "B3"]).rename("NDTI")
vis_ndti = dict(
    min=-0.2, max=0.4,
    palette=["313695","4575b4","74add1","abd9e9","e0f3f8","fee090","fdae61","f46d43","d73027"]
)

ee.batch.Export.image.toDrive(
    image=ndti,
    description="Mapa7_NDTI_Turbidez",
    folder=EXPORT_FOLDER, region=area, scale=10, maxPixels=1e13
).start()
print("✓ Mapa 7: NDTI — Turbidez / Sedimentos en Suspensión exportado")

# ────────────────────────────────────────────────────────────────────────────
# MAPA 8: Análisis Multitemporal — Migración Lateral del Cauce (Landsat 8/9)
# Composición RGB falso color: R=2024, G=2019, B=2014
# Canales que cambiaron de color → zonas de migración lateral
# ────────────────────────────────────────────────────────────────────────────
def get_landsat_sr(year, point, aoi):
    col = (ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
           .merge(ee.ImageCollection("LANDSAT/LC08/C02/T1_L2"))
           .filterBounds(point)
           .filterDate(str(year) + "-06-01", str(year) + "-09-30")
           .filter(ee.Filter.lt("CLOUD_COVER", 20))
           .map(lambda img: img.multiply(0.0000275).add(-0.2).clip(aoi))
           .median()
           .select(["SR_B4","SR_B3","SR_B2"], ["R","G","B"]))
    return col

sc_2014 = get_landsat_sr(2014, punto, area)
sc_2019 = get_landsat_sr(2019, punto, area)
sc_2024 = get_landsat_sr(2024, punto, area)

rgb_cambio = ee.Image.cat([
    sc_2024.select("R").rename("R"),
    sc_2019.select("R").rename("G"),
    sc_2014.select("R").rename("B"),
])

ee.batch.Export.image.toDrive(
    image=rgb_cambio,
    description="Mapa8_Migracion_Lateral_2014_2019_2024",
    folder=EXPORT_FOLDER, region=area, scale=30, maxPixels=1e13
).start()
print("✓ Mapa 8: Migración Lateral Multitemporal (RGB: 2024/2019/2014) exportado")

# ────────────────────────────────────────────────────────────────────────────
# MAPA 9: Zonas de Riesgo de Movilidad de Sedimentos
# Índice compuesto (0–100):
#   Pendiente   40%  (>2% = riesgo creciente)
#   NDWI        30%  (zonas húmedas/riparias)
#   NDVI inv.   30%  (suelo expuesto = sin cubierta vegetal)
# ────────────────────────────────────────────────────────────────────────────
r_slope = pct_slope.subtract(2).divide(18).clamp(0, 1).multiply(40)
r_ndwi  = ndwi.add(0.3).divide(0.8).clamp(0, 1).multiply(30)
r_ndvi  = ndvi.multiply(-1).add(0.2).divide(0.5).clamp(0, 1).multiply(30)

indice_riesgo = (r_slope.add(r_ndwi).add(r_ndvi)
                 .clamp(0, 100)
                 .rename("riesgo_movilidad"))

vis_riesgo = dict(
    min=0, max=100,
    palette=["1a9850","91cf60","d9ef8b","fee08b","fc8d59","d73027","7f0000"]
)

ee.batch.Export.image.toDrive(
    image=indice_riesgo,
    description="Mapa9_Riesgo_Movilidad_Sedimentos",
    folder=EXPORT_FOLDER, region=area, scale=SCALE, maxPixels=1e13
).start()
print("✓ Mapa 9: Zonas de Riesgo de Movilidad de Sedimentos exportado")

# ────────────────────────────────────────────────────────────────────────────
print()
print("=" * 60)
print("  EXPORTACIÓN COMPLETADA — 9 MAPAS EN COLA")
print(f"  Carpeta destino: Google Drive / {EXPORT_FOLDER}")
print("  Revisar estado: https://code.earthengine.google.com/tasks")
print("=" * 60)
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
