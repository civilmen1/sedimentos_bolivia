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
                "transport": {
                    "meyer_peter_muller": round(qb_mpm, 8),
                    "engelund_hansen": round(qt_eh, 8),
                    "van_rijn": round(qb_vr, 8),
                },
            },
        )
    except Exception as e:
        return str(e), 400


if __name__ == "__main__":
    app.run(debug=True, port=5000)
