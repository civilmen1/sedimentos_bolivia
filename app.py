from flask import Flask, render_template, request, jsonify, make_response
import numpy as np
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
from utils.pdf_utils import create_pdf

app = Flask(__name__)

# Initialize GEE at startup
initialize_gee()

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/calculate', methods=['POST'])
def calculate():
    data = request.json
    try:
        lat = float(data.get('lat'))
        lon = float(data.get('lon'))
        d50_mm = float(data.get('d50'))
        d90_mm = float(data.get('d90', d50_mm * 2))
        rho_s = float(data.get('rho_s', 2650))
        temp = float(data.get('temp', 20))

        # Hydraulic variables (optional/default)
        depth = float(data.get('depth', 1.0))
        velocity = float(data.get('velocity', 1.0))

        # 1. Fetch data from GEE
        slope = get_slope_from_dem(lat, lon)
        landcover = get_landcover_at_point(lat, lon)
        ndti = get_ndti_turbidity(lat, lon)

        # 2. Calculate physical properties
        rho_w = calculate_water_density(temp)
        nu = calculate_kinematic_viscosity(temp)
        s = calculate_specific_gravity(rho_s, rho_w)

        # 3. Particle parameters
        ws = calculate_fall_velocity(d50_mm, s, nu)
        dstar = calculate_dimensionless_particle_parameter(d50_mm, s, nu)
        tau_c, theta_c = calculate_critical_shear_stress_shields(dstar, s, rho_w, d50_mm)

        # 4. Sediment Transport Models
        qb_mpm = meyer_peter_muller(s, d50_mm, slope, depth, rho_w)
        qt_eh = engelund_hansen(velocity, depth, slope, d50_mm, s, rho_w)
        qb_vr = van_rijn_bedload(velocity, depth, d50_mm, dstar, s, rho_w, nu)

        # 5. Map URLs
        maps = {
            'slope': get_map_url(lat, lon, 'slope'),
            'landcover': get_map_url(lat, lon, 'landcover')
        }

        results = {
            'physical': {
                'rho_w': round(rho_w, 2),
                'nu': f"{nu:.3e}",
                's': round(s, 3)
            },
            'particle': {
                'ws': round(ws, 5),
                'dstar': round(dstar, 2),
                'tau_c': round(tau_c, 3),
                'theta_c': round(theta_c, 4)
            },
            'gee': {
                'slope': round(slope, 6),
                'landcover': landcover,
                'ndti': round(ndti, 4)
            },
            'transport': {
                'meyer_peter_muller': round(qb_mpm, 4),
                'engelund_hansen': round(qt_eh, 4),
                'van_rijn': round(qb_vr, 4)
            },
            'maps': maps
        }

        return jsonify(results)

    except Exception as e:
        return jsonify({'error': str(e)}), 400

def get_report_results(args):
    lat = float(args.get('lat', -12.0464))
    lon = float(args.get('lon', -77.0428))
    d50 = float(args.get('d50', 0.45))
    depth = float(args.get('depth', 1.0))
    velocity = float(args.get('velocity', 1.0))
    temp = float(args.get('temp', 20))

    # Calculate for report
    rho_w = calculate_water_density(temp)
    nu = calculate_kinematic_viscosity(temp)
    s = 2.65
    ws = calculate_fall_velocity(d50, s, nu)
    dstar = calculate_dimensionless_particle_parameter(d50, s, nu)

    # Try to get GEE data for report (simplified for now)
    try:
        slope = get_slope_from_dem(lat, lon)
        landcover = get_landcover_at_point(lat, lon)
        ndti = get_ndti_turbidity(lat, lon)
    except:
        slope, landcover, ndti = 0.005, "Unknown", 0.0

    qb_mpm = meyer_peter_muller(s, d50, slope, depth, rho_w)
    qt_eh = engelund_hansen(velocity, depth, slope, d50, s, rho_w)
    qb_vr = van_rijn_bedload(velocity, depth, d50, dstar, s, rho_w, nu)

    return {
        'lat': lat, 'lon': lon, 'd50': d50, 'temp': temp,
        'depth': depth, 'velocity': velocity,
        'rho_w': round(rho_w, 2),
        'nu': f"{nu:.3e}",
        'ws': round(ws, 5),
        'dstar': round(dstar, 2),
        'slope': round(slope, 6),
        'landcover': landcover,
        'ndti': round(ndti, 4),
        'transport': {
            'meyer_peter_muller': round(qb_mpm, 4),
            'engelund_hansen': round(qt_eh, 4),
            'van_rijn': round(qb_vr, 4)
        }
    }

@app.route('/report')
def report():
    try:
        results = get_report_results(request.args)
        return render_template('report.html', results=results, is_pdf=False)
    except Exception as e:
        return str(e), 400

@app.route('/download_report')
def download_report():
    try:
        results = get_report_results(request.args)

        html = render_template('report.html', results=results, is_pdf=True)
        pdf = create_pdf(html)

        if pdf:
            response = make_response(pdf)
            response.headers['Content-Type'] = 'application/pdf'
            response.headers['Content-Disposition'] = f"attachment; filename=informe_sedimentos_{results['lat']}_{results['lon']}.pdf"
            return response
        else:
            return "Error generating PDF", 500

    except Exception as e:
        return str(e), 400

if __name__ == '__main__':
    app.run(debug=True, port=5000)
