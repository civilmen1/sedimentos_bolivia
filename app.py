from flask import Flask, render_template, request, jsonify
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
    van_rijn_bedload,
    calculate_hydraulic_downscaling
)
from utils.gee_handler import (
    initialize_gee,
    get_slope_from_dem,
    get_map_url,
    get_landcover_at_point,
    get_manning_n_from_lc,
    get_ndti_turbidity
)

app = Flask(__name__)

# Initialize GEE at startup
gee_available = initialize_gee()

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

        # 1. Fetch data from GEE
        if gee_available:
            try:
                slope = get_slope_from_dem(lat, lon)
                landcover = get_landcover_at_point(lat, lon)
                ndti = get_ndti_turbidity(lat, lon)
            except Exception as gee_err:
                print(f"GEE Fetch error (using fallbacks): {gee_err}")
                slope = 0.005
                landcover = 80 # Water
                ndti = 0.0
        else:
            print("GEE not available, using fallbacks")
            slope = 0.005
            landcover = 80 # Water
            ndti = 0.0

        # Hydraulic variables
        hydraulic_results = {}
        if data.get('use_downscaling'):
            q = float(data.get('q', 10.0))
            b = float(data.get('b', 5.0))
            n_manning = get_manning_n_from_lc(landcover)
            hydraulic_results = calculate_hydraulic_downscaling(q, b, slope, n_manning)
            depth = hydraulic_results['depth']
            velocity = hydraulic_results['velocity']
        else:
            depth = float(data.get('depth', 1.0))
            velocity = float(data.get('velocity', 1.0))

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
        if gee_available:
            try:
                maps = {
                    'slope': get_map_url(lat, lon, 'slope'),
                    'landcover': get_map_url(lat, lon, 'landcover')
                }
            except:
                maps = {'slope': '', 'landcover': ''}
        else:
            maps = {'slope': '', 'landcover': ''}

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
            'hydraulic': {
                'depth': round(depth, 3),
                'velocity': round(velocity, 3),
                'critical_depth': round(hydraulic_results.get('critical_depth', 0), 3),
                'froude': round(hydraulic_results.get('froude', 0), 3),
                'regime': hydraulic_results.get('regime', 'N/A')
            },
            'maps': maps
        }

        return jsonify(results)

    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'error': str(e)}), 400

@app.route('/report')
def report():
    try:
        lat = float(request.args.get('lat', -12.0464))
        lon = float(request.args.get('lon', -77.0428))
        d50 = float(request.args.get('d50', 0.45))
        depth = float(request.args.get('depth', 1.0))
        velocity = float(request.args.get('velocity', 1.0))
        temp = float(request.args.get('temp', 20))

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

        results = {
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
        return render_template('report.html', results=results)
    except Exception as e:
        return str(e), 400

if __name__ == '__main__':
    app.run(debug=True, port=5000)
