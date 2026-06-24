import ee

def initialize_gee():
    """
    Initializes Google Earth Engine.
    Requires service account credentials or user authentication in a real environment.
    For this sandbox, we assume ee is already authenticated or we mock it for tests.
    """
    try:
        ee.Initialize()
        return True
    except Exception as e:
        print(f"GEE Initialization failed: {e}")
        return False

def get_slope_from_dem(lat, lon, buffer_m=100):
    """
    Fetches the average slope at a given location using SRTM DEM.
    """
    point = ee.Geometry.Point([lon, lat])
    dem = ee.Image("USGS/SRTMGL1_003")
    slope_img = ee.Terrain.slope(dem)

    # Reduce region to get average slope in a small buffer
    stats = slope_img.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=point.buffer(buffer_m),
        scale=30
    ).getInfo()

    # Slope is in degrees, convert to m/m
    slope_deg = stats.get('slope', 0)
    slope_m_m = (slope_deg * 3.14159 / 180.0) # tan(theta) approx theta for small angles
    return slope_m_m

def get_map_url(lat, lon, layer_type='slope', zoom=15):
    """
    Generates a thumbnail or map URL for the given location and layer.
    """
    point = ee.Geometry.Point([lon, lat])
    region = point.buffer(2000).bounds()

    if layer_type == 'slope':
        dem = ee.Image("USGS/SRTMGL1_003")
        img = ee.Terrain.slope(dem)
        vis_params = {'min': 0, 'max': 30, 'palette': ['blue', 'green', 'red']}
    elif layer_type == 'landcover':
        img = ee.Image("ESA/WorldCover/v100/2020").select('Map')
        vis_params = {} # Default palette
    else:
        img = ee.Image("USGS/SRTMGL1_003")
        vis_params = {}

    # Using getThumbURL for static map previews
    params = {
        'region': region,
        'dimensions': 512,
        'format': 'png',
        'min': vis_params.get('min', 0),
        'max': vis_params.get('max', 30),
        'palette': ','.join(vis_params.get('palette', ['blue', 'green', 'red']))
    }
    return img.getThumbURL(params)

def get_landcover_at_point(lat, lon):
    """
    Gets the ESA WorldCover class at a point.
    """
    point = ee.Geometry.Point([lon, lat])
    lc = ee.Image("ESA/WorldCover/v100/2020").select('Map')
    val = lc.reduceRegion(
        reducer=ee.Reducer.first(),
        geometry=point,
        scale=10
    ).getInfo()
    return val.get('Map', 'Unknown')

def get_manning_n_from_lc(lc_class):
    """
    Mapping ESA WorldCover classes to estimated Manning's n.
    Source: Chow (1959) estimates for similar land covers.
    """
    # 10: Trees, 20: Shrubland, 30: Grassland, 40: Cropland, 50: Built-up,
    # 60: Bare/sparse vegetation, 70: Snow/ice, 80: Permanent water,
    # 90: Herbaceous wetland, 95: Mangroves, 100: Moss and lichen
    mapping = {
        10: 0.080, # Densely wooded
        20: 0.050, # Shrubland
        30: 0.035, # Grassland
        40: 0.040, # Cropland
        50: 0.015, # Concrete/Urban
        60: 0.030, # Bare soil
        80: 0.025, # Clean water channel
        90: 0.060, # Wetland
    }
    return mapping.get(lc_class, 0.035) # Default to 0.035

def get_ndti_turbidity(lat, lon):
    """
    Estimates NDTI (Normalized Difference Turbidity Index) using Sentinel-2.
    NDTI = (Red - Green) / (Red + Green)
    """
    point = ee.Geometry.Point([lon, lat])
    s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED") \
        .filterBounds(point) \
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20)) \
        .sort('system:time_start', False) \
        .first()

    if not s2:
        return 0.0

    ndti = s2.normalizedDifference(['B4', 'B3']) # Red: B4, Green: B3
    val = ndti.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=point.buffer(50),
        scale=10
    ).getInfo()
    return val.get('nd', 0.0)
