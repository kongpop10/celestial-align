import ephem
import pandas as pd
import math
from datetime import timedelta
from math import radians, degrees, sin, cos, asin, atan2, sqrt

# Define Ceres orbital elements (approximate values as of 2023)
# These are simplified and should be updated for more accurate calculations
CERES_ORBITAL_ELEMENTS = {
    'a': 2.7691,  # semi-major axis in AU
    'e': 0.0758,  # eccentricity
    'i': 10.593,  # inclination in degrees
    'O': 80.393,  # longitude of ascending node in degrees
    'w': 73.597,  # argument of perihelion in degrees
    'M': 77.37,   # mean anomaly at epoch in degrees
    'epoch': ephem.date('2023/1/1')  # epoch for the orbital elements
}

def calculate_ceres_position(dt):
    """Calculate the position of Ceres for a given datetime."""
    # Convert datetime to Julian date
    jd = ephem.julian_date(dt)

    # Get orbital elements
    a = CERES_ORBITAL_ELEMENTS['a']
    e = CERES_ORBITAL_ELEMENTS['e']
    i = radians(CERES_ORBITAL_ELEMENTS['i'])
    O = radians(CERES_ORBITAL_ELEMENTS['O'])
    w = radians(CERES_ORBITAL_ELEMENTS['w'])

    # Calculate days since epoch
    days_since_epoch = jd - ephem.julian_date(CERES_ORBITAL_ELEMENTS['epoch'])

    # Mean motion in degrees per day (approximate)
    n = 0.214  # This is an approximation for Ceres

    # Calculate current mean anomaly
    M = radians(CERES_ORBITAL_ELEMENTS['M'] + n * days_since_epoch)

    # Solve Kepler's equation for eccentric anomaly (E)
    # Using a simple iterative approach
    E = M
    for _ in range(10):  # Usually converges quickly
        E = M + e * sin(E)

    # Calculate true anomaly (v)
    v = 2 * atan2(sqrt(1 + e) * sin(E/2), sqrt(1 - e) * cos(E/2))

    # Calculate distance from Sun (r)
    r = a * (1 - e * cos(E))

    # Calculate heliocentric coordinates
    x = r * (cos(O) * cos(v + w) - sin(O) * sin(v + w) * cos(i))
    y = r * (sin(O) * cos(v + w) + cos(O) * sin(v + w) * cos(i))
    z = r * sin(v + w) * sin(i)

    # Convert to geocentric coordinates (simplified)
    # This is a very simplified approach and doesn't account for Earth's position
    # For a more accurate calculation, we would need Earth's position as well
    observer = ephem.Observer()
    observer.date = dt
    sun = ephem.Sun(observer)
    sun.compute(dt)

    # Get Earth's position (from the Sun's geocentric position)
    earth_x = -r * math.cos(float(sun.ra)) * math.cos(float(sun.dec))
    earth_y = -r * math.sin(float(sun.ra)) * math.cos(float(sun.dec))
    earth_z = -r * math.sin(float(sun.dec))

    # Calculate geocentric coordinates
    geo_x = x - earth_x
    geo_y = y - earth_y
    geo_z = z - earth_z

    # Convert to RA and Dec
    geo_r = sqrt(geo_x**2 + geo_y**2 + geo_z**2)
    geo_ra = atan2(geo_y, geo_x)
    geo_dec = asin(geo_z / geo_r)

    # Convert to degrees
    ra_deg = degrees(geo_ra)
    dec_deg = degrees(geo_dec)

    # Calculate Earth distance
    earth_distance = sqrt((x - earth_x)**2 + (y - earth_y)**2 + (z - earth_z)**2)

    return ra_deg, dec_deg, earth_distance

def get_planet_position(dt, body_name):
    observer = ephem.Observer()
    observer.date = dt

    # Special case for the Sun
    if body_name == "Sun":
        sun = ephem.Sun(observer)
        sun.compute(dt)
        lat = math.degrees(float(sun.dec))
        ra = math.degrees(float(sun.ra))
        lst = math.degrees(float(observer.sidereal_time()))
        lon = ra - lst
        lon = ((lon + 180) % 360) - 180
        return lat, lon
    # Special case for Ceres
    elif body_name == "Ceres":
        ra_deg, dec_deg, _ = calculate_ceres_position(dt)
        lat = dec_deg
        lst = math.degrees(float(observer.sidereal_time()))
        lon = ra_deg - lst
        lon = ((lon + 180) % 360) - 180
        return lat, lon
    else:
        # Regular planets and Moon
        body = getattr(ephem, body_name)()
        body.compute(observer)
        lat = math.degrees(float(body.dec))
        ra = math.degrees(float(body.ra))
        lst = math.degrees(float(observer.sidereal_time()))
        lon = ra - lst
        lon = ((lon + 180) % 360) - 180
        return lat, lon

def calculate_planet_path(base_dt, body_name):
    path_data = []
    for hour in range(25):
        dt = base_dt + timedelta(hours=hour)
        lat, lon = get_planet_position(dt, body_name)
        path_data.append({
            'lat': lat,
            'lon': lon,
            'time': dt.strftime('%Y-%m-%d %H:%M UTC'),
            'hour': hour
        })
    return pd.DataFrame(path_data)
# --- ESI Calculation Utilities

PLANET_MASS = {
    "Sun": 1.9885e30,
    "Moon": 7.3477e22,
    "Mercury": 3.3011e23,
    "Venus": 4.8675e24,
    "Mars": 6.4171e23,
    "Jupiter": 1.8982e27,
    "Saturn": 5.6834e26,
    "Uranus": 8.6810e25,
    "Neptune": 1.0243e26,
    "Pluto": 1.303e22,
    "Ceres": 9.3835e20  # Mass of Ceres in kg
}

EARTH_RADIUS = 6371.0  # km
AU_IN_KM = 149597870.7  # km

def great_circle_distance(lat1, lon1, lat2, lon2):
    # Compute great-circle distance in kilometers
    from math import radians, sin, cos, acos
    phi1, lam1, phi2, lam2 = map(radians, (lat1, lon1, lat2, lon2))
    # Normalize distance by Earth's radius: return central angle in radians
    return acos(
        sin(phi1)*sin(phi2) +
        cos(phi1)*cos(phi2)*cos(lam1 - lam2)
    )

def compute_esi(planet_positions):
    """
    Compute the Event Significance Index (ESI) for a set of planets.
    planet_positions: dict of planet -> { 'lat', 'lon', 'distance_au' }
    Returns: float
    """
    names = list(planet_positions.keys())
    N = len(names)
    if N < 2:
        return 0.0

    # 1) Cluster Proximity (P)
    inv_d_sum = 0.0
    eps = 1e-3
    for i in range(N):
        for j in range(i+1, N):
            p_i, p_j = names[i], names[j]
            lat1, lon1 = planet_positions[p_i]['lat'], planet_positions[p_i]['lon']
            lat2, lon2 = planet_positions[p_j]['lat'], planet_positions[p_j]['lon']
            d = great_circle_distance(lat1, lon1, lat2, lon2)
            inv_d_sum += 1.0 / (d + eps)
    P = (2.0 / (N*(N-1))) * inv_d_sum

    # 2) Mass-Distance Weight (W)
    W = 0.0
    for name in names:
        m = PLANET_MASS.get(name, 0.0)
        D_km = planet_positions[name]['distance_au'] * AU_IN_KM
        W += m / (D_km + eps)

    # 3) Combined ESI
    return P * W


# --- Ground-Based Event Significance Index (GESI)
def compute_gesi(lat0, lon0, dt, bodies):
    # Compute Ground-Based Event Significance Index (GESI) for observation from (lat0, lon0) at datetime dt
    from ephem import Observer
    from itertools import combinations
    from math import radians, sin, cos, acos, asin, pi

    # Set up observer
    obs = Observer()
    obs.lat = radians(lat0)
    obs.lon = radians(lon0)
    obs.date = dt

    # Gather azimuth, altitude, and distance for each body
    positions = {}
    for name in bodies:
        # Special case for Ceres
        if name == "Ceres":
            # Get Ceres position using our custom function
            ra_deg, dec_deg, earth_distance = calculate_ceres_position(dt)

            # Convert to radians for calculations
            ra = radians(ra_deg)
            dec = radians(dec_deg)

            # Calculate local sidereal time
            lst = float(obs.sidereal_time())

            # Calculate hour angle
            ha = lst - ra

            # Calculate altitude and azimuth
            sin_alt = sin(dec) * sin(obs.lat) + cos(dec) * cos(obs.lat) * cos(ha)
            alt = asin(sin_alt)

            cos_az = (sin(dec) - sin(obs.lat) * sin_alt) / (cos(obs.lat) * cos(alt))
            # Handle potential numerical errors
            cos_az = max(min(cos_az, 1.0), -1.0)
            az = acos(cos_az)

            # Adjust azimuth based on hour angle
            if sin(ha) > 0:
                az = 2 * pi - az

            positions[name] = {
                'az': az,
                'alt': alt,
                'distance_au': earth_distance
            }
        # Special case for Pluto
        elif name == "Pluto":
            try:
                # Try to use ephem.Pluto if available
                pluto = ephem.Pluto(obs)
                pluto.compute(obs)
                positions[name] = {
                    'az': float(pluto.az),
                    'alt': float(pluto.alt),
                    'distance_au': pluto.earth_distance
                }
            except (AttributeError, TypeError):
                # If ephem.Pluto is not available, skip it
                continue
        else:
            # Regular planets and other bodies
            try:
                body = getattr(ephem, name)(obs)
                body.compute(obs)
                positions[name] = {
                    'az': float(body.az),
                    'alt': float(body.alt),
                    'distance_au': body.earth_distance
                }
            except AttributeError:
                # Skip bodies that aren't available in ephem
                continue

    # If we have fewer than 2 valid bodies, return 0
    if len(positions) < 2:
        return 0.0

    # Calculate angular separations
    theta = []
    valid_bodies = list(positions.keys())
    for i, j in combinations(valid_bodies, 2):
        ai = positions[i]['alt']
        aj = positions[j]['alt']
        Ai = positions[i]['az']
        Aj = positions[j]['az']
        theta.append(acos(sin(ai)*sin(aj) + cos(ai)*cos(aj)*cos(Ai - Aj)))

    # Cluster proximity P_loc
    N = len(valid_bodies)
    eps = 1e-3
    P_loc = (2.0 / (N*(N-1))) * sum(1.0 / (t + eps) for t in theta) if N > 1 else 0.0

    # Mass-distance weight W
    W = 0.0
    for name in valid_bodies:
        m = PLANET_MASS.get(name, 0.0)
        D_km = positions[name]['distance_au'] * AU_IN_KM
        W += m / (D_km + eps)

    return P_loc * W

# --- Local Alignment Significance Index (LASI)
def compute_lasi(lat0, lon0, dt, bodies):
    """
    Compute Local Alignment Significance Index (LASI) for an observer at (lat0, lon0) at datetime dt.
    LASI = P_loc * W_loc, where:
      P_loc = (2 / [N(N-1)]) * Σ_{i<j} 1 / (θᵢⱼ + ε), θᵢⱼ pairwise angular separation
      W_loc = Σ_{k=1..N} mₖ / (Dₖ + ε)
    Parameters:
      lat0: observer latitude in degrees
      lon0: observer longitude in degrees
      dt: datetime for calculation
      bodies: list of body names
    Returns:
      float LASI value
    """
    return compute_gesi(lat0, lon0, dt, bodies)

# --- Gravitational Force Index (GFI)
def compute_gfi(lat0, lon0, dt, bodies):
    """
    Compute Gravitational Force Index (GFI) for a location on Earth at (lat0, lon0) at datetime dt.
    GFI represents the net gravitational force exerted by celestial bodies on a specific location.

    The formula uses Newton's law of universal gravitation: F = G * (m1 * m2) / r^2
    where:
    - G is the gravitational constant (not needed for relative comparisons)
    - m1 is the mass of the celestial body
    - m2 is a unit mass at the Earth location (constant, so not included in calculation)
    - r is the distance from the celestial body to the specific Earth location

    Parameters:
      lat0: location latitude in degrees
      lon0: location longitude in degrees
      dt: datetime for calculation
      bodies: list of body names
    Returns:
      float GFI value (relative, not in physical units)
    """
    import ephem
    from math import radians, sin, cos, sqrt

    # Set up observer
    obs = ephem.Observer()
    obs.lat = radians(lat0)
    obs.lon = radians(lon0)
    obs.date = dt

    # Earth radius in AU for distance calculations
    earth_radius_au = EARTH_RADIUS / AU_IN_KM

    # Gather positions and calculate force vectors
    force_x, force_y, force_z = 0.0, 0.0, 0.0
    eps = 1e-10  # Small constant to avoid division by zero

    for name in bodies:
        # Skip Earth if it's in the list
        if name.lower() == "earth":
            continue

        # Get body position
        if name == "Ceres":
            # Get Ceres position using our custom function
            ra_deg, dec_deg, earth_distance = calculate_ceres_position(dt)
            ra = radians(ra_deg)
            dec = radians(dec_deg)
        elif name == "Pluto":
            try:
                pluto = ephem.Pluto(obs)
                pluto.compute(dt)
                ra = float(pluto.ra)
                dec = float(pluto.dec)
                earth_distance = pluto.earth_distance
            except (AttributeError, TypeError):
                continue
        else:
            try:
                body = getattr(ephem, name)()
                body.compute(dt)
                ra = float(body.ra)
                dec = float(body.dec)
                earth_distance = body.earth_distance
            except AttributeError:
                continue

        # Convert celestial coordinates to cartesian coordinates
        # First, get the geocentric cartesian coordinates of the celestial body
        body_x = earth_distance * cos(dec) * cos(ra)
        body_y = earth_distance * cos(dec) * sin(ra)
        body_z = earth_distance * sin(dec)

        # Now get the cartesian coordinates of the observer location on Earth's surface
        # We need to convert from observer's local coordinates to geocentric coordinates
        # First, get the local sidereal time to rotate the coordinates properly
        lst = float(obs.sidereal_time())

        # Calculate the geocentric position of the observer
        observer_x = earth_radius_au * cos(obs.lat) * cos(lst)
        observer_y = earth_radius_au * cos(obs.lat) * sin(lst)
        observer_z = earth_radius_au * sin(obs.lat)

        # Calculate the vector from the observer to the celestial body
        dx = body_x - observer_x
        dy = body_y - observer_y
        dz = body_z - observer_z

        # Calculate the distance from the observer to the celestial body
        distance = sqrt(dx*dx + dy*dy + dz*dz)

        # Get the mass of the celestial body
        mass = PLANET_MASS.get(name, 0.0)

        # Calculate the gravitational force magnitude (F = G * m1 * m2 / r^2)
        # We omit G and m2 (unit mass) as we're calculating relative forces
        force_magnitude = mass / ((distance * AU_IN_KM) ** 2 + eps)

        # Calculate the unit vector components
        distance_eps = distance + eps  # Avoid division by zero
        unit_x = dx / distance_eps
        unit_y = dy / distance_eps
        unit_z = dz / distance_eps

        # Add the force vector components to the net force
        force_x += force_magnitude * unit_x
        force_y += force_magnitude * unit_y
        force_z += force_magnitude * unit_z

    # Calculate the magnitude of the net force vector
    net_force = sqrt(force_x*force_x + force_y*force_y + force_z*force_z)

    # Apply a logarithmic scale to make the values more differentiable
    # Add 1 before taking log to ensure positive values
    return net_force
