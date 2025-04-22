# Event Significance Index (ESI)

## Overview

The Event Significance Index (ESI) quantifies how "special" a grouping of planets is based on:

- Physical Proximity on Earth's surface
- Planetary Mass
- Geocentric Distance from Earth

A higher ESI indicates a tighter cluster of massive bodies that are nearer to Earth.

## Definitions

- **N**: Number of selected planets  
- **(φᵢ, λᵢ)**: Latitude and longitude of planet *i* (degrees)  
- **dᵢⱼ**: Great‑circle distance between planets *i* and *j* (km)  
- **Dᵢ**: Geocentric distance from Earth to planet *i* (AU)  
- **mᵢ**: Mass of planet *i* (kg)  

## Calculations

1. **Cluster Proximity (P)**  
   ```text
   P = (2 / (N*(N-1))) * ∑_{i<j} 1 / (dᵢⱼ + ε)
   ```
   - ε is a small constant (e.g., 1e-3 km) to avoid division by zero.

2. **Mass‑Distance Weight (W)**  
   ```text
   W = ∑_{i=1..N} (mᵢ / (Dᵢ * AU_IN_KM + ε))
   ```
   - Convert geocentric distance from AU to km using `AU_IN_KM`.

3. **Combined Index**  
   ```text
   ESI = P * W
   ```

## Implementation

Add the following to **astro_utils.py**:

```python
# Constants
PLANET_MASS = {
    'Sun': 1.9885e30,
    'Moon': 7.3477e22,
    'Mercury': 3.3011e23,
    # … other bodies
}
EARTH_RADIUS = 6371.0  # km
AU_IN_KM = 149597870.7  # km

def great_circle_distance(lat1, lon1, lat2, lon2):
    # Returns distance in km between two lat/lon points
    from math import radians, sin, cos, acos
    φ1, λ1, φ2, λ2 = map(radians, (lat1, lon1, lat2, lon2))
    return EARTH_RADIUS * acos(
        sin(φ1)*sin(φ2) +
        cos(φ1)*cos(φ2)*cos(λ1 - λ2)
    )

def compute_esi(planet_positions):
    names = list(planet_positions)
    N = len(names)
    # 1) Proximity P
    inv_d_sum = 0.0
    eps = 1e-3
    for i in range(N):
        for j in range(i+1, N):
            pi, pj = names[i], names[j]
            lat1, lon1 = planet_positions[pi]['lat'], planet_positions[pi]['lon']
            lat2, lon2 = planet_positions[pj]['lat'], planet_positions[pj]['lon']
            d = great_circle_distance(lat1, lon1, lat2, lon2)
            inv_d_sum += 1.0 / (d + eps)
    P = (2.0 / (N*(N-1))) * inv_d_sum

    # 2) Weight W
    W = 0.0
    for name in names:
        m = PLANET_MASS[name]
        D_km = planet_positions[name]['distance_au'] * AU_IN_KM
        W += m / (D_km + eps)

    # 3) Combined index
    return P * W
```

In **planets_app.py**, after computing `planet_data`:

```python
from astro_utils import compute_esi

# Extend planet_data with distance_au
for planet in selected_planets:
    body = getattr(ephem, planet)()
    body.compute(selected_datetime)
    planet_data[planet]['distance_au'] = body.earth_distance

# Compute Event Significance Index
esi_value = compute_esi({
    p: {
        'lat': planet_data[p]['current_lat'],
        'lon': planet_data[p]['current_lon'],
        'distance_au': planet_data[p]['distance_au']
    } for p in selected_planets
})

st.sidebar.metric("Event Significance Index (ESI)", f"{esi_value:.2e}")
```

## Usage

- Select planets and date/time as usual.  
- The sidebar will display the ESI.  
- Higher values indicate more significant clusterings.  
- For normalized values (0–1), divide by a historical maximum or user‑configured scale.

## Tuning

- Adjust `eps` to change numerical stability.  
- Use `1/d**2` in proximity for steeper sensitivity.  
- Experiment with geometric mean: `ESI = sqrt(P * W)`.

---