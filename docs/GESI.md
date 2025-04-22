# Ground-Based Event Significance Index (GESI)

## Overview
The Ground-Based Event Significance Index (GESI) extends the geocentric ESI by evaluating how “significant” a clustering of celestial bodies appears from a specific point on Earth’s surface at a given date/time.

## Definitions
- **Observer Position (φ₀, λ₀)**: Geographic latitude (φ₀) and longitude (λ₀) in degrees.  
- **Date/Time (t)**: UTC datetime of observation.  
- **N**: Number of selected celestial bodies.  
- **(Aᵢ, aᵢ)**: Azimuth (Aᵢ) and altitude (aᵢ) of body *i* from observer.  
- **Dᵢ**: Geocentric distance of body *i* in astronomical units (AU).  
- **mᵢ**: Mass of body *i* in kilograms (kg).  
- **ε**: Small constant to avoid division by zero (e.g., 1e-3).

## Mathematical Formulation

1. **Angular Separation (θᵢⱼ)**  
   θᵢⱼ = arccos [ sin(aᵢ)·sin(aⱼ) + cos(aᵢ)·cos(aⱼ)·cos(Aᵢ − Aⱼ) ]

2. **Cluster Proximity (P_loc)**  
   P_loc = (2 / (N·(N−1))) · ∑₍ᵢ<ⱼ₎ 1 / (θᵢⱼ + ε)

3. **Mass‑Distance Weight (W)**  
   W = ∑ᵢ mᵢ / (Dᵢ · AU_IN_KM + ε)

4. **Ground‑Based Index (GESI)**  
   GESI(φ₀,λ₀,t) = P_loc × W

## Computation (Pseudo‑Code)
```python
from ephem import Observer
from math import radians, sin, cos, acos
from itertools import combinations

def compute_gesi(lat0, lon0, dt, bodies):
    obs = Observer()
    obs.lat, obs.lon, obs.date = radians(lat0), radians(lon0), dt

    positions = {}
    for name in bodies:
        body = getattr(ephem, name)(obs)
        body.compute(obs)
        positions[name] = {
            'az': float(body.az),
            'alt': float(body.alt),
            'distance_au': body.earth_distance
        }

    # 1) Angular separations
    theta = {}
    for i, j in combinations(bodies, 2):
        ai, aj = positions[i]['alt'], positions[j]['alt']
        Ai, Aj = positions[i]['az'], positions[j]['az']
        theta[(i, j)] = acos(sin(ai)*sin(aj) + cos(ai)*cos(aj)*cos(Ai - Aj))

    # 2) Cluster proximity
    N = len(bodies)
    P_loc = (2.0 / (N*(N-1))) * sum(1.0/(t+1e-3) for t in theta.values())

    # 3) Mass‑distance weight
    W = sum(MASS[name] / (pos['distance_au']*AU_IN_KM + 1e-3)
            for name, pos in positions.items())

    return P_loc * W
```

## UI Visualization
- **World Map Heatmap**: Render GESI as a heatmap overlay on a global map (e.g., Plotly, Folium).  
- **Grid Sampling**: Compute GESI over a regular latitude/longitude grid (e.g., 1°×1°).  
- **Interactive Slider**: Allow users to change date/time and update the heatmap live.  
- **Best Location Marker**: Highlight the point with maximum GESI on each frame.

### Example (Streamlit + PyDeck)
```python
import streamlit as st
import pydeck as pdk

st.title("Ground‑Based Event Significance Index")
dt = st.slider("Select Date/Time", ...)
latitudes = [...]
longitudes = [...]
gesi_grid = compute_grid_gesi(latitudes, longitudes, dt, bodies)

layer = pdk.Layer(
    "HeatmapLayer",
    data=gesi_grid,
    get_position=["lon", "lat"],
    get_weight="gesi",
    radiusPixels=50,
)

st.pydeck_chart(pdk.Deck(layers=[layer], initial_view_state=...))
```

## Next Steps
- Optimize sampling (e.g., adaptive or higher resolution).  
- Support dynamic body selection and time ranges.  
- Combine GESI with normalized ESI for comparative analysis.

## References
- See `ESI_documentation.md` for geocentric ESI details.  
- Implementation in `astro_utils.py` and `planets_app.py`.