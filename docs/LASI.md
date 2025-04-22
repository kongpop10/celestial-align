# Local Alignment Significance Index (LASI)

## Overview
The Local Alignment Significance Index (LASI) measures how significant a celestial alignment appears from a specific location on Earth at a given date/time.

## Definitions
- **Observer Position (φ₀, λ₀)**: latitude and longitude in degrees.
- **Date/Time (t)**: UTC datetime.
- **N**: number of selected celestial bodies.
- **(Aᵢ, aᵢ)**: azimuth (Aᵢ) and altitude (aᵢ) of body i from observer.
- **Dᵢ**: geocentric distance of body i in astronomical units (AU).
- **mᵢ**: mass of body i in kg.
- **ε**: small constant to avoid divide-by-zero (e.g. 1e-3).

## Mathematical Formulation
1. Angular separation:  
   θᵢⱼ = arccos[ sin(aᵢ) sin(aⱼ) + cos(aᵢ) cos(aⱼ) cos(Aᵢ − Aⱼ) ]

2. Cluster Proximity P_loc:  
   P_loc = (2 / [N(N−1)]) · Σ₍ᵢ<ⱼ₎ 1 / (θᵢⱼ + ε)

3. Mass–Distance Weight W_loc:  
   W_loc = Σₖ mₖ / (Dₖ * AU_IN_KM + ε)

4. LASI calculation:  
   LASI(φ₀,λ₀,t) = P_loc × W_loc

5. Normalization (optional):  
   LASI_norm = (LASI − LASI_min) / (LASI_max − LASI_min)

## Pseudocode
```python
from astro_utils import compute_lasi
lat0, lon0, dt = 37.7749, -122.4194, datetime.utcnow()
bodies = ['Sun','Moon','Mars']
lasi = compute_lasi(lat0, lon0, dt, bodies)
```

## References
- Implementation in `astro_utils.compute_lasi` and `planets_appV2.py`.