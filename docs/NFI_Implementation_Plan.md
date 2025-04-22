# Net Force Index (NFI) Implementation Plan

**Overview**  
We aim to compute the magnitude of the net gravitational force exerted by selected planets on a 1kg test mass at specified locations on Earth’s surface and visualize it on the globe in `planets_appV3.py`.

## Steps

1. **Define NFI in astro_utils.py**
   - Add constants: Gravitational constant `G`, planet masses, Earth radius.
   - Implement `compute_nfi(lat, lon, dt, selected_planets)`.

2. **Grid Calculation in planets_appV3.py**
   - Create `calculate_nfi_grid(dt, selected_planets, grid_resolution)`.

3. **UI Integration in planets_appV3.py**
   - Remove ESI code.
   - Add sidebar controls:
     ```python
     st.sidebar.markdown("### NFI Visualization")
     show_nfi = st.sidebar.checkbox("Show NFI Heatmap", value=True)
     if show_nfi:
         nfi_resolution = st.sidebar.slider(...)
         marker_size = st.sidebar.slider(...)
     ```
   - Compute grid when enabled, display max value:
     ```python
     with st.spinner("Calculating NFI values..."):
         nfi_df = calculate_nfi_grid(...)
     max_nfi_point = nfi_df.loc[nfi_df['nfi'].idxmax()]
     st.sidebar.metric("Max NFI", f"{max_nfi:.2e} N/kg")
     st.sidebar.write(f"Location: {max_nfi_lat:.1f}°, {max_nfi_lon:.1f}°")
     ```
   - Add heatmap scattergeo trace:
     ```python
     fig.add_trace(go.Scattergeo(
         lon=nfi_df['lon'],
         lat=nfi_df['lat'],
         mode='markers',
         marker=dict(
             color=nfi_df['nfi'],
             colorscale='Viridis',
             size=marker_size,
             opacity=0.7,
             colorbar=dict(title='NFI (N/kg)')
         ),
         name='NFI Heatmap',
         hoverinfo='text',
         text=[...]
     ))
     ```

```mermaid
graph TD
    A[compute_nfi in astro_utils.py] --> B[calculate_nfi_grid in planets_appV3.py]
    B --> C[Remove ESI, add NFI UI controls]
    C --> D[Compute NFI grid, display max, add heatmap trace]
    D --> E[Test and refine]