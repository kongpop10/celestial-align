import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import datetime
import ephem
import pytz
from math import degrees, radians, sin, cos, acos
from astro_utils import get_planet_position, calculate_planet_path, compute_esi, calculate_ceres_position, compute_lasi


st.set_page_config(
    page_title="Celestial Alignment Analyzer",
    page_icon="🪐",
    layout="wide",
    initial_sidebar_state="expanded"  # Start with sidebar expanded
)

# Title and description
st.title("🪐 Celestial Alignment Analyzer")
st.markdown("""
Track planetary positions, analyze significant celestial alignments, and visualize the Event Significance Index (ESI) over time.
Select multiple planets to display, adjust the date and time, and discover when planetary configurations reach their peak significance.
""")

# Function to calculate Moon's position for a given time

# Function to calculate positions for the Moon's path over 24 hours

# Function to calculate LASI values across a global grid
def calculate_lasi_grid(dt, selected_planets, grid_resolution=10, exclude_pole_degrees=10):
    """Calculate LASI values across a global grid.

    Args:
        dt: Datetime for the calculation
        selected_planets: List of planet names to include in LASI calculation
        grid_resolution: Number of degrees between grid points (lower = more detailed but slower)
        exclude_pole_degrees: Degrees from poles to exclude (90-exclude_pole_degrees to 90)

    Returns:
        DataFrame with lat, lon, and lasi values
    """
    lasi_data = []

    # Create a grid of points across the globe, excluding poles
    lat_range = np.arange(-90 + exclude_pole_degrees, 90 - exclude_pole_degrees + 1, grid_resolution)
    lon_range = np.arange(-180, 181, grid_resolution)

    # Calculate LASI for each point
    for lat in lat_range:
        for lon in lon_range:
            # Calculate LASI for this location
            lasi_value = compute_lasi(lat, lon, dt, selected_planets)

            lasi_data.append({
                'lat': lat,
                'lon': lon,
                'lasi': lasi_value
            })

    # Convert to DataFrame
    lasi_df = pd.DataFrame(lasi_data)

    return lasi_df

# Function to calculate day/night terminator
def calculate_terminator(dt):
    # Get the Sun's position
    observer = ephem.Observer()
    observer.date = dt
    sun = ephem.Sun(observer)
    sun.compute(dt)

    # Get Sun's declination (latitude)
    sun_dec = degrees(float(sun.dec))

    # Convert right ascension to longitude
    ra_deg = degrees(float(sun.ra))
    gst = degrees(observer.sidereal_time())
    sun_ra = (ra_deg - gst) % 360
    if sun_ra > 180:
        sun_ra -= 360

    # Create terminator points
    terminator_lats = []
    terminator_lons = []

    # Create a more efficient grid for day/night visualization
    day_lats = []
    day_lons = []
    night_lats = []
    night_lons = []

    # Convert sun position to radians for calculations
    sun_lat_rad = radians(sun_dec)
    sun_lon_rad = radians(sun_ra)

    # Create a grid of points (reduced density for performance)
    lat_step = 5
    lon_step = 5

    for lat in np.arange(-90, 91, lat_step):
        for lon in np.arange(-180, 181, lon_step):
            lat_rad = radians(lat)
            lon_rad = radians(lon)

            # Calculate the angle between the sun and the point
            angle = acos(sin(sun_lat_rad) * sin(lat_rad) +
                         cos(sun_lat_rad) * cos(lat_rad) * cos(lon_rad - sun_lon_rad))

            # Determine if the point is in day or night
            # The terminator is approximately at 90 degrees from the sun
            if abs(degrees(angle) - 90) < 2:  # Points near the terminator
                terminator_lats.append(lat)
                terminator_lons.append(lon)
            elif degrees(angle) < 90:  # Day side
                day_lats.append(lat)
                day_lons.append(lon)
            else:  # Night side
                night_lats.append(lat)
                night_lons.append(lon)

    return {
        'terminator': (terminator_lats, terminator_lons),
        'day': (day_lats, day_lons),
        'night': (night_lats, night_lons)
    }

# Sidebar for date and time selection
st.sidebar.header("Date and Time Settings")

# Get current UTC time as default
initial_default_time = datetime.datetime.now(pytz.UTC)

# Initialize session state for date and time if they don't exist
if 'selected_date' not in st.session_state:
    st.session_state.selected_date = initial_default_time.date()
if 'selected_time' not in st.session_state:
    st.session_state.selected_time = initial_default_time.time()

# Define time navigation callback functions
def fast_backward():
    # Create a datetime object from the current date and time
    current_dt = datetime.datetime.combine(
        st.session_state.selected_date,
        st.session_state.selected_time,
        tzinfo=pytz.UTC
    )
    # Calculate new datetime (1 hour back)
    new_dt = current_dt - datetime.timedelta(hours=1)
    # Update session state for the next render
    st.session_state.selected_date = new_dt.date()
    st.session_state.selected_time = new_dt.time()

def backward():
    # Create a datetime object from the current date and time
    current_dt = datetime.datetime.combine(
        st.session_state.selected_date,
        st.session_state.selected_time,
        tzinfo=pytz.UTC
    )
    # Calculate new datetime (15 minutes back)
    new_dt = current_dt - datetime.timedelta(minutes=15)
    # Update session state for the next render
    st.session_state.selected_date = new_dt.date()
    st.session_state.selected_time = new_dt.time()

def forward():
    # Create a datetime object from the current date and time
    current_dt = datetime.datetime.combine(
        st.session_state.selected_date,
        st.session_state.selected_time,
        tzinfo=pytz.UTC
    )
    # Calculate new datetime (15 minutes forward)
    new_dt = current_dt + datetime.timedelta(minutes=15)
    # Update session state for the next render
    st.session_state.selected_date = new_dt.date()
    st.session_state.selected_time = new_dt.time()

def fast_forward():
    # Create a datetime object from the current date and time
    current_dt = datetime.datetime.combine(
        st.session_state.selected_date,
        st.session_state.selected_time,
        tzinfo=pytz.UTC
    )
    # Calculate new datetime (1 hour forward)
    new_dt = current_dt + datetime.timedelta(hours=1)
    # Update session state for the next render
    st.session_state.selected_date = new_dt.date()
    st.session_state.selected_time = new_dt.time()

# Date picker
selected_date = st.sidebar.date_input(
    "Select Date",
    key='selected_date' # Use key to automatically update session state
)

# Time picker
selected_time = st.sidebar.time_input(
    "Select Time (UTC)",
    key='selected_time' # Use key to automatically update session state
)

# Combine date and time
selected_datetime = datetime.datetime.combine(
    st.session_state.selected_date, # Use session state values
    st.session_state.selected_time, # Use session state values
    tzinfo=pytz.UTC
)

# Display selected date and time
st.sidebar.write(f"Selected: {selected_datetime.strftime('%Y-%m-%d %H:%M UTC')}")

# Add custom CSS for button styling
st.markdown("""
<style>
    /* Custom button style for navigation buttons */
    section[data-testid="stSidebar"] div[data-testid="stHorizontalBlock"] button {
        background-color: #00A6ED !important;
        color: white !important;
        border-color: #00A6ED !important;
    }
</style>
""", unsafe_allow_html=True)

# Time navigation buttons
st.sidebar.markdown("### Time Navigation")
# Create a container for the buttons
button_container = st.sidebar.container()

# Create columns for the buttons
col1, col2, col3, col4 = button_container.columns(4)

# Fast backward button (-1 hour)
col1.button("⏪", help="Go back 1 hour", key="fast_backward", on_click=fast_backward)

# Backward button (-15 minutes)
col2.button("◀️", help="Go back 15 minutes", key="backward", on_click=backward)

# Forward button (+15 minutes)
col3.button("▶️", help="Go forward 15 minutes", key="forward", on_click=forward)

# Fast forward button (+1 hour)
col4.button("⏩", help="Go forward 1 hour", key="fast_forward", on_click=fast_forward)

# Planet selection
st.sidebar.markdown("### Planet Selection")
planet_options = ["Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto", "Ceres"]

# Initialize session state for selected planets if it doesn't exist
if 'selected_planets' not in st.session_state:
    st.session_state.selected_planets = ["Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn"]

# Multi-select for planets
selected_planets = st.sidebar.multiselect(
    "Select Planets to Display",
    options=planet_options,
    default=st.session_state.selected_planets,
    key="selected_planets"
)

# Ensure at least one planet is selected
if not selected_planets:
    st.sidebar.warning("Please select at least one planet to display.")
    selected_planets = ["Moon"]  # Default to Moon if nothing selected

# Option to select which planet to center the map on
if len(selected_planets) > 1:
    center_planet = st.sidebar.selectbox(
        "Center Map On",
        options=selected_planets,
        index=0
    )
else:
    center_planet = selected_planets[0]  # If only one planet is selected, center on it

# Calculate positions and paths for all selected planets
planet_data = {}
for planet in selected_planets:
    current_lat, current_lon = get_planet_position(selected_datetime, planet)
    path_df = calculate_planet_path(selected_datetime, planet)
    planet_data[planet] = {
        'current_lat': current_lat,
        'current_lon': current_lon,
        'path_df': path_df
    }

# Extend planet_data with geocentric distance (AU)
for planet in selected_planets:
    if planet == "Ceres":
        # For Ceres, use our custom function
        _, _, earth_distance = calculate_ceres_position(selected_datetime)
        planet_data[planet]['distance_au'] = earth_distance
    else:
        # For other planets, use ephem
        body = getattr(ephem, planet)()
        body.compute(selected_datetime)
        planet_data[planet]['distance_au'] = body.earth_distance

# Compute Event Significance Index (ESI)
esi_value = compute_esi({
    p: {
        'lat': planet_data[p]['current_lat'],
        'lon': planet_data[p]['current_lon'],
        'distance_au': planet_data[p]['distance_au']
    } for p in selected_planets
})
# Display ESI in the sidebar
st.sidebar.metric("Event Significance Index (ESI)", f"{esi_value:.2e}")

# LASI Visualization Section
st.sidebar.markdown("### LASI Visualization")
show_lasi = st.sidebar.checkbox("Show LASI Heatmap", value=False,
                              help="Display Local Alignment Significance Index (LASI) as a heatmap on the globe")

# Only show these options if LASI is enabled
if show_lasi:
    st.sidebar.markdown("#### LASI Map Settings")
    lasi_resolution = st.sidebar.slider("Grid Resolution", 5, 30, 15, 5,
                                     help="Lower values give more detailed but slower heatmap (degrees between points)")
    exclude_poles = st.sidebar.slider("Exclude Pole Regions (degrees)", 5, 30, 10, 5,
                                   help="Exclude regions within this many degrees of the poles")
    marker_size = st.sidebar.slider("Marker Size", 4, 15, 8, 1,
                                 help="Size of the dots on the LASI heatmap")

    # Calculate LASI grid if enabled
    with st.sidebar.status("Calculating LASI values...", expanded=False) as status:
        lasi_df = calculate_lasi_grid(selected_datetime, selected_planets,
                                    grid_resolution=lasi_resolution,
                                    exclude_pole_degrees=exclude_poles)
        # Find the point with maximum LASI value
        max_lasi_point = lasi_df.loc[lasi_df['lasi'].idxmax()]
        status.update(label="LASI calculation complete!", state="complete", expanded=False)

    # Display the maximum LASI value and its location
    st.sidebar.metric("Maximum LASI Value", f"{max_lasi_point['lasi']:.2e}")
    st.sidebar.write(f"Max LASI Location: {max_lasi_point['lat']:.1f}°, {max_lasi_point['lon']:.1f}°")


# Calculate day/night terminator
terminator_data = calculate_terminator(selected_datetime)

# Create a two-column layout for large screens
left_col, right_col = st.columns([2.2, 2.8])  # 2.2:2.8 ratio for a more balanced layout

# Create the map visualization in the left column
with left_col:
    st.subheader(f"Planetary Paths and Day/Night Visualization")

    # Create a figure without using a container to maximize space
    fig = go.Figure()

    # Add night side as a scatter plot with dark blue color
    fig.add_trace(go.Scattergeo(
        lon=terminator_data['night'][1],
        lat=terminator_data['night'][0],
        mode='markers',
        marker=dict(
            size=4,
            color='rgba(25, 25, 112, 0.7)',  # Dark blue with transparency
            opacity=0.7
        ),
        name='Night',
        showlegend=True
    ))

    # Add day side as a scatter plot with light blue color
    fig.add_trace(go.Scattergeo(
        lon=terminator_data['day'][1],
        lat=terminator_data['day'][0],
        mode='markers',
        marker=dict(
            size=4,
            color='rgba(135, 206, 250, 0.7)',  # Light blue with transparency
            opacity=0.7
        ),
        name='Day',
        showlegend=True
    ))

    # Add terminator line
    fig.add_trace(go.Scattergeo(
        lon=terminator_data['terminator'][1],
        lat=terminator_data['terminator'][0],
        mode='markers',
        marker=dict(
            size=3,
            color='rgba(255, 165, 0, 0.8)',  # Orange with transparency
            opacity=0.8
        ),
        name='Terminator',
        showlegend=True
    ))

    # Define planet colors and symbols
    planet_colors = {
        "Sun": "#FFD700",  # Bright gold for the Sun
        "Moon": "#E6E6E6",  # Light gray for the Moon
        "Mercury": "#A9A9A9",  # Gray for Mercury
        "Venus": "#FFA500",  # Orange for Venus
        "Mars": "#FF4500",  # Red-orange for Mars
        "Jupiter": "#F0E68C",  # Khaki for Jupiter
        "Saturn": "#DAA520",   # Golden rod for Saturn
        "Uranus": "#00FFFF",  # Cyan for Uranus
        "Neptune": "#0000FF",  # Blue for Neptune
        "Pluto": "#8B4513",   # SaddleBrown for Pluto
        "Ceres": "#C0C0C0"    # Silver for Ceres (dwarf planet)
    }

    # Define more representative symbols for planets
    planet_symbols = {
        "Sun": "circle",       # Circle for the Sun
        "Moon": "circle",      # Circle for the Moon
        "Mercury": "circle",   # Circle for Mercury
        "Venus": "circle",     # Circle for Venus
        "Mars": "circle",      # Circle for Mars
        "Jupiter": "circle",   # Circle for Jupiter
        "Saturn": "circle",    # Circle for Saturn
        "Uranus": "circle",    # Circle for Uranus
        "Neptune": "circle",   # Circle for Neptune
        "Pluto": "circle",     # Circle for Pluto
        "Ceres": "circle"      # Circle for Ceres
    }

    # Define planet sizes to better represent their relative sizes
    planet_sizes = {
        "Sun": 40,        # Largest
        "Moon": 15,       # Small
        "Mercury": 12,    # Very small
        "Venus": 18,      # Medium
        "Mars": 14,       # Small-medium
        "Jupiter": 30,    # Large
        "Saturn": 28,     # Large
        "Uranus": 24,     # Medium-large
        "Neptune": 23,    # Medium-large
        "Pluto": 8,       # Smallest
        "Ceres": 10       # Dwarf planet, slightly larger than Pluto
    }

    # Add paths and current positions for all selected planets
    for planet in selected_planets:
        # Get planet data
        planet_path_df = planet_data[planet]['path_df']
        current_planet_lat = planet_data[planet]['current_lat']
        current_planet_lon = planet_data[planet]['current_lon']

        # Add planet path
        fig.add_trace(go.Scattergeo(
            lon=planet_path_df['lon'],
            lat=planet_path_df['lat'],
            mode='lines+markers',
            line=dict(
                width=2,
                color=planet_colors[planet]
            ),
            marker=dict(
                size=4,
                color=planet_colors[planet],
                symbol=planet_symbols[planet]
            ),
            text=planet_path_df['time'],
            hoverinfo='text+lon+lat',
            name=f'{planet} Path'
        ))

        # Add current planet position with a marker
        fig.add_trace(go.Scattergeo(
            lon=[current_planet_lon],
            lat=[current_planet_lat],
            mode='markers',
            marker=dict(
                size=planet_sizes[planet],
                color=planet_colors[planet],
                line=dict(width=2, color='black'),
                symbol=planet_symbols[planet]
            ),
            name=f'Current {planet} Position',
            text=[f"{planet} at {selected_datetime.strftime('%Y-%m-%d %H:%M UTC')}"],
            hoverinfo='text+lon+lat'
        ))


    # Add only the maximum LASI point marker to the main map if LASI is enabled
    if show_lasi:
        # Add marker for maximum LASI point
        fig.add_trace(go.Scattergeo(
            lon=[max_lasi_point['lon']],
            lat=[max_lasi_point['lat']],
            mode='markers',
            marker=dict(
                size=15,  # Slightly larger for better visibility
                color='red',
                symbol='star',
                line=dict(width=2, color='white')
            ),
            name='Maximum LASI',
            text=[f"Maximum LASI: {max_lasi_point['lasi']:.2e}"],
            hoverinfo='text+lon+lat'
        ))

    # Configure the map to maximize the globe view
    # Center the map on the selected center planet
    fig.update_geos(
        projection_type="orthographic",
        projection_rotation=dict(lon=planet_data[center_planet]['current_lon'], lat=planet_data[center_planet]['current_lat'], roll=0), # Rotate view to center on selected planet
        landcolor='rgb(217, 217, 217)',
        oceancolor='rgb(55, 97, 164)',
        showland=True,
        showocean=True,
        showcoastlines=True,
        showcountries=True,
        countrycolor='white',
        coastlinecolor='white',
        coastlinewidth=0.5,
        countrywidth=0.5,
        bgcolor='rgba(0,0,0,0)',
        # projection_scale=1.1, # Removed to allow auto-scaling
        fitbounds=False,  # Don't auto-fit bounds
        visible=True  # Ensure the map is visible
    )

    # Update layout to maximize space and remove blank areas
    fig.update_layout(
        title="",
        height=600,  # Reduced height to make it smaller than the ESI section
        autosize=True,  # Allow auto-sizing based on container
        margin=dict(l=0, r=0, t=30, b=0, pad=0),  # Remove all margins and padding
        legend=dict(
            x=0,
            y=1,
            bgcolor='rgba(255, 255, 255, 0.7)',
            orientation='h'  # Horizontal legend to save space
        ),
        paper_bgcolor='rgba(0,0,0,0)',  # Transparent background
        plot_bgcolor='rgba(0,0,0,0)'  # Transparent plot area
    )

    # Display the map with full width in the left column
    st.plotly_chart(fig, use_container_width=True, config={
        'responsive': True,
        'displayModeBar': False,
        'scrollZoom': True,  # Enable scroll zoom for better interaction
        'modeBarButtonsToRemove': ['select2d', 'lasso2d']  # Remove unnecessary buttons
    })

    # Add brief explanation about the red star if LASI is enabled
    if show_lasi:
        st.markdown(f"""
        **Note:** The **red star** marks the location with the maximum LASI value ({max_lasi_point['lat']:.1f}°, {max_lasi_point['lon']:.1f}°).
        See the dedicated LASI map below for more details.
        """)

# Add a dedicated LASI heatmap section if enabled
if show_lasi:
    st.subheader("Local Alignment Significance Index (LASI) Heatmap")

    # Create a dedicated figure for LASI visualization
    lasi_fig = go.Figure()

    # Add LASI heatmap with larger markers
    lasi_fig.add_trace(go.Scattergeo(
        lon=lasi_df['lon'],
        lat=lasi_df['lat'],
        mode='markers',
        marker=dict(
            size=marker_size,  # Use the user-selected marker size
            color=lasi_df['lasi'],
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(
                title='LASI Value',
                thickness=15,
                len=0.5
            ),
            opacity=0.8,  # Slightly more opaque
            cmin=lasi_df['lasi'].min(),
            cmax=lasi_df['lasi'].max(),
        ),
        name='LASI Values',
        text=[f"LASI: {lasi:.2e}" for lasi in lasi_df['lasi']],
        hoverinfo='text+lon+lat'
    ))

    # Add marker for maximum LASI point
    lasi_fig.add_trace(go.Scattergeo(
        lon=[max_lasi_point['lon']],
        lat=[max_lasi_point['lat']],
        mode='markers+text',
        marker=dict(
            size=20,  # Even larger for emphasis
            color='red',
            symbol='star',
            line=dict(width=2, color='white')
        ),
        name='Maximum LASI',
        text=["Maximum LASI"],
        textposition="top center",
        hovertext=[f"Maximum LASI: {max_lasi_point['lasi']:.2e}"],
        hoverinfo='text+lon+lat'
    ))

    # Add coastlines and country borders for reference
    lasi_fig.update_geos(
        projection_type="natural earth",  # Use a flat map projection for better global view
        showland=True,
        landcolor='rgb(217, 217, 217)',
        oceancolor='rgb(55, 97, 164)',
        showocean=True,
        showcoastlines=True,
        showcountries=True,
        countrycolor='white',
        coastlinecolor='white',
        coastlinewidth=0.5,
        countrywidth=0.5,
        bgcolor='rgba(0,0,0,0)',
        lataxis_showgrid=True,
        lonaxis_showgrid=True
    )

    # Update layout
    lasi_fig.update_layout(
        title=f"LASI Values at {selected_datetime.strftime('%Y-%m-%d %H:%M UTC')}",
        height=600,
        margin=dict(l=0, r=0, t=50, b=0),
        legend=dict(
            x=0,
            y=1,
            bgcolor='rgba(255, 255, 255, 0.7)',
            orientation='h'
        )
    )

    # Display the LASI map
    st.plotly_chart(lasi_fig, use_container_width=True, config={
        'responsive': True,
        'displayModeBar': True,  # Show mode bar for this map to allow zooming
        'scrollZoom': True
    })

    # Add detailed explanation for LASI heatmap
    st.markdown(f"""
    **About the LASI Heatmap:**

    The colored dots represent the Local Alignment Significance Index (LASI) calculated for each location on Earth at {selected_datetime.strftime('%Y-%m-%d %H:%M UTC')}.

    - **Higher values (yellow/green)** indicate locations where the selected celestial bodies appear more closely aligned in the sky
    - **Lower values (purple/blue)** indicate locations where the alignments are less significant
    - The **red star** marks the location with the maximum LASI value ({max_lasi_point['lat']:.1f}°, {max_lasi_point['lon']:.1f}°)
    - Pole regions (within {exclude_poles}° of the poles) are excluded as requested

    LASI considers both the angular separation between celestial bodies as seen from each location and the mass-distance relationship of each body.

    **How to use this map:**
    - Hover over points to see exact LASI values
    - Use the toolbar in the top-right to zoom, pan, or download the visualization
    - Try different dates/times and planet selections to see how LASI values change globally
    """)

# Function to calculate ESI values for each hour of the selected date
def calculate_esi_over_time(base_date, selected_planets):
    esi_data = []
    # Start at midnight and calculate for each hour of the day
    for hour in range(24):
        # Create datetime for this hour
        dt = datetime.datetime.combine(base_date, datetime.time(hour=hour), tzinfo=pytz.UTC)

        # Calculate positions for all selected planets at this time
        hour_planet_data = {}
        for planet in selected_planets:
            lat, lon = get_planet_position(dt, planet)

            # Handle Ceres separately
            if planet == "Ceres":
                _, _, earth_distance = calculate_ceres_position(dt)
                hour_planet_data[planet] = {
                    'lat': lat,
                    'lon': lon,
                    'distance_au': earth_distance
                }
            else:
                body = getattr(ephem, planet)()
                body.compute(dt)
                hour_planet_data[planet] = {
                    'lat': lat,
                    'lon': lon,
                    'distance_au': body.earth_distance
                }

        # Calculate ESI for this hour
        esi_value = compute_esi(hour_planet_data)

        esi_data.append({
            'hour': hour,
            'time': dt.strftime('%H:%M'),
            'datetime': dt,
            'esi': esi_value
        })

    return pd.DataFrame(esi_data)

# Function to calculate ESI values for a custom date range
def calculate_esi_over_range(base_date, days_before, days_after, selected_planets):
    esi_data = []
    # Calculate start and end dates based on the selected range
    start_date = base_date - datetime.timedelta(days=abs(days_before))
    end_date = base_date + datetime.timedelta(days=days_after)

    # Calculate for each day in the range
    current_day = start_date
    day_index = 0
    while current_day <= end_date:
        # Create datetime for noon on this day (representative time)
        dt = datetime.datetime.combine(current_day, datetime.time(hour=12), tzinfo=pytz.UTC)

        # Calculate positions for all selected planets at this time
        day_planet_data = {}
        for planet in selected_planets:
            lat, lon = get_planet_position(dt, planet)

            # Handle Ceres separately
            if planet == "Ceres":
                _, _, earth_distance = calculate_ceres_position(dt)
                day_planet_data[planet] = {
                    'lat': lat,
                    'lon': lon,
                    'distance_au': earth_distance
                }
            else:
                body = getattr(ephem, planet)()
                body.compute(dt)
                day_planet_data[planet] = {
                    'lat': lat,
                    'lon': lon,
                    'distance_au': body.earth_distance
                }

        # Calculate ESI for this day
        esi_value = compute_esi(day_planet_data)

        esi_data.append({
            'day': day_index,
            'date': current_day,
            'display_date': current_day.strftime('%d'),
            'full_date': current_day.strftime('%Y-%m-%d'),
            'datetime': dt,
            'esi': esi_value
        })

        # Move to next day
        current_day += datetime.timedelta(days=1)
        day_index += 1

    return pd.DataFrame(esi_data)


# Calculate ESI values over time for the selected date
# Display in the right column
with right_col:
    st.subheader("Event Significance Index (ESI) Over Time")
    esi_df = calculate_esi_over_time(selected_date, selected_planets)

    # Find min and max ESI values
    min_esi = esi_df['esi'].min()
    max_esi = esi_df['esi'].max()
    min_esi_time = esi_df.loc[esi_df['esi'] == min_esi, 'time'].iloc[0]
    max_esi_time = esi_df.loc[esi_df['esi'] == max_esi, 'time'].iloc[0]

    # Get the selected hour for the vertical line
    selected_hour = selected_datetime.hour + selected_datetime.minute/60

    # Create ESI chart
    esi_fig = go.Figure()

    # Add ESI line
    esi_fig.add_trace(go.Scatter(
        x=esi_df['hour'],
        y=esi_df['esi'],
        mode='lines',
        name='ESI',
        line=dict(color='#1f77b4', width=2)
    ))

    # Add minimum ESI marker
    esi_fig.add_trace(go.Scatter(
        x=[esi_df.loc[esi_df['esi'] == min_esi, 'hour'].iloc[0]],
        y=[min_esi],
        mode='markers+text',
        marker=dict(color='blue', size=10, symbol='circle'),
        text=[f'Min: {min_esi:.2e}'],
        textposition='top center',
        name='Minimum ESI'
    ))

    # Add maximum ESI marker
    esi_fig.add_trace(go.Scatter(
        x=[esi_df.loc[esi_df['esi'] == max_esi, 'hour'].iloc[0]],
        y=[max_esi],
        mode='markers+text',
        marker=dict(color='red', size=10, symbol='circle'),
        text=[f'Max: {max_esi:.2e}'],
        textposition='top center',
        name='Maximum ESI'
    ))

    # Add vertical line for selected time
    esi_fig.add_vline(
        x=selected_hour,
        line_width=2,
        line_dash='dash',
        line_color='green',
        annotation_text=f'Selected Time: {selected_datetime.strftime("%H:%M")}',
        annotation_position='top right'
    )

    # Update layout
    esi_fig.update_layout(
        title='',
        xaxis_title='Hour of Day (UTC)',
        yaxis_title='Event Significance Index (ESI)',
        xaxis=dict(
            tickmode='array',
            tickvals=list(range(0, 24, 2)),
            ticktext=[f'{h:02d}:00' for h in range(0, 24, 2)]
        ),
        yaxis=dict(type='log'),  # Log scale for better visualization of ESI values
        height=650,  # Slightly larger than the map visualization
        margin=dict(l=0, r=0, t=50, b=0),
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='right', x=1)
    )

    # Display the ESI chart in the right column
    st.plotly_chart(esi_fig, use_container_width=True, config={
        'responsive': True,
        'displayModeBar': False
    })

    # No explanation here - moved to consolidated section below

# Add the Date Range ESI Variation section below the daily chart
st.subheader("ESI Variation in the Selected Date Range")

# Add sliders for adjusting the date range
col1, col2 = st.columns(2)
with col1:
    days_before = -st.slider("Days before selected date", 0, 60, 30, 1,
                          help="Select how many days before the selected date to include")
with col2:
    days_after = st.slider("Days after selected date", 0, 60, 15, 1,
                         help="Select how many days after the selected date to include")

# Calculate ESI values for the custom date range
range_esi_df = calculate_esi_over_range(selected_date, days_before, days_after, selected_planets)

# Find min and max ESI values for the date range
range_min_esi = range_esi_df['esi'].min()
range_max_esi = range_esi_df['esi'].max()
range_min_esi_date = range_esi_df.loc[range_esi_df['esi'] == range_min_esi, 'full_date'].iloc[0]
range_max_esi_date = range_esi_df.loc[range_esi_df['esi'] == range_max_esi, 'full_date'].iloc[0]

# Get the selected day for the vertical line
# Check if the selected date is in the range
if selected_date in range_esi_df['date'].values:
    selected_day_index = range_esi_df[range_esi_df['date'] == selected_date]['day'].iloc[0]
    show_selected_line = True
else:
    # If selected date is outside the range, don't show the line
    selected_day_index = 0
    show_selected_line = False

# Create date range ESI chart
range_esi_fig = go.Figure()

# Add ESI line for the date range
range_esi_fig.add_trace(go.Scatter(
    x=range_esi_df['day'],
    y=range_esi_df['esi'],
    mode='lines',
    name='ESI',
    line=dict(color='#1f77b4', width=2)
))

# Add minimum ESI marker
range_esi_fig.add_trace(go.Scatter(
    x=[range_esi_df.loc[range_esi_df['esi'] == range_min_esi, 'day'].iloc[0]],
    y=[range_min_esi],
    mode='markers+text',
    marker=dict(color='blue', size=10, symbol='circle'),
    text=[f'Min: {range_min_esi:.2e}'],
    textposition='top center',
    name='Minimum ESI'
))

# Add maximum ESI marker
range_esi_fig.add_trace(go.Scatter(
    x=[range_esi_df.loc[range_esi_df['esi'] == range_max_esi, 'day'].iloc[0]],
    y=[range_max_esi],
    mode='markers+text',
    marker=dict(color='red', size=10, symbol='circle'),
    text=[f'Max: {range_max_esi:.2e}'],
    textposition='top center',
    name='Maximum ESI'
))

# Add vertical line for selected date only if it's in the range
if show_selected_line:
    range_esi_fig.add_vline(
        x=selected_day_index,
        line_width=2,
        line_dash='dash',
        line_color='green',
        annotation_text=f'Selected Date: {selected_date.strftime("%Y-%m-%d")}',
        annotation_position='top right'
    )

# Calculate the date range for the title
start_date = selected_date + datetime.timedelta(days=days_before)
end_date = selected_date + datetime.timedelta(days=days_after)
date_range_str = f'{start_date.strftime("%Y-%m-%d")} to {end_date.strftime("%Y-%m-%d")}'

# Update layout for date range chart
range_esi_fig.update_layout(
    title=f'ESI Variation from {date_range_str}',
    xaxis_title='Day in Range',
    yaxis_title='Event Significance Index (ESI)',
    xaxis=dict(
        tickmode='array',
        tickvals=list(range(0, len(range_esi_df), max(1, len(range_esi_df) // 10))),  # Adjust tick frequency based on range length
        ticktext=[range_esi_df.iloc[i]['full_date'] if i < len(range_esi_df) else '' for i in range(0, len(range_esi_df), max(1, len(range_esi_df) // 10))]
    ),
    yaxis=dict(type='log'),  # Log scale for better visualization of ESI values
    height=650,  # Same height as the daily chart for consistency
    margin=dict(l=0, r=0, t=50, b=0),
    legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='right', x=1)
)

# Display the date range ESI chart with full width
st.plotly_chart(range_esi_fig, use_container_width=True, config={
    'responsive': True,
    'displayModeBar': False
})

# Add consolidated explanation for both ESI charts
st.markdown("""
**About the ESI Charts:**

*Daily ESI Chart:*
- Shows how the Event Significance Index (ESI) changes throughout the selected date
- Each point represents the ESI value at that specific hour
- The green dashed line marks your selected time
- Blue and red markers show the minimum and maximum ESI values for the day

*Date Range ESI Chart:*
- Shows how the ESI changes throughout your selected date range (adjustable from 0-60 days before to 0-60 days after)
- Use the sliders to customize the date range you want to analyze
- Each point represents the ESI value at noon (UTC) for that day
- The green dashed line marks your selected date
- Blue and red markers show the minimum and maximum ESI values for the selected range

*General Information:*
- Higher ESI values indicate more significant planetary alignments
- The y-axis uses a logarithmic scale to better visualize the variations
- ESI considers the proximity of celestial bodies, their masses, and distances from Earth
""")

# How to Use and About section moved to sidebar

# Add a note about the calculations
st.sidebar.markdown("""
### Note on Calculations
This visualization provides an approximation of:
- The selected celestial bodies' positions (where they appear directly overhead)
- The day/night terminator
- Each body's path over 24 hours

The Sun's position indicates where it appears directly overhead (at the zenith) at the selected time.

For precise astronomical calculations, specialized tools should be used.
""")

# Add the How to Use and About section to the sidebar
st.sidebar.markdown("---") # Add a separator

# Create a collapsible section for How to Use
with st.sidebar.expander("How to Use", expanded=False):
    st.markdown("""
    - Select multiple celestial bodies (including the Sun, planets, and dwarf planets like Ceres and Pluto) to visualize simultaneously using the multi-select dropdown.
    - Use the date and time selectors in the sidebar to choose when to view the positions.
    - The colored lines show each body's path over 24 hours from the selected time.
    - The larger markers show each body's position at the selected time.
    - Blue areas are in daylight, while dark blue areas are in night.
    - The orange line represents the day/night terminator (dawn/dusk line).
    - You can choose which body to center the map on when multiple are selected.
    - The Sun's position shows where it appears directly overhead on Earth.
    - Enable the LASI heatmap to visualize the Local Alignment Significance Index across the globe.
    - The LASI heatmap shows where celestial alignments appear most significant from different Earth locations.
    - A red star marks the location with the maximum LASI value.
    - Adjust the grid resolution and pole exclusion settings to customize the LASI visualization.
    """)

# Create a collapsible section for About
with st.sidebar.expander("About", expanded=False):
    st.markdown("""
    This visualization uses astronomical calculations to determine:
    - The position of the Sun and planets relative to Earth (where they appear directly overhead)
    - The day/night terminator based on the Sun's position
    - The projected path of celestial bodies over a 24-hour period
    - The Event Significance Index (ESI) for global planetary alignments
    - The Local Alignment Significance Index (LASI) for location-specific alignments

    **LASI (Local Alignment Significance Index)** measures how significant a celestial alignment appears from a specific location on Earth. It considers:
    - The angular separation between celestial bodies as seen from a specific location
    - The mass and distance of each celestial body
    - The combined effect of all selected bodies

    Higher LASI values indicate locations where the selected celestial bodies appear more closely aligned in the sky.
    """)
