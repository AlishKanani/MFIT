"""
Shared colorscale definitions and utilities for Plotly visualizations.
"""

# Standard hot_r colorscale (white=cool to dark red=hot)
HOT_R_COLORSCALE = [
    [0.0, 'rgb(255, 255, 255)'],  # white (cool)
    [0.2, 'rgb(255, 255, 0)'],     # yellow
    [0.4, 'rgb(255, 200, 0)'],     # orange-yellow
    [0.6, 'rgb(255, 100, 0)'],     # orange
    [0.8, 'rgb(255, 0, 0)'],       # red
    [1.0, 'rgb(128, 0, 0)']        # dark red (hot)
]

# Use Plotly's built-in Viridis for flux magnitude
VIRIDIS_COLORSCALE = 'Viridis'


def temp_to_color(temp, temp_min, temp_max, colorscale=None):
    """
    Convert temperature to RGB color using specified colorscale.
    
    Args:
        temp: Temperature value
        temp_min: Minimum temperature for normalization
        temp_max: Maximum temperature for normalization
        colorscale: List of [position, 'rgb(r,g,b)'] pairs (default: HOT_R_COLORSCALE)
    
    Returns:
        RGB color string like 'rgb(255,128,0)'
    """
    if colorscale is None:
        colorscale = HOT_R_COLORSCALE
    
    # Normalize temperature to [0, 1]
    norm = (temp - temp_min) / (temp_max - temp_min) if temp_max > temp_min else 0.0
    norm = max(0, min(1, norm))
    
    # Find the two colors to interpolate between
    for i in range(len(colorscale) - 1):
        if norm <= colorscale[i+1][0]:
            t1, c1 = colorscale[i]
            t2, c2 = colorscale[i+1]
            # Linear interpolation
            alpha = (norm - t1) / (t2 - t1) if t2 > t1 else 0.0
            
            # Parse RGB values
            rgb1 = [int(x) for x in c1.replace('rgb(', '').replace(')', '').split(',')]
            rgb2 = [int(x) for x in c2.replace('rgb(', '').replace(')', '').split(',')]
            
            # Interpolate
            r = int(rgb1[0] + alpha * (rgb2[0] - rgb1[0]))
            g = int(rgb1[1] + alpha * (rgb2[1] - rgb1[1]))
            b = int(rgb1[2] + alpha * (rgb2[2] - rgb1[2]))
            
            return f'rgb({r},{g},{b})'
    
    return colorscale[-1][1]


def create_colorbar_dict(title, position='right', x_pos=1.02, length=0.8):
    """
    Generate standardized colorbar configuration for Plotly traces.
    
    Args:
        title: Colorbar title text
        position: Side for title ('right' or 'left')
        x_pos: X position of colorbar (default 1.02 for right side)
        length: Fractional length of colorbar (default 0.8)
    
    Returns:
        Dictionary with colorbar configuration
    """
    return dict(
        title=dict(
            text=title,
            side=position
        ),
        x=x_pos,
        len=length,
        thickness=20,
        yanchor="middle",
        y=0.5
    )


def matplotlib_to_plotly_colorscale(mpl_colorscale_name='hot_r', n_samples=256):
    """
    Convert matplotlib colorscale to Plotly format.
    
    Args:
        mpl_colorscale_name: Name of matplotlib colormap
        n_samples: Number of samples to use
    
    Returns:
        List of [position, 'rgb(r,g,b)'] pairs
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return HOT_R_COLORSCALE
    
    cmap = plt.cm.get_cmap(mpl_colorscale_name)
    positions = np.linspace(0, 1, n_samples)
    colorscale = []
    
    for pos in positions:
        rgba = cmap(pos)
        r, g, b = int(rgba[0] * 255), int(rgba[1] * 255), int(rgba[2] * 255)
        colorscale.append([float(pos), f'rgb({r},{g},{b})'])
    
    return colorscale

