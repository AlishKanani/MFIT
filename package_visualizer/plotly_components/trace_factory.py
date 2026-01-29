"""
Factory functions for creating standardized Plotly trace objects.
"""

try:
    import plotly.graph_objects as go
except ImportError:
    go = None

from .colorscales import HOT_R_COLORSCALE, create_colorbar_dict
import numpy as np


def create_mesh_trace(vertices_x, vertices_y, vertices_z, intensities,
                      i_indices, j_indices, k_indices,
                      temp_min, temp_max, body_name,
                      colorscale='Hot_r', show_colorbar=True,
                      show_legend=True, legendgroup=None,
                      hover_texts=None, hover_template=None):
    """
    Create a Mesh3d trace for temperature visualization.
    
    Args:
        vertices_x, vertices_y, vertices_z: Vertex coordinate lists
        intensities: Temperature values for each vertex
        i_indices, j_indices, k_indices: Triangle face indices
        temp_min, temp_max: Temperature range for color normalization
        body_name: Name of the body for legend
        colorscale: Plotly colorscale name
        show_colorbar: Whether to show the colorbar
        hover_texts: Optional per-vertex hover text list
        hover_template: Optional hovertemplate override
    
    Returns:
        Plotly Mesh3d trace
    """
    if go is None:
        raise ImportError("plotly is required for mesh trace creation")
    text_values = hover_texts if hover_texts else [body_name] * len(vertices_x)
    default_hover = (
        '<b>%{text}</b><br>Temp: %{intensity:.2f}°C<extra></extra>'
        if hover_texts is None else
        '%{text}<extra></extra>'
    )
    
    return go.Mesh3d(
        x=vertices_x,
        y=vertices_y,
        z=vertices_z,
        i=i_indices,
        j=j_indices,
        k=k_indices,
        intensity=intensities,
        colorscale=colorscale,
        cmin=temp_min,
        cmax=temp_max,
        name=body_name,
        showscale=show_colorbar,
        colorbar=dict(
            title=dict(
                text="Temp (°C)",
                side="right"
            ),
            x=1.02,
            len=0.7,
            thickness=20,
            yanchor="top",
            y=0.95
        ) if show_colorbar else None,
        flatshading=True,
        lighting=dict(
            ambient=0.8,
            diffuse=0.8,
            specular=0.2,
            roughness=0.9
        ),
        hovertemplate=hover_template or default_hover,
        text=text_values,
        showlegend=show_legend,
        legendgroup=legendgroup or body_name
    )


def create_edge_trace(edge_x, edge_y, edge_z, body_name, 
                     color='lightgray', width=1):
    """
    Create a Scatter3d trace for wireframe edges.
    
    Args:
        edge_x, edge_y, edge_z: Edge coordinate lists (with None separators)
        body_name: Name of the body for legend grouping
        color: Line color
        width: Line width
    
    Returns:
        Plotly Scatter3d trace
    """
    if go is None:
        raise ImportError("plotly is required for edge trace creation")
    
    return go.Scatter3d(
        x=edge_x,
        y=edge_y,
        z=edge_z,
        mode='lines',
        line=dict(color=color, width=width),
        showlegend=False,
        hoverinfo='skip',
        legendgroup=body_name
    )


def create_cone_trace(x_centers, y_centers, z_centers,
                     u_scaled, v_scaled, w_scaled,
                     magnitudes, body_name,
                     colorscale='Viridis', visible=False):
    """
    Create a Cone trace for flux vector visualization.
    
    Args:
        x_centers, y_centers, z_centers: Cone anchor positions
        u_scaled, v_scaled, w_scaled: Scaled vector components
        magnitudes: Flux magnitudes for coloring
        body_name: Name of the body for legend
        colorscale: Plotly colorscale name
        visible: Initial visibility state
    
    Returns:
        Plotly Cone trace
    """
    if go is None:
        raise ImportError("plotly is required for cone trace creation")
    
    if not x_centers:
        # Return empty trace
        return go.Cone(
            x=[], y=[], z=[],
            u=[], v=[], w=[],
            colorscale=colorscale,
            sizemode="absolute",
            sizeref=1.0,
            showscale=False,
            name=f"{body_name} flux",
            visible=visible,
            anchor="tail",
            legendgroup=body_name,
        )
    
    return go.Cone(
        x=x_centers,
        y=y_centers,
        z=z_centers,
        u=u_scaled,
        v=v_scaled,
        w=w_scaled,
        colorscale=colorscale,
        cmin=float(np.min(magnitudes)),
        cmax=float(np.max(magnitudes)),
        sizemode="absolute",
        sizeref=1.0,
        showscale=False,
        name=f"{body_name} flux",
        visible=visible,
        anchor="tail",
        legendgroup=body_name,
    )


def create_scatter_trace(x, y, z, mode='markers', marker_size=8,
                        marker_color=None, colorscale='Viridis',
                        colorbar_title=None, hover_text=None,
                        name='', show_legend=False):
    """
    Create a Scatter3d trace for markers or lines.
    
    Args:
        x, y, z: Coordinate arrays
        mode: 'markers', 'lines', or 'markers+lines'
        marker_size: Size of markers
        marker_color: Color values (can be array for colormap)
        colorscale: Plotly colorscale name
        colorbar_title: Title for colorbar (if using color array)
        hover_text: List of hover text strings
        name: Trace name for legend
        show_legend: Whether to show in legend
    
    Returns:
        Plotly Scatter3d trace
    """
    if go is None:
        raise ImportError("plotly is required for scatter trace creation")
    
    marker_dict = dict(size=marker_size)
    
    if marker_color is not None:
        if isinstance(marker_color, (list, np.ndarray)):
            marker_dict['color'] = marker_color
            marker_dict['colorscale'] = colorscale
            if colorbar_title:
                marker_dict['colorbar'] = create_colorbar_dict(colorbar_title)
                marker_dict['showscale'] = True
        else:
            marker_dict['color'] = marker_color
    
    return go.Scatter3d(
        x=x,
        y=y,
        z=z,
        mode=mode,
        marker=marker_dict if 'markers' in mode else None,
        hovertext=hover_text if hover_text else [],
        hoverinfo='text' if hover_text else 'skip',
        name=name,
        showlegend=show_legend
    )

