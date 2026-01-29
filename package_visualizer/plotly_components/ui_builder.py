"""
UI component builders for Plotly interactive controls.
"""

import numpy as np


def create_axis_clip_slider(
    axis_name,
    clip_values,
    body_names,
    body_node_data,
    filter_func,
    mesh_builder_func,
    edge_builder_func,
    flux_builder_func,
    trace_indices,
    position=(0.02, -0.05),
    package_bounds=None
):
    """
    Create a generic axis clipping slider for X/Y/Z planes.
    
    This function eliminates 120+ lines × 3 axes of repetition by providing
    a unified slider builder that works for any axis.
    
    Args:
        axis_name: 'x', 'y', or 'z'
        clip_values: Array of clipping plane positions
        body_names: List of body names in order
        body_node_data: Dictionary of body_name -> node list
        filter_func: Function(nodes, clip_val, axis) -> filtered_nodes
        mesh_builder_func: Function(nodes) -> mesh data
        edge_builder_func: Function(nodes) -> edge data
        flux_builder_func: Function(nodes) -> flux data
        trace_indices: List of trace indices to update
        position: (x, y) position for slider
        package_bounds: Dict with 'x', 'y', 'z' tuples of (min, max) to fix axis ranges
    
    Returns:
        Plotly slider dictionary
    """
    steps = []
    
    for clip_val in clip_values:
        # Build separate lists for each property across all traces
        x_updates = []
        y_updates = []
        z_updates = []
        intensity_updates = []
        i_updates = []
        j_updates = []
        k_updates = []
        text_updates = []
        u_updates = []
        v_updates = []
        w_updates = []
        
        for body_name in body_names:
            filtered_nodes = filter_func(body_node_data[body_name], clip_val, axis_name)
            
            if filtered_nodes:
                # Rebuild mesh for filtered nodes
                vx, vy, vz, intensities, i_i, j_i, k_i, texts = mesh_builder_func(filtered_nodes)
                
                # Mesh trace data
                x_updates.append(vx)
                y_updates.append(vy)
                z_updates.append(vz)
                intensity_updates.append(intensities)
                i_updates.append(i_i)
                j_updates.append(j_i)
                k_updates.append(k_i)
                text_updates.append(texts if texts else [])
                u_updates.append([])
                v_updates.append([])
                w_updates.append([])
                
                # Edge trace data
                edge_x, edge_y, edge_z = edge_builder_func(filtered_nodes)
                x_updates.append(edge_x)
                y_updates.append(edge_y)
                z_updates.append(edge_z)
                intensity_updates.append([])
                i_updates.append([])
                j_updates.append([])
                k_updates.append([])
                text_updates.append([])
                u_updates.append([])
                v_updates.append([])
                w_updates.append([])
                
                # Flux cones trace data
                fx, fy, fz, fu, fv, fw, _ = flux_builder_func(filtered_nodes)
                x_updates.append(fx)
                y_updates.append(fy)
                z_updates.append(fz)
                intensity_updates.append([])
                i_updates.append([])
                j_updates.append([])
                k_updates.append([])
                text_updates.append([])
                u_updates.append(fu)
                v_updates.append(fv)
                w_updates.append(fw)
            else:
                # Empty mesh data
                for _ in range(3):  # 3 traces per body (mesh, edge, flux)
                    x_updates.append([])
                    y_updates.append([])
                    z_updates.append([])
                    intensity_updates.append([])
                    i_updates.append([])
                    j_updates.append([])
                    k_updates.append([])
                    text_updates.append([])
                    u_updates.append([])
                    v_updates.append([])
                    w_updates.append([])
        
        # Build step with both trace updates and layout updates to fix axis ranges
        # For method="update": args=[{trace_dict}, {layout_dict}, trace_selector]
        trace_update = {
                'x': x_updates,
                'y': y_updates,
                'z': z_updates,
                'intensity': intensity_updates,
                'i': i_updates,
                'j': j_updates,
                'k': k_updates,
                'text': text_updates,
                'u': u_updates,
                'v': v_updates,
                'w': w_updates
        }
        
        # Add layout update to preserve fixed axis ranges (prevents auto-scaling)
        layout_update = {}
        if package_bounds is not None:
            layout_update = {
                'scene.xaxis.range': list(package_bounds['x']),
                'scene.yaxis.range': list(package_bounds['y']),
                'scene.zaxis.range': list(package_bounds['z']),
                'scene.xaxis.autorange': False,
                'scene.yaxis.autorange': False,
                'scene.zaxis.autorange': False,
                'scene.aspectmode': 'manual'
            }
        
        steps.append(dict(
            method="update",
            args=[trace_update, layout_update, trace_indices],
            label=f"{clip_val:.1f}"
        ))
    
    return dict(
        active=0,
        yanchor="bottom",
        y=position[1],
        xanchor="left",
        x=position[0],
        currentvalue=dict(
            prefix=f"{axis_name.upper()} Clip: ",
            suffix=" mm",
            visible=True,
            xanchor="left"
        ),
        pad=dict(b=10, t=10),
        len=0.25,
        steps=steps
    )


def create_z_scale_slider(
    package_z_len,
    scale_values=None,
    position=(0.02, -0.15),
    for_3d=True,
    package_bounds=None
):
    """
    Create a Z-axis scaling slider for vertical stretching.
    
    Args:
        package_z_len: Package Z dimension in mm
        scale_values: List of scale multipliers (default: [1.0, 1.5, ..., 10.0])
        position: (x, y) position for slider
        for_3d: If True, use scene.aspectratio.z; if False, use yaxis.range
        package_bounds: Dict with 'x', 'y', 'z' tuples of (min, max) to fix axis ranges
    
    Returns:
        Plotly slider dictionary
    """
    if scale_values is None:
        scale_values = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10.0]
    
    steps = []
    
    if for_3d:
        # 3D scene aspect ratio adjustment
        for scale in scale_values:
            relayout_args = {'scene.aspectratio.z': package_z_len * scale}
            
            # Also preserve fixed axis ranges to prevent auto-scaling
            if package_bounds is not None:
                relayout_args.update({
                    'scene.xaxis.range': list(package_bounds['x']),
                    'scene.yaxis.range': list(package_bounds['y']),
                    'scene.zaxis.range': list(package_bounds['z']),
                    'scene.xaxis.autorange': False,
                    'scene.yaxis.autorange': False,
                    'scene.zaxis.autorange': False,
                    'scene.aspectmode': 'manual'
                })
            
            steps.append(dict(
                method="relayout",
                args=[relayout_args],
                label=f"{scale:.1f}x"
            ))
    else:
        # 2D vertical plot y-axis range adjustment
        base_ymin = -0.5
        base_ymax = package_z_len + 0.5
        y_center = 0.5 * (base_ymin + base_ymax)
        y_half = 0.5 * (base_ymax - base_ymin)
        
        for scale in scale_values:
            if scale <= 0:
                continue
            # Larger scale => smaller y-range (visual stretching of Z)
            new_half = y_half / scale
            new_min = y_center - new_half
            new_max = y_center + new_half
            steps.append(dict(
                method="relayout",
                args=[{"yaxis.range": [new_min, new_max]}],
                label=f"{scale:.1f}x"
            ))
    
    return dict(
        active=0,
        yanchor="bottom",
        y=position[1],
        xanchor="left",
        x=position[0],
        currentvalue=dict(
            prefix="Z Scale: ",
            suffix="",
            visible=True,
            xanchor="left"
        ),
        pad=dict(b=10, t=10),
        len=0.25,
        steps=steps
    )


def create_mode_buttons(n_traces, mesh_indices, edge_indices, flux_indices):
    """
    Create Temperature/Flux/Both mode toggle buttons.
    
    Args:
        n_traces: Total number of traces in figure
        mesh_indices: List of mesh trace indices
        edge_indices: List of edge trace indices
        flux_indices: List of flux trace indices
    
    Returns:
        List of button dictionaries
    """
    temp_visible = [False] * n_traces
    flux_visible = [False] * n_traces
    both_visible = [False] * n_traces
    
    mesh_set = set(mesh_indices)
    edge_set = set(edge_indices)
    flux_set = set(flux_indices)
    
    for idx in range(n_traces):
        # Temperature mode shows meshes + edges
        if idx in mesh_set or idx in edge_set:
            temp_visible[idx] = True
            both_visible[idx] = True
        # Flux mode shows cones + wireframe
        if idx in flux_set or idx in edge_set:
            flux_visible[idx] = True
        # Both mode shows everything
        if idx in flux_set:
            both_visible[idx] = True
    
    return [
        dict(
            label="Temperature",
            method="restyle",
            args=["visible", temp_visible],
        ),
        dict(
            label="Flux",
            method="restyle",
            args=["visible", flux_visible],
        ),
        dict(
            label="Both",
            method="restyle",
            args=["visible", both_visible],
        ),
    ]


def create_visibility_buttons(n_traces):
    """
    Create Show All / Hide All body visibility buttons.
    
    Args:
        n_traces: Total number of traces in figure
    
    Returns:
        List of button dictionaries
    """
    return [
        dict(
            label="Show All Bodies",
            method="restyle",
            args=["visible", [True] * n_traces]
        ),
        dict(
            label="Hide All Bodies",
            method="restyle",
            args=["visible", [False] * n_traces]
        ),
    ]


def create_3d_layout(
    title,
    package_bounds,
    sliders,
    buttons_visibility,
    buttons_mode,
    width=1200,
    height=900
):
    """
    Create standardized 3D layout configuration.
    
    Args:
        title: Plot title
        package_bounds: Dict with 'x', 'y', 'z' tuples of (min, max)
        sliders: List of slider dictionaries
        buttons_visibility: List of visibility button dicts
        buttons_mode: List of mode button dicts
        width: Figure width
        height: Figure height
    
    Returns:
        Dictionary of layout parameters
    """
    x_range = package_bounds['x']
    y_range = package_bounds['y']
    z_range = package_bounds['z']
    
    x_len = x_range[1] - x_range[0]
    y_len = y_range[1] - y_range[0]
    z_len = z_range[1] - z_range[0]
    
    return dict(
        title=dict(
            text=title,
            x=0.5,
            xanchor='center',
            y=0.97,
            yanchor='top'
        ),
        scene=dict(
            xaxis=dict(
                title='X (mm)',
                gridcolor='lightgray',
                gridwidth=0.25,
                range=list(x_range),
                autorange=False
            ),
            yaxis=dict(
                title='Y (mm)',
                gridcolor='lightgray',
                gridwidth=0.25,
                range=list(y_range),
                autorange=False
            ),
            zaxis=dict(
                title='Z (mm)',
                gridcolor='lightgray',
                gridwidth=0.25,
                range=list(z_range),
                autorange=False
            ),
            aspectmode='manual',
            aspectratio=dict(
                x=x_len,
                y=y_len,
                z=z_len
            )
        ),
        updatemenus=[
            dict(
                type="buttons",
                direction="right",
                x=0.01,
                y=1.02,
                xanchor='left',
                yanchor='bottom',
                buttons=buttons_visibility,
                showactive=True,
                bgcolor='rgba(255, 255, 255, 0.9)',
                bordercolor='gray',
                borderwidth=1,
                pad=dict(t=5, b=5, l=5, r=5)
            ),
            dict(
                type="buttons",
                direction="right",
                x=0.5,
                y=1.02,
                xanchor='center',
                yanchor='bottom',
                buttons=buttons_mode,
                showactive=True,
                bgcolor='rgba(255, 255, 255, 0.9)',
                bordercolor='gray',
                borderwidth=1,
                pad=dict(t=5, b=5, l=5, r=5)
            ),
        ],
        legend=dict(
            x=0.01,
            y=0.99,
            xanchor='left',
            yanchor='top',
            bgcolor='rgba(255, 255, 255, 0.9)',
            bordercolor='gray',
            borderwidth=1,
            itemsizing='constant',
            itemclick='toggle',
            itemdoubleclick='toggleothers'
        ),
        sliders=sliders,
        width=width,
        height=height,
        hovermode='closest',
        margin=dict(l=0, r=100, t=100, b=50)
    )

