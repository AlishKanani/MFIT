"""
Vertical cross-section heatmap visualization for XZ and YZ planes.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

try:
    import plotly.graph_objects as go
except ImportError:
    go = None

from ..base import BaseVisualizer
from ..geometry.flux_utils import compute_flux_field
from ..geometry.flux_3d_vectors import compute_vector_lengths, FluxVectorConfig
from ..plotly_components.colorscales import (
    HOT_R_COLORSCALE,
    temp_to_color,
    create_colorbar_dict
)
from ..plotly_components.ui_builder import create_z_scale_slider, create_mode_buttons


class VerticalHeatmapVisualizer(BaseVisualizer):
    """Handles vertical cross-section (XZ/YZ plane) heatmap visualizations."""
    
    def _collect_cut_data(self, temperature_all_map, cut_value, plane_type, index_heatmap):
        """
        Collect all nodes across all layers at a specific cut plane.
        
        Args:
            temperature_all_map: Temperature array (nodes x timesteps) in Celsius
            cut_value: X or Y coordinate where the cut is made
            plane_type: 'XZ' or 'YZ'
            index_heatmap: Timestep index to visualize
            
        Returns:
            List of dictionaries containing node data at the cut
        """
        all_cut_data = []
        
        num_nodes = 0
        for layer in self.layers:
            layer_start = num_nodes
            num_nodes += layer.layer_total_nodes()
            
            # Get nodes at this cut for this layer
            selected_nodes = layer.get_nodes_at_cut(cut_value, plane_type)
            
            for node_idx, node, temp_idx in selected_nodes:
                # Get temperature for this node at the specified timestep
                global_temp_idx = layer_start + temp_idx
                temp_value = temperature_all_map[global_temp_idx, index_heatmap]
                
                # Determine horizontal position and length based on plane type
                if plane_type == 'YZ':
                    horizontal_pos = node.y
                    horizontal_length = node.y_length
                elif plane_type == 'XZ':
                    horizontal_pos = node.x
                    horizontal_length = node.x_length
                
                all_cut_data.append({
                    'layer': layer,
                    'node': node,
                    'temp': temp_value,
                    'z_pos': layer.start_point['z'],
                    'thickness': layer.thickness,
                    'horizontal_pos': horizontal_pos,
                    'horizontal_length': horizontal_length
                })
        
        return all_cut_data
    
    def plot_vertical_heatmap(self, temperature_all_map, cut_value, plane_type, index_heatmap):
        """
        Generate a vertical cross-section heatmap showing all layers at a specific cut plane.
        
        Args:
            temperature_all_map: Temperature array (nodes x timesteps) in Celsius
            cut_value: X or Y coordinate where the cut is made
            plane_type: 'XZ' or 'YZ'
            index_heatmap: Timestep index to visualize
        """
        # Collect all nodes across all layers at this cut
        all_cut_data = self._collect_cut_data(temperature_all_map, cut_value, plane_type, index_heatmap)
        
        if not all_cut_data:
            print(f"Warning: No nodes found at {plane_type} cut at {cut_value}")
            return
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Get temperature range for color normalization
        temps = [d['temp'] for d in all_cut_data]
        temps_arr, ambient_mask, temp_min, temp_max = self._prepare_temperature_scale(temps, pad=1.0)
        norm = plt.Normalize(temp_min, temp_max)
        cmap = plt.cm.hot_r
        
        # Draw rectangles for each node
        for idx, data in enumerate(all_cut_data):
            color = self._matplotlib_color(temps_arr[idx], ambient_mask[idx], cmap, norm)
            rect = Rectangle(
                (data['horizontal_pos'], data['z_pos']),
                data['horizontal_length'],
                data['thickness'],
                linewidth=1,
                edgecolor='black',
                facecolor=color
            )
            ax.add_patch(rect)
        
        # Set axis limits
        bounds = self._get_package_bounds()
        if plane_type == 'YZ':
            ax.set_xlim(bounds['y'])
            ax.set_xlabel('Y dimension (mm)')
        elif plane_type == 'XZ':
            ax.set_xlim(bounds['x'])
            ax.set_xlabel('X dimension (mm)')
        
        ax.set_ylim(bounds['z'])
        ax.set_ylabel('Z dimension (mm, scaled 3x for visualization)')
        
        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array(self._active_temperature_values(temps_arr, ambient_mask))
        plt.colorbar(sm, ax=ax).set_label('Temperature (C)')
        
        # Set equal aspect ratio
        ax.set_aspect('equal', adjustable='box')
        
        # Set title
        title = f'{plane_type} Plane Cross-Section at '
        title += f'X = {cut_value:.2f} mm' if plane_type == 'YZ' else f'Y = {cut_value:.2f} mm'
        plt.title(title)
        
        # Save figure
        subdir = f'heatmaps/{plane_type}_heatmaps'
        filename = f'{plane_type.lower()}_cut_at_'
        filename += f'x{cut_value:.2f}mm.png' if plane_type == 'YZ' else f'y{cut_value:.2f}mm.png'
        
        output_dir = self._get_output_dir(subdir)
        plt.savefig(os.path.join(output_dir, filename), dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Generated {plane_type} heatmap (PNG) at {cut_value:.2f} mm")
        
        # Generate interactive HTML version if requested
        if self.args.interactive_heatmaps:
            self.plot_vertical_combined_interactive(temperature_all_map, cut_value, plane_type, index_heatmap)
    
    def plot_vertical_combined_interactive(self, temperature_all_map, cut_value, plane_type, index_heatmap):
        """
        Generate a combined interactive HTML visualization with temperature and flux for vertical cross-sections.
        
        Args:
            temperature_all_map: Temperature array (nodes x timesteps) in Celsius
            cut_value: X or Y coordinate where the cut is made
            plane_type: 'XZ' or 'YZ'
            index_heatmap: Timestep index to visualize
        """
        if go is None:
            print("Warning: plotly not installed. Skipping interactive heatmap generation.")
            print("Install with: pip install plotly")
            return
        
        # ==================== TEMPERATURE DATA ====================
        # Collect all nodes across all layers at this cut
        rectangles_data = []
        
        num_nodes = 0
        for layer in self.layers:
            layer_start = num_nodes
            num_nodes += layer.layer_total_nodes()
            
            # Get nodes at this cut for this layer
            selected_nodes = layer.get_nodes_at_cut(cut_value, plane_type)
            
            for node_idx, node, temp_idx in selected_nodes:
                # Get temperature for this node at the specified timestep
                global_temp_idx = layer_start + temp_idx
                temp_value = temperature_all_map[global_temp_idx, index_heatmap]
                
                # Determine horizontal position and length based on plane type
                if plane_type == 'YZ':
                    horizontal_pos = node.y
                    horizontal_length = node.y_length
                    h_label = 'Y'
                elif plane_type == 'XZ':
                    horizontal_pos = node.x
                    horizontal_length = node.x_length
                    h_label = 'X'
                
                rectangles_data.append({
                    'h_min': horizontal_pos,
                    'h_max': horizontal_pos + horizontal_length,
                    'z_min': layer.start_point['z'],
                    'z_max': layer.start_point['z'] + layer.thickness,
                    'temp': temp_value,
                    'layer_name': layer.layer_name,
                    'h_label': h_label
                })
        
        if not rectangles_data:
            print(f"Warning: No nodes found at {plane_type} cut at {cut_value}")
            return
        
        # Create shapes for temperature rectangles
        temp_shapes = []
        temps = [r['temp'] for r in rectangles_data]
        temps_arr, ambient_mask, temp_min, temp_max = self._prepare_temperature_scale(temps, pad=0.0)
        
        for idx, rect in enumerate(rectangles_data):
            color = 'rgb(255,255,255)' if ambient_mask[idx] else temp_to_color(
                rect['temp'], temp_min, temp_max, HOT_R_COLORSCALE
            )
            temp_shapes.append(dict(
                type='rect',
                x0=rect['h_min'],
                x1=rect['h_max'],
                y0=rect['z_min'],
                y1=rect['z_max'],
                fillcolor=color,
                line=dict(color='black', width=1.25),  # Thicker by 25%
                layer='below'
            ))
        
        # ==================== FLUX DATA ====================
        # Compute full 3D flux field (in Kelvin)
        temperature_all_map_K = temperature_all_map + 273.15
        qx, qy, qz = compute_flux_field(self.layers, self.package, temperature_all_map_K, index_heatmap)
        
        if qx is None:
            print("Warning: No layers available for flux visualization")
            return
        
        # Extract 2D flux components and geometry on the requested cut
        (h, z, qh, qv, mag, layer_names, h_min, h_max, z_min, z_max) = \
            self._collect_flux_data_for_cut(qx, qy, qz, cut_value, plane_type)
        
        if h is None:
            print(f"Warning: No flux data found at {plane_type} cut at {cut_value}")
            return
        
        # Determine geometric scale for arrows
        bounds = self._get_package_bounds()
        if plane_type == "YZ":
            span_h = bounds['y'][1] - bounds['y'][0]
        else:
            span_h = bounds['x'][1] - bounds['x'][0]
        span_z = bounds['z'][1] - bounds['z'][0]
        span = max(span_h, span_z)
        
        # Use shared flux vector scaling
        flux_config = FluxVectorConfig(
            max_vec_frac=0.05,
            min_vec_frac=0.01,
            percentile_ref=90.0,
            min_visible_factor=0.15
        )
        length_factor = compute_vector_lengths(mag, span, flux_config) / span
        base_scale = span
        
        # Create outline shapes for flux view (thin dotted grey rectangles)
        flux_shapes = []
        num_nodes = 0
        for layer in self.layers:
            layer_start = num_nodes
            num_nodes += layer.layer_total_nodes()
            selected_nodes = layer.get_nodes_at_cut(cut_value, plane_type)
            for _, node, _ in selected_nodes:
                z0 = layer.start_point["z"]
                z1 = layer.start_point["z"] + layer.thickness
                if plane_type == "YZ":
                    h0 = node.y
                    h1 = node.y + node.y_length
                else:
                    h0 = node.x
                    h1 = node.x + node.x_length
                flux_shapes.append(
                    dict(
                        type="rect",
                        x0=h0,
                        x1=h1,
                        y0=z0,
                        y1=z1,
                        line=dict(
                            color="rgba(120,120,120,0.7)",
                            width=0.5,
                            dash="dot",
                        ),
                        fillcolor="rgba(0,0,0,0)",
                        layer="below",
                    )
                )
        
        # ==================== CREATE FIGURE WITH TRACES ====================
        fig = go.Figure()
        
        # --- Temperature Traces (indices 0-1) ---
        # Trace 0: Temperature hover info
        hover_x_temp = []
        hover_y_temp = []
        hover_text_temp = []
        for rect in rectangles_data:
            hover_x_temp.append((rect['h_min'] + rect['h_max']) / 2)
            hover_y_temp.append((rect['z_min'] + rect['z_max']) / 2)
            hover_text_temp.append(
                f"Layer: {rect['layer_name']}<br>" +
                f"{rect['h_label']}: {rect['h_min']:.3f} - {rect['h_max']:.3f} mm<br>" +
                f"Z: {rect['z_min']:.3f} - {rect['z_max']:.3f} mm<br>" +
                f"Temperature: {rect['temp']:.2f} °C"
            )
        
        fig.add_trace(go.Scatter(
            x=hover_x_temp,
            y=hover_y_temp,
            mode='markers',
            marker=dict(size=0.1, opacity=0),
            hovertext=hover_text_temp,
            hoverinfo='text',
            showlegend=False,
            visible=True
        ))
        
        # Trace 1: Temperature colorbar (positioned on the right)
        plotly_colorscale = [[c[0], c[1]] for c in HOT_R_COLORSCALE]
        temp_colorbar = create_colorbar_dict("Temp (°C)")
        temp_colorbar['x'] = 1.02  # Position to the right
        temp_colorbar['xanchor'] = 'left'
        
        fig.add_trace(go.Scatter(
            x=[None],
            y=[None],
            mode='markers',
            marker=dict(
                size=0.1,
                opacity=0,
                colorscale=plotly_colorscale,
                cmin=temp_min,
                cmax=temp_max,
                colorbar=temp_colorbar,
                showscale=True
            ),
            showlegend=False,
            hoverinfo='skip',
            visible=True
        ))
        
        # --- Flux Traces (indices 2-3) ---
        # Trace 2: Flux magnitude markers with colorbar (positioned further right)
        h_label = "Y" if plane_type == "YZ" else "X"
        hover_text_flux = []
        for name, hm0, hm1, zm0, zm1, m_val, uh, uv in zip(
            layer_names, h_min, h_max, z_min, z_max, mag, qh, qv
        ):
            hover_text_flux.append(
                f"Layer: {name}<br>"
                f"{h_label}: {hm0:.3f} - {hm1:.3f} mm<br>"
                f"Z: {zm0:.3f} - {zm1:.3f} mm<br>"
                f"|q|={m_val:.3e} W/m²<br>"
                f"qh={uh:.3e} W/m²<br>"
                f"qv={uv:.3e} W/m²"
            )
        
        flux_colorbar = create_colorbar_dict("|q| (W/m²)")
        flux_colorbar['x'] = 1.15  # Position further right to avoid overlap with temp colorbar
        flux_colorbar['xanchor'] = 'left'
        
        fig.add_trace(go.Scatter(
            x=h,
            y=z,
            mode="markers",
            marker=dict(
                size=8,
                color=mag,
                colorscale="Viridis",
                colorbar=flux_colorbar,
                showscale=True,
            ),
            hovertext=hover_text_flux,
            hoverinfo="text",
            showlegend=False,
            visible=False
        ))
        
        # Trace 3: Dummy trace for flux (needed for consistent indexing)
        fig.add_trace(go.Scatter(
            x=[None],
            y=[None],
            mode='markers',
            marker=dict(size=0.1, opacity=0),
            showlegend=False,
            hoverinfo='skip',
            visible=False
        ))
        
        # ==================== CREATE ANNOTATIONS FOR FLUX ARROWS ====================
        flux_annotations = []
        for xi, zi, ui, vi, lf, m in zip(h, z, qh, qv, length_factor, mag):
            if lf <= 0.0 or m <= 0.0:
                continue
            ux = ui / m
            uz = vi / m
            dx = base_scale * lf * ux
            dz = base_scale * lf * uz
            flux_annotations.append(
                dict(
                    x=xi + dx,
                    y=zi + dz,
                    ax=xi,
                    ay=zi,
                    xref="x",
                    yref="y",
                    axref="x",
                    ayref="y",
                    showarrow=True,
                    arrowhead=2,
                    arrowsize=1.0,
                    arrowwidth=1,
                    arrowcolor="black",
                )
            )
        
        # ==================== MODE BUTTONS ====================
        # Create mode buttons similar to 3D version
        # Temperature mode: show traces 0-1, hide 2-3, use temp_shapes, no annotations
        # Flux mode: show traces 2-3, hide 0-1, use flux_shapes, show annotations
        # Both mode: show all traces, use temp_shapes, show annotations
        
        buttons_mode = [
            dict(
                label="Temperature",
                method="update",
                args=[
                    {"visible": [True, True, False, False]},
                    {
                        "shapes": temp_shapes,
                        "annotations": []
                    }
                ]
            ),
            dict(
                label="Heat Flux",
                method="update",
                args=[
                    {"visible": [False, False, True, False]},
                    {
                        "shapes": flux_shapes,
                        "annotations": flux_annotations
                    }
                ]
            ),
            dict(
                label="Both",
                method="update",
                args=[
                    {"visible": [True, True, True, False]},
                    {
                        "shapes": temp_shapes,
                        "annotations": flux_annotations
                    }
                ]
            ),
        ]
        
        # ==================== LAYOUT ====================
        if plane_type == 'YZ':
            x_label = 'Y dimension (mm)'
            title = f'YZ Plane Cross-Section at X = {cut_value:.2f} mm (Interactive)'
            x_range = list(bounds['y'])
        else:
            x_label = 'X dimension (mm)'
            title = f'XZ Plane Cross-Section at Y = {cut_value:.2f} mm (Interactive)'
            x_range = list(bounds['x'])
        
        # Create Z-scale slider (2D mode, adjusts y-axis range)
        z_scale_slider = create_z_scale_slider(
            self.common_utils.package_z_len,
            position=(0.02, -0.1),
            for_3d=False
        )
        
        fig.update_layout(
            title=title,
            xaxis=dict(
                title=x_label,
                range=x_range,
                constrain='domain',
                showgrid=False,
            ),
            yaxis=dict(
                title='Z dimension (mm)',
                # Start from full package Z range; slider will adjust this
                range=[-0.5, self.common_utils.package_z_len + 0.5],
                showgrid=False,
            ),
            width=1100,
            height=600,
            hovermode='closest',
            margin=dict(l=60, r=180, t=80, b=100),  # Increased right margin for two colorbars
            shapes=temp_shapes,  # Start with temperature view
            annotations=[],  # Start without flux arrows
            sliders=[z_scale_slider],
            updatemenus=[
                dict(
                    type="buttons",
                    direction="left",
                    buttons=buttons_mode,
                    pad={"r": 10, "t": 10},
                    showactive=True,
                    x=0.5,
                    xanchor="center",
                    y=1.08,
                    yanchor="top"
                )
            ]
        )
        
        # Save to HTML with scroll zoom enabled
        subdir = f'heatmaps/{plane_type}_combined'
        filename = f'{plane_type.lower()}_cut_at_'
        filename += f'x{cut_value:.2f}mm.html' if plane_type == 'YZ' else f'y{cut_value:.2f}mm.html'
        
        config = {'scrollZoom': True}
        self._save_plotly_figure(fig, subdir, filename, config=config)
        print(f"Generated combined {plane_type} visualization at {cut_value:.2f} mm")
    
    def generate_vertical_heatmaps(self, temperature_all_map, index_heatmap):
        """
        Generate all requested vertical heatmaps based on command-line arguments.
        
        Args:
            temperature_all_map: Temperature array (nodes x timesteps) in Celsius
            index_heatmap: Timestep index to visualize
        """
        if not self.args.vertical_planes:
            return
        
        # Parse which planes to generate
        planes = [p.strip().upper() for p in self.args.vertical_planes.split(',') if p.strip()]
        
        # Generate YZ plane heatmaps (cuts perpendicular to X-axis)
        if 'YZ' in planes and self.args.yz_cuts:
            yz_cut_values = [float(v.strip()) for v in self.args.yz_cuts.split(',') if v.strip()]
            for cut_value in yz_cut_values:
                self.plot_vertical_heatmap(temperature_all_map, cut_value, 'YZ', index_heatmap)
        
        # Generate XZ plane heatmaps (cuts perpendicular to Y-axis)
        if 'XZ' in planes and self.args.xz_cuts:
            xz_cut_values = [float(v.strip()) for v in self.args.xz_cuts.split(',') if v.strip()]
            for cut_value in xz_cut_values:
                self.plot_vertical_heatmap(temperature_all_map, cut_value, 'XZ', index_heatmap)

    # ----------------------- Vertical heat-flux vector plots -----------------------
    
    def _collect_flux_data_for_cut(self, qx, qy, qz, cut_value, plane_type):
        """
        Collect node-center positions and 2D flux components for a given cut.

        Returns:
            Tuple of (h_centers, z_centers, qh, qv, magnitudes, layer_names,
                     h_min, h_max, z_min, z_max)
        """
        h_list = []
        z_list = []
        qh_list = []
        qv_list = []
        mag_list = []
        layer_names = []
        h_min_list = []
        h_max_list = []
        z_min_list = []
        z_max_list = []

        num_nodes = 0
        for layer in self.layers:
            layer_start = num_nodes
            num_nodes += layer.layer_total_nodes()

            selected_nodes = layer.get_nodes_at_cut(cut_value, plane_type)

            for _, node, temp_idx in selected_nodes:
                global_idx = layer_start + temp_idx

                # Node center and extents in Z
                z0 = layer.start_point['z']
                z1 = layer.start_point['z'] + layer.thickness
                z_center = z0 + 0.5 * layer.thickness

                if plane_type == 'YZ':
                    h_center = node.y + 0.5 * node.y_length
                    qh_val = qy[global_idx]
                    qv_val = qz[global_idx]
                    h0 = node.y
                    h1 = node.y + node.y_length
                else:  # 'XZ'
                    h_center = node.x + 0.5 * node.x_length
                    qh_val = qx[global_idx]
                    qv_val = qz[global_idx]
                    h0 = node.x
                    h1 = node.x + node.x_length

                mag = float(np.sqrt(qh_val ** 2 + qv_val ** 2))

                h_list.append(h_center)
                z_list.append(z_center)
                qh_list.append(qh_val)
                qv_list.append(qv_val)
                mag_list.append(mag)
                layer_names.append(layer.layer_name)
                h_min_list.append(h0)
                h_max_list.append(h1)
                z_min_list.append(z0)
                z_max_list.append(z1)

        if not h_list:
            return None, None, None, None, None, None, None, None, None

        return (
            np.array(h_list),
            np.array(z_list),
            np.array(qh_list),
            np.array(qv_list),
            np.array(mag_list),
            np.array(layer_names),
            np.array(h_min_list),
            np.array(h_max_list),
            np.array(z_min_list),
            np.array(z_max_list),
        )

    def generate_vertical_flux_maps(self, temperature_all_map, index_heatmap):
        """
        DEPRECATED: Flux visualization is now combined with temperature in generate_vertical_heatmaps.
        
        This method is kept for backward compatibility but does nothing.
        Temperature and flux are now shown together in a single interactive visualization
        with toggle buttons to switch between views.

        Args:
            temperature_all_map: Temperature array (nodes x timesteps) in Kelvin
            index_heatmap: Timestep index to visualize
        """
        # Flux is now generated as part of the combined interactive visualization
        # in generate_vertical_heatmaps(), so this method no longer needs to do anything
        pass

