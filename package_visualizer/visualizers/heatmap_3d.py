"""
3D interactive heatmap visualization (Plotly only).
"""

import numpy as np
try:
    import plotly.graph_objects as go
except ImportError:
    go = None

from ..base import BaseVisualizer
from ..geometry.flux_utils import compute_flux_field
from ..geometry.flux_3d_vectors import build_flux_cones_for_body, FluxVectorConfig
from ..geometry.mesh_builder import (
    collect_body_data,
    build_mesh_from_nodes,
    build_edges_from_nodes
)
from ..plotly_components.trace_factory import (
    create_mesh_trace,
    create_edge_trace,
    create_cone_trace
)
from ..plotly_components.ui_builder import (
    create_z_scale_slider,
    create_mode_buttons,
    create_visibility_buttons,
    create_3d_layout
)


class Heatmap3D(BaseVisualizer):
    """Handles 3D interactive heatmap visualizations (Plotly only)."""
    
    def plot_3d_heatmap(self, temperature_all_map, index_heatmap):
        """
        Generate an interactive 3D visualization with temperature coloring and body filtering.
        
        Args:
            temperature_all_map: Temperature array (nodes x timesteps) in Celsius
            index_heatmap: Timestep index to visualize
        """
        if go is None:
            print("Warning: plotly not installed. Skipping 3D heatmap generation.")
            print("Install with: pip install plotly")
            return
        
        if self.layers is None:
            print("Warning: No layers available for 3D visualization")
            return
        
        # Collect all nodes organized by body with ambient tagging
        body_data = collect_body_data(
            self.layers,
            self.package,
            temperature_all_map,
            index_heatmap,
            ambient_temp_c=self.ambient_temp_C,
            ambient_tol=self.ambient_tol
        )
        
        if not body_data:
            print("Warning: No node data collected for 3D visualization")
            return
        
        # Get global temperature range for consistent coloring (ignore ambient nodes)
        active_temps = [
            node['temp']
            for body_nodes in body_data.values()
            for node in body_nodes
            if not node.get('is_ambient', False)
        ]
        if active_temps:
            active_arr = np.asarray(active_temps, dtype=float)
            temp_min = float(active_arr.min())
            temp_max = float(active_arr.max())
            if temp_max <= temp_min:
                temp_max = temp_min + 1.0
        else:
            temp_min = self.ambient_temp_C
            temp_max = self.ambient_temp_C + 1.0
        
        # Compute global flux field (qx, qy, qz) for this timestep (Kelvin input)
        temperature_all_map_K = temperature_all_map + 273.15
        qx_all, qy_all, qz_all = compute_flux_field(
            self.layers, self.package, temperature_all_map_K, index_heatmap
        )
        
        # Get package span for flux vector scaling
        span = self._get_package_span()
        
        # Configure flux vector appearance
        flux_config = FluxVectorConfig(
            max_vec_frac=0.08,  # 8% of span for largest vectors
            min_vec_frac=0.01   # 1% of span for smallest vectors
        )
        
        # Create figure and add traces for each body
        fig = go.Figure()
        
        temp_traces = []
        flux_traces = []
        mesh_indices = []
        edge_indices = []
        body_render_order = []
        processed_body_nodes = {}
        show_colorbar = True
        
        for body_name, nodes in body_data.items():
            if not nodes:
                continue
            body_render_order.append(body_name)

            prepped_nodes = []
            has_active = False
            for node in nodes:
                node_copy = node.copy()
                if node.get('is_ambient', False):
                    node_copy['display_temp'] = temp_min
                else:
                    has_active = True
                    node_copy['display_temp'] = node['temp']
                prepped_nodes.append(node_copy)

            # Keep detailed hover text as requested
            (
                vertices_x,
                vertices_y,
                vertices_z,
                intensities,
                i_i,
                j_i,
                k_i,
                hover_texts
            ) = build_mesh_from_nodes(
                prepped_nodes,
                hover_formatter=self._format_hover_text
            )
            mesh_trace = create_mesh_trace(
                vertices_x, vertices_y, vertices_z, intensities,
                i_i, j_i, k_i, temp_min, temp_max, body_name,
                colorscale='Hot_r',
                show_colorbar=(show_colorbar and has_active),
                hover_texts=hover_texts
            )
            fig.add_trace(mesh_trace)
            temp_traces.append(len(fig.data) - 1)
            mesh_indices.append(len(fig.data) - 1)
            if show_colorbar and has_active:
                show_colorbar = False

            # Build edge wireframe (keep as requested)
            edge_x, edge_y, edge_z = build_edges_from_nodes(nodes)
            edge_trace = create_edge_trace(edge_x, edge_y, edge_z, body_name)
            fig.add_trace(edge_trace)
            temp_traces.append(len(fig.data) - 1)
            edge_indices.append(len(fig.data) - 1)

            # Build flux cones (optional - very expensive for file size)
            # NOTE: 3D flux disabled per user request - only 2D flux is needed
            show_flux = False  # Set to True if you need 3D flux visualization (greatly increases file size)
            if show_flux:
                (x_centers, y_centers, z_centers, u_scaled, v_scaled, w_scaled, mag_arr) = \
                    build_flux_cones_for_body(nodes, qx_all, qy_all, qz_all, span, flux_config)
                cone_trace = create_cone_trace(
                    x_centers, y_centers, z_centers,
                    u_scaled, v_scaled, w_scaled,
                    mag_arr, body_name,
                    visible=False
                )
                fig.add_trace(cone_trace)
                flux_traces.append(len(fig.data) - 1)
            processed_body_nodes[body_name] = prepped_nodes
        
        # Total number of traces
        n_traces = len(fig.data)
        
        # Create UI control buttons
        buttons_visibility = create_visibility_buttons(n_traces)
        buttons_mode = create_mode_buttons(n_traces, mesh_indices, edge_indices, flux_traces)
        
        bounds = self._get_package_bounds()
        
        # Create Z-scale slider
        z_scale_slider = create_z_scale_slider(
            self.common_utils.package_z_len,
            position=(0.02, -0.15),
            for_3d=True,
            package_bounds=bounds
        )
        
        sliders = [z_scale_slider]
        
        # Configure layout
        title = '3D Thermal Visualization (Interactive)'
        layout_config = create_3d_layout(
            title=title,
            package_bounds=bounds,
            sliders=sliders,
            buttons_visibility=buttons_visibility,
            buttons_mode=buttons_mode
        )
        
        fig.update_layout(**layout_config)
        
        # Save to HTML (in 3d/plotly subfolder)
        self._save_plotly_figure(fig, 'heatmaps/3d/plotly', '3d_thermal_view.html')

    @staticmethod
    def _format_hover_text(node):
        """
        Build a descriptive hover label for a node using its geometry.
        """
        layer_name = node.get('layer_name', 'Unknown layer')
        body_name = node.get('body_name', 'Unknown body')
        x_min = float(node.get('x_min', 0.0))
        x_max = float(node.get('x_max', 0.0))
        y_min = float(node.get('y_min', 0.0))
        y_max = float(node.get('y_max', 0.0))
        z_min = float(node.get('z_min', 0.0))
        z_max = float(node.get('z_max', 0.0))
        x_len = x_max - x_min
        y_len = y_max - y_min
        z_len = z_max - z_min
        temp = float(node.get('display_temp', node.get('temp', 0.0)))
        
        return (
            f"<b>Layer:</b> {layer_name}<br>"
            f"<b>Body:</b> {body_name}<br>"
            f"<b>X:</b> {x_min:.3f} – {x_max:.3f} mm (Δ={x_len:.3f} mm)<br>"
            f"<b>Y:</b> {y_min:.3f} – {y_max:.3f} mm (Δ={y_len:.3f} mm)<br>"
            f"<b>Z:</b> {z_min:.3f} – {z_max:.3f} mm (Δ={z_len:.3f} mm)<br>"
            f"<b>Temp:</b> {temp:.2f} °C"
        )
