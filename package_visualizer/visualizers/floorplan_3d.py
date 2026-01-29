"""
3D interactive floorplan visualization (Plotly only).
"""

try:
    import plotly.graph_objects as go
except ImportError:
    go = None

from ..base import BaseVisualizer
from ..geometry.mesh_builder import collect_body_data, build_mesh_from_nodes, build_edges_from_nodes
from ..plotly_components.trace_factory import create_solid_mesh_trace, create_edge_trace
from ..plotly_components.ui_builder import (
    create_z_scale_slider,
    create_visibility_buttons,
    create_floorplan_mode_buttons,
    create_3d_layout
)


class Floorplan3D(BaseVisualizer):
    """Handles 3D interactive floorplan visualizations (Plotly only)."""

    _DEFAULT_COLORS = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
        "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173",
        "#3182bd", "#31a354", "#756bb1", "#636363", "#e6550d"
    ]

    def plot_3d_floorplan(self):
        """
        Generate an interactive 3D floorplan visualization with body coloring.
        """
        if go is None:
            print("Warning: plotly not installed. Skipping 3D floorplan generation.")
            print("Install with: pip install plotly")
            return

        if self.layers is None:
            print("Warning: No layers available for 3D visualization")
            return

        body_data = collect_body_data(
            self.layers,
            self.package,
            temperature_all_map=None,
            index_heatmap=0,
            ambient_temp_c=self.ambient_temp_C,
            ambient_tol=self.ambient_tol,
            split_by_region=True
        )

        if not body_data:
            print("Warning: No node data collected for 3D visualization")
            return

        fig = go.Figure()

        mesh_indices = []
        edge_indices = []

        body_names = list(body_data.keys())
        for idx, body_name in enumerate(body_names):
            nodes = body_data.get(body_name, [])
            if not nodes:
                continue

            color = self._DEFAULT_COLORS[idx % len(self._DEFAULT_COLORS)]

            (
                vertices_x,
                vertices_y,
                vertices_z,
                _intensities,
                i_i,
                j_i,
                k_i,
                hover_texts
            ) = build_mesh_from_nodes(
                nodes,
                hover_formatter=self._format_hover_text
            )

            mesh_trace = create_solid_mesh_trace(
                vertices_x, vertices_y, vertices_z,
                i_i, j_i, k_i,
                body_name,
                color=color,
                hover_texts=hover_texts
            )
            fig.add_trace(mesh_trace)
            mesh_indices.append(len(fig.data) - 1)

            edge_x, edge_y, edge_z = build_edges_from_nodes(nodes)
            edge_trace = create_edge_trace(edge_x, edge_y, edge_z, body_name)
            fig.add_trace(edge_trace)
            edge_indices.append(len(fig.data) - 1)

        n_traces = len(fig.data)

        buttons_visibility = create_visibility_buttons(n_traces)
        buttons_mode = create_floorplan_mode_buttons(n_traces, mesh_indices, edge_indices)

        bounds = self._get_package_bounds()

        z_scale_slider = create_z_scale_slider(
            self.common_utils.package_z_len,
            position=(0.02, -0.15),
            for_3d=True,
            package_bounds=bounds
        )

        layout_config = create_3d_layout(
            title='3D Floorplan Visualization (Interactive)',
            package_bounds=bounds,
            sliders=[z_scale_slider],
            buttons_visibility=buttons_visibility,
            buttons_mode=buttons_mode
        )

        fig.update_layout(**layout_config)

        self._save_plotly_figure(fig, 'floorplan/3d', '3d_floorplan_view.html')

    @staticmethod
    def _format_hover_text(node):
        layer_name = node.get('layer_name', 'Unknown layer')
        body_name = node.get('body_name', 'Unknown body')
        material_name = node.get('material_name', None)
        material_region = node.get('material_region', None)
        x_min = float(node.get('x_min', 0.0))
        x_max = float(node.get('x_max', 0.0))
        y_min = float(node.get('y_min', 0.0))
        y_max = float(node.get('y_max', 0.0))
        z_min = float(node.get('z_min', 0.0))
        z_max = float(node.get('z_max', 0.0))
        x_len = x_max - x_min
        y_len = y_max - y_min
        z_len = z_max - z_min

        material_line = f"<b>Material:</b> {material_name}<br>" if material_name else ""
        region_line = f"<b>Region:</b> {material_region}<br>" if material_region else ""

        return (
            f"<b>Layer:</b> {layer_name}<br>"
            f"<b>Body:</b> {body_name}<br>"
            f"{material_line}"
            f"{region_line}"
            f"<b>X:</b> {x_min:.3f} – {x_max:.3f} mm (Δ={x_len:.3f} mm)<br>"
            f"<b>Y:</b> {y_min:.3f} – {y_max:.3f} mm (Δ={y_len:.3f} mm)<br>"
            f"<b>Z:</b> {z_min:.3f} – {z_max:.3f} mm (Δ={z_len:.3f} mm)"
        )
