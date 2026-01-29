"""
2D heatmap visualization for XY plane heatmaps and floorplans.
"""

import os

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from ..base import BaseVisualizer


class Heatmap2D(BaseVisualizer):
    """Handles 2D (XY plane) heatmap and floorplan visualizations."""
    
    def plot_layer_heatmap(self, layer, temperature_all_layer):
        """
        Generate a horizontal (XY plane) heatmap for a single layer.
        
        Args:
            layer: Layer_chiplet instance
            temperature_all_layer: Temperature array for this layer (in Celsius)
        """
        fig, ax = plt.subplots()

        temps_arr, ambient_mask, temp_min, temp_max = self._prepare_temperature_scale(
            temperature_all_layer, pad=1.0
        )
        norm = plt.Normalize(temp_min, temp_max)
        cmap = plt.cm.hot_r

        if layer.is_layer_under_chiplet() and not layer.args.is_homogeneous:
            for i in range(layer.total_nodes):
                colour = self._matplotlib_color(
                    temps_arr[i], ambient_mask[i], cmap, norm
                )
                rect = Rectangle((layer.nodes[i].x, layer.nodes[i].y), 
                                    layer.nodes[i].x_length, layer.nodes[i].y_length, 
                                    linewidth=1, edgecolor='black',
                                    facecolor=colour)
                ax.add_patch(rect)
        
        else:
            for i in range(layer.total_x_nodes):
                for j in range(layer.total_y_nodes):
                    idx = i*layer.total_y_nodes + j
                    colour = self._matplotlib_color(
                        temps_arr[idx], ambient_mask[idx], cmap, norm
                    )
                    rect = Rectangle((layer.nodes[i][j].x, layer.nodes[i][j].y), 
                                        layer.nodes[i][j].x_length, layer.nodes[i][j].y_length,
                                        linewidth=1, edgecolor='black',
                                        facecolor=colour) 
                    ax.add_patch(rect)

        ax.set_xlim(-0.5, self.common_utils.package_x_len + 0.5 + 1)
        ax.set_ylim(-0.5, self.common_utils.package_y_len + 0.5 + 1)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array(self._active_temperature_values(temps_arr, ambient_mask))
        plt.colorbar(sm, ax=ax).set_label('Temperature (C)')

        plt.gca().set_aspect('equal', adjustable='box')
        plt.title(layer.layer_name + ' Heatmap')
        plt.xlabel('X dimension (mm)')
        plt.ylabel('Y dimension (mm)')

        heatmap_dir = self._get_output_dir('heatmaps/XY_heatmaps')
        
        plt.savefig(os.path.join(heatmap_dir, layer.layer_name + '_heatmap.png'), dpi=300, bbox_inches='tight')

        # close the plot
        plt.close(fig)
    
    def generate_floorplan_visual(self):
        """Generate the floorplan of the package and chiplet layers."""
        num_nodes = 0
        for layer in self.layers:
            layer.plot_layer(utils=self.common_utils, layer_start=num_nodes)
            num_nodes += layer.layer_total_nodes()
