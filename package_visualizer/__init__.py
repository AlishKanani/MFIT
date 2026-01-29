"""
Visualization module for thermal package analysis.
Handles generation of floorplans, heatmaps, and interactive plots.

This module provides a facade pattern where PackageVisualizer delegates
to specialized visualization classes for different visualization types.
"""

from .base import BaseVisualizer
from .visualizers.heatmap_2d import Heatmap2D
from .visualizers.heatmap_vertical import VerticalHeatmapVisualizer
from .visualizers.heatmap_3d import Heatmap3D
from .visualizers.floorplan_3d import Floorplan3D


class PackageVisualizer:
    """
    Main visualizer facade that delegates to specialized visualization classes.
    
    This class maintains backward compatibility with the original API while
    organizing code into smaller, more maintainable modules.
    """
    
    def __init__(self, package):
        """
        Initialize visualizer with a package instance.
        
        Args:
            package: Chiplet_package instance containing layers and simulation data
        """
        self.package = package
        self.args = package.args
        self.common_utils = package.common_utils
        self.layers = getattr(package, 'layers', None)
        
        # Initialize specialized visualizers
        self._heatmap_2d = Heatmap2D(package)
        self._heatmap_vertical = VerticalHeatmapVisualizer(package)
        self._heatmap_3d = Heatmap3D(package)
        self._floorplan_3d = Floorplan3D(package)
    
    # Delegate 2D heatmap methods
    def plot_layer_heatmap(self, layer, temperature_all_layer):
        """Generate a horizontal (XY plane) heatmap for a single layer."""
        return self._heatmap_2d.plot_layer_heatmap(layer, temperature_all_layer)
    
    def generate_floorplan_visual(self):
        """Generate the floorplan of the package and chiplet layers (2D)."""
        return self._heatmap_2d.generate_floorplan_visual()

    def generate_floorplan_2d(self):
        """Generate the 2D floorplan of the package and chiplet layers."""
        return self._heatmap_2d.generate_floorplan_visual()
    
    # Delegate vertical heatmap methods
    def plot_vertical_heatmap(self, temperature_all_map, cut_value, plane_type, index_heatmap):
        """Generate a vertical cross-section heatmap showing all layers at a specific cut plane."""
        return self._heatmap_vertical.plot_vertical_heatmap(
            temperature_all_map, cut_value, plane_type, index_heatmap
        )
    
    def plot_vertical_heatmap_interactive(self, temperature_all_map, cut_value, plane_type, index_heatmap):
        """Generate an interactive HTML heatmap using Plotly for vertical cross-sections."""
        return self._heatmap_vertical.plot_vertical_heatmap_interactive(
            temperature_all_map, cut_value, plane_type, index_heatmap
        )

    def generate_vertical_heatmaps(self, temperature_all_map, index_heatmap):
        """Generate all requested vertical heatmaps based on command-line arguments."""
        return self._heatmap_vertical.generate_vertical_heatmaps(temperature_all_map, index_heatmap)

    # Delegate vertical flux vector plot methods
    def plot_vertical_flux_interactive(self, temperature_all_map, cut_value, plane_type, index_heatmap):
        """Generate an interactive 2D heat-flux vector plot for a vertical cross-section."""
        return self._heatmap_vertical.plot_vertical_flux_interactive(
            temperature_all_map, cut_value, plane_type, index_heatmap
        )

    def generate_vertical_flux_maps(self, temperature_all_map, index_heatmap):
        """Generate all requested vertical heat-flux maps based on command-line arguments."""
        return self._heatmap_vertical.generate_vertical_flux_maps(temperature_all_map, index_heatmap)
    
    # Delegate 3D heatmap methods
    def plot_3d_heatmap(self, temperature_all_map, index_heatmap):
        """Generate an interactive 3D visualization with temperature coloring and body filtering."""
        return self._heatmap_3d.plot_3d_heatmap(temperature_all_map, index_heatmap)

    def plot_3d_floorplan(self):
        """Generate an interactive 3D floorplan visualization with body coloring."""
        return self._floorplan_3d.plot_3d_floorplan()
