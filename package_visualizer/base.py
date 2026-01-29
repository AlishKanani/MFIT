"""
Base visualizer class with common utilities and properties.
"""

import os
import numpy as np


class BaseVisualizer:
    """Base class for all visualization components with shared utilities."""
    
    def __init__(self, package):
        """
        Initialize base visualizer with a package instance.
        
        Args:
            package: Chiplet_package instance containing layers and simulation data
        """
        self.package = package
        self.args = package.args
        self.common_utils = package.common_utils
        self.layers = getattr(package, 'layers', None)
        self.ambient_temp_K = float(getattr(self.args, 'ambient_temp_K', 300.0))
        self.ambient_temp_C = float(getattr(self.args, 'ambient_temp_C', self.ambient_temp_K - 273.15))
        self.ambient_tol = float(getattr(self.args, 'ambient_temp_tolerance', 1e-6))
    
    def _get_output_dir(self, subdir):
        """
        Helper to get output directory path and create it if needed.
        
        Args:
            subdir: Subdirectory path relative to output_dir
            
        Returns:
            Full path to the output subdirectory
        """
        base = getattr(self.args, "output_rc_dir", self.args.output_dir)
        full_path = os.path.join(base, subdir)
        if not os.path.exists(full_path):
            os.makedirs(full_path)
        return full_path
    
    def _save_plotly_figure(self, fig, subdir, filename, config=None):
        """
        Save Plotly figure to HTML with consistent patterns.
        
        Args:
            fig: Plotly Figure object
            subdir: Subdirectory relative to output_dir
            filename: Output filename
            config: Optional Plotly config dict (e.g., {'scrollZoom': True})
        
        Returns:
            Full path to saved file
        """
        output_dir = self._get_output_dir(subdir)
        output_path = os.path.join(output_dir, filename)
        
        if config is None:
            config = {}
        
        fig.write_html(output_path, config=config)
        print(f"Generated {filename}: {output_path}")
        return output_path
    
    def _get_package_bounds(self):
        """
        Return package dimensions for axis ranges.
        
        Returns:
            Dictionary with 'x', 'y', 'z' keys containing (min, max) tuples
        """
        return {
            'x': (-0.5, self.common_utils.package_x_len + 0.5),
            'y': (-0.5, self.common_utils.package_y_len + 0.5),
            'z': (-0.5, self.common_utils.package_z_len + 0.5),
        }
    
    def _get_package_span(self):
        """
        Get the maximum package dimension for vector scaling.
        
        Returns:
            Maximum dimension (mm)
        """
        return max(
            self.common_utils.package_x_len + 1.0,
            self.common_utils.package_y_len + 1.0,
            self.common_utils.package_z_len + 1.0
        )

    def _prepare_temperature_scale(self, temps, pad=0.0):
        """
        Normalize temperature arrays by separating ambient nodes and computing a scale that ignores them.

        Args:
            temps: Iterable of temperature values (°C)
            pad:   Optional absolute padding (°C) added to min/max for plotting

        Returns:
            tuple -> (temps_array, ambient_mask, temp_min, temp_max)
        """
        temps_arr = np.asarray(temps, dtype=float).reshape(-1)
        ambient_mask = np.isclose(temps_arr, self.ambient_temp_C, atol=self.ambient_tol)
        active = temps_arr[~ambient_mask]
        if active.size == 0:
            active = np.array([self.ambient_temp_C], dtype=float)
        temp_min = float(active.min())
        temp_max = float(active.max())
        if temp_max <= temp_min:
            temp_max = temp_min + 1.0
        if pad:
            temp_min -= pad
            temp_max += pad
        return temps_arr, ambient_mask, temp_min, temp_max

    def _active_temperature_values(self, temps_arr, ambient_mask):
        active = temps_arr[~ambient_mask]
        if active.size == 0:
            active = np.array([self.ambient_temp_C], dtype=float)
        return active

    @staticmethod
    def _matplotlib_color(temp, is_ambient, cmap, norm):
        if is_ambient:
            return (1.0, 1.0, 1.0, 1.0)
        return cmap(norm(temp))
