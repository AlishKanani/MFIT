"""
3D flux vector computation and scaling utilities.
"""

from dataclasses import dataclass
import numpy as np


@dataclass
class FluxVectorConfig:
    """Configuration for flux vector visualization scaling."""
    max_vec_frac: float = 0.08  # Fraction of domain span for largest vectors
    min_vec_frac: float = 0.01  # Fraction of domain span for smallest vectors
    percentile_ref: float = 90.0  # Percentile to use as magnitude reference
    min_visible_factor: float = 0.0  # Minimum visibility factor (0.0 = no floor)


def compute_vector_lengths(magnitudes, span, config=None):
    """
    Compute scaled vector lengths using log-magnitude mapping.
    
    This function maps flux magnitudes to physical vector lengths using a
    log-based normalization that compresses the dynamic range and makes
    small but significant vectors visible.
    
    Args:
        magnitudes: Array of flux magnitudes (W/m²)
        span: Domain span in mm (for scaling)
        config: FluxVectorConfig instance (uses defaults if None)
    
    Returns:
        Array of physical lengths (mm) for each vector
    """
    if config is None:
        config = FluxVectorConfig()
    
    mag_arr = np.array(magnitudes, dtype=float)
    max_vec_len = config.max_vec_frac * span
    min_vec_len = config.min_vec_frac * span
    
    # Handle empty or zero-magnitude arrays
    if mag_arr.size == 0:
        return np.array([], dtype=float)
    
    max_mag = float(np.max(mag_arr)) if mag_arr.size > 0 else 0.0
    if max_mag <= 0.0:
        return np.zeros_like(mag_arr)
    
    # Use high percentile as reference to prevent outliers from dominating
    mag_pos = mag_arr[mag_arr > 0]
    if mag_pos.size > 0:
        mag_ref = float(np.percentile(mag_pos, config.percentile_ref))
        if mag_ref <= 0.0:
            mag_ref = max_mag
    else:
        mag_ref = max_mag
    
    # Log-magnitude scaling compresses dynamic range
    with np.errstate(divide="ignore", invalid="ignore"):
        log_mag = np.log10(1.0 + mag_arr / mag_ref)
    
    log_max = float(np.max(log_mag))
    if log_max > 0.0:
        mag_norm = log_mag / log_max
    else:
        mag_norm = np.zeros_like(mag_arr)
    
    # Apply minimum visibility floor if configured
    if config.min_visible_factor > 0.0:
        mag_norm = np.where(mag_arr > 0, 
                           np.maximum(mag_norm, config.min_visible_factor), 
                           0.0)
    
    # Map normalized magnitudes to [min_vec_len, max_vec_len]
    return min_vec_len + mag_norm * (max_vec_len - min_vec_len)


def build_flux_cones_for_body(nodes, qx_all, qy_all, qz_all, span, config=None):
    """
    Build 3D flux cone vectors for a given body's node list.
    
    Args:
        nodes: List of node dictionaries with 'global_idx' and geometry
        qx_all: Global array of x-direction flux (W/m²)
        qy_all: Global array of y-direction flux (W/m²)
        qz_all: Global array of z-direction flux (W/m²)
        span: Domain span for vector scaling (mm)
        config: FluxVectorConfig instance (uses defaults if None)
    
    Returns:
        Tuple of (x_centers, y_centers, z_centers, u_scaled, v_scaled, w_scaled, magnitudes)
        All as lists (empty if no valid vectors)
    """
    if config is None:
        config = FluxVectorConfig()
    
    if qx_all is None or qy_all is None or qz_all is None:
        return [], [], [], [], [], [], np.array([], dtype=float)
    
    # First pass: collect valid nodes and their magnitudes
    valid_nodes = []
    valid_indices = []
    mags = []
    
    for node in nodes:
        idx = node.get("global_idx", None)
        if idx is None:
            continue
        qx = float(qx_all[idx])
        qy = float(qy_all[idx])
        qz = float(qz_all[idx])
        mag = float(np.sqrt(qx ** 2 + qy ** 2 + qz ** 2))
        if mag <= 0.0 or not np.isfinite(mag):
            continue
        valid_nodes.append(node)
        valid_indices.append(idx)
        mags.append(mag)
    
    if not mags:
        return [], [], [], [], [], [], np.array([], dtype=float)
    
    # Compute scaled lengths for all vectors in this body
    mag_arr = np.array(mags, dtype=float)
    length_mm = compute_vector_lengths(mag_arr, span, config)
    
    # Second pass: build scaled vectors
    x_centers = []
    y_centers = []
    z_centers = []
    u_scaled = []
    v_scaled = []
    w_scaled = []
    
    for node, idx, mag, L in zip(valid_nodes, valid_indices, mag_arr, length_mm):
        if mag <= 0.0:
            continue
        
        # Anchor at node center (use original center for clipped nodes)
        if 'original_x_center' in node:
            x0 = node['original_x_center']
            y0 = node['original_y_center']
            z0 = node['original_z_center']
        else:
            x0 = 0.5 * (node["x_min"] + node["x_max"])
            y0 = 0.5 * (node["y_min"] + node["y_max"])
            z0 = 0.5 * (node["z_min"] + node["z_max"])
        
        # Unit direction scaled by the precomputed length L
        qx = float(qx_all[idx])
        qy = float(qy_all[idx])
        qz = float(qz_all[idx])
        ux = qx / mag
        uy = qy / mag
        uz = qz / mag
        
        x_centers.append(x0)
        y_centers.append(y0)
        z_centers.append(z0)
        u_scaled.append(ux * L)
        v_scaled.append(uy * L)
        w_scaled.append(uz * L)
    
    if not x_centers:
        return [], [], [], [], [], [], np.array([], dtype=float)
    
    return (
        x_centers,
        y_centers,
        z_centers,
        u_scaled,
        v_scaled,
        w_scaled,
        mag_arr,
    )

