"""
Geometric computation utilities for thermal visualization.
"""

from .flux_utils import compute_flux_field
from .flux_3d_vectors import (
    FluxVectorConfig,
    compute_vector_lengths,
    build_flux_cones_for_body
)
from .mesh_builder import (
    collect_body_data,
    filter_nodes_by_clip,
    build_mesh_from_nodes,
    build_edges_from_nodes
)

__all__ = [
    'compute_flux_field',
    'FluxVectorConfig',
    'compute_vector_lengths',
    'build_flux_cones_for_body',
    'collect_body_data',
    'filter_nodes_by_clip',
    'build_mesh_from_nodes',
    'build_edges_from_nodes'
]
