"""
Case stats writer for MFIT runs.

This is invoked by the simulator (thermal_RC.py) to write a compact summary of
what was run (inputs/args), what was generated (node counts), and how long it took.
"""

from __future__ import annotations

import os
import platform
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any, Dict, Optional

import yaml


@dataclass(frozen=True)
class CaseStatsInputs:
    elapsed_wall_s: float
    simulated_steps: int
    dt_s: float
    timing: Optional[Dict[str, float]] = None


def _iso_utc_now() -> str:
    return datetime.now(tz=timezone.utc).isoformat()


def _layer_nodes_summary(layer: Any) -> Dict[str, Any]:
    nodes_dict = getattr(layer, "nodes_dict", {}) or {}
    return {
        "layer": getattr(layer, "layer_name", None),
        "total_nodes": int(layer.layer_total_nodes()),
        "under_chiplet": bool(nodes_dict.get("under_chiplet", False)),
        "x_nodes": nodes_dict.get("x_nodes", None),
        "y_nodes": nodes_dict.get("y_nodes", None),
        "material": getattr(layer, "material", None),
        "thickness_mm": getattr(layer, "thickness", None),
    }


def write_case_stats(
    *,
    output_dir: str,
    args: Any,
    package: Any,
    geometry_dict: Optional[Dict[str, Any]],
    inputs: CaseStatsInputs,
) -> str:
    """
    Write case_stats.yml into output_dir and return the written path.

    This function should be best-effort: call it from the simulator and catch exceptions
    so failures here never fail the simulation.
    """
    os.makedirs(output_dir, exist_ok=True)

    total_nodes = int(package.package_total_nodes())
    layers = list(getattr(package, "layers", []) or [])
    layer_summaries = [_layer_nodes_summary(L) for L in layers]

    power_steps = None
    try:
        power_steps = int(getattr(package, "power").shape[1])
    except Exception:
        power_steps = None

    saved_steps = None
    try:
        saved_steps = int(len(getattr(package, "temperature_all_save")))
    except Exception:
        saved_steps = None

    stats: Dict[str, Any] = {
        "case": {
            "written_utc": _iso_utc_now(),
            "output_dir": os.path.abspath(output_dir),
            "inputs": {
                "material_prop_file": getattr(args, "material_prop_file", None),
                "geometry_file": getattr(args, "geometry_file", None),
                "power_config_file": getattr(args, "power_config_file", None),
                "power_seq_file": getattr(args, "power_seq_file", None),
            },
        },
        "execution": {
            "elapsed_wall_s": float(inputs.elapsed_wall_s),
            "timing_s": inputs.timing or {},
            "python": {
                "version": sys.version.split()[0],
                "executable": sys.executable,
            },
            "platform": {
                "system": platform.system(),
                "release": platform.release(),
                "machine": platform.machine(),
            },
        },
        "simulation": {
            "simulation_type": getattr(args, "simulation_type", None),
            "total_duration_s": float(getattr(args, "total_duration", 0.0)),
            "time_step_s": float(getattr(args, "time_step", 0.0)),
            "power_interval_s": float(getattr(args, "power_interval", 0.0)),
            "requested_steps": int(inputs.simulated_steps),
            "simulated_time_s": float(inputs.simulated_steps * inputs.dt_s),
            "saved_states": saved_steps,  # typically requested_steps + 1 (includes t=0)
            "power_sequence_points": power_steps,  # columns in power sequence after resampling
        },
        "solver": {
            "sparse_assembly": True,  # MFIT always uses sparse
            "generate_DSS": bool(getattr(args, "generate_DSS", False)),
            "use_tuned_C": bool(getattr(args, "use_tuned_C", False)),
        },
        "outputs": {
            "generate_floorplan": bool(getattr(args, "generate_floorplan", False)),
            "generate_2d_heatmap": bool(getattr(args, "generate_2d_heatmap", False)),
            "time_heatmap_s": float(getattr(args, "time_heatmap", 0.0)),
            "vertical_planes": getattr(args, "vertical_planes", ""),
            "xz_cuts": getattr(args, "xz_cuts", ""),
            "yz_cuts": getattr(args, "yz_cuts", ""),
            "interactive_heatmaps": bool(getattr(args, "interactive_heatmaps", False)),
            "generate_3d_heatmap": bool(getattr(args, "generate_3d_heatmap", False)),
        },
        "model": {
            "total_nodes": total_nodes,
            "n_layers": int(len(layers)),
            "layers": layer_summaries,
        },
    }

    out_path = os.path.join(output_dir, "case_stats.yml")
    with open(out_path, "w") as f:
        yaml.safe_dump(stats, f, sort_keys=False)
    return out_path
