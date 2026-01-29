#!/usr/bin/env python3
"""
MFIT system.yml generator

Reads a compact system.yml describing a multi-chip package and generates standard MFIT input files:
- geometry.generated.yml
- power_cfg.generated.yml
- power_seq.generated.csv

Features:
- Arbitrary chip placement (use MFIT --is_homogeneous false)
- Uniform meshes (integer x_nodes/y_nodes)
- Material regions for heterogeneous substrates (e.g., embedded silicon links)
"""

from __future__ import annotations

import argparse
import csv
import os
from dataclasses import dataclass
from typing import Any, Dict, List, Tuple

import yaml


def _load_yaml(path: str) -> Dict[str, Any]:
    with open(path, "r") as f:
        data = yaml.safe_load(f)
    if not isinstance(data, dict):
        raise ValueError(f"Expected YAML mapping at top-level: {path}")
    return data


def _ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def _req(d: Dict[str, Any], key: str, where: str) -> Any:
    if key not in d:
        raise KeyError(f"Missing required key '{key}' in {where}")
    return d[key]


def _as_float(x: Any, where: str) -> float:
    try:
        return float(x)
    except Exception as e:
        raise ValueError(f"Expected float for {where}, got {x!r}") from e


def _as_int(x: Any, where: str) -> int:
    try:
        v = int(x)
    except Exception as e:
        raise ValueError(f"Expected int for {where}, got {x!r}") from e
    if v <= 0:
        raise ValueError(f"Expected positive int for {where}, got {v}")
    return v


def _region_rect(region: Dict[str, Any], where: str) -> Tuple[float, float, float, float]:
    sx = _as_float(region.get("start_x", region.get("x", 0.0)), f"{where}.start_x")
    sy = _as_float(region.get("start_y", region.get("y", 0.0)), f"{where}.start_y")
    lx_raw = region.get("length_x", region.get("lx", None))
    ly_raw = region.get("length_y", region.get("ly", None))
    if lx_raw is None or ly_raw is None:
        raise KeyError(f"{where} must define length_x/length_y (or lx/ly)")
    lx = _as_float(lx_raw, f"{where}.length_x")
    ly = _as_float(ly_raw, f"{where}.length_y")
    if lx <= 0 or ly <= 0:
        raise ValueError(f"{where} length_x/length_y must be positive")
    return sx, sy, lx, ly


def _rects_overlap(a: Tuple[float, float, float, float], b: Tuple[float, float, float, float]) -> bool:
    ax, ay, alx, aly = a
    bx, by, blx, bly = b
    return (ax < bx + blx) and (ax + alx > bx) and (ay < by + bly) and (ay + aly > by)


def _validate_regions_no_overlap(regions: List[Dict[str, Any]], where: str) -> None:
    rects = []
    for i, region in enumerate(regions):
        if not isinstance(region, dict):
            raise ValueError(f"{where}[{i}] must be a mapping")
        rects.append(_region_rect(region, f"{where}[{i}]"))
    for i in range(len(rects)):
        for j in range(i + 1, len(rects)):
            if _rects_overlap(rects[i], rects[j]):
                name_i = regions[i].get("name", f"{where}[{i}]")
                name_j = regions[j].get("name", f"{where}[{j}]")
                raise ValueError(f"Overlapping regions are not allowed: {name_i} overlaps {name_j} in {where}")


@dataclass(frozen=True)
class Rect:
    x: float
    y: float
    lx: float
    ly: float

    @staticmethod
    def from_dict(d: Dict[str, Any], where: str) -> "Rect":
        return Rect(
            x=_as_float(_req(d, "x", where), f"{where}.x"),
            y=_as_float(_req(d, "y", where), f"{where}.y"),
            lx=_as_float(_req(d, "lx", where), f"{where}.lx"),
            ly=_as_float(_req(d, "ly", where), f"{where}.ly"),
        )


@dataclass(frozen=True)
class Chip:
    name: str
    origin: Tuple[float, float]
    size: Tuple[float, float]


def _parse_chips(system: Dict[str, Any]) -> List[Chip]:
    chips_raw = _req(system, "chips", "system")
    if not isinstance(chips_raw, list) or not chips_raw:
        raise ValueError("system.chips must be a non-empty list")

    seen = set()
    chips: List[Chip] = []
    for i, c in enumerate(chips_raw):
        where = f"system.chips[{i}]"
        if not isinstance(c, dict):
            raise ValueError(f"{where} must be a mapping")
        name = str(_req(c, "name", where))
        if name in seen:
            raise ValueError(f"Duplicate chip name: {name}")
        seen.add(name)
        origin_d = _req(c, "origin", where)
        size_d = _req(c, "size", where)
        if not isinstance(origin_d, dict) or not isinstance(size_d, dict):
            raise ValueError(f"{where}.origin and {where}.size must be mappings")
        ox = _as_float(_req(origin_d, "x", f"{where}.origin"), f"{where}.origin.x")
        oy = _as_float(_req(origin_d, "y", f"{where}.origin"), f"{where}.origin.y")
        sx = _as_float(_req(size_d, "x", f"{where}.size"), f"{where}.size.x")
        sy = _as_float(_req(size_d, "y", f"{where}.size"), f"{where}.size.y")
        chips.append(Chip(name=name, origin=(ox, oy), size=(sx, sy)))
    return chips


def _parse_power_traces(system: Dict[str, Any]) -> Dict[str, List[float]]:
    traces_raw = system.get("power_traces", {}) or {}
    if not isinstance(traces_raw, dict):
        raise ValueError("system.power_traces must be a mapping trace_name -> [percent,...]")
    traces: Dict[str, List[float]] = {}
    for name, vals in traces_raw.items():
        if not isinstance(vals, list) or not vals:
            raise ValueError(f"power_traces.{name} must be a non-empty list")
        traces[str(name)] = [_as_float(v, f"power_traces.{name}[]") for v in vals]
    return traces


@dataclass(frozen=True)
class PowerBlock:
    chip: str
    layer: str
    name: str
    rect: Rect
    max_power_w: float
    trace: str


def _parse_power_blocks(system: Dict[str, Any]) -> List[PowerBlock]:
    blocks_raw = system.get("power_blocks", []) or []
    if not isinstance(blocks_raw, list):
        raise ValueError("system.power_blocks must be a list")
    blocks: List[PowerBlock] = []
    for i, b in enumerate(blocks_raw):
        where = f"system.power_blocks[{i}]"
        if not isinstance(b, dict):
            raise ValueError(f"{where} must be a mapping")
        chip = str(_req(b, "chip", where))
        layer = str(_req(b, "layer", where))
        name = str(_req(b, "name", where))
        rect = Rect.from_dict(_req(b, "rect", where), f"{where}.rect")
        mp = _as_float(_req(b, "max_power_w", where), f"{where}.max_power_w")
        trace = str(_req(b, "trace", where))
        blocks.append(PowerBlock(chip=chip, layer=layer, name=name, rect=rect, max_power_w=mp, trace=trace))
    return blocks


def _stack_thickness(layers: List[Dict[str, Any]]) -> float:
    return sum(_as_float(_req(L, "thickness", "stack layer"), "stack layer thickness") for L in layers)


def _emit_power_seq_csv(path: str, blocks: List[PowerBlock], traces: Dict[str, List[float]]) -> None:
    # Ensure at least one trace exists
    if not traces:
        traces = {"steady_100": [100.0]}

    max_len = max(len(v) for v in traces.values())

    def trace_values(name: str) -> List[float]:
        if name not in traces:
            raise KeyError(f"Power block references unknown trace: {name}")
        vals = list(traces[name])
        if len(vals) < max_len:
            vals.extend([vals[-1]] * (max_len - len(vals)))
        return vals

    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        for blk in blocks:
            row_name = f"{blk.layer}_{blk.chip}_{blk.name}"
            w.writerow([row_name] + trace_values(blk.trace))


def build(system: Dict[str, Any], outdir: str) -> Tuple[str, str, str]:
    version = int(system.get("version", 1))
    if version != 1:
        raise ValueError(f"Unsupported system.yml version: {version}")

    _ensure_dir(outdir)

    package = _req(system, "package", "system")
    if not isinstance(package, dict):
        raise ValueError("system.package must be a mapping")

    stack = _req(system, "stack", "system")
    if not isinstance(stack, dict):
        raise ValueError("system.stack must be a mapping")

    chips = _parse_chips(system)
    traces = _parse_power_traces(system)
    pblocks = _parse_power_blocks(system)

    chip_by_name = {c.name: c for c in chips}
    for b in pblocks:
        if b.chip not in chip_by_name:
            raise ValueError(f"power_blocks references unknown chip: {b.chip}")

    # --- Stack parsing ---
    substrate_cfg = _req(stack, "substrate", "stack")
    chip_cfg = _req(stack, "chip", "stack")
    tim_cfg = _req(stack, "tim", "stack")
    lid_cfg = _req(stack, "lid", "stack")
    ubump_cfg = stack.get("ubump", None)

    # substrate (global) - support multi-layer substrates with regions
    substrate_layers: List[Dict[str, Any]] = []
    if isinstance(substrate_cfg, dict) and substrate_cfg.get("layers", None) is not None:
        sub_layers_raw = substrate_cfg.get("layers", [])
        if not isinstance(sub_layers_raw, list) or not sub_layers_raw:
            raise ValueError("stack.substrate.layers must be a non-empty list when present")
        for i, L in enumerate(sub_layers_raw):
            where = f"stack.substrate.layers[{i}]"
            if not isinstance(L, dict):
                raise ValueError(f"{where} must be a mapping")
            lname = str(_req(L, "name", where))
            th = _as_float(_req(L, "thickness", where), f"{where}.thickness")
            mat = str(_req(L, "material", where))
            nd = _req(L, "nodes", where)
            xn = _as_int(_req(nd, "x_nodes", f"{where}.nodes"), f"{where}.nodes.x_nodes")
            yn = _as_int(_req(nd, "y_nodes", f"{where}.nodes"), f"{where}.nodes.y_nodes")
            regions = L.get("regions", []) or []
            if regions is not None and not isinstance(regions, list):
                raise ValueError(f"{where}.regions must be a list when present")
            if regions:
                _validate_regions_no_overlap(regions, f"{where}.regions")
            substrate_layers.append(
                {
                    "name": lname,
                    "thickness": float(th),
                    "material": mat,
                    "x_nodes": int(xn),
                    "y_nodes": int(yn),
                    "regions": regions,
                }
            )
    else:
        # Legacy single-layer substrate
        sub_th = _as_float(_req(substrate_cfg, "thickness", "stack.substrate"), "stack.substrate.thickness")
        sub_mat = str(_req(substrate_cfg, "material", "stack.substrate"))
        sub_nodes = _req(substrate_cfg, "nodes", "stack.substrate")
        sub_xn = _as_int(_req(sub_nodes, "x_nodes", "stack.substrate.nodes"), "stack.substrate.nodes.x_nodes")
        sub_yn = _as_int(_req(sub_nodes, "y_nodes", "stack.substrate.nodes"), "stack.substrate.nodes.y_nodes")
        regions = substrate_cfg.get("regions", []) or []
        if regions is not None and not isinstance(regions, list):
            raise ValueError("stack.substrate.regions must be a list when present")
        if regions:
            _validate_regions_no_overlap(regions, "stack.substrate.regions")
        substrate_layers = [
            {
                "name": "substrate",
                "thickness": float(sub_th),
                "material": sub_mat,
                "x_nodes": int(sub_xn),
                "y_nodes": int(sub_yn),
                "regions": regions,
            }
        ]

    # chip stack (per chip)
    chip_layers = _req(chip_cfg, "layers", "stack.chip")
    if not isinstance(chip_layers, dict):
        raise ValueError("stack.chip.layers must be a mapping")
    chip_count = _as_int(_req(chip_layers, "count", "stack.chip.layers"), "stack.chip.layers.count")
    chip_prefix = str(_req(chip_layers, "layer_name_prefix", "stack.chip.layers"))
    chip_th_each = _as_float(_req(chip_layers, "thickness_each", "stack.chip.layers"), "stack.chip.layers.thickness_each")
    chip_mat = str(_req(chip_layers, "material", "stack.chip.layers"))
    chip_nodes = _req(chip_layers, "nodes", "stack.chip.layers")
    chip_xn = _as_int(_req(chip_nodes, "x_nodes", "stack.chip.layers.nodes"), "stack.chip.layers.nodes.x_nodes")
    chip_yn = _as_int(_req(chip_nodes, "y_nodes", "stack.chip.layers.nodes"), "stack.chip.layers.nodes.y_nodes")
    power_src_layers = _req(chip_layers, "power_src_layers", "stack.chip.layers")
    if not isinstance(power_src_layers, list) or not power_src_layers:
        raise ValueError("stack.chip.layers.power_src_layers must be a non-empty list of 1-based indices")
    power_src_layers_i = set(int(x) for x in power_src_layers)

    # tim layers (per chip)
    tim_layers = _req(tim_cfg, "layers", "stack.tim")
    if not isinstance(tim_layers, list) or not tim_layers:
        raise ValueError("stack.tim.layers must be a non-empty list")
    for L in tim_layers:
        if not isinstance(L, dict):
            raise ValueError("stack.tim.layers entries must be mappings")

    # ubump layers (per chip) - optional
    ubump_layers = []
    if ubump_cfg is not None:
        if not isinstance(ubump_cfg, dict):
            raise ValueError("stack.ubump must be a mapping when present")
        ubump_layers = ubump_cfg.get("layers", []) or []
        if not isinstance(ubump_layers, list) or not ubump_layers:
            raise ValueError("stack.ubump.layers must be a non-empty list when stack.ubump is present")
        for L in ubump_layers:
            if not isinstance(L, dict):
                raise ValueError("stack.ubump.layers entries must be mappings")

    # lid layers (global)
    lid_layers = _req(lid_cfg, "layers", "stack.lid")
    if not isinstance(lid_layers, list) or not lid_layers:
        raise ValueError("stack.lid.layers must be a non-empty list")
    for L in lid_layers:
        if not isinstance(L, dict):
            raise ValueError("stack.lid.layers entries must be mappings")

    # --- Derive package Z ---
    chip_total_th = chip_count * chip_th_each
    sub_th = sum(float(L["thickness"]) for L in substrate_layers)
    ubump_total_th = _stack_thickness(ubump_layers) if ubump_layers else 0.0
    tim_total_th = _stack_thickness(tim_layers)
    lid_total_th = _stack_thickness(lid_layers)
    z_length = sub_th + ubump_total_th + chip_total_th + tim_total_th + lid_total_th

    # --- Geometry generation ---
    geom: Dict[str, Any] = {}
    geom["common"] = {
        "x_length": float(_req(package, "x_length", "package")),
        "y_length": float(_req(package, "y_length", "package")),
        "z_length": float(z_length),
        # Required by MFIT's Utils even in non-homogeneous mode; set to first chip size.
        "chiplet_x": float(chips[0].size[0]),
        "chiplet_y": float(chips[0].size[1]),
        "chiplet_spacing": 0.0,
        "n_chiplet_x": 1,
        "n_chiplet_y": 1,
        "bc_top_htc": float(_req(package, "bc_top_htc", "package")),
        "bc_bottom_htc": float(_req(package, "bc_bottom_htc", "package")),
        "ambient_temp": float(package.get("ambient_temp", 300.0)),
    }

    layers_dict: Dict[str, Any] = {}

    # Substrate layers (global, not under chiplet)
    z_sub = 0.0
    for L in substrate_layers:
        lname = str(L["name"])
        layer_entry: Dict[str, Any] = {
            "thickness": float(L["thickness"]),
            "nodes": {
                "uniform": True,
                "under_chiplet": False,
                "x_nodes": int(L["x_nodes"]),
                "y_nodes": int(L["y_nodes"])
            },
            "start_point": {"x": 0.0, "y": 0.0, "z": float(z_sub)},
            "power_src": False,
            "material": str(L["material"]),
        }
        regions = L.get("regions", []) or []
        if regions:
            layer_entry["regions"] = regions
        layers_dict[lname] = layer_entry
        z_sub += float(L["thickness"])

    z = sub_th

    # Ubump layers (per-chip footprint)
    if ubump_layers:
        for L in ubump_layers:
            lname = str(_req(L, "name", "stack.ubump.layers[]"))
            th = _as_float(_req(L, "thickness", f"stack.ubump.layers[{lname}]"), f"stack.ubump.layers[{lname}].thickness")
            mat = str(_req(L, "material", f"stack.ubump.layers[{lname}]"))
            nd = _req(L, "nodes", f"stack.ubump.layers[{lname}]")
            xn = _as_int(_req(nd, "x_nodes", f"stack.ubump.layers[{lname}].nodes"), f"stack.ubump.layers[{lname}].nodes.x_nodes")
            yn = _as_int(_req(nd, "y_nodes", f"stack.ubump.layers[{lname}].nodes"), f"stack.ubump.layers[{lname}].nodes.y_nodes")
            layers_dict[lname] = {
                "thickness": float(th),
                "nodes": {"uniform": True, "under_chiplet": True, "x_nodes": int(xn), "y_nodes": int(yn)},
                "start_point": {"x": 0.0, "y": 0.0, "z": float(z)},
                "power_src": False,
                "material": mat,
            }
            z += th

    # Chiplet layers (per-chip footprint)
    for i in range(1, chip_count + 1):
        lname = f"{chip_prefix}{i}"
        layers_dict[lname] = {
            "thickness": float(chip_th_each),
            "nodes": {"uniform": True, "under_chiplet": True, "x_nodes": int(chip_xn), "y_nodes": int(chip_yn)},
            "start_point": {"x": 0.0, "y": 0.0, "z": float(z)},
            "power_src": (i in power_src_layers_i),
            "material": chip_mat,
        }
        z += chip_th_each

    # TIM layers (per-chip footprint)
    for L in tim_layers:
        lname = str(_req(L, "name", "stack.tim.layers[]"))
        th = _as_float(_req(L, "thickness", f"stack.tim.layers[{lname}]"), f"stack.tim.layers[{lname}].thickness")
        mat = str(_req(L, "material", f"stack.tim.layers[{lname}]"))
        nd = _req(L, "nodes", f"stack.tim.layers[{lname}]")
        xn = _as_int(_req(nd, "x_nodes", f"stack.tim.layers[{lname}].nodes"), f"stack.tim.layers[{lname}].nodes.x_nodes")
        yn = _as_int(_req(nd, "y_nodes", f"stack.tim.layers[{lname}].nodes"), f"stack.tim.layers[{lname}].nodes.y_nodes")
        layers_dict[lname] = {
            "thickness": float(th),
            "nodes": {"uniform": True, "under_chiplet": True, "x_nodes": int(xn), "y_nodes": int(yn)},
            "start_point": {"x": 0.0, "y": 0.0, "z": float(z)},
            "power_src": False,
            "material": mat,
        }
        z += th

    # Lid layers (global, not under chiplet)
    for L in lid_layers:
        lname = str(_req(L, "name", "stack.lid.layers[]"))
        th = _as_float(_req(L, "thickness", f"stack.lid.layers[{lname}]"), f"stack.lid.layers[{lname}].thickness")
        mat = str(_req(L, "material", f"stack.lid.layers[{lname}]"))
        nd = _req(L, "nodes", f"stack.lid.layers[{lname}]")
        xn = _as_int(_req(nd, "x_nodes", f"stack.lid.layers[{lname}].nodes"), f"stack.lid.layers[{lname}].nodes.x_nodes")
        yn = _as_int(_req(nd, "y_nodes", f"stack.lid.layers[{lname}].nodes"), f"stack.lid.layers[{lname}].nodes.y_nodes")
        layers_dict[lname] = {
            "thickness": float(th),
            "nodes": {"uniform": True, "under_chiplet": False, "x_nodes": int(xn), "y_nodes": int(yn)},
            "start_point": {"x": 0.0, "y": 0.0, "z": float(z)},
            "power_src": False,
            "material": mat,
        }
        z += th

    geom["layers"] = layers_dict

    geom_path = os.path.join(outdir, "geometry.generated.yml")
    with open(geom_path, "w") as f:
        yaml.safe_dump(geom, f, sort_keys=False)

    # --- Power config generation ---
    pwr: Dict[str, Any] = {}

    # Under-chiplet layers that must be present in power_cfg for non-homogeneous node generation
    under_layers: List[Tuple[str, int, int, bool]] = []
    if ubump_layers:
        for L in ubump_layers:
            lname = str(_req(L, "name", "stack.ubump.layers[]"))
            nd = _req(L, "nodes", f"stack.ubump.layers[{lname}]")
            xn = _as_int(_req(nd, "x_nodes", f"stack.ubump.layers[{lname}].nodes"), f"stack.ubump.layers[{lname}].nodes.x_nodes")
            yn = _as_int(_req(nd, "y_nodes", f"stack.ubump.layers[{lname}].nodes"), f"stack.ubump.layers[{lname}].nodes.y_nodes")
            under_layers.append((lname, xn, yn, False))
    for i in range(1, chip_count + 1):
        under_layers.append((f"{chip_prefix}{i}", chip_xn, chip_yn, (i in power_src_layers_i)))
    for L in tim_layers:
        lname = str(_req(L, "name", "stack.tim.layers[]"))
        nd = _req(L, "nodes", f"stack.tim.layers[{lname}]")
        xn = _as_int(_req(nd, "x_nodes", f"stack.tim.layers[{lname}].nodes"), f"stack.tim.layers[{lname}].nodes.x_nodes")
        yn = _as_int(_req(nd, "y_nodes", f"stack.tim.layers[{lname}].nodes"), f"stack.tim.layers[{lname}].nodes.y_nodes")
        under_layers.append((lname, xn, yn, False))

    # Index power blocks by (layer, chip)
    blocks_by_layer_chip: Dict[Tuple[str, str], List[PowerBlock]] = {}
    for b in pblocks:
        blocks_by_layer_chip.setdefault((b.layer, b.chip), []).append(b)

    for (layer_name, nodes_x, nodes_y, _) in under_layers:
        pwr[layer_name] = {}
        for chip in chips:
            entry: Dict[str, Any] = {
                "start_chiplet_x": float(chip.origin[0]),
                "start_chiplet_y": float(chip.origin[1]),
                "length_chiplet_x": float(chip.size[0]),
                "length_chiplet_y": float(chip.size[1]),
                "nodes_x": int(nodes_x),
                "nodes_y": int(nodes_y),
            }
            blks = blocks_by_layer_chip.get((layer_name, chip.name), [])
            if blks:
                entry["layout_blocks"] = {
                    blk.name: {
                        "start_point_x": float(blk.rect.x),
                        "start_point_y": float(blk.rect.y),
                        "length_x": float(blk.rect.lx),
                        "length_y": float(blk.rect.ly),
                        "max_power": float(blk.max_power_w),
                    }
                    for blk in blks
                }
            pwr[layer_name][chip.name] = entry

    pwr_path = os.path.join(outdir, "power_cfg.generated.yml")
    with open(pwr_path, "w") as f:
        yaml.safe_dump(pwr, f, sort_keys=False)

    # --- Power sequence generation ---
    pseq_path = os.path.join(outdir, "power_seq.generated.csv")
    _emit_power_seq_csv(pseq_path, pblocks, traces)

    return geom_path, pwr_path, pseq_path


def main() -> None:
    ap = argparse.ArgumentParser(description="Generate MFIT inputs from system.yml.")
    ap.add_argument("--system", required=True, help="Path to system.yml")
    ap.add_argument("--outdir", required=True, help="Output directory for generated MFIT inputs")
    args = ap.parse_args()

    system = _load_yaml(args.system)
    geom_path, pwr_path, pseq_path = build(system, args.outdir)
    print("Generated:")
    print(f"  - {geom_path}")
    print(f"  - {pwr_path}")
    print(f"  - {pseq_path}")
    print("")
    print("Run MFIT with:")
    print("  python3 thermal_RC.py --is_homogeneous false \\")
    print(f"    --geometry_file {geom_path} --power_config_file {pwr_path} --power_seq_file {pseq_path} ...")


if __name__ == "__main__":
    main()
