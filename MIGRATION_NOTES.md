# MFIT Output Directory Structure Migration

## Summary

Ported the cleaner output directory structure and case statistics feature from CTEC/Thermoelectric-Modeling to MFIT.

## Changes Made

### 1. New Files Created

- **`tools/case_stats.py`**: Simplified version of CTEC's case stats writer
  - Generates comprehensive YAML summary of simulation run
  - Includes: case metadata, execution timing, simulation parameters, solver config, output flags, model summary
  - Removed CTEC-specific features: TEC metrics, body membership details
  
- **`tools/__init__.py`**: Module initialization for tools package

### 2. Modified Files

#### `thermal_RC.py`
- Added import for `case_stats` module
- Creates new output directory structure on startup:
  ```python
  output/
  ├── RC/     # Thermal RC simulation results
  └── DSS/    # Discrete state-space matrices
  ```
- Collects detailed timing breakdown for all phases
- Calls `write_case_stats()` at end of simulation to generate `case_stats.yml`

#### `package_class.py`
- Updated `generate_DSS()` to save matrices to `output/DSS/` directory
- Updated `write_temperature_to_file()` to save results to `output/RC/` directory
- Updated `set_initial_conditions()` to save power data to `output/DSS/` directory
- All methods use `getattr()` with fallback for backward compatibility

#### `layer_class.py`
- Updated `write_temperature_to_file()` to save per-layer results to `output/RC/` directory
- Updated `generate_floorplan_visual()` to save floorplan images to `output/RC/floorplan/` directory

### 3. New Output Directory Structure

**Before:**
```
<output_dir>/
├── output/
│   ├── temperature_all_*.csv
│   ├── temperature_<layer>_*.csv
│   ├── disc_A_matrix.csv
│   ├── disc_B_matrix.csv
│   └── power_all.csv
├── floorplan/
└── heatmaps/
```

**After:**
```
<output_dir>/
├── output/
│   ├── RC/                         # ← NEW: RC simulation results
│   │   ├── case_stats.yml          # ← NEW: Comprehensive case summary
│   │   ├── temperature_all_*.csv
│   │   ├── temperature_<layer>_*.csv
│   │   └── floorplan/              # ← MOVED: Floorplan images
│   │       ├── <layer_name>.png
│   │       └── <power_source>_power_.png
│   └── DSS/                        # ← NEW: DSS matrices
│       ├── disc_A_matrix.csv
│       ├── disc_B_matrix.csv
│       └── power_all.csv
└── heatmaps/
```

### 4. Case Stats Contents

The `case_stats.yml` file includes:

- **case**: Timestamp, output directory path, input file paths
- **execution**: Wall time, detailed timing breakdown, Python version, platform info
- **simulation**: Type, duration, time step, power interval, step counts
- **solver**: Sparse assembly flag, DSS generation flag, tuned C flag
- **outputs**: All visualization flags (floorplan, heatmaps, vertical planes, 3D)
- **model**: Total nodes, layer count, per-layer details (nodes, material, thickness)

## Benefits

1. **Better Organization**: Separates RC and DSS outputs into logical subdirectories
2. **Reproducibility**: Case stats file captures all simulation parameters and timing
3. **Debugging**: Comprehensive timing breakdown helps identify performance bottlenecks
4. **Documentation**: Self-documenting runs with all metadata in one YAML file
5. **Consistency**: Matches CTEC simulator's output structure

#### `power_class.py`
- Updated `generate_floorplan_visual()` to save power floorplan images to `output/RC/floorplan/` directory

### 3. CLI Argument Changes (Breaking Change)

**Renamed for clarity:**
- `--generate_heatmap` → `--generate_2d_heatmap` (makes it clear this is for 2D layer heatmaps only)

This is **NOT backward compatible** - old scripts using `--generate_heatmap` must be updated.

**Updated files:**
- All example run scripts (`run.sh`, `run3.sh`, etc.)
- README.md documentation

## Backward Compatibility

- Code uses `getattr()` with fallbacks to handle older code that doesn't set `output_rc_dir`/`output_dss_dir`
- If those attributes are missing, falls back to old `output_dir/output/RC` path
- Heatmap outputs remain in their original `heatmaps/` location
- **BREAKING:** CLI argument `--generate_heatmap` renamed to `--generate_2d_heatmap` (must update run scripts)

## Testing

The case stats functionality has been tested and verified to:
- ✓ Generate valid YAML output
- ✓ Capture all timing metrics
- ✓ Include model summary with layer details
- ✓ Record all input/output parameters

## Usage

No changes required to run scripts. The new directory structure and case stats are created automatically on every run.

Example case stats location:
```
./example_3_heterogeneous_chiplets/output/RC/case_stats.yml
```
