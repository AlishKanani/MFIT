#!/bin/bash
# AMD-style multi-chiplet example with embedded silicon links
# Uses build_system.py to generate configuration from compact system.yml

set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
mfroot="$(cd "$here/.." && pwd)"

outdir="$here/generated"
mkdir -p "$outdir"

# # Generate configuration files from system.yml
# python3 "$mfroot/tools/build_system.py" \
#   --system "$here/system.yml" \
#   --outdir "$outdir"

# Run simulation
python3 "$mfroot/thermal_RC.py" \
  --material_prop_file "$mfroot/material_prop.yml" \
  --geometry_file "$outdir/geometry.generated.yml" \
  --power_config_file "$outdir/power_cfg.generated.yml" \
  --power_seq_file "$outdir/power_seq.generated.csv" \
  --output_dir "$here" \
  --is_homogeneous false \
  --total_duration 100 \
  --power_interval 1 \
  --time_step 1 \
  --time_heatmap 100 \
  --generate_2d_floorplan false \
  --generate_3d_floorplan false \
  --generate_2d_heatmap false \
  --vertical_planes "XZ" \
  --xz_cuts "15.0" \
  --interactive_heatmaps true \
  --generate_3d_heatmap false
