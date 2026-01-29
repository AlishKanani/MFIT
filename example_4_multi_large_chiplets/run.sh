#!/bin/bash
# High-fidelity ubump modeling example with 4 chiplets in 2x2 array
# Features fine-grained substrate meshing near ubumps for better thermal accuracy

cd ../
python3 thermal_RC.py \
  --material_prop_file material_prop.yml \
  --geometry_file example_4_high_fidelity_ubumps/geometry.yml \
  --power_config_file example_4_high_fidelity_ubumps/power_cfg.yml \
  --power_seq_file example_4_high_fidelity_ubumps/power_seq.csv \
  --output_dir example_4_high_fidelity_ubumps/ \
  --is_homogeneous false \
  --total_duration 100 \
  --power_interval 1 \
  --time_step 1 \
  --time_heatmap 100 \
  --generate_floorplan true \
  --generate_2d_heatmap false \
  --generate_3d_heatmap true
cd -
