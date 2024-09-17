#   In this example, we have 3 chiplets, one bigger and two smaller ones. 
  
#   When `is_homogeneous` is set to False, the model will ignore chiplet-specific parameters from `geometry_file` and instead use the parameters from `power_config_file` for each chiplet. 
  
#   For each layer marked with `under_chiplet: True` in the `geometry_file`, individual blocks need to be defined in the `power_config_file`. In this case, the power sequence should be defined only for blocks with the `layout_blocks:` dictionary.

#   Also, the material properties can be overridden for each block using the `material:` key in the `power_config_file`.

#   Check the geometry file chiplet_geometry_3_chiplets_uniform_nodes.yml and power configuration file power_dist_config_heterogeneous.yml for more details.


cd ../
python thermal_RC.py --material_prop_file material_prop.yml \
                       --geometry_file example_3_heterogeneous_chiplets/chiplet_geometry_3_chiplets_uniform_nodes.yml \
                       --power_config_file example_3_heterogeneous_chiplets/power_dist_config_heterogeneous.yml \
                       --power_seq_file example_3_heterogeneous_chiplets/power_seq_random_3.csv \
                       --output_dir example_3_heterogeneous_chiplets/ \
                       --is_homogeneous false
cd -