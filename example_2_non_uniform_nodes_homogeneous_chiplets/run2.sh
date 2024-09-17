#   This configuration is similar to the previous one but nodes are defined as x and y coordinates for some layers. You can also mix the uniform and non-uniform for different layers, checkout the geometry file chiplet_geometry_4_chiplets_non_uniform.yml. Also notice the power configuration file example_2_non_uniform_nodes_homogeneous_chiplets that power sources are defined as small blocks in the chiplets.
cd ../  
python thermal_RC.py --material_prop_file material_prop.yml \
                       --geometry_file example_2_non_uniform_nodes_homogeneous_chiplets/chiplet_geometry_4_chiplets_non_uniform.yml \
                       --power_config_file example_2_non_uniform_nodes_homogeneous_chiplets/power_dist_config_tiles_tx_rx.yml \
                       --power_seq_file example_2_non_uniform_nodes_homogeneous_chiplets/power_seq_30s_tiles_tx_rs.csv \
                       --output_dir example_2_non_uniform_nodes_homogeneous_chiplets/ \
                       --total_duration 30
cd -

#   Other parameters are set to default values. 

#   The average temperature of each chiplet is not calculated in this case as the nodes are not uniformly distributed. But read node temperature from the `temperature_all_<time_step>.csv` file. Node numberings are provided in the floorplan images.