
# This is most basic configuration where all the chiplets are identical. chiplet_geometry_4_chiplets_uniform_nodes.yml that nodes are defined as number of nodes in x and y direction for each layer.

cd ../
python thermal_RC.py --material_prop_file material_prop.yml \
                       --geometry_file example_1_uniform_nodes_homogeneous_chiplets/chiplet_geometry_4_chiplets_uniform_nodes.yml \
                       --power_config_file example_1_uniform_nodes_homogeneous_chiplets/power_dist_config_homogeneous.yml \
                       --power_seq_file example_1_uniform_nodes_homogeneous_chiplets/power_seq_random_4.csv \
                       --output_dir example_1_uniform_nodes_homogeneous_chiplets/
cd -
# Other parameters are set to default values.
# From output floorplan or heatmaps notice that the nodes are uniformly distributed in each layer. Still the nodes can have different granularity by changing `x_nodes`, `y_nodes` in the geometry file.
