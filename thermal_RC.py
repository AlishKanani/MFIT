import common
from power_class import Power_grid
from package_class import Chiplet_package
import time
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--power_config_file', type=str, default='example_3_heterogeneous_chiplets/power_dist_config_heterogeneous.yml', help='Power distribution configuration file')
    parser.add_argument('--power_seq_file', type=str, default='example_3_heterogeneous_chiplets/power_seq_random_3.csv', help='Power sequence file')
    parser.add_argument('--material_prop_file', type=str, default='material_prop.yml', help='Material properties file')
    parser.add_argument('--geometry_file', type=str, default='example_3_heterogeneous_chiplets/chiplet_geometry_3_chiplets_uniform_nodes.yml', help='Geometry properties file')
    parser.add_argument('--output_dir', type=str, default='./example_3_heterogeneous_chiplets/', help='Output directory')

    parser.add_argument('--simulation_type', type=str, default='transient', choices=['transient', 'steady'], help='transient or steady state simulation')
    parser.add_argument('--generate_DSS', type=lambda x: (str(x).lower() in ['true','1', 'yes']), default=True, help='Generate A and B for DSS')

    parser.add_argument('--is_homogeneous', type=lambda x: (str(x).lower() in ['true','1', 'yes']), default=True, help='Are chiplet placement homogeneous?')
    parser.add_argument('--time_step', type=float, default=0.1, help='Time step for transient simulation in sec')
    parser.add_argument('--power_interval', type=float, default=1, help='Power interval for transient simulation in sec')
    parser.add_argument('--total_duration', type=float, default=50, help='Total time for transient simulation in sec')

    parser.add_argument('--use_tuned_C', type=lambda x: (str(x).lower() in ['true','1', 'yes']), default=True, help='Use tuned C matrix for simulation')

    parser.add_argument('--generate_heatmap', type=lambda x: (str(x).lower() in ['true','1', 'yes']), default=True, help='Generate heatmap of final temperature')
    parser.add_argument('--time_heatmap', type=float, default=4, help='Time for heatmap generation in sec')

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    # load material properties from yaml file
    start = time.time()

    args = parse_args()

    material_properties = common.load_dict_yaml(args.material_prop_file)

    # load geometry properties from yaml file
    geometry_dict = common.load_dict_yaml(args.geometry_file)

    # power_config_dict = common.load_dict_yaml('power_dist_config.yml')
    power_config_dict = common.load_dict_yaml(args.power_config_file)

    common_utils = common.Utils(geometry_dict['common'])

    # create power grid object
    power_grid_class = Power_grid(power_config_dict, args)
    power_grid_class.create_power_seq_grid(utils=common_utils)

    package = Chiplet_package(material_properties, geometry_dict, power_grid_class, args, common_utils)

    package.create_layers()

    package.connect_nodes()

    package.generate_floorplan()

    package.run_simulation_c_superlu()
    print('Time taken for simulation: ', time.time()-start)
