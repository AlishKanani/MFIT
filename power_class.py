import csv
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import os
import common

class Power_block:
    def __init__(self, power_node_dict, name, chiplet_x, chiplet_y, args):
        self.name = name
        self.x = power_node_dict['start_point_x'] + chiplet_x
        self.y = power_node_dict['start_point_y'] + chiplet_y
        self.length_x = power_node_dict['length_x']
        self.length_y = power_node_dict['length_y']
        self.max_power = float(power_node_dict['max_power'])
        self.args = args

    def create_power_seq(self, layer_name, chiplet_name, power_seq):
        found_flag = False
        for i in range(0, len(power_seq)):
            power_sq_name =  layer_name + '_' + chiplet_name + '_' + self.name
            if power_seq[i][0] == power_sq_name:
                self.power_seq = power_seq[i][1:]
                found_flag = True
                break
        
        if not found_flag:
            print('Power sequence not found for block: ', layer_name + '_' + chiplet_name + '_' + self.name)
            return
        
        # power seq is in percentage, convert to float
        if self.args.simulation_type == 'steady':
            self.power_seq = np.array([1])
        else:
            self.power_seq = np.array([float(i)/100.0 for i in self.power_seq])

        # print
        # print('Power sequence for node: ', layer , '_', self.name, ' is: ', self.power_seq)

    def get_area(self):
        return self.length_x * self.length_y

class Power_chiplet:
    def __init__(self, power_chiplet_dict, name, args):
        self.name = name
        self.power_blocks = []
        self.args = args
        self.chiplet_x = power_chiplet_dict['start_chiplet_x']
        self.chiplet_y = power_chiplet_dict['start_chiplet_y']
        self.is_power_src = True
        if not self.args.is_homogeneous:
            self.length_chiplet_x = power_chiplet_dict['length_chiplet_x']
            self.length_chiplet_y = power_chiplet_dict['length_chiplet_y']
            self.nodes_x = power_chiplet_dict['nodes_x']
            self.nodes_y = power_chiplet_dict['nodes_y']

            if 'material' in power_chiplet_dict:
                self.material = power_chiplet_dict['material']
            
        if 'layout_blocks' in power_chiplet_dict:
            for block in power_chiplet_dict['layout_blocks']:
                power_block = Power_block(power_chiplet_dict['layout_blocks'][block], block, self.chiplet_x, self.chiplet_y, self.args)
                self.power_blocks.append(power_block)
        else:
            self.is_power_src = False
    
    def create_power_seq_chiplet(self, layer_name, power_seq):
        for block in self.power_blocks:
            block.create_power_seq(layer_name, self.name, power_seq)

class Power_layer:
    def __init__(self, power_layer_dict, layer, args):
        self.name = layer
        self.power_chiplets = []
        self.args = args
        for chiplet in power_layer_dict[layer]:
            power_chiplet = Power_chiplet(power_layer_dict[layer][chiplet], chiplet, self.args)
            self.power_chiplets.append(power_chiplet)
        
        # sort the chiplets based on x and y
        sorted_by_y = sorted(self.power_chiplets, key=lambda chiplet: chiplet.chiplet_y)
        self.power_chiplets = sorted(sorted_by_y, key=lambda chiplet: chiplet.chiplet_x)

        # for chiplet in self.power_chiplets:
        #     print(chiplet.chiplet_x, chiplet.chiplet_y)

    def create_power_seq_layer(self, power_seq):
        for chiplet in self.power_chiplets:
            chiplet.create_power_seq_chiplet(self.name, power_seq)
    
    def plot_layer(self, utils):
        fig, ax = plt.subplots(figsize=(6,6))
        
        for power_chiplet in self.power_chiplets:
            for power_block in power_chiplet.power_blocks:
                rect = Rectangle((power_block.x , power_block.y), 
                                    power_block.length_x, power_block.length_y, 
                                    fc="none", ec="black", lw=0.5)
                ax.add_patch(rect)
                name = power_chiplet.name + '-' + power_block.name
                # name = power_block.name
                plt.plot(power_block.x, power_block.y, 'ro', markersize=0.5)
                plt.text(power_block.x, power_block.y, name, fontsize=5)

        ax.set_xlim(-0.5, utils.package_x_len + 0.5)
        ax.set_ylim(-0.5, utils.package_y_len + 0.5)

        plt.title('Power grid for ' + self.name)
        plt.xlabel('X dimension (mm)')
        plt.ylabel('Y dimension (mm)')

        # check if floorplan directory exists, if not create it

        if not os.path.exists(self.args.output_dir + '/floorplan'):
            os.makedirs(self.args.output_dir + '/floorplan')

        fig.savefig(self.args.output_dir + '/floorplan/' + self.name + '_power_' + '.png', dpi=300, bbox_inches='tight')
        plt.close()


class Power_grid:
    def __init__(self, power_dict, args):
        self.power_dict = power_dict
        self.args = args

        with open(self.args.power_seq_file, 'r') as f:
            reader = csv.reader(f)
            self.power_seq = list(reader)
        
        self.power_layers = []
        for layer in self.power_dict:
            power_layer = Power_layer(self.power_dict, layer, self.args)
            self.power_layers.append(power_layer)
    
    def create_power_seq_grid(self, utils):
        for layer in self.power_layers:
            layer.plot_layer(utils)
            layer.create_power_seq_layer(self.power_seq)

if __name__ == '__main__':
    power_config_dict = common.load_dict_yaml('power_dist_config.yml')
    grid = Power_grid(power_config_dict, 'power_seq.csv')

    grid.create_power_seq_grid()

    for layer in grid.power_layers:
        if layer.name == 'chiplet_1':
            for node in layer.power_nodes:
                print(node.name, node.power_seq, node.x, node.y, node.length_x, node.length_y)