import common
from node_class import Node
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import os, sys

# layer class for the chiplet thermal model, it takes a dictionary of material properties and geometry properties
class Layer_chiplet:
    def __init__(self, name, layer_dict, material_properties, power_grid_class, args):
        self.layer_name = name
        self.material = layer_dict['material']
        self.thickness = layer_dict['thickness']
        self.power_src = layer_dict['power_src']
        self.nodes_dict = layer_dict['nodes']
        self.start_point = layer_dict['start_point']
        self.args = args

        self.nodes = None
        self.layer_material_properties = material_properties[self.material]
        self.power_grid_class = power_grid_class

        self.total_x_nodes = None
        self.total_y_nodes = None

        self.total_nodes = None

    def is_power_src(self):
        return self.power_src
        
    def is_uniform(self):
        return self.nodes_dict['uniform']
    
    def is_layer_under_chiplet(self):
        return self.nodes_dict['under_chiplet']

    def layer_total_nodes(self):
        if self.is_layer_under_chiplet() and not self.args.is_homogeneous:
            return self.total_nodes
        else:
            return self.total_x_nodes*self.total_y_nodes
    
    def calculate_length_for_non_uniform_nodes(self, iter, nodes, start_cordinates, end_cordinates):
        if iter == 0:
            length = nodes[iter] + ((nodes[iter] + nodes[iter+1])/2.0 - nodes[iter])
            start_point = start_cordinates
        elif iter == len(nodes)-1 :
            length = (end_cordinates - (nodes[iter] + start_cordinates)) + (nodes[iter] - (nodes[iter] + nodes[iter-1])/2.0)
            start_point = (nodes[iter] + nodes[iter-1])/2.0 + start_cordinates
        else:
            length = (nodes[iter] + nodes[iter+1])/2.0 - (nodes[iter] + nodes[iter-1])/2.0
            start_point = (nodes[iter] + nodes[iter-1])/2.0 + start_cordinates

        return length, start_point
    
    def create_nodes(self, utils, material_properties):
        # creates nodes for the layer and updates thermal properties - RC 

        z_dim = self.start_point['z'] 

        bc = ''
        if self.start_point['z'] == 0.0:
            bc = 'z-'
        elif self.start_point['z'] + self.thickness == utils.package_z_len:
            bc = 'z+'

        if not self.is_layer_under_chiplet():

            if self.is_uniform():
                # uniform nodes, not under chiplet
                self.total_x_nodes = self.nodes_dict['x_nodes']
                self.total_y_nodes = self.nodes_dict['y_nodes']

                # check if x_nodes is not a list, error out if it is
                if isinstance(self.total_x_nodes, list):
                    print('Error: x_nodes should be an integer')
                    sys.exit(1)

                x_length = (utils.package_x_len - 2*self.start_point['x'])/self.total_x_nodes
                y_length = (utils.package_y_len - 2*self.start_point['y'])/self.total_y_nodes

                # create a 2D array of nodes
                self.nodes = [[None for i in range(self.total_y_nodes)] for j in range(self.total_x_nodes)]


                for i in range(self.total_x_nodes):
                    x_dim = i*float(x_length) + self.start_point['x']

                    for j in range(self.total_y_nodes):

                        y_dim = j*float(y_length) + self.start_point['y']

                        self.nodes[i][j] = Node(x=x_dim, y=y_dim, z=z_dim,
                                                x_length=x_length,
                                                y_length=y_length,
                                                thickness=self.thickness)
        
            else:
                # non uniform nodes, not under chiplet
                x_nodes = self.nodes_dict['x_nodes']
                y_nodes = self.nodes_dict['y_nodes']

                # check if x_nodes is a list, error out if it is not
                if not isinstance(x_nodes, list):
                    print('Error: x_nodes should be a list')
                    sys.exit(1)

                self.total_x_nodes = len(x_nodes)
                self.total_y_nodes = len(y_nodes)

                self.nodes = [[None for i in range(self.total_y_nodes)] for j in range(self.total_x_nodes)]

                for i in range(self.total_x_nodes):
                    # assume package is symmetric in x direction, empty space is same on both sides
                    x_length, x_dim = self.calculate_length_for_non_uniform_nodes(iter=i,
                                                    nodes=x_nodes,
                                                    start_cordinates=self.start_point['x'],
                                                    end_cordinates=(utils.package_x_len - self.start_point['x'])) 
                    for j in range(self.total_y_nodes):
                        # assume package is symmetric in y direction, empty space is same on both sides
                        y_length, y_dim = self.calculate_length_for_non_uniform_nodes(iter=j,
                                                        nodes=y_nodes,
                                                        start_cordinates=self.start_point['y'],
                                                        end_cordinates=(utils.package_y_len- self.start_point['y']))

                        self.nodes[i][j] = Node(x=x_dim, y=y_dim, z=z_dim,
                                                x_length=x_length,
                                                y_length=y_length,
                                                thickness=self.thickness)

        elif self.is_layer_under_chiplet() and self.args.is_homogeneous:
            # if layer is under chiplet, chiplets can be homogeneous or non homogeneous
            # for homogeneous chiples, with non uniform nodes - read from geometry file
            if not self.is_uniform():
                x_nodes = self.nodes_dict['x_nodes']
                y_nodes = self.nodes_dict['y_nodes']

                # check if x_nodes is a list, error out if it is not
                if not isinstance(x_nodes, list):
                    print('Error: x_nodes should be a list')
                    sys.exit(1)

                self.total_x_nodes = len(x_nodes)*utils.num_chiplet_x
                self.total_y_nodes = len(y_nodes)*utils.num_chiplet_y

                self.nodes = [[None for i in range(self.total_y_nodes)] for j in range(self.total_x_nodes)]

                for i_chiplet in range(utils.num_chiplet_x):
                    if i_chiplet == 0:
                        x_start = self.start_point['x']
                        x_end = x_start + utils.chiplet_x_len
                    else:
                        x_start = x_end + utils.chiplet_spacing 
                        x_end = x_start + utils.chiplet_x_len

                    for j_chiplet in range(utils.num_chiplet_y):
                        if j_chiplet == 0:
                            y_start = self.start_point['y']
                            y_end = y_start + utils.chiplet_y_len
                        else:
                            y_start = y_end + utils.chiplet_spacing 
                            y_end = y_start + utils.chiplet_y_len

                        for i in range(len(x_nodes)):
                            x_length, x_dim = self.calculate_length_for_non_uniform_nodes(iter=i,
                                                            nodes=x_nodes,
                                                            start_cordinates=x_start,
                                                            end_cordinates=x_end) 
                            for j in range(len(y_nodes)):
                                y_length, y_dim = self.calculate_length_for_non_uniform_nodes(iter=j,
                                                                nodes=y_nodes,
                                                                start_cordinates=y_start,
                                                                end_cordinates=y_end) 

                                self.nodes[i_chiplet*len(x_nodes) + i][j_chiplet*len(y_nodes) + j] = Node(x=x_dim, y=y_dim, z=z_dim,
                                                                                                        x_length=x_length,
                                                                                                        y_length=y_length,
                                                                                                        thickness=self.thickness)
            
            elif self.is_uniform():
                # for homogeneous chiplets, with uniform nodes add logic to create nodes - read from geometry file
                x_nodes_chiplet = self.nodes_dict['x_nodes']
                y_nodes_chiplet = self.nodes_dict['y_nodes']

                # check if x_nodes_chiplet is not a list, error out if it is
                if isinstance(x_nodes_chiplet, list):
                    print('Error: x_nodes_chiplet should be an integer')
                    sys.exit(1)


                self.total_x_nodes = x_nodes_chiplet*utils.num_chiplet_x
                self.total_y_nodes = y_nodes_chiplet*utils.num_chiplet_y

                self.nodes = [[None for i in range(self.total_y_nodes)] for j in range(self.total_x_nodes)]

                for i_chiplet in range(utils.num_chiplet_x):
                    if i_chiplet == 0:
                        x_start = self.start_point['x']
                        x_end = x_start + utils.chiplet_x_len
                    else:
                        x_start = x_end + utils.chiplet_spacing 
                        x_end = x_start + utils.chiplet_x_len

                    for j_chiplet in range(utils.num_chiplet_y):
                        if j_chiplet == 0:
                            y_start = self.start_point['y']
                            y_end = y_start + utils.chiplet_y_len
                        else:
                            y_start = y_end + utils.chiplet_spacing 
                            y_end = y_start + utils.chiplet_y_len

                        for i in range(x_nodes_chiplet):
                            x_dim = i*float(x_end - x_start)/x_nodes_chiplet + x_start
                            x_length = float(x_end - x_start)/x_nodes_chiplet

                            for j in range(y_nodes_chiplet):
                                y_dim = j*float(y_end - y_start)/y_nodes_chiplet + y_start
                                y_length = float(y_end - y_start)/y_nodes_chiplet

                                self.nodes[i_chiplet*x_nodes_chiplet + i][j_chiplet*y_nodes_chiplet + j] = Node(x=x_dim, y=y_dim, z=z_dim,
                                                                                                        x_length=x_length,
                                                                                                        y_length=y_length,
                                                                                                        thickness=self.thickness)
            
        elif self.is_layer_under_chiplet() and not self.args.is_homogeneous:
            # for non homogeneous chiplets, only uniform nodes are supported 
            # it would be read from power grid file, so all the layers would have same nodes
            # if layer is chiplet, use the nodes from the power grid class
            # else use the nodes from the first layer if layer is ubump 
            #  use the nodes from the last layer if layer is tim

            geometry_layer = None

            for power_layer in self.power_grid_class.power_layers:
                if power_layer.name == self.layer_name:
                    geometry_layer = power_layer
                    break
            
            if geometry_layer is None:
                # if 'ubump' in self.layer_name:
                #     geometry_layer = self.power_grid_class.power_layers[0]
                # elif 'tim' in self.layer_name:
                #     geometry_layer = self.power_grid_class.power_layers[-1]
                # else:
                print(f'Error: Layer {self.layer_name} not found')
                sys.exit(1)
            
            self.total_nodes = 0
            for chiplets in geometry_layer.power_chiplets:
                self.total_nodes = self.total_nodes + chiplets.nodes_x*chiplets.nodes_y

            self.nodes = [None for i in range(self.total_nodes)]

            prev_nodes = 0
            for chiplets in geometry_layer.power_chiplets:
                x_nodes_c = chiplets.nodes_x
                y_nodes_c = chiplets.nodes_y

                x_length = chiplets.length_chiplet_x/x_nodes_c
                y_length = chiplets.length_chiplet_y/y_nodes_c
                
                if hasattr(chiplets, 'material'):
                    chiplet_material = material_properties[chiplets.material]
                else:
                    chiplet_material = self.layer_material_properties

                for i in range(chiplets.nodes_x):
                    for j in range(chiplets.nodes_y):
                        self.nodes[prev_nodes + i*chiplets.nodes_y + j] = Node(x=chiplets.chiplet_x + i*x_length,
                                                                y=chiplets.chiplet_y + j*y_length,
                                                                z=self.start_point['z'],
                                                                x_length=x_length,
                                                                y_length=y_length,
                                                                thickness=self.thickness,
                                                                material_properties=chiplet_material,
                                                                is_power_src=chiplets.is_power_src)
                        
                prev_nodes = prev_nodes + x_nodes_c*y_nodes_c

        # set thermal capacitance and resistance for each node
        if self.is_layer_under_chiplet() and not self.args.is_homogeneous:
            for i in range(self.total_nodes):
                self.nodes[i].set_thermal_properties(self.nodes[i].material_properties)
                self.nodes[i].set_boundary_condition(boundary=bc, 
                                                    utils=utils)
        else:
            for i in range(self.total_x_nodes):
                for j in range(self.total_y_nodes):
                    self.nodes[i][j].set_thermal_properties(self.layer_material_properties)
                    self.nodes[i][j].set_boundary_condition(boundary=bc, 
                                                            utils=utils)

    def connect_nodes(self):
        # 1. each layer would create 1D and 2D array for capacitance and resistance. 
            # layer would return 1D array of capacitanace, 2D array of resistance, 
            # 1D array of Z resistance, length_x, length_y and x,y coordinates of the nodes
        
        self.capacitance = np.zeros(self.layer_total_nodes())
        self.xy_conductance = np.zeros((self.layer_total_nodes(), self.layer_total_nodes()))
        self.z_conductance = np.zeros(self.layer_total_nodes())
        self.x_lengths = np.zeros(self.layer_total_nodes())
        self.y_lengths = np.zeros(self.layer_total_nodes())
        self.xy_area = np.zeros(self.layer_total_nodes())
        self.x_cordinates = np.zeros(self.layer_total_nodes())
        self.y_cordinates = np.zeros(self.layer_total_nodes())

        # convection resistance
        self.z_plus_conductance = np.zeros(self.layer_total_nodes())
        self.z_minus_conductance = np.zeros(self.layer_total_nodes())


        if self.is_layer_under_chiplet() and not self.args.is_homogeneous:
            for ij in range(self.layer_total_nodes()):
                self.capacitance[ij] = self.nodes[ij].thermal_capacitance
                self.z_conductance[ij] = self.nodes[ij].thermal_conductance_z
                self.x_lengths[ij] = self.nodes[ij].x_length
                self.y_lengths[ij] = self.nodes[ij].y_length
                self.xy_area[ij] = self.nodes[ij].get_area_xy()
                self.x_cordinates[ij] = self.nodes[ij].x
                self.y_cordinates[ij] = self.nodes[ij].y

                if self.nodes[ij].boundary_condition['z+']:
                    self.z_plus_conductance[ij] = self.nodes[ij].convection_conductance*self.nodes[ij].thermal_conductance_z/(self.nodes[ij].convection_conductance + self.nodes[ij].thermal_conductance_z)
                if self.nodes[ij].boundary_condition['z-']:
                    self.z_minus_conductance[ij] = self.nodes[ij].convection_conductance*self.nodes[ij].thermal_conductance_z/(self.nodes[ij].convection_conductance + self.nodes[ij].thermal_conductance_z)

                for kl in range(self.layer_total_nodes()):
                    if ij == kl:
                        self.xy_conductance[ij][kl] = 0
                    else:
                        self.xy_conductance[ij][kl] = self.nodes[ij].get_thermal_conductance_bw_nodes_of_same_layer(self.nodes[kl])

        else:
            for i in range(self.total_x_nodes):
                for j in range(self.total_y_nodes):
                    self.capacitance[i*self.total_y_nodes + j] = self.nodes[i][j].thermal_capacitance
                    self.z_conductance[i*self.total_y_nodes + j] = self.nodes[i][j].thermal_conductance_z
                    self.x_lengths[i*self.total_y_nodes + j] = self.nodes[i][j].x_length
                    self.y_lengths[i*self.total_y_nodes + j] = self.nodes[i][j].y_length
                    self.xy_area[i*self.total_y_nodes + j] = self.nodes[i][j].get_area_xy()
                    self.x_cordinates[i*self.total_y_nodes + j] = self.nodes[i][j].x
                    self.y_cordinates[i*self.total_y_nodes + j] = self.nodes[i][j].y

                    if self.nodes[i][j].boundary_condition['z+']:
                        self.z_plus_conductance[i*self.total_y_nodes + j] = self.nodes[i][j].convection_conductance*self.nodes[i][j].thermal_conductance_z/(self.nodes[i][j].convection_conductance + self.nodes[i][j].thermal_conductance_z)
                        # self.z_plus_conductance[i*self.total_y_nodes + j] = self.nodes[i][j].convection_resistance + self.nodes[i][j].thermal_resistance_z
                    if self.nodes[i][j].boundary_condition['z-']:
                        self.z_minus_conductance[i*self.total_y_nodes + j] = self.nodes[i][j].convection_conductance*self.nodes[i][j].thermal_conductance_z/(self.nodes[i][j].convection_conductance + self.nodes[i][j].thermal_conductance_z)
                        # self.z_minus_conductance[i*self.total_y_nodes + j] = self.nodes[i][j].convection_resistance + self.nodes[i][j].thermal_resistance_z

            # add logic to calculate the resistance between the nodes
                    for k in range(self.total_x_nodes):
                        for l in range(self.total_y_nodes):
                            if i == k and j == l:
                                self.xy_conductance[i*self.total_y_nodes + j][k*self.total_y_nodes + l] = 0
                            else:
                                self.xy_conductance[i*self.total_y_nodes + j][k*self.total_y_nodes + l] = self.nodes[i][j].get_thermal_conductance_bw_nodes_of_same_layer(self.nodes[k][l])

    def get_power(self, total_steps):
        
        self.power = np.zeros((self.layer_total_nodes(), total_steps))
        
        if self.is_power_src():
            # get corresponding power layer from power grid class and if it doesn't exist, error out

            found_power_layer = False
            now_power_layer = None
            for power_layer in self.power_grid_class.power_layers:
                if power_layer.name == self.layer_name:
                    found_power_layer = True
                    now_power_layer = power_layer

            if not found_power_layer:
                print(f'Error: Power layer: {self.layer_name} not found')
                sys.exit(1)

            # get area overlap between power block and nodes and assign power to the nodes
            if self.is_layer_under_chiplet() and not self.args.is_homogeneous:
                for node_iter in range(self.layer_total_nodes()):
                    node = self.nodes[node_iter]
                    node_area = node.get_area_xy()
                    node_power = np.zeros(total_steps)
                    if node.is_power_src:
                        for power_chiplet in now_power_layer.power_chiplets:
                            for power_block in power_chiplet.power_blocks:
                                overlap_area = common.calculate_overlapping_area(x1 = node.x,
                                                                                y1 = node.y,
                                                                                x2 = power_block.x,
                                                                                y2 = power_block.y,
                                                                                x_len1 = node.x_length,
                                                                                y_len1 = node.y_length,
                                                                                x_len2 = power_block.length_x,
                                                                                y_len2 = power_block.length_y)
                                if overlap_area > 0:
                                    node_max_power =  power_block.max_power*overlap_area/(power_block.get_area())
                                    node_area = node_area - overlap_area
                                    node_power = node_power + node_max_power*power_block.power_seq

                        self.power[node_iter] = node_power
                        if (node_area > 0.0001) or (node_area < -0.0001):
                            print(f'node power {self.layer_name}, {node_area} not found')
                            sys.exit(1)
            else:
                node_iter = 0
                for node_x in range(self.total_x_nodes):
                    for node_y in range(self.total_y_nodes):
                        node = self.nodes[node_x][node_y]

                        # get area of the node to check if it is fully covered by power blocks
                        node_area = node.get_area_xy()
                        node_power = np.zeros(total_steps)
                        for power_chiplet in now_power_layer.power_chiplets:
                            for power_block in power_chiplet.power_blocks:
                                overlap_area = common.calculate_overlapping_area(x1 = node.x,
                                                                                y1 = node.y,
                                                                                x2 = power_block.x,
                                                                                y2 = power_block.y,
                                                                                x_len1 = node.x_length,
                                                                                y_len1 = node.y_length,
                                                                                x_len2 = power_block.length_x,
                                                                                y_len2 = power_block.length_y)
                                if overlap_area > 0:
                                    node_max_power =  power_block.max_power*overlap_area/(power_block.get_area())
                                    node_area = node_area - overlap_area

                                    node_power = node_power + node_max_power*power_block.power_seq

                        self.power[node_iter] = self.power[node_iter] + node_power
                        if (node_area > 0.0001) or (node_area < -0.0001):
                            print(f'node power {self.layer_name}, {node_area} not found')
                            sys.exit(1)
                        node_iter = node_iter + 1
                
            
            # for node_iter in range(self.layer_total_nodes()):
            #     self.power[node_iter] = 3/8.0
        return self.power
    
    def get_capacitance(self):
        return self.capacitance
    
    def get_conductance(self):
        return self.xy_conductance
    
    def get_convective_conductance(self):
        return self.z_plus_conductance + self.z_minus_conductance
    
    def write_temperature_to_file(self, block_temperatures, ts):
        file_name = f'{self.args.output_dir}/output/temperature_{self.layer_name}_{ts}.csv'
        with open(file_name, 'w') as f:
            for block_temperature in block_temperatures:
                f.write(','.join([str(i) for i in block_temperature]) + '\n')
    

    def map_temperature_to_blk(self, temperature_all, utils, ts):
        # only for chiplet layers - is power source true

        # its list, each element is again a list, where 1st element is the name of the chiplet
        block_temperatures = []
        geometry_layer = None

        for power_layer in self.power_grid_class.power_layers:
            if power_layer.name == self.layer_name:
                geometry_layer = power_layer
                break
        if geometry_layer is None:
            # not found, error out
            print(f'Error: Layer {self.layer_name} not found')
            sys.exit(1)

        if self.is_layer_under_chiplet():
            if self.args.is_homogeneous:
                # create a new numpy array to group the temperature of the chiplets
                sorted_temperature = np.zeros((self.total_x_nodes*self.total_y_nodes, temperature_all.shape[1]))
                if self.is_uniform():
                    x_nodes_chiplet = self.nodes_dict['x_nodes']
                    y_nodes_chiplet = self.nodes_dict['y_nodes']

                    groupped_index = 0
                    for i_chiplet in range(utils.num_chiplet_x):
                        for j_chiplet in range(utils.num_chiplet_y):
                            for i in range(x_nodes_chiplet):
                                for j in range(y_nodes_chiplet):
                                    sorted_temperature[groupped_index] = temperature_all[i_chiplet*utils.num_chiplet_y*y_nodes_chiplet*x_nodes_chiplet + 
                                                                                         j_chiplet*y_nodes_chiplet + 
                                                                                         i*y_nodes_chiplet*utils.num_chiplet_y + j]
                                    groupped_index = groupped_index + 1
                
                    start_node = 0
                    for chiplet in geometry_layer.power_chiplets:
                        total_nodes = x_nodes_chiplet*y_nodes_chiplet
                        name = chiplet.name
                        temperature_this_chiplet = sorted_temperature[start_node:start_node + total_nodes]
                        start_node = start_node + total_nodes
                        temperature_this_chiplet = (np.mean(temperature_this_chiplet, axis=0)).tolist()
                        block_temperatures.append([name] + temperature_this_chiplet)

                else:
                    # not supported, print an error
                    print('Mapping temperature to block not supported for non uniform chiplet layers')
                    sys.exit(1)
            else:
                sorted_temperature = temperature_all
                start_node = 0
                for chiplet in geometry_layer.power_chiplets:
                    total_nodes = chiplet.nodes_x*chiplet.nodes_y
                    name = chiplet.name
                    temperature_this_chiplet = sorted_temperature[start_node:start_node + total_nodes]
                    start_node = start_node + total_nodes

                    # temperature_this_chiplet is total_nodes x total_steps, average it over rows
                    temperature_this_chiplet = (np.mean(temperature_this_chiplet, axis=0)).tolist()
                    block_temperatures.append([name] + temperature_this_chiplet)

        else:
            # not supported, print an error
            print('Error: Mapping temperature to block not supported for non chiplet layers')
            sys.exit(1)
        
        # save the temperature to a file
        self.write_temperature_to_file(block_temperatures, ts)

    def plot_layer(self, utils, layer_start):
        fig, ax = plt.subplots()

        if self.is_layer_under_chiplet() and not self.args.is_homogeneous:
            for i in range(self.total_nodes):
                rect = Rectangle((self.nodes[i].x , self.nodes[i].y), 
                                    self.nodes[i].x_length, self.nodes[i].y_length, 
                                    fc="none", ec="black")
                ax.add_patch(rect)
                plt.plot(self.nodes[i].x, self.nodes[i].y, 'ro')
                plt.text(self.nodes[i].x, self.nodes[i].y, str(layer_start+i), fontsize=5)
        else:
            for i in range(self.total_x_nodes):
                for j in range(self.total_y_nodes):
                    rect = Rectangle((self.nodes[i][j].x , self.nodes[i][j].y), 
                                        self.nodes[i][j].x_length, self.nodes[i][j].y_length, 
                                        fc="none", ec="black")
                    ax.add_patch(rect)
                    plt.plot(self.nodes[i][j].x, self.nodes[i][j].y, 'ro')
                    plt.text(self.nodes[i][j].x, self.nodes[i][j].y, str(layer_start + i*self.total_y_nodes + j), fontsize=5)

        ax.set_xlim(-0.5, utils.package_x_len + 0.5)
        ax.set_ylim(-0.5, utils.package_y_len + 0.5)

        plt.title(self.layer_name + ' Floorplan')
        plt.xlabel('X dimension (mm)')
        plt.ylabel('Y dimension (mm)')

        # check if floorplan directory exists, if not create it

        if not os.path.exists(self.args.output_dir + '/floorplan'):
            os.makedirs(self.args.output_dir + '/floorplan')

        fig.savefig(self.args.output_dir + '/floorplan/' + self.layer_name + '.png', dpi=300, bbox_inches='tight')
        plt.close(fig)

    def plot_heatmap(self, temperature_all_layer, utils):
        fig, ax = plt.subplots()

        norm = plt.Normalize(temperature_all_layer.min()-1, temperature_all_layer.max()+1)
        cmap = plt.cm.hot_r

        if self.is_layer_under_chiplet() and not self.args.is_homogeneous:
            for i in range(self.total_nodes):
                colour = cmap(norm(temperature_all_layer[i]))
                rect = Rectangle((self.nodes[i].x , self.nodes[i].y), 
                                    self.nodes[i].x_length, self.nodes[i].y_length, 
                                    linewidth=1, edgecolor='black',
                                    facecolor=colour)
                ax.add_patch(rect)
        
        else:
            for i in range(self.total_x_nodes):
                for j in range(self.total_y_nodes):
                    colour = cmap(norm(temperature_all_layer[i*self.total_y_nodes + j]))
                    rect = Rectangle((self.nodes[i][j].x , self.nodes[i][j].y), 
                                        self.nodes[i][j].x_length, self.nodes[i][j].y_length,
                                        linewidth=1, edgecolor='black',
                                        facecolor=colour) 
                    ax.add_patch(rect)

        ax.set_xlim(-0.5, utils.package_x_len + 0.5 +1)
        ax.set_ylim(-0.5, utils.package_y_len + 0.5 +1)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array(temperature_all_layer)
        plt.colorbar(sm, ax=ax).set_label('Temperature (C)')

        plt.gca().set_aspect('equal', adjustable='box')
        plt.title(self.layer_name + ' Heatmap')
        plt.xlabel('X dimension (mm)')
        plt.ylabel('Y dimension (mm)')

        if not os.path.exists(self.args.output_dir + '/heatmaps'):
            os.makedirs(self.args.output_dir + '/heatmaps')
        
        plt.savefig(self.args.output_dir + '/heatmaps/' + self.layer_name + '_heatmap.png', dpi=300, bbox_inches='tight')

        # close the plot
        plt.close(fig)

