import common 
from layer_class import Layer_chiplet
import numpy as np
import os
from scipy.signal import lti
from ctypes import *

# load the shared C library, for the entire c-based solver
lib_solver = CDLL('c_files/chiplet_ode.so')  

lib_solver.chiplet_ode.argtypes = [
    POINTER(POINTER(c_double)),  # output temperature
    POINTER(POINTER(c_double)),  # power input
    POINTER(POINTER(c_double)),  # G_all
    POINTER(POINTER(c_int)),   # non_zero_indexes
    np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='C_CONTIGUOUS'),   # C
    c_int,  # total_nodes
    c_int,  # non_zero_columns
    c_double, # total_duration
    c_double, # time_step
    c_double, # power_interval
]

lib_solver.chiplet_ode.restype = None

def chiplet_ode_c(output_temperature, power_input, G_all, non_zero_indexes, C, total_nodes, non_zero_columns, total_duration, time_step, power_interval):

    lib_solver.chiplet_ode(output_temperature, power_input, G_all, non_zero_indexes, C, total_nodes, non_zero_columns, total_duration, time_step, power_interval)

class Chiplet_package:
    def __init__(self, material_properties, geometry_dict, power_grid_class, args, utils=common.Utils):
        self.material_properties = material_properties
        self.geometry_dict = geometry_dict
        self.layers = []
        self.common_utils = utils
        self.power_grid_class = power_grid_class
        self.args = args
    
    def create_layers(self):
        for layer in self.geometry_dict['layers']:
            self.layers.append(Layer_chiplet(name=layer, 
                                             layer_dict=self.geometry_dict['layers'][layer], 
                                             material_properties=self.material_properties,
                                             power_grid_class=self.power_grid_class,
                                             args=self.args))
            
        for layer in self.layers:
            layer.create_nodes(utils=self.common_utils, 
                               material_properties=self.material_properties)
    
    def connect_nodes(self):
        # create the thermal network, create a 1D numpy array of nodes with thermal capacitance.
        # create a 2D array of thermal resistance between the nodes
        
        # 1. each layer would create 1D and 2D array for capacitance and resistance. 
            # layer would return 1D array of capacitanace, 2D array of resistance, 
            # 1D array of Z resistance, 1D array of area for Z direction and x,y coordinates of the nodes
        # 2. connect the layers for Z direction, using the Z resistance and area overlap
        # 3. Add a ground node to connect top and bottom layer with the convection resistance

        # Three resistances: layer's internal, layer's with top and bottom layer and convection resistance

        # connects nodes and calculate the RC between the nodes, and convection resistance
        for layer in self.layers:
            layer.connect_nodes()

        self.shared_conductance_list = []
        # connect the layers for Z direction
        for i in range(len(self.layers)-1):
            bottom_layer = self.layers[i]
            top_layer = self.layers[i+1]
            shared_conductance = self.connect_layers(top_layer, bottom_layer)
            self.shared_conductance_list.append(shared_conductance)

        self.capacitance_all = []
        for layer in self.layers:
            self.capacitance_all.append(layer.get_capacitance())
        
        self.capacitance_all = np.concatenate(self.capacitance_all, axis=0)
        
        # invert the capacitance
        self.capacitance_all = 1/self.capacitance_all

        self.conductance_all = np.zeros((self.package_total_nodes(), self.package_total_nodes() + 1))
        layer_global_iter = 0
        for layer_iter in range(len(self.layers)):
            current_layer = self.layers[layer_iter]
            convective_cond = current_layer.get_convective_conductance()
            xy_conductance = current_layer.get_conductance()
            current_layer_length = current_layer.layer_total_nodes()

            if layer_iter == 0: 
                bottom_layer = 0
                bottom_layer_length = 0
                bottom_conductance = np.zeros(1)

                top_layer = self.layers[layer_iter+1]
                top_layer_length = top_layer.layer_total_nodes()
                top_conductance = self.shared_conductance_list[layer_iter]
            
            elif layer_iter == len(self.layers)-1:
                top_layer = 0
                top_layer_length = 0
                top_conductance = np.zeros(1)

                bottom_layer = self.layers[layer_iter-1]
                bottom_layer_length = bottom_layer.layer_total_nodes()
                bottom_conductance = np.transpose(self.shared_conductance_list[layer_iter-1])

            else:
                top_layer = self.layers[layer_iter+1]
                top_layer_length = top_layer.layer_total_nodes()
                top_conductance = self.shared_conductance_list[layer_iter]

                bottom_layer = self.layers[layer_iter-1]
                bottom_layer_length = bottom_layer.layer_total_nodes()
                bottom_conductance = np.transpose(self.shared_conductance_list[layer_iter-1])
            
            for layer_node_iter in range(current_layer_length):
                node_iter = layer_global_iter + layer_node_iter
                
                # convective resistance
                self.conductance_all[node_iter, -1] = convective_cond[layer_node_iter]
                
                # layer resistance
                self.conductance_all[node_iter, layer_global_iter:layer_global_iter+current_layer_length] = xy_conductance[layer_node_iter]

                # bottom layer resistance
                if layer_iter != 0:
                    self.conductance_all[node_iter, layer_global_iter - bottom_layer_length:layer_global_iter] = bottom_conductance[layer_node_iter]
                
                # top layer resistance
                if layer_iter != len(self.layers)-1:
                    self.conductance_all[node_iter, layer_global_iter + current_layer_length:layer_global_iter + current_layer_length + top_layer_length] = top_conductance[layer_node_iter]

            layer_global_iter += current_layer_length

        conductance_diag = np.sum(self.conductance_all, axis=1)
        conductance_diag = np.diag(conductance_diag)

        self.conductance_all = self.conductance_all[:,:-1]
        self.conductance_all = conductance_diag - self.conductance_all
        
        if self.args.use_tuned_C:
            self.apply_tuned_C()
        else:
            pass

        if self.args.generate_DSS:
            self.generate_DSS()

        non_zero_conductance = [row[row != 0] for row in self.conductance_all]
        non_zero_index = [list(np.nonzero(row)[0]) for row in self.conductance_all]
        max_len = max(len(row) for row in non_zero_conductance)

        self.non_zero_index = np.array([row + [row[-1]] * (max_len - len(row)) for row in non_zero_index])
        self.conductance_all = np.zeros((self.conductance_all.shape[0], max_len))

        for i, row in enumerate(non_zero_conductance):
            self.conductance_all[i, :len(row)] = row

    def connect_layers(self, top_layer, bottom_layer):
        shared_conductance = np.zeros((bottom_layer.layer_total_nodes(), top_layer.layer_total_nodes()))
        
        for i in range(bottom_layer.layer_total_nodes()):
            for j in range(top_layer.layer_total_nodes()):
                # calculate overlapping area
                common_area = common.calculate_overlapping_area(x1=bottom_layer.x_cordinates[i], y1=bottom_layer.y_cordinates[i], 
                                                                x2=top_layer.x_cordinates[j], y2=top_layer.y_cordinates[j], 
                                                                x_len1=bottom_layer.x_lengths[i], y_len1=bottom_layer.y_lengths[i], 
                                                                x_len2=top_layer.x_lengths[j], y_len2=top_layer.y_lengths[j])
                
                if common_area > 0:
                    shared_conductance[i,j] = (common_area*bottom_layer.z_conductance[i]*top_layer.z_conductance[j])/(bottom_layer.z_conductance[i]*top_layer.xy_area[j] + top_layer.z_conductance[j]*bottom_layer.xy_area[i])
                    # shared_conductance[i,j] = (top_layer.z_resistance[i]*top_layer.xy_area[i] + bottom_layer.z_resistance[j]*bottom_layer.xy_area[j])/(common_area)

        return shared_conductance

    def apply_tuned_C(self):
        intial_C_guess = np.array([1.35551358, 1.3345646, 0.46572207, 0.85322922, 2.0361129, 1.77131198, 2.0619255, 0.94317305, 0.6672266])

        num_nodes = 0

        for layer in self.layers:
            if 'substrate_1' in layer.layer_name:
                layer_start = num_nodes
                num_nodes += layer.layer_total_nodes()
                layer_end = num_nodes
                self.capacitance_all[layer_start:layer_end] = intial_C_guess[0]*self.capacitance_all[layer_start:layer_end]
            
            elif 'substrate_2' in layer.layer_name:
                layer_start = num_nodes
                num_nodes += layer.layer_total_nodes()
                layer_end = num_nodes
                self.capacitance_all[layer_start:layer_end] = intial_C_guess[1]*self.capacitance_all[layer_start:layer_end] 
            
            elif 'c4' in layer.layer_name:
                layer_start = num_nodes
                num_nodes += layer.layer_total_nodes()
                layer_end = num_nodes
                self.capacitance_all[layer_start:layer_end] = intial_C_guess[2]*self.capacitance_all[layer_start:layer_end]

            elif 'interposer' in layer.layer_name:
                layer_start = num_nodes
                num_nodes += layer.layer_total_nodes()
                layer_end = num_nodes
                self.capacitance_all[layer_start:layer_end] = intial_C_guess[3]*self.capacitance_all[layer_start:layer_end]
            
            elif 'ubump' in layer.layer_name:
                layer_start = num_nodes
                num_nodes += layer.layer_total_nodes()
                layer_end = num_nodes
                self.capacitance_all[layer_start:layer_end] = intial_C_guess[4]*self.capacitance_all[layer_start:layer_end]
            
            elif 'chiplet' in layer.layer_name:
                layer_start = num_nodes
                num_nodes += layer.layer_total_nodes()
                layer_end = num_nodes
                self.capacitance_all[layer_start:layer_end] = intial_C_guess[5]*self.capacitance_all[layer_start:layer_end]
            
            elif 'tim' in layer.layer_name:
                layer_start = num_nodes
                num_nodes += layer.layer_total_nodes()
                layer_end = num_nodes
                self.capacitance_all[layer_start:layer_end] = intial_C_guess[6]*self.capacitance_all[layer_start:layer_end]
            
            elif 'lid1' in layer.layer_name:
                layer_start = num_nodes
                num_nodes += layer.layer_total_nodes()
                layer_end = num_nodes
                self.capacitance_all[layer_start:layer_end] = intial_C_guess[7]*self.capacitance_all[layer_start:layer_end]
            
            elif 'lid2' in layer.layer_name:
                layer_start = num_nodes
                num_nodes += layer.layer_total_nodes()
                layer_end = num_nodes
                self.capacitance_all[layer_start:layer_end] = intial_C_guess[8]*self.capacitance_all[layer_start:layer_end]
    
    def generate_floorplan(self):
        # generate the floorplan of the package, and chiplet
        num_nodes = 0
        for layer in self.layers:
            layer.plot_layer(utils=self.common_utils, layer_start=num_nodes)
            num_nodes += layer.layer_total_nodes()

    def package_total_nodes(self):
        total_nodes = 0
        for layer in self.layers:
            total_nodes += layer.layer_total_nodes()
        return total_nodes
    
    def generate_DSS(self):
        # generate the A and B matrix for DSS
        # A = -C^-1*G
        # B = C^-1
        # C = I
        # D = 0
        capacitance_matrix = np.diag(self.capacitance_all)
        A = -(capacitance_matrix @ self.conductance_all)
        B = capacitance_matrix
        C = np.eye(self.package_total_nodes())
        D = np.zeros((self.package_total_nodes(), self.package_total_nodes()))

        l_sys = lti(A, B, C, D)
        d_sys = l_sys.to_discrete(self.args.time_step, method='zoh')

        discrete_A = d_sys.A
        discrete_B = d_sys.B
        
        if not os.path.exists(self.args.output_dir + '/output'):
            os.makedirs(self.args.output_dir + '/output')

        np.savetxt(self.args.output_dir + '/output/disc_A_matrix.csv', discrete_A, delimiter=',')
        np.savetxt(self.args.output_dir + '/output/disc_B_matrix.csv', discrete_B, delimiter=',')

    def write_temperature_to_file(self, ts):
        self.temperature_all_save = np.array(self.temperature_all_save) + 300.0

        # save the temperature to a file
        file_name = f'{self.args.output_dir}/output/temperature_all_{ts}.csv'
        np.savetxt(file_name, self.temperature_all_save, delimiter=',')

        
        temperature_all_map = self.temperature_all_save.T

        if self.args.generate_heatmap:
            index_heatmap = int(self.args.time_heatmap/self.args.time_step)
            plot_temperature = temperature_all_map[:, index_heatmap] - 300.0
            num_nodes = 0
            for layer in self.layers:
                layer_start = num_nodes
                num_nodes += layer.layer_total_nodes()
                layer_end = num_nodes
                layer.plot_heatmap(plot_temperature[layer_start:layer_end], utils=self.common_utils)

        num_nodes = 0
        for layer in self.layers:
            layer_start = num_nodes
            num_nodes += layer.layer_total_nodes()
            layer_end = num_nodes
            if layer.is_power_src():
                layer.map_temperature_to_blk(temperature_all_map[layer_start:layer_end], utils=self.common_utils, ts=ts)

    def convert_to_np_array(self, pointer):
        def dereference_pointer(pointer, length):
            double_array = cast(pointer, POINTER(c_double * length))
            double_list = list(double_array.contents)
            return double_list
        
        self.temperature_all_save = [dereference_pointer(pointers,self.package_total_nodes()) for pointers in pointer]

    def set_initial_conditions(self):
        if self.args.simulation_type == 'steady':
            power_steps = 1
        else:
            power_steps = int(self.args.total_duration/float(self.args.power_interval))
        self.temperature_all_save = []
        self.temperature_all = np.zeros(self.package_total_nodes())
        self.power = np.zeros((self.package_total_nodes(), power_steps))
        
        # set power for chiplet nodes
        global_iter = 0
        for layer in self.layers:
            self.power[global_iter:global_iter+layer.layer_total_nodes()] = layer.get_power(power_steps)
            global_iter += layer.layer_total_nodes()

        # export power to csv file
        np.savetxt(self.args.output_dir + '/output/power_all.csv', self.power.T, delimiter=',')

    def run_simulation_c_lsoda(self):

        total_duration = self.args.total_duration
        self.set_initial_conditions()
        power = self.power.T
        
        if self.args.simulation_type == 'steady':
            dt = float(total_duration)
            power_interval = float(total_duration)

            np_temperature_all = np.zeros((2, self.package_total_nodes()))
            c_temperature_all = (POINTER(c_double) * 2)()
            c_temperature_all[0] = (c_double * self.package_total_nodes())(*np_temperature_all[0])
            c_temperature_all[1] = (c_double * self.package_total_nodes())(*np_temperature_all[1])

            c_power = (POINTER(c_double) * 1)()
            c_power[0] = (c_double * power.shape[1])(*power[0])

        else:
            dt = self.args.time_step
            power_interval = self.args.power_interval

            np_temperature_all = np.zeros((int(total_duration/dt) + 1, self.package_total_nodes()))
            c_temperature_all = (POINTER(c_double) * (int(total_duration/dt) + 1))()
            for i in range(int(total_duration/dt) + 1):
                c_temperature_all[i] = (c_double * self.package_total_nodes())(*np_temperature_all[i])

            c_power = (POINTER(c_double) * power.shape[0])()
            for i in range(power.shape[0]):
                c_power[i] = (c_double * power.shape[1])(*power[i])

        G_all = self.conductance_all
        non_zero_index = self.non_zero_index

        c_G_all = (POINTER(c_double) * G_all.shape[0])()
        for i in range(G_all.shape[0]):
            c_G_all[i] = (c_double * G_all.shape[1])(*G_all[i])

        c_non_zero_index = (POINTER(c_int) * non_zero_index.shape[0])()
        for i in range(non_zero_index.shape[0]):
            c_non_zero_index[i] = (c_int * non_zero_index.shape[1])(*non_zero_index[i])

        chiplet_ode_c(c_temperature_all, c_power, 
                      c_G_all, 
                      c_non_zero_index, 
                      self.capacitance_all, 
                      self.package_total_nodes(), 
                      self.non_zero_index.shape[1], 
                      total_duration, dt, power_interval)
        

        self.convert_to_np_array(c_temperature_all)
        self.write_temperature_to_file(dt)
        
