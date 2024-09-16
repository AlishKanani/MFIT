import common

class Node:
    def __init__(self,x,y,z, x_length, y_length, thickness, material_properties=None, is_power_src=False):
        self.x = x
        self.y = y
        self.z = z
        self.x_length = x_length
        self.y_length = y_length
        self.thickness = thickness
        self.material_properties = material_properties
        self.is_power_src = is_power_src

        self.thermal_conductance_x = None
        self.thermal_conductance_y = None
        self.thermal_conductance_z = None
        self.thermal_capacitance = None

        self.convection_conductance = None

        # list to store the boundary condition, if nodes are on the boundary. add convective resistance to the node
        # x+ x- y+ y- z+ z-
        self.boundary_condition = {'z+': False, 'z-': False}

    def get_volume(self):
        return self.x_length * self.y_length * self.thickness
    
    def set_boundary_condition(self, boundary, utils):
        # boundary is string containing z+, z- or empty list
        if boundary == 'z+':
            self.boundary_condition['z+'] = True
            self.convection_conductance = common.cal_convection_conductance(h=utils.top_htc, 
                                                                          area=self.get_area_xy())
        elif boundary == 'z-':
            self.boundary_condition['z-'] = True
            self.convection_conductance = common.cal_convection_conductance(h=utils.bottom_htc, 
                                                                          area=self.get_area_xy())
    
    def get_area_xy(self):
        return self.x_length * self.y_length

    def get_area_xz(self):
        return self.x_length * self.thickness
    
    def get_area_yz(self):
        return self.y_length * self.thickness
    
    def set_thermal_properties(self, material_properties):
        self.set_thermal_capacitance(material_properties)
        self.set_thermal_conductance(material_properties)

    def set_thermal_capacitance(self, material_properties):
        density = material_properties['density']
        specific_heat = material_properties['specific_heat']
        self.thermal_capacitance = common.cal_thermal_capacitance(density=density, 
                                                                  specific_heat=specific_heat, 
                                                                  volume=self.get_volume())
        
    def set_thermal_conductance(self, material_properties):
        thermal_conductivity = material_properties['thermal_conductivity']

        if type(thermal_conductivity) is dict:
            tc_x = thermal_conductivity['kx']
            tc_y = thermal_conductivity['ky']
            tc_z = thermal_conductivity['kz']
        else:
            tc_x = thermal_conductivity
            tc_y = thermal_conductivity
            tc_z = thermal_conductivity
        self.thermal_conductance_x = common.cal_conductance(thermal_conductivity=tc_x, 
                                                          area=self.get_area_yz(), 
                                                          thickness=self.x_length)
        self.thermal_conductance_y = common.cal_conductance(thermal_conductivity=tc_y, 
                                                          area=self.get_area_xz(), 
                                                          thickness=self.y_length)
        self.thermal_conductance_z = common.cal_conductance(thermal_conductivity=tc_z, 
                                                          area=self.get_area_xy(), 
                                                          thickness=self.thickness)
    
    def get_thermal_conductance(self, direction):
        if direction == 'x':
            return self.thermal_conductance_x
        elif direction == 'y':
            return self.thermal_conductance_y
        elif direction == 'z':
            return self.thermal_conductance_z
        else:
            # raise an error
            pass

    def get_thermal_capacitance(self):
        return self.thermal_capacitance
    
    def get_thermal_conductance_bw_nodes_of_same_layer(self, node):
        # check if nodes corditnates can be in any of 4 neighboring directions
        # if yes, calculate overlap area (length, since thickness would be same) and return the conductance

        if abs((node.y + node.y_length ) - self.y) < 0.0001:
            # node can be below the current node y-
            common_len = common.calculate_overlapping_length(x1=self.x, x2=node.x, x_len1=self.x_length, x_len2=node.x_length)
            return (common_len*self.thermal_conductance_y*node.thermal_conductance_y)/(self.thermal_conductance_y*node.x_length + node.thermal_conductance_y*self.x_length)
        
        elif abs((node.y - self.y_length) - self.y) < 0.0001:
            # node can be above the current node y+
            common_len = common.calculate_overlapping_length(x1=self.x, x2=node.x, x_len1=self.x_length, x_len2=node.x_length)
            return (common_len*self.thermal_conductance_y*node.thermal_conductance_y)/(self.thermal_conductance_y*node.x_length + node.thermal_conductance_y*self.x_length)
        
        elif abs((node.x + node.x_length) - self.x) < 0.0001:
            # node can be to the left of the current node x-
            common_len = common.calculate_overlapping_length(x1=self.y, x2=node.y, x_len1=self.y_length, x_len2=node.y_length)
            return (common_len*self.thermal_conductance_x*node.thermal_conductance_x)/(self.thermal_conductance_x*node.y_length + node.thermal_conductance_x*self.y_length)
        
        elif abs((node.x - self.x_length) - self.x) < 0.0001:
            # node can be to the right of the current node x+
            common_len = common.calculate_overlapping_length(x1=self.y, x2=node.y, x_len1=self.y_length, x_len2=node.y_length)
            return (common_len*self.thermal_conductance_x*node.thermal_conductance_x)/(self.thermal_conductance_x*node.y_length + node.thermal_conductance_x*self.y_length)

        else:
            return 0.0
