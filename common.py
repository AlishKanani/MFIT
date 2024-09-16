import yaml

# registance function
# thickness, area are in mm, give half resistance because of 6 faces
def cal_conductance(thermal_conductivity, area, thickness):
    return (thermal_conductivity * area*2)/(thickness*1000.0)

def cal_convection_conductance(h, area):
    return (h * area)/(1*1000000.0) 

# thermal capacitance function, with density, specific heat, and volume
def cal_thermal_capacitance(specific_heat, volume, density):
    return specific_heat * volume * density/1000000000.0

# load material properties from yaml file
def load_dict_yaml(yaml_file):
    with open(yaml_file, 'r') as f:
        dict = yaml.safe_load(f)
    return dict

# return data from dictionary with key as input
def get_data_from_dict(data, key):
    return data[key]

def calculate_overlapping_area(x1, y1, x2, y2, x_len1, y_len1, x_len2, y_len2):
    x_overlap = max(0, min(x1 + x_len1, x2 + x_len2) - max(x1, x2))
    y_overlap = max(0, min(y1 + y_len1, y2 + y_len2) - max(y1, y2))
    return x_overlap * y_overlap

def calculate_overlapping_length(x1, x2, x_len1, x_len2):
    x_overlap = max(0, min(x1 + x_len1, x2 + x_len2) - max(x1, x2))
    return x_overlap

# class to store common properties of chiplet

class Utils:
    def __init__(self, common_dict):
        self.package_x_len = common_dict['x_length']
        self.package_y_len = common_dict['y_length']
        self.package_z_len = common_dict['z_length']

        # self.is_homogeneous_chiplets = common_dict['homogeneous_chiplets']

        self.chiplet_x_len = common_dict['chiplet_x']
        self.chiplet_y_len = common_dict['chiplet_y']
        self.chiplet_spacing = common_dict['chiplet_spacing']

        self.num_chiplet_x = common_dict['n_chiplet_x']
        self.num_chiplet_y = common_dict['n_chiplet_y']

        self.top_htc = common_dict['bc_top_htc']
        self.bottom_htc = common_dict['bc_bottom_htc']
