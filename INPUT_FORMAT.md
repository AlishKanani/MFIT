# Input format for all the configuration files

This document describes the input format for each configuration file required to run the thermal RC model. 
There are three main configuration files which are in `yaml` format. `yaml` is dictionary format which is easy to read and write and automate using python dictionary.

1. `material_prop_file` : yaml file containing the material properties of the layers in the 2.5D/3D chiplets.
2. `geometry_file` : yaml file containing the geometry of the 2.5D/3D chiplets.
3. `power_config_file` : yaml file containing information about power sources in the 2.5D/3D chiplets.
4. `power_sequence_file` : csv file with power traces for each power source defined in `power_config_file`.

## 1. Material Properties File

This file contains the material properties of the materials used in the 2.5D/3D chiplets. The file should be in `yaml` format. The file should be formatted as follows:

```yaml
<material_name>:
  density: <density in kg/m^3>
  specific_heat: <specific heat in J/m^3>
  thermal_conductivity: <thermal conductivity in W/m-K> or
  thermal_conductivity: 
    kx: <thermal conductivity in W/m-K along x direction>
    ky: <thermal conductivity in W/m-K along y direction>
    kz: <thermal conductivity in W/m-K along z direction>
```

Thermal conductivity can be defined as scalar or as a dictionary with `kx`, `ky`, and `kz` values. If the material is isotropic, then scalar value can be used. If the material is anisotropic, then the dictionary format should be used.

## 2. Geometry File

The geometry file contains two dictionaries, `common` and `layers`. Common dictionary contains the common properties of package, while layer dictionary contains the properties of each layer in the 2.5D/3D chiplets. The file should be in `yaml` format. The file should be formatted as follows:

#### Common Properties
```yaml
common:
  x_length: <length of the package in x direction in mm>
  y_length: <length of the package in y direction in mm>
  z_length: <height of the package in z direction in mm>
  chiplet_x: <length of the chiplet in x direction in mm> # used only for homogeneous chiplets
  chiplet_y: <length of the chiplet in y direction in mm> # used only for homogeneous chiplets
  chiplet_spacing: <spacing between chiplets in mm> # used only for homogeneous chiplets
  n_chiplet_x: <number of chiplets in x direction> # used only for homogeneous chiplets
  n_chiplet_y: <number of chiplets in y direction> # used only for homogeneous chiplets
  bc_top_htc : <heat transfer coefficient at top boundary in W/m^2-K>
  bc_bottom_htc : <heat transfer coefficient at bottom boundary in W/m^2-K>
```

#### Layer Properties
```yaml

layers:
  <layer_name>:
    thickness: <thickness of the layer in mm>
    nodes: 
      uniform: <True/False> # True if nodes are uniformly spaced and defined as number of nodes in x and y direction. False if nodes are non-uniformly spaced and defined as list of x and y coordinates.
      under_chiplet: <True/False> # True if layer's geometry is similar to chiplet's geometry. For heterogeneous chiplets, individual block's geometry should be defined in the power_config_file.
      x_nodes: <scalar or list of x coordinates of nodes> 
      y_nodes: <scalar or list of y coordinates of nodes>
    start_point: 
      x: <x coordinate of the start point in mm>
      y: <y coordinate of the start point in mm>
      z: <z coordinate of the start point in mm>
    power_src: <True/False> # True if power source is present in the layer. If true, power source should be defined in power_config_file.
    material: <material_name> # material name defined in material_prop_file
```

## 3. Power Configuration File

This file contains the information about power sources in the 2.5D/3D chiplets. 

If `is_homogeneous` is True, then only layers with `power_src` set to True will have to be defined. If `is_homogeneous` is False, then individual blocks should be defined for each layer with `under_chiplet` set to True in the geometry file. The file should be in `yaml` format. The file should be formatted as follows: 

#### Homogeneous Chiplets
```yaml

<layer_name>: # should be same as defined in geometry file
  <chiplet_name>:
    start_chiplet_x: <x coordinate of the start point of the chiplet in mm>
    start_chiplet_y: <y coordinate of the start point of the chiplet in mm>
    
    layout_blocks:
      <block_name>:
        start_point_x: <x coordinate of the start point of the block w.r.t. start_chiplet_x in mm>
        start_point_y: <y coordinate of the start point of the block w.r.t. start_chiplet_y in mm>
        length_x: <length of the block in x direction in mm>
        length_y: <length of the block in y direction in mm>
        max_power: <maximum power in W> # this number is used to normalize the power trace to 100%.
```

#### Heterogeneous Chiplets
```yaml

<layer_name>: # should be same as defined in geometry file
  <chiplet_name>:
    start_chiplet_x: <x coordinate of the start point of the chiplet in mm>
    start_chiplet_y: <y coordinate of the start point of the chiplet in mm>
    length_chiplet_x: <length of the chiplet in x direction in mm> # extra parameter for heterogeneous chiplets
    length_chiplet_y: <length of the chiplet in y direction in mm> # extra parameter for heterogeneous chiplets
    nodes_x: <number of nodes in this chiplet for x direction> # extra parameter for heterogeneous chiplets
    nodes_y: <number of nodes in this chiplet for y direction> # extra parameter for heterogeneous chiplets
    material: <material_name> # optional argument to override the material properties for this chiplet

    layout_blocks: # only for blocks with power sources, will have to provide power traces in power_sequence_file
      <block_name>:
        start_point_x: <x coordinate of the start point of the block w.r.t. start_chiplet_x in mm>
        start_point_y: <y coordinate of the start point of the block w.r.t. start_chiplet_y in mm>
        length_x: <length of the block in x direction in mm>
        length_y: <length of the block in y direction in mm>
        max_power: <maximum power in W> # this number is used to normalize the power trace to 100%.
```

## 4. Power Sequence File

This file contains the power traces for each power source defined in the `power_config_file`. The file should be in `csv` format. Each row should have comma separated values in % of the `max_power`. The tool expects total `--total_duration/--power_interval` columns, ex. with `--total_duration 50 --power_interval 1` tool expects 50 values in each row, each value should be running average of power for that interval. 

The file should be formatted as follows:

```csv
<layer_name from `geometry_file`>_<chiplet_name from `power_config_file`>_<block_name from `power_config_file`>, <power values in % of max_power>, <power values in % of max_power>, ...
```

