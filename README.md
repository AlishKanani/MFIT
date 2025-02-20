# Thermal RC and DSS models for 2.5D and 3D chiplets

This repository contains the code for the thermal RC and Discrete State Space (DSS) models for 2.5D and 3D chiplets. These models are part of our multi-fidelity thermal modeling for 2.5D and 3D multi-chiplet architectures.

## Installation
The thermal RC model is implemented using a Python wrapper that uses the Backward Euler Solver with [[SuperLU](https://github.com/xiaoyeli/superlu)]. However, the executable is already included in the repository, you'd need to install superlu as dynamic shared library. Follow [INSTALL_SuperLU.md](INSTALL_SuperLU.md) for installation instructions.

To install other dependencies, run the following command:
```bash
pip install -r requirements.txt
```

## Major features
- Thermal RC and DSS models are validated against FEM simulations (using [ANSYS Fluent](https://www.ansys.com/products/fluids/ansys-fluent)) for a range of system sizes and power traces for 2.5D and 3D chiplets.
- Each layer and block can have different granularity of nodes. This allows a dense grid for active power sources and a coarse grid for passive layers. Which makes the simulation faster.
- The thermal RC model uses an Backward Euler Solver with SuperLU, to solve the system of ODEs, which is faster than traditional explicit solvers.
- In 2.5D and 3D chiplets systems, the layers are expected to have different material properties in the x-y-z directions. The model supports anisotropic material properties.
- Each layer can have different material blocks with different material properties, allowing for the modeling of heterogeneous architectures.

## Examples
Thermal RC model supports broadly three types of geometries:

1. Homogeneous 2.5D/3D chiplets regular grid of nodes for each layer.
2. Homogeneous 2.5D/3D chiplets with non uniform grid of nodes for each layer.
3. Heterogeneous 2.5D/3D chiplets.

We've provided examples for each of these configurations for 2x2 - 2.5D chiplets.

<details>
  <summary>1. Homogeneous chiplets, uniform nodes</summary>

  This is most basic configuration where all the chiplets are identical. Check the geometry file [here](example_1_uniform_nodes_homogeneous_chiplets/chiplet_geometry_4_chiplets_uniform_nodes.yml) that nodes are defined as number of nodes in x and y direction for each layer.
  
  ```bash
  python thermal_RC.py --material_prop_file material_prop.yml \
                       --geometry_file example_1_uniform_nodes_homogeneous_chiplets/chiplet_geometry_4_chiplets_uniform_nodes.yml \
                       --power_config_file example_1_uniform_nodes_homogeneous_chiplets/power_dist_config_homogeneous.yml \
                       --power_seq_file example_1_uniform_nodes_homogeneous_chiplets/power_seq_random_4.csv \
                       --output_dir example_1_uniform_nodes_homogeneous_chiplets/
  ``` 

  Other parameters are set to default values.

  From output floorplan or heatmaps notice that the nodes are uniformly distributed in each layer. Still the nodes can have different granularity by changing `x_nodes`, `y_nodes` in the geometry file.
</details>

<details>
  <summary>2. Homogeneous chiplets, non-uniform nodes</summary>

  This configuration is similar to the previous one but nodes are defined as x and y coordinates for some layers. You can also mix the uniform and non-uniform for different layers, checkout the geometry file [here](example_2_non_uniform_nodes_homogeneous_chiplets/chiplet_geometry_4_chiplets_non_uniform.yml). Also notice the power configuration file [here](example_2_non_uniform_nodes_homogeneous_chiplets/power_dist_config_tiles_tx_rx.yml) that power sources are defined as small blocks in the chiplets.
  
  ```bash
  python thermal_RC.py --material_prop_file material_prop.yml \
                       --geometry_file example_2_non_uniform_nodes_homogeneous_chiplets/chiplet_geometry_4_chiplets_non_uniform.yml \
                       --power_config_file example_2_non_uniform_nodes_homogeneous_chiplets/power_dist_config_tiles_tx_rx.yml \
                       --power_seq_file example_2_non_uniform_nodes_homogeneous_chiplets/power_seq_30s_tiles_tx_rs.csv \
                       --output_dir example_2_non_uniform_nodes_homogeneous_chiplets/ \
                       --total_duration 30
  ```
  Other parameters are set to default values. 

  The average temperature of each chiplet is not calculated in this case as the nodes are not uniformly distributed. But read node temperature from the `temperature_all_<time_step>.csv` file. Node numberings are provided in the floorplan images.
</details>

<details>
  <summary>3. Heterogeneous chiplets</summary>
  In this example, we have 3 chiplets, one bigger and two smaller ones. 
  
  When `is_homogeneous` is set to False, the model will ignore chiplet-specific parameters from `geometry_file` and instead use the parameters from `power_config_file` for each chiplet. 
  
  For each layer marked with `under_chiplet: True` in the `geometry_file`, individual blocks need to be defined in the `power_config_file`. In this case, the power sequence should be defined only for blocks with the `layout_blocks:` dictionary.

  Also, the material properties can be overridden for each block using the `material:` key in the `power_config_file`.

  Check the geometry file [here](example_3_heterogeneous_chiplets/chiplet_geometry_3_chiplets_uniform_nodes.yml) and power configuration file [here](example_3_heterogeneous_chiplets/power_dist_config_heterogeneous.yml) for more details.


  ```bash
  python thermal_RC.py --material_prop_file material_prop.yml \
                       --geometry_file example_3_heterogeneous_chiplets/chiplet_geometry_3_chiplets_uniform_nodes.yml \
                       --power_config_file example_3_heterogeneous_chiplets/power_dist_config_heterogeneous.yml \
                       --power_seq_file example_3_heterogeneous_chiplets/power_seq_random_3.csv \
                       --output_dir example_3_heterogeneous_chiplets/ \
                       --is_homogeneous false
  ```

</details>

## Usage

Detailed usage of the thermal RC model is as follows:

```bash
python thermal_RC.py --material_prop_file <path to material property yaml file> \
                     --geometry_file <path to geometry yaml file> \
                     --power_config_file <path to power config yaml file> \
                     --power_seq_file <path to power sequence csv file> \
                     --output_dir <path to output directory> \

                     --simulation_type {transient,steady} \
                     --generate_DSS {True,False} \ # generates A and B for DSS model
                     --is_homogeneous {True,False} \ # True if all chiplets are identical

                     --time_step <Time step for transient simulation in sec> \
                     --power_interval <Power interval to read power_seq_file in sec>\
                     --total_duration <Total duration for transient simulation in sec> \ # should match power sequence length
                     --use_tuned_C {True,False} \ # True if using tuned Capacitance values for each layer

                     --generate_heatmap {True,False} \
                     --time_heatmap <Time for heatmap generation in sec> 

```
Note: For timing analysis, turn off the `generate_heatmap` and `generate_DSS` flags.

### Inputs
There are four main input configuration files which are used to run the thermal RC model:

1. `material_prop_file` : yaml file containing the material properties of the layers in the 2.5D/3D chiplets.
2. `geometry_file` : yaml file containing the geometry of the 2.5D/3D chiplets.
3. `power_config_file` : yaml file containing information about power sources in the 2.5D/3D chiplets.
4. `power_sequence_file` : csv file with power traces for each power source defined in `power_config_file`.

more details about the format of these files is available in the [input_format.md](INPUT_FORMAT.md) file.

### Outputs
The model generates the following outputs:
1. Floorplan of each layer in the 2.5D/3D chiplet with the node numbers in `<output_dir>/floorplan/` directory.
2. Floorplan of power sources in the 2.5D/3D chiplet in `<output_dir>/floorplan/` directory.
3. Temperature of each node in the 2.5D/3D chiplet at each time step. The output is saved in `<output_dir>/output/temperature_all_<time_step>.csv` file.
4. Average temperature of each chiplet at each time step. The output is saved in `<output_dir>/output/temperature_<chiplet_layer>_<time_step>.csv` file.
5. Output file with the A and B matrices and mapped power for all nodes for the DSS model in `<output_dir>/output/` directory.
6. Heatmap of the temperature of each layer at a specific time step in `<output_dir>/heatmaps/` directory.


### DSS Model
If `--generate_DSS` is set to True, the code will generate the A and B matrices for the Discrete State Space (DSS) model, which can be used as the following equation to predict the temperature of all the nodes:

```
T(k+1) = A * T(k) + B * P(k)
```

where `T(k)` is the temperature of a node at step k, `P(k)` is the power input at step k, and `A` and `B` are the matrices generated by the code. More details about the DSS model can be found in our paper:

## Reference
If you use the library in your research, please refer to the following paper:
```bibtex
@article{MFIT_2024,
      title={MFIT: Multi-Fidelity Thermal Modeling for 2.5D and 3D Multi-Chiplet Architectures}, 
      author={Pfromm, Lukas and Kanani, Alish and Sharma, Harsh and Solanki, Parth and Tervo, Eric and Park, Jaehyun and Doppa, Janardhan Rao and Pande, Partha Pratim and Ogras, Umit Y},
      journal={arXiv preprint arXiv:2410.09188},
      year={2024},
      url={https://arxiv.org/abs/2410.09188}, 
}
```
