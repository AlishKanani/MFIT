import common
from power_class import Power_grid
from package_class import Chiplet_package
from tools.case_stats import write_case_stats, CaseStatsInputs
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu
import time
import argparse
import os

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

    parser.add_argument('--use_tuned_C', type=lambda x: (str(x).lower() in ['true','1', 'yes']), default=False, help='Use tuned C matrix for simulation')

    # Visualization arguments
    parser.add_argument('--generate_floorplan', type=lambda x: (str(x).lower() in ['true','1', 'yes']), default=False, help='Generate floorplan images (slow for large node counts)')
    parser.add_argument('--generate_2d_heatmap', type=lambda x: (str(x).lower() in ['true','1', 'yes']), default=True, help='Generate 2D layer heatmaps of final temperature')
    parser.add_argument('--time_heatmap', type=float, default=4, help='Time for heatmap generation in sec')
    parser.add_argument('--vertical_planes', type=str, default='', help='Comma-separated list of vertical planes to generate (XZ, YZ, or XZ,YZ)')
    parser.add_argument('--xz_cuts', type=str, default='', help='Comma-separated Y values for XZ plane cuts (e.g., 1.5,3.0,5.5)')
    parser.add_argument('--yz_cuts', type=str, default='', help='Comma-separated X values for YZ plane cuts (e.g., 2.0,4.0)')
    parser.add_argument('--interactive_heatmaps', type=lambda x: (str(x).lower() in ['true','1', 'yes']), default=True, help='Generate interactive HTML heatmaps for vertical cuts')
    parser.add_argument('--generate_3d_heatmap', type=lambda x: (str(x).lower() in ['true','1', 'yes']), default=False, help='Generate interactive 3D Plotly visualization (large file size)')

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    # load material properties from yaml file
    start = time.time()

    args = parse_args()
    
    # ---- Output layout ----
    # Top-level: <output_dir>/{floorplan, heatmaps, ...}
    # Simulation outputs: <output_dir>/output/{RC,DSS}/...
    args.output_root = os.path.join(args.output_dir, "output")
    args.output_rc_dir = os.path.join(args.output_root, "RC")
    args.output_dss_dir = os.path.join(args.output_root, "DSS")
    os.makedirs(args.output_rc_dir, exist_ok=True)
    os.makedirs(args.output_dss_dir, exist_ok=True)

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

    t0 = time.time()
    package.create_layers()
    create_layers_time = time.time() - t0

    t0 = time.time()
    package.connect_nodes()
    connect_nodes_time = time.time() - t0

    if args.generate_floorplan:
        t0 = time.time()
        package.generate_floorplan()
        floorplan_time = time.time() - t0
    else:
        floorplan_time = 0.0

    # Python SuperLU solver (fast, sparse-aware)
    t0 = time.time()
    package.set_initial_conditions()
    init_time = time.time() - t0
    
    total_duration = args.total_duration
    dt = args.time_step
    n_steps = int(total_duration/dt)
    
    # Build A and B (already computed by connect_nodes())
    A = package.cont_A.copy()
    B = package.cont_B.copy()
    
    # Factorize A once (constant across timesteps)
    print(f"Factorizing system matrix ({package.package_total_nodes()} nodes)...")
    factor_start = time.time()
    A_csc = csc_matrix(A)
    slu = splu(A_csc)
    factor_time = time.time() - factor_start
    
    def solve_A(rhs):
        return slu.solve(rhs)
    
    # Base power resampling
    power_steps = package.power.shape[1]
    P_base = package.power  # shape (N, power_steps)
    
    def power_col_for_step(k):
        idx = int(((k+1)*dt)/float(args.power_interval))
        if idx < 0:
            idx = 0
        if idx >= power_steps:
            idx = power_steps - 1
        return idx
    
    # State (deltaT) vector
    T = package.temperature_all.copy()
    package.temperature_all_save.append(T.copy())
    
    # Time-stepping loop
    print(f"Running {n_steps} timesteps...")
    solve_start = time.time()
    for k in range(n_steps):
        base_idx = power_col_for_step(k)
        Pk = P_base[:, base_idx].copy()
        
        rhs = T + B * Pk
        T = solve_A(rhs)
        package.temperature_all_save.append(T.copy())
        
        if (k+1) % 10 == 0 or k == n_steps-1:
            print(f"  Step {k+1}/{n_steps}")
    
    solve_time = time.time() - solve_start
    
    # Save results
    t0 = time.time()
    package.write_temperature_to_file(dt)
    save_time = time.time() - t0
    
    total_time = time.time() - start
    print(f'\nTotal simulation time: {total_time:.3f}s')
    
    # ---- Case stats (written by simulator, not run scripts) ----
    try:
        rc_timing = getattr(package, "rc_timing", {}) or {}
        timing_summary = {
            "elapsed_total_s": float(total_time),
            "create_layers_s": float(create_layers_time),
            "xy_conductance_s": float(rc_timing.get("xy_conductance_s", 0.0)),
            "z_conductance_s": float(rc_timing.get("z_conductance_s", 0.0)),
            "connect_nodes_s": float(connect_nodes_time),
            "generate_floorplan_s": float(floorplan_time),
            "set_initial_conditions_s": float(init_time),
            "factorization_s": float(factor_time),
            "solver_loop_s": float(solve_time),
            "save_outputs_s": float(save_time),
        }
        stats_path = write_case_stats(
            output_dir=args.output_rc_dir,
            args=args,
            package=package,
            geometry_dict=geometry_dict,
            inputs=CaseStatsInputs(
                elapsed_wall_s=total_time,
                simulated_steps=n_steps,
                dt_s=dt,
                timing=timing_summary,
            ),
        )
        print(f"\nWrote case stats: {stats_path}")
    except Exception as e:
        print(f"\nWarning: failed to write case stats: {e}")
        stats_path = None
