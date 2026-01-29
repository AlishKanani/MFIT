import common 
from layer_class import Layer_chiplet
import numpy as np
import os
import time
from scipy.signal import lti
from scipy.sparse import coo_matrix, diags, identity
from package_visualizer import PackageVisualizer

class Chiplet_package:
    def __init__(self, material_properties, geometry_dict, power_grid_class, args, utils=common.Utils):
        self.material_properties = material_properties
        self.geometry_dict = geometry_dict
        self.layers = []
        self.common_utils = utils
        self.power_grid_class = power_grid_class
        self.args = args
        self.visualizer = None  # Will be initialized after layers are created
    
    def create_layers(self):
        for layer in self.geometry_dict['layers']:
            self.layers.append(Layer_chiplet(name=layer, 
                                             layer_dict=self.geometry_dict['layers'][layer], 
                                             material_properties=self.material_properties,
                                             power_grid_class=self.power_grid_class,
                                             args=self.args))
        
        # Create nodes for each layer
        for layer in self.layers:
            layer.create_nodes(utils=self.common_utils, 
                             material_properties=self.material_properties)
    
    def connect_nodes(self):
        # Sparse assembly: build conductance matrix from edge lists
        # This avoids dense N×N allocation and enables scaling to 100k+ nodes
        
        # Coarse timing for case stats (seconds)
        rc_timing = {
            "xy_conductance_s": 0.0,
            "z_conductance_s": 0.0,
        }

        # connects nodes and calculate the RC between the nodes, and convection resistance
        for layer in self.layers:
            layer.connect_nodes()

        self.shared_conductance_list = []
        self.shared_overlap_area_list = []
        # connect the layers for Z direction
        t0 = time.time()
        for i in range(len(self.layers)-1):
            bottom_layer = self.layers[i]
            top_layer = self.layers[i+1]
            # Fast path: aligned grids -> diagonal coupling (vector form)
            shared_conductance, overlap_area, meta = self.connect_layers(top_layer, bottom_layer)
            self.shared_conductance_list.append(shared_conductance)
            self.shared_overlap_area_list.append(overlap_area)
        rc_timing["z_conductance_s"] = time.time() - t0

        self.capacitance_all = []
        for layer in self.layers:
            self.capacitance_all.append(layer.get_capacitance())
        
        self.capacitance_all = np.concatenate(self.capacitance_all, axis=0)
        
        # invert the capacitance
        self.capacitance_all = 1/self.capacitance_all

        # ===== Sparse assembly =====
        N = self.package_total_nodes()
        diag_accum = np.zeros(N, dtype=float)
        rows = []
        cols = []
        data = []

        # Layer offsets
        layer_offsets = []
        off = 0
        for L in self.layers:
            layer_offsets.append(off)
            off += L.layer_total_nodes()

        # XY edges + convection (diagonal)
        t0 = time.time()
        for li, L in enumerate(self.layers):
            offL = layer_offsets[li]
            nL = L.layer_total_nodes()
            conv = L.get_convective_conductance()
            if conv is not None and len(conv) == nL:
                diag_accum[offL:offL + nL] += conv

            edges = L.get_xy_edges()
            if edges is None:
                # Dense fallback for layers that didn't produce edges (keeps correctness, may be slow)
                xy_conductance = L.get_conductance()
                if xy_conductance is not None:
                    for i in range(nL):
                        gi = xy_conductance[i]
                        nz = np.nonzero(gi)[0]
                        for j in nz:
                            g = float(gi[j])
                            if g <= 0:
                                continue
                            u = offL + i
                            v = offL + int(j)
                            rows.append(u); cols.append(v); data.append(-g)
                            diag_accum[u] += g
            else:
                u_loc, v_loc, g_loc = edges
                for k in range(len(g_loc)):
                    g = float(g_loc[k])
                    if g <= 0:
                        continue
                    u = offL + int(u_loc[k])
                    v = offL + int(v_loc[k])
                    rows.append(u); cols.append(v); data.append(-g)
                    rows.append(v); cols.append(u); data.append(-g)
                    diag_accum[u] += g
                    diag_accum[v] += g
        rc_timing["xy_conductance_s"] = time.time() - t0

        # Z edges between adjacent layers
        for pi in range(len(self.layers) - 1):
            offB = layer_offsets[pi]
            offT = layer_offsets[pi + 1]
            bottom = self.layers[pi]
            top = self.layers[pi + 1]
            Gb = self.shared_conductance_list[pi]
            if Gb is None:
                continue

            # Diagonal coupling vector (aligned grids)
            if isinstance(Gb, np.ndarray) and Gb.ndim == 1:
                n = min(bottom.layer_total_nodes(), top.layer_total_nodes(), Gb.shape[0])
                for i in range(n):
                    g = float(Gb[i])
                    if g <= 0:
                        continue
                    u = offB + i
                    v = offT + i
                    rows.append(u); cols.append(v); data.append(-g)
                    rows.append(v); cols.append(u); data.append(-g)
                    diag_accum[u] += g
                    diag_accum[v] += g
            elif isinstance(Gb, dict) and Gb.get("type") == "edge_list":
                u_loc = Gb.get("u")
                v_loc = Gb.get("v")
                g_loc = Gb.get("g")
                if u_loc is None or v_loc is None or g_loc is None:
                    continue
                for k in range(len(g_loc)):
                    g = float(g_loc[k])
                    if g <= 0:
                        continue
                    u = offB + int(u_loc[k])
                    v = offT + int(v_loc[k])
                    rows.append(u); cols.append(v); data.append(-g)
                    rows.append(v); cols.append(u); data.append(-g)
                    diag_accum[u] += g
                    diag_accum[v] += g
            else:
                # Dense matrix fallback
                conductance = Gb
                if hasattr(conductance, "shape") and len(conductance.shape) == 2:
                    idxs = np.argwhere(conductance > 0)
                    for ij in idxs:
                        i = int(ij[0]); j = int(ij[1])
                        g = float(conductance[i, j])
                        if g <= 0:
                            continue
                        u = offB + i
                        v = offT + j
                        rows.append(u); cols.append(v); data.append(-g)
                        rows.append(v); cols.append(u); data.append(-g)
                        diag_accum[u] += g
                        diag_accum[v] += g

        G_off = coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()
        self.conductance_all = (G_off + diags(diag_accum, offsets=0, shape=(N, N), format='csr'))

        if self.args.use_tuned_C:
            self.apply_tuned_C()

        Cinv = self.capacitance_all
        A = (identity(N, format='csc') + self.args.time_step * (self.conductance_all @ diags(Cinv, format='csc'))).T
        self.cont_A = A.tocsc()
        self.cont_B = self.args.time_step * Cinv
        self.rc_timing = rc_timing


    def connect_layers(self, top_layer, bottom_layer):
        # Fast path: aligned discretization in x/y => diagonal coupling vectors.
        # This avoids O(n^2) pairwise overlap checks and avoids storing dense nb×nt matrices.
        method_meta = {"method": "dense"}
        
        try:
            tol = 1e-9
            nb = bottom_layer.layer_total_nodes()
            nt = top_layer.layer_total_nodes()
            if nb == nt:
                # Check if node footprints match (x/y coords + lengths). If so, overlap is diagonal.
                if (np.max(np.abs(bottom_layer.x_cordinates - top_layer.x_cordinates)) < tol and
                    np.max(np.abs(bottom_layer.y_cordinates - top_layer.y_cordinates)) < tol and
                    np.max(np.abs(bottom_layer.x_lengths - top_layer.x_lengths)) < tol and
                    np.max(np.abs(bottom_layer.y_lengths - top_layer.y_lengths)) < tol):
                    common_area = np.minimum(bottom_layer.xy_area, top_layer.xy_area)
                    gz_bottom = bottom_layer.z_conductance
                    gz_top = top_layer.z_conductance
                    denom = gz_bottom*top_layer.xy_area + gz_top*bottom_layer.xy_area
                    g = np.zeros(nb, dtype=float)
                    mask = (common_area > 0) & (gz_bottom > 0) & (gz_top > 0) & (denom > 0)
                    g[mask] = (common_area[mask]*gz_bottom[mask]*gz_top[mask]) / denom[mask]
                    method_meta = {"method": "diagonal", "edges": int(np.count_nonzero(g))}
                    return g, common_area, method_meta
        except Exception:
            pass

        # Fast path: rectilinear grids (possibly different resolutions), build sparse edge list.
        try:
            tol = 1e-9
            t0 = time.time()
            def _safe_int(val):
                if val is None:
                    return 0
                try:
                    return int(val)
                except Exception:
                    return 0
            def _node_at(layer, start, ny, i, j):
                if layer.is_layer_under_chiplet() and not layer.args.is_homogeneous:
                    return layer.nodes[start + i * ny + j]
                return layer.nodes[i][j]

            def _grid_edges(layer, start, nx, ny):
                x_starts = np.zeros(nx, dtype=float)
                x_ends = np.zeros(nx, dtype=float)
                y_starts = np.zeros(ny, dtype=float)
                y_ends = np.zeros(ny, dtype=float)
                for i in range(nx):
                    n = _node_at(layer, start, ny, i, 0)
                    x_starts[i] = n.x
                    x_ends[i] = n.x + n.x_length
                for j in range(ny):
                    n = _node_at(layer, start, ny, 0, j)
                    y_starts[j] = n.y
                    y_ends[j] = n.y + n.y_length
                if (np.any(np.diff(x_starts) < -tol) or np.any(np.diff(y_starts) < -tol)):
                    return None
                return x_starts, x_ends, y_starts, y_ends

            def _interval_overlaps(a0, a1, b0, b1):
                i = 0
                j = 0
                overlaps = []
                while i < len(a0) and j < len(b0):
                    left = max(a0[i], b0[j])
                    right = min(a1[i], b1[j])
                    if right > left:
                        overlaps.append((i, j, right - left))
                    if a1[i] <= b1[j]:
                        i += 1
                    else:
                        j += 1
                return overlaps

            def _add_edges_for_range(start_b, nx_b, ny_b, start_t, nx_t, ny_t, u_list, v_list, g_list):
                gb = _grid_edges(bottom_layer, start_b, nx_b, ny_b)
                gt = _grid_edges(top_layer, start_t, nx_t, ny_t)
                if gb is None or gt is None:
                    return False
                xb0, xb1, yb0, yb1 = gb
                xt0, xt1, yt0, yt1 = gt

                t_xy = time.time()
                x_over = _interval_overlaps(xb0, xb1, xt0, xt1)
                y_over = _interval_overlaps(yb0, yb1, yt0, yt1)
                method_meta["interval_overlap_s"] = method_meta.get("interval_overlap_s", 0.0) + (time.time() - t_xy)
                if not x_over or not y_over:
                    return True

                t_edges = time.time()
                for (ixb, ixt, xlen) in x_over:
                    for (iyb, iyt, ylen) in y_over:
                        area = xlen * ylen
                        if area <= 0:
                            continue
                        u = start_b + ixb * ny_b + iyb
                        v = start_t + ixt * ny_t + iyt
                        gz_bottom = bottom_layer.z_conductance[u]
                        gz_top = top_layer.z_conductance[v]
                        if (gz_bottom <= 0.0) or (gz_top <= 0.0):
                            continue
                        denom = gz_bottom * top_layer.xy_area[v] + gz_top * bottom_layer.xy_area[u]
                        if denom <= 0.0:
                            continue
                        g = (area * gz_bottom * gz_top) / denom
                        if g > 0:
                            u_list.append(u)
                            v_list.append(v)
                            g_list.append(g)
                method_meta["edge_build_s"] = method_meta.get("edge_build_s", 0.0) + (time.time() - t_edges)
                return True

            u_list = []
            v_list = []
            g_list = []

            if bottom_layer.is_layer_under_chiplet() and top_layer.is_layer_under_chiplet() and (not self.args.is_homogeneous):
                ranges_b = getattr(bottom_layer, "_chiplet_ranges", None)
                ranges_t = getattr(top_layer, "_chiplet_ranges", None)
                if ranges_b and ranges_t and len(ranges_b) == len(ranges_t):
                    ok = True
                    for (sb, nxb, nyb), (st, nxt, nyt) in zip(ranges_b, ranges_t):
                        if not _add_edges_for_range(sb, int(nxb), int(nyb), st, int(nxt), int(nyt), u_list, v_list, g_list):
                            ok = False
                            break
                    if ok and u_list:
                        method_meta["method"] = "rect_grid"
                        method_meta["edges"] = int(len(u_list))
                        method_meta["total_s"] = time.time() - t0
                        return {
                            "type": "edge_list",
                            "u": np.array(u_list, dtype=np.int64),
                            "v": np.array(v_list, dtype=np.int64),
                            "g": np.array(g_list, dtype=float),
                        }, None, method_meta
            else:
                ranges_b = getattr(bottom_layer, "_chiplet_ranges", None)
                ranges_t = getattr(top_layer, "_chiplet_ranges", None)
                nx_b = _safe_int(getattr(bottom_layer, "total_x_nodes", None))
                ny_b = _safe_int(getattr(bottom_layer, "total_y_nodes", None))
                nx_t = _safe_int(getattr(top_layer, "total_x_nodes", None))
                ny_t = _safe_int(getattr(top_layer, "total_y_nodes", None))

                if ranges_b and not ranges_t:
                    ok = True
                    for (sb, nxb, nyb) in ranges_b:
                        if not _add_edges_for_range(int(sb), int(nxb), int(nyb), 0, nx_t, ny_t, u_list, v_list, g_list):
                            ok = False
                            break
                    if ok and u_list:
                        method_meta["method"] = "rect_grid"
                        method_meta["edges"] = int(len(u_list))
                        method_meta["total_s"] = time.time() - t0
                        return {
                            "type": "edge_list",
                            "u": np.array(u_list, dtype=np.int64),
                            "v": np.array(v_list, dtype=np.int64),
                            "g": np.array(g_list, dtype=float),
                        }, None, method_meta

                if ranges_t and not ranges_b:
                    ok = True
                    for (st, nxt, nyt) in ranges_t:
                        if not _add_edges_for_range(0, nx_b, ny_b, int(st), int(nxt), int(nyt), u_list, v_list, g_list):
                            ok = False
                            break
                    if ok and u_list:
                        method_meta["method"] = "rect_grid"
                        method_meta["edges"] = int(len(u_list))
                        method_meta["total_s"] = time.time() - t0
                        return {
                            "type": "edge_list",
                            "u": np.array(u_list, dtype=np.int64),
                            "v": np.array(v_list, dtype=np.int64),
                            "g": np.array(g_list, dtype=float),
                        }, None, method_meta
            if not ranges_b and not ranges_t:
                nx_b = int(bottom_layer.total_x_nodes)
                ny_b = int(bottom_layer.total_y_nodes)
                nx_t = int(top_layer.total_x_nodes)
                ny_t = int(top_layer.total_y_nodes)
                if _add_edges_for_range(0, nx_b, ny_b, 0, nx_t, ny_t, u_list, v_list, g_list) and u_list:
                    method_meta["method"] = "rect_grid"
                    method_meta["edges"] = int(len(u_list))
                    method_meta["total_s"] = time.time() - t0
                    return {
                        "type": "edge_list",
                        "u": np.array(u_list, dtype=np.int64),
                        "v": np.array(v_list, dtype=np.int64),
                        "g": np.array(g_list, dtype=float),
                    }, None, method_meta
        except Exception:
            pass

        # Dense fallback
        shared_conductance = np.zeros((bottom_layer.layer_total_nodes(), top_layer.layer_total_nodes()))
        overlap_area = np.zeros((bottom_layer.layer_total_nodes(), top_layer.layer_total_nodes()))
        
        for i in range(bottom_layer.layer_total_nodes()):
            for j in range(top_layer.layer_total_nodes()):
                # calculate overlapping area (mm^2)
                common_area = common.calculate_overlapping_area(x1=bottom_layer.x_cordinates[i], y1=bottom_layer.y_cordinates[i], 
                                                                x2=top_layer.x_cordinates[j], y2=top_layer.y_cordinates[j], 
                                                                x_len1=bottom_layer.x_lengths[i], y_len1=bottom_layer.y_lengths[i], 
                                                                x_len2=top_layer.x_lengths[j], y_len2=top_layer.y_lengths[j])
                overlap_area[i, j] = common_area
                if common_area > 0:
                    gz_bottom = bottom_layer.z_conductance[i]
                    gz_top = top_layer.z_conductance[j]
                    if (gz_bottom <= 0.0) or (gz_top <= 0.0):
                        shared_conductance[i, j] = 0.0
                    else:
                        denom = gz_bottom*top_layer.xy_area[j] + gz_top*bottom_layer.xy_area[i]
                        if denom <= 0.0:
                            shared_conductance[i, j] = 0.0
                        else:
                            shared_conductance[i,j] = (common_area*gz_bottom*gz_top)/denom
        return shared_conductance, overlap_area, method_meta
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
    
    def generate_floorplan_visual(self):
        """Generate floorplan using the new visualizer."""
        if self.visualizer is None:
            self.visualizer = PackageVisualizer(self)
        self.visualizer.generate_floorplan_visual()

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
        
        # conductance_all is sparse; densify only for small N.
        if hasattr(self.conductance_all, "toarray"):
            G = self.conductance_all.toarray()
        else:
            G = np.asarray(self.conductance_all, dtype=float)
        
        capacitance_matrix = np.diag(self.capacitance_all)
        A = -(capacitance_matrix @ G)
        B = capacitance_matrix
        C = np.eye(self.package_total_nodes())
        D = np.zeros((self.package_total_nodes(), self.package_total_nodes()))

        l_sys = lti(A, B, C, D)
        d_sys = l_sys.to_discrete(self.args.time_step, method='zoh')

        discrete_A = d_sys.A
        discrete_B = d_sys.B
        
        # Use new output/DSS/ directory structure
        dss_dir = getattr(self.args, "output_dss_dir", os.path.join(self.args.output_dir, "output", "DSS"))
        os.makedirs(dss_dir, exist_ok=True)

        np.savetxt(os.path.join(dss_dir, 'disc_A_matrix.csv'), discrete_A, delimiter=',')
        np.savetxt(os.path.join(dss_dir, 'disc_B_matrix.csv'), discrete_B, delimiter=',')

    def write_temperature_to_file(self, ts):
        self.temperature_all_save = np.array(self.temperature_all_save) + self.common_utils.ambient_temp

        # Use new output/RC/ directory structure
        rc_dir = getattr(self.args, "output_rc_dir", os.path.join(self.args.output_dir, "output", "RC"))
        os.makedirs(rc_dir, exist_ok=True)
        
        # save the temperature to a file
        file_name = os.path.join(rc_dir, f'temperature_all_{ts}.csv')
        np.savetxt(file_name, self.temperature_all_save, delimiter=',')

        
        temperature_all_map = self.temperature_all_save.T

        if self.args.generate_2d_heatmap:
            if self.visualizer is None:
                self.visualizer = PackageVisualizer(self)
            
            index_heatmap = int(self.args.time_heatmap/self.args.time_step)
            temperature_all_map_celsius = temperature_all_map - 273.15
            
            # Generate 2D layer heatmaps using new visualizer
            num_nodes = 0
            for layer in self.layers:
                layer_start = num_nodes
                num_nodes += layer.layer_total_nodes()
                layer_end = num_nodes
                layer_temps = temperature_all_map_celsius[layer_start:layer_end, index_heatmap]
                self.visualizer.plot_layer_heatmap(layer, layer_temps)
            
            # Generate vertical heatmaps if requested
            if hasattr(self.args, 'vertical_planes') and self.args.vertical_planes:
                self.visualizer.generate_vertical_heatmaps(temperature_all_map_celsius, index_heatmap)
        
        # Generate 3D visualization if requested
        if hasattr(self.args, 'generate_3d_heatmap') and self.args.generate_3d_heatmap:
            if self.visualizer is None:
                self.visualizer = PackageVisualizer(self)
            index_heatmap = int(self.args.time_heatmap/self.args.time_step)
            temperature_all_map_celsius = temperature_all_map - 273.15
            self.visualizer.plot_3d_heatmap(temperature_all_map_celsius, index_heatmap)

        num_nodes = 0
        for layer in self.layers:
            layer_start = num_nodes
            num_nodes += layer.layer_total_nodes()
            layer_end = num_nodes
            if layer.is_power_src():
                layer.map_temperature_to_blk(temperature_all_map[layer_start:layer_end], utils=self.common_utils, ts=ts)

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

        # if dss
        if self.args.generate_DSS:
            # export power to csv file - use DSS directory
            dss_dir = getattr(self.args, "output_dss_dir", os.path.join(self.args.output_dir, "output", "DSS"))
            os.makedirs(dss_dir, exist_ok=True)
            np.savetxt(os.path.join(dss_dir, 'power_all.csv'), self.power.T, delimiter=',')

