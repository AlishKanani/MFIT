"""
Shared utilities for computing node-level 3D heat-flux vectors from the RC network.
"""

import numpy as np


def compute_flux_field(layers, package, temperature_all_map, index_heatmap, T_ambient=300.0):
    """
    Compute per-node 3D heat-flux components (qx, qy, qz) in W/m^2.

    This uses the same RC network as the solver:
      - In-layer conductances (layer.xy_conductance)
      - Cross-layer conductances (package.shared_conductance_list)
      - Convection branches (layer.z_plus_conductance/z_minus_conductance)

    Args:
        layers: sequence of Layer_chiplet instances in package order.
        package: Chiplet_package instance (provides shared_conductance_list).
        temperature_all_map: ndarray, shape (nodes, timesteps), absolute temperature in Kelvin.
        index_heatmap: timestep index to sample.
        T_ambient: ambient temperature in Kelvin (default 300.0).

    Returns:
        (qx, qy, qz): three 1D arrays of length N_nodes with flux components in W/m^2.
    """
    if layers is None:
        return None, None, None

    T_global = np.asarray(temperature_all_map[:, index_heatmap], dtype=float)
    n_nodes = T_global.shape[0]

    # Heat flow (W) accumulated per node in physical axes
    Qx = np.zeros(n_nodes, dtype=float)
    Qy = np.zeros(n_nodes, dtype=float)
    Qz = np.zeros(n_nodes, dtype=float)

    # Precompute layer offsets in global indexing
    layer_offsets = []
    off = 0
    for layer in layers:
        layer_offsets.append(off)
        off += layer.layer_total_nodes()

    # ---- In-layer conduction (XY plane) ----
    for layer_idx, layer in enumerate(layers):
        nL = layer.layer_total_nodes()
        if nL == 0:
            continue
        g_off = layer_offsets[layer_idx]
        T_layer = T_global[g_off:g_off + nL]

        xy_cond = getattr(layer, "xy_conductance", None)
        if xy_cond is None:
            continue

        # Local helper for mapping flat index to Node
        def get_node(local_idx):
            if layer.is_layer_under_chiplet() and not layer.args.is_homogeneous:
                return layer.nodes[local_idx]
            i = local_idx // layer.total_y_nodes
            j = local_idx % layer.total_y_nodes
            return layer.nodes[i][j]

        # Scan upper triangle; each non-zero entry represents a contacting neighbor pair
        for i in range(nL):
            for j in range(i + 1, nL):
                G_ij = xy_cond[i, j]
                if G_ij == 0:
                    continue
                Ti = T_layer[i]
                Tj = T_layer[j]
                Q_ij = G_ij * (Ti - Tj)  # W, positive from i -> j

                node_i = get_node(i)
                node_j = get_node(j)

                dx = node_j.x - node_i.x
                dy = node_j.y - node_i.y
                dist = np.hypot(dx, dy)
                if dist == 0.0:
                    continue
                ux = dx / dist
                uy = dy / dist

                # Assign signed components to both nodes so local flux direction
                # matches heat transport (hot -> cold).
                Qx[g_off + i] += ux * Q_ij
                Qx[g_off + j] += ux * Q_ij
                Qy[g_off + i] += uy * Q_ij
                Qy[g_off + j] += uy * Q_ij

    # ---- Cross-layer conduction (Z direction) ----
    shared_list = getattr(package, "shared_conductance_list", [])
    if shared_list:
        for k, shared_cond in enumerate(shared_list):
            bottom_layer = layers[k]
            top_layer = layers[k + 1]

            bottom_off = layer_offsets[k]
            top_off = layer_offsets[k + 1]

            nB = bottom_layer.layer_total_nodes()
            nT = top_layer.layer_total_nodes()

            # shared_cond is either:
            # - dense matrix (nB, nT)
            # - diagonal vector (aligned grids): shared_cond[i] couples bottom i to top i
            if shared_cond is None:
                continue
            if isinstance(shared_cond, np.ndarray) and shared_cond.ndim == 1:
                n = min(nB, nT, shared_cond.shape[0])
                if n <= 0:
                    continue
                Tb = T_global[bottom_off:bottom_off + n]
                Tt = T_global[top_off:top_off + n]
                Q_bt = shared_cond[:n] * (Tb - Tt)  # W, positive from bottom -> top
                Qz[bottom_off:bottom_off + n] += Q_bt
                Qz[top_off:top_off + n] += Q_bt
            else:
                for i in range(nB):
                    for j in range(nT):
                        G_bt = shared_cond[i, j]
                        if G_bt == 0:
                            continue
                        Tb = T_global[bottom_off + i]
                        Tt = T_global[top_off + j]
                        Q_bt = G_bt * (Tb - Tt)  # W, positive from bottom -> top (upward)

                        # Assign the same signed flow to both nodes so direction
                        # matches physical heat transport (upward if Tb > Tt).
                        Qz[bottom_off + i] += Q_bt
                        Qz[top_off + j] += Q_bt

    # ---- Convection branches (Z direction) ----
    for layer_idx, layer in enumerate(layers):
        nL = layer.layer_total_nodes()
        if nL == 0:
            continue
        g_off = layer_offsets[layer_idx]
        T_layer = T_global[g_off:g_off + nL]

        z_plus = getattr(layer, "z_plus_conductance", None)
        z_minus = getattr(layer, "z_minus_conductance", None)

        if z_plus is not None:
            for i in range(nL):
                Gp = z_plus[i]
                if Gp <= 0.0:
                    continue
                Ti = T_layer[i]
                Q_conv = Gp * (Ti - T_ambient)  # W, positive upward to ambient
                Qz[g_off + i] += Q_conv

        if z_minus is not None:
            for i in range(nL):
                Gm = z_minus[i]
                if Gm <= 0.0:
                    continue
                Ti = T_layer[i]
                Q_conv = Gm * (Ti - T_ambient)  # W, positive downward to ambient
                Qz[g_off + i] -= Q_conv

    # ---- Convert heat flow (W) to flux density (W/m^2) per node ----
    qx = np.zeros_like(Qx)
    qy = np.zeros_like(Qy)
    qz = np.zeros_like(Qz)

    for layer_idx, layer in enumerate(layers):
        nL = layer.layer_total_nodes()
        if nL == 0:
            continue
        g_off = layer_offsets[layer_idx]

        def get_node(local_idx):
            if layer.is_layer_under_chiplet() and not layer.args.is_homogeneous:
                return layer.nodes[local_idx]
            i = local_idx // layer.total_y_nodes
            j = local_idx % layer.total_y_nodes
            return layer.nodes[i][j]

        for i in range(nL):
            node = get_node(i)
            idx = g_off + i
            # Side areas for lateral flux (mm^2)
            area_yz = node.get_area_yz()  # x-normal faces
            area_xz = node.get_area_xz()  # y-normal faces
            area_xy = node.get_area_xy()  # z-normal faces

            # Avoid divide-by-zero; if area is zero, leave flux as 0
            if abs(area_yz) > 0.0:
                # Qx is W, area_yz in mm^2 -> convert to W/m^2 via 1e6
                qx[idx] = Qx[idx] / (area_yz * 1e-6)
            if abs(area_xz) > 0.0:
                qy[idx] = Qy[idx] / (area_xz * 1e-6)
            if abs(area_xy) > 0.0:
                qz[idx] = Qz[idx] / (area_xy * 1e-6)

    return qx, qy, qz


