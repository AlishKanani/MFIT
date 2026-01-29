"""
Mesh geometry utilities for building 3D visualizations from node data.
"""


def collect_body_data(layers, package, temperature_all_map, index_heatmap,
                      ambient_temp_c=26.85, ambient_tol=1e-6):
    """
    Collect all nodes organized by body membership.
    
    Args:
        layers: List of Layer_chiplet instances
        package: Chiplet_package instance
        temperature_all_map: Temperature array (nodes x timesteps) in Celsius
        index_heatmap: Timestep index to visualize
    
    Returns:
        Dictionary mapping body_name -> list of node data dicts
    """
    # Get body membership from geometry dict
    body_membership = package.geometry_dict.get('body_membership', {})
    
    # If no body membership defined, use layer names as body names
    # This gives each layer its own toggle in the 3D view
    if not body_membership:
        body_membership = {layer.layer_name: layer.layer_name for layer in layers}
    
    # Collect all nodes organized by body
    body_data = {}  # body_name -> list of node data
    
    num_nodes = 0
    for layer in layers:
        layer_start = num_nodes
        num_nodes += layer.layer_total_nodes()
        
        # Determine which body this layer belongs to
        body_name = body_membership.get(layer.layer_name, layer.layer_name)
        
        if body_name not in body_data:
            body_data[body_name] = []
        
        # Collect node information
        if layer.is_layer_under_chiplet() and not layer.args.is_homogeneous:
            for i in range(layer.total_nodes):
                node = layer.nodes[i]
                global_temp_idx = layer_start + i
                temp = temperature_all_map[global_temp_idx, index_heatmap]
                is_ambient = abs(temp - ambient_temp_c) <= ambient_tol
                
                body_data[body_name].append({
                    'x_min': node.x,
                    'x_max': node.x + node.x_length,
                    'y_min': node.y,
                    'y_max': node.y + node.y_length,
                    'z_min': layer.start_point['z'],
                    'z_max': layer.start_point['z'] + layer.thickness,
                    'temp': temp,
                    'display_temp': temp,
                    'is_ambient': is_ambient,
                    'global_idx': global_temp_idx,
                    'layer_name': layer.layer_name,
                    'body_name': body_name,
                })
        else:
            for i in range(layer.total_x_nodes):
                for j in range(layer.total_y_nodes):
                    node = layer.nodes[i][j]
                    global_temp_idx = layer_start + i * layer.total_y_nodes + j
                    temp = temperature_all_map[global_temp_idx, index_heatmap]
                    is_ambient = abs(temp - ambient_temp_c) <= ambient_tol
                    
                    body_data[body_name].append({
                        'x_min': node.x,
                        'x_max': node.x + node.x_length,
                        'y_min': node.y,
                        'y_max': node.y + node.y_length,
                        'z_min': layer.start_point['z'],
                        'z_max': layer.start_point['z'] + layer.thickness,
                        'temp': temp,
                        'display_temp': temp,
                        'is_ambient': is_ambient,
                        'global_idx': global_temp_idx,
                        'layer_name': layer.layer_name,
                        'body_name': body_name,
                    })
    
    return body_data


def filter_nodes_by_clip(nodes, clip_val, clip_axis):
    """
    Filter nodes and clip them at the exact cut position.
    
    The clipping plane acts as a cut through the geometry:
    - Nodes completely behind (max <= clip_val) are removed
    - Nodes completely in front (min >= clip_val) are kept unchanged  
    - Nodes intersecting the plane are clipped at the plane position
    
    Args:
        nodes: List of node dictionaries with x_min/max, y_min/max, z_min/max
        clip_val: Clipping plane coordinate value (mm)
        clip_axis: 'x', 'y', or 'z'
    
    Returns:
        Filtered list of node dictionaries with clipped coordinates
    """
    filtered = []
    for node in nodes:
        if clip_axis == 'x':
            if node['x_max'] > clip_val:  # Node extends past clip plane
                clipped_node = node.copy()
                # Store original center for flux vectors before clipping
                clipped_node['original_x_center'] = 0.5 * (node['x_min'] + node['x_max'])
                clipped_node['original_y_center'] = 0.5 * (node['y_min'] + node['y_max'])
                clipped_node['original_z_center'] = 0.5 * (node['z_min'] + node['z_max'])
                # Clip the node at the plane
                if node['x_min'] < clip_val:
                    clipped_node['x_min'] = clip_val  # Truncate at plane
                filtered.append(clipped_node)
        elif clip_axis == 'y':
            if node['y_max'] > clip_val:
                clipped_node = node.copy()
                clipped_node['original_x_center'] = 0.5 * (node['x_min'] + node['x_max'])
                clipped_node['original_y_center'] = 0.5 * (node['y_min'] + node['y_max'])
                clipped_node['original_z_center'] = 0.5 * (node['z_min'] + node['z_max'])
                if node['y_min'] < clip_val:
                    clipped_node['y_min'] = clip_val
                filtered.append(clipped_node)
        elif clip_axis == 'z':
            if node['z_max'] > clip_val:
                clipped_node = node.copy()
                clipped_node['original_x_center'] = 0.5 * (node['x_min'] + node['x_max'])
                clipped_node['original_y_center'] = 0.5 * (node['y_min'] + node['y_max'])
                clipped_node['original_z_center'] = 0.5 * (node['z_min'] + node['z_max'])
                if node['z_min'] < clip_val:
                    clipped_node['z_min'] = clip_val
                filtered.append(clipped_node)
    return filtered


def build_mesh_from_nodes(nodes, hover_formatter=None):
    """
    Build Plotly Mesh3d data from filtered node list with clipping support.
    
    For clipped nodes, adds a flat face at the clipping plane.
    
    Args:
        nodes: List of node dictionaries (may include 'clipped' flag)
    
    Returns:
        Tuple of (vertices_x, vertices_y, vertices_z, intensities, 
                 i_indices, j_indices, k_indices)
    """
    vertices_x = []
    vertices_y = []
    vertices_z = []
    intensities = []
    i_indices = []
    j_indices = []
    k_indices = []
    vertex_texts = [] if hover_formatter is not None else None
    
    vertex_count = 0
    for node in nodes:
        x0, x1 = node['x_min'], node['x_max']
        y0, y1 = node['y_min'], node['y_max']
        z0, z1 = node['z_min'], node['z_max']
        
        display_temp = node.get('display_temp', node['temp'])

        # Define 8 vertices of the (possibly clipped) box
        box_vertices_x = [x0, x1, x1, x0, x0, x1, x1, x0]
        box_vertices_y = [y0, y0, y1, y1, y0, y0, y1, y1]
        box_vertices_z = [z0, z0, z0, z0, z1, z1, z1, z1]
        
        vertices_x.extend(box_vertices_x)
        vertices_y.extend(box_vertices_y)
        vertices_z.extend(box_vertices_z)
        intensities.extend([display_temp] * 8)

        if hover_formatter is not None:
            hover_text = hover_formatter(node)
            vertex_texts.extend([hover_text] * 8)
        
        # Define all 12 triangular faces for the (possibly clipped) box
        # The clipped box naturally has faces at the clipping plane
        # Vertices: 0-3 bottom (z0), 4-7 top (z1)
        # 0:(x0,y0,z0), 1:(x1,y0,z0), 2:(x1,y1,z0), 3:(x0,y1,z0)
        # 4:(x0,y0,z1), 5:(x1,y0,z1), 6:(x1,y1,z1), 7:(x0,y1,z1)
        
        # Bottom face (z=z0)
        i_indices.extend([vertex_count+0, vertex_count+0])
        j_indices.extend([vertex_count+1, vertex_count+2])
        k_indices.extend([vertex_count+2, vertex_count+3])
        
        # Top face (z=z1)
        i_indices.extend([vertex_count+4, vertex_count+4])
        j_indices.extend([vertex_count+6, vertex_count+7])
        k_indices.extend([vertex_count+5, vertex_count+6])
        
        # Front face (y=y0) - includes clipping face if clipped on Y
        i_indices.extend([vertex_count+0, vertex_count+0])
        j_indices.extend([vertex_count+1, vertex_count+5])
        k_indices.extend([vertex_count+5, vertex_count+4])
        
        # Back face (y=y1)
        i_indices.extend([vertex_count+2, vertex_count+2])
        j_indices.extend([vertex_count+3, vertex_count+7])
        k_indices.extend([vertex_count+7, vertex_count+6])
        
        # Left face (x=x0) - includes clipping face if clipped on X
        i_indices.extend([vertex_count+0, vertex_count+0])
        j_indices.extend([vertex_count+3, vertex_count+7])
        k_indices.extend([vertex_count+7, vertex_count+4])
        
        # Right face (x=x1)
        i_indices.extend([vertex_count+1, vertex_count+1])
        j_indices.extend([vertex_count+2, vertex_count+6])
        k_indices.extend([vertex_count+6, vertex_count+5])
        
        vertex_count += 8
    
    if vertex_texts is None:
        vertex_texts = []
    return vertices_x, vertices_y, vertices_z, intensities, i_indices, j_indices, k_indices, vertex_texts


def build_edges_from_nodes(nodes):
    """
    Build edge line data from node list for wireframe visualization with clipping support.
    
    For clipped nodes, shows edges of the clipped box including edges on the flat clipping face.
    
    Args:
        nodes: List of node dictionaries (may include 'clipped' flag)
    
    Returns:
        Tuple of (edge_x, edge_y, edge_z) lists with None separators
    """
    edge_x = []
    edge_y = []
    edge_z = []
    
    for node in nodes:
        x0, x1 = node['x_min'], node['x_max']
        y0, y1 = node['y_min'], node['y_max']
        z0, z1 = node['z_min'], node['z_max']
        
        # The clipped node already has adjusted bounds (x0/y0/z0 set to clip_val)
        # So we just draw all 12 edges of the (possibly clipped) box
        edges = [
            # Bottom face edges
            [(x0,y0,z0), (x1,y0,z0)], [(x1,y0,z0), (x1,y1,z0)],
            [(x1,y1,z0), (x0,y1,z0)], [(x0,y1,z0), (x0,y0,z0)],
            # Top face edges
            [(x0,y0,z1), (x1,y0,z1)], [(x1,y0,z1), (x1,y1,z1)],
            [(x1,y1,z1), (x0,y1,z1)], [(x0,y1,z1), (x0,y0,z1)],
            # Vertical edges
            [(x0,y0,z0), (x0,y0,z1)], [(x1,y0,z0), (x1,y0,z1)],
            [(x1,y1,z0), (x1,y1,z1)], [(x0,y1,z0), (x0,y1,z1)]
        ]
        
        for edge in edges:
            edge_x.extend([edge[0][0], edge[1][0], None])
            edge_y.extend([edge[0][1], edge[1][1], None])
            edge_z.extend([edge[0][2], edge[1][2], None])
    
    return edge_x, edge_y, edge_z

