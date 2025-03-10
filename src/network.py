from compas_fd.solvers import fd_constrained_numpy
from compas.geometry import distance_point_point
from compas.datastructures import Graph
import numpy as np
import warnings as Warning
from src.shapeoptimizer import ShapeOptimizer 
from collections import Counter

class Network_custom(object):
    def __init__(self):
        self.vertices = []              # Initial vertices of the network before updating shape or the final vertices after updating the shape
        self.vertices_2d = []           # Projection of vertices onto a single plane
        self.vertices_optimized = []    # Vertices of the network after optimization
        self.vertices_scaled = []       # Scaled vertices of the network
        self.edges = []                 # Edges of the network
        self.q = []                     # Force densities of the edges
        self.f = []                     # Forces in the edges
        self.fixed = []                 # Fixed vertices
        self.loads = []                 # Loads on the vertices
        self.l0 = []                    # Initial lengths of the edges
        self.l1 = []                    # Stressed lengths of the edges
        self.l1_scaled = []             # Stressed lengths of the edges after scaling
        self.R  = []                    # Radius of the arcs
        self.th = []                    # Angle of the arcs
        self.xyz_vec= []                # Points on the arcs
        self.path = []                  # Printable paths of connected edges
        self.paths_xyz = []             # Points on the paths
        self.constraints = []           # vertices that are constraint to a line (not tested yet)
        self.dir = []                   # Direction of the arcs either 1 or -1
        self.n_split = 5                # Number of points to split the arc into
        self.intersections = []         # list of list of vertices that are also in one or more previous paths
        self.crossings = {}             # Dictionary of crossing edges. The key is the edge number that crosses another edge. The value is the edges that it crosses.
        self.mat_dict = {}              # Dictionary of material properties
        self.ShapeOptimizer = None      # ShapeOptimizer object
        self.result = None              # Result of the optimization
        self.leaf_edges = []            # Leaf edges of the network
        self.leafs = []                 # Leaf vertices of the network

    def get_geometric_edge_keys(self):
        """Get the keys of the geometric edges."""
        self.g_keys = {tuple(edge): idx for idx, edge in enumerate(self.edges)}

    def initialize_shape_optimizer(self, function_type = 'standard',  method = 'L-BFGS-B', params = None, options = {}):
        """Initialize the ShapeOptimizer object.
        parameters:
        function_type: str
            Type of error function ('standard', 'sigmoid', or 'no optimization').
        method: str
            Optimization method to use.
        params: dict
            Parameters for the error function (e.g., a, b for sigmoid).
        """
        self.ShapeOptimizer = ShapeOptimizer(self.vertices_2d, self.edges, self.l0)
        self.ShapeOptimizer.set_optimization_method(method)
        self.ShapeOptimizer.set_options(options)
        if function_type is not None:
            if params is None:
                self.ShapeOptimizer.set_error_function(function_type)
            else:
                self.ShapeOptimizer.set_error_function(function_type, **params)
        return

    def _set_leafs(self):
        """Set the leafs of the network."""
        vertex_list = [edge[0] for edge in self.edges] + [edge[1] for edge in self.edges]
        counter = Counter(vertex_list)
        leafs = [vertex for vertex, count in counter.items() if count == 1]
        leafs = sorted(leafs)
        self.leafs = leafs
        self.leaf_edges = []
        for edge_i, edge in enumerate(self.edges):
            if edge[0] in leafs or edge[1] in leafs:
                self.leaf_edges.append(edge_i)
        
        # for path in self.path:
        #     if not self.is_closed(path):
        #         self.leaf_edges.append(path[0])
        #         self.leaf_edges.append(path[-1])
        #         self.leafs.append(self.edges[path[0]][0])
        #         self.leafs.append(self.edges[path[-1]][1])
        return

    def get_net_projection(self, plane_normal):
    # Normalize the plane normal to ensure it's a unit vector
        plane_normal = plane_normal / np.linalg.norm(plane_normal)
        
        # Define an arbitrary vector to construct the local basis
        arbitrary_vector = np.array([1, 0, 0]) if not np.allclose(plane_normal, [1, 0, 0]) else np.array([0, 1, 0])
        
        # Compute the local x' axis (orthogonal to the plane normal)
        x_prime = np.cross(arbitrary_vector, plane_normal)
        x_prime /= np.linalg.norm(x_prime)
        
        # Compute the local y' axis (orthogonal to both the plane normal and x')
        y_prime = np.cross(plane_normal, x_prime)
        y_prime /= np.linalg.norm(y_prime)
        
        # Project vertices into the local coordinate system
        x_y_coordinates = []
        for vertex in self.vertices:
            # Subtract any offset (e.g., if plane passes through a point other than origin)
            # For simplicity, assuming the plane passes through the origin
            x_prime_coord = np.dot(vertex, x_prime)
            y_prime_coord = np.dot(vertex, y_prime)
            x_y_coordinates.append([x_prime_coord, y_prime_coord])
        self.vertices_2d = np.array(x_y_coordinates)
        return self.vertices_2d

    def set_material(self, material_dict):
        """Set the material properties of the network.
        parameters:
        material_dict: dict
            Dictionary of material properties
        """
        self.mat_dict = material_dict
        return
    
    def num_vertices(self):
        return len(self.vertices)
    
    def num_edges(self):
        return len(self.edges)
    
    def update_shape(self, q, edges = None):
        """Update the vertices of the network based on the force densities.
        parameters:
        q: list of force densities
        edges: list of edges. If not provided, the function will assume that a new q is provided for all edges in the same order as the previous q.
        """
        if edges is None:
            if len(q) != len(self.edges):
                raise ValueError("Number of force densities does not match the number of edges.")
            self.q = np.array(q)
        else:
            for edge in edges:
                self.q[edge] = q[edge]
        
        result = fd_constrained_numpy(
            vertices=self.vertices,
            fixed=self.fixed,
            edges=self.edges,
            forcedensities=q,
            loads=self.loads,
            constraints=self.constraints,
        )

        self.vertices = np.array(result.vertices).reshape(-1, 3)
        self.vertices_2d = self.vertices[..., :2]
        self.f = np.array(result.forces).reshape(-1)
        self.l1 = np.array(result.lengths).reshape(-1)

    def update_shape_nlf(self, q, loads_0, edges = None):
        """Update the vertices of the network based on the force densities.
        parameters:
        q: list of force densities
        edges: list of edges. If not provided, the function will assume that a new q is provided for all edges in the same order as the previous q.
        """
        if edges is None:
            if len(q) != len(self.edges):
                raise ValueError("Number of force densities does not match the number of edges.")
            self.q = np.array(q)
        else:
            for edge in edges:
                self.q[edge] = q[edge]
        loads = np.copy(loads_0)
        while True:
            result = fd_constrained_numpy(
                vertices=self.vertices,
                fixed=self.fixed,
                edges=self.edges,
                forcedensities=q,
                loads=loads,
                constraints=self.constraints,
            )
            loads_new = ... # a function that redefines the loads such that they are perpendicular to the surface of the cylindar
            if np.linalg.norm(loads - loads_new) < 1e-3:
                break
            loads = loads_new

        self.vertices = np.array(result.vertices).reshape(-1, 3)
        self.vertices_2d = self.vertices[..., :2]
        self.f = np.array(result.forces).reshape(-1)
        self.l1 = np.array(result.lengths).reshape(-1)

    def is_closed(self, path):
        """Check if path is closed."""
        return self.edges[path[0]][0] == self.edges[path[-1]][1]

    @classmethod
    def from_fd(cls, vertices, edges, q, fixed, loads = None, constraints = None, paths = None, dir = None):
        """Construct a Network object from the input parameters for the force density solver.
        parameters:
        vertices: list of vertices coordinates
        edges: list of edges as pairs of vertex indices
        q: list of forcedensities of the edges
        fixed: list of indices of fixed vertices
        loads: list of loads on the vertices
        constraints: list of vertex constraints. An element in the list can be constraint to a line defined by a compas.geometry type, e.g. line = Line([5, 0, 0], [5, 10, 0])
        """
        net = cls()
        net.edges = np.array(edges)
        net.q = np.array(q)
        net.dir = np.array(dir)

        if paths is None:
            Warning.warn("No paths provided. Make sure to set them later.")
        else:
            net.path = paths

        
        net.fixed = np.array(fixed)
        if loads is None:
            loads = [[0, 0, 0] for _ in range(len(vertices))]
        net.loads = np.array(loads)
        if constraints is None:
            constraints = [None] * len(vertices)
        net.constraints = constraints
        
        result = fd_constrained_numpy(
            vertices=vertices,
            fixed=fixed,
            edges=edges,
            forcedensities=q,
            loads=loads,
            constraints=constraints,
        )

        net.vertices = np.array(result.vertices).reshape(-1, 3)
        net.vertices_2d = net.vertices[..., :2]
        net.f = np.array(result.forces).reshape(-1)
        net.l1 = np.array(result.lengths).reshape(-1)

        net.get_geometric_edge_keys()
        net._set_leafs()
        return net
    
    def materialize(self, E, A):
        """Materialize the network by assigning material properties to the edges.
        parameters:
        E: list
            Young's modulus
        A: list
            Cross sectional area
        
        Returns
        -------
        l0 : array
            Initial lengths of the edges
        l_scalar : float
            minimum ratio of the initial length to the stressed length
        """
        from scipy.sparse.linalg import spsolve
        from scipy.sparse import diags
        from numpy import ones
        self.E = np.array(E)
        self.A = np.array(A)
        ea = self.E * self.A

        v = len(self.q)
        Q = diags([self.q], [0])
        L = diags([self.l1.flatten()], [0])
        EA = diags([ea.flatten()], [0])
        I = diags(ones(v))
        A = spsolve(EA, Q.dot(L)) + I
        self.l0 = spsolve(A,self.l1)
        self.l_scalar = np.min(self.l0/self.l1)
        return self.l0, self.l_scalar
    
    def materialize_nonlinear(self, A, stress_data, strain_data, interpolation_kind = 'cubic'):
        """Materialize the network by assigning material properties to the edges.
        parameters:
        A: list
            Cross sectional area
        strain_data: np.array
            strain data from tensile test
        strain: np.array

        
        Returns
        -------
        l0 : array
            Initial lengths of the edges
        l_scalar : float
            minimum ratio of the initial length to the stressed length
        """
        from scipy.interpolate import interp1d
        self.A = np.array(A)

        #ensure stress values are ordered, positive, and unique:
        stress_mask = stress_data >= 0
        order = np.argsort(stress_data[stress_mask])
        stress_sorted = stress_data[stress_mask][order]
        strain_sorted = strain_data[stress_mask][order]

        # Find unique stress values and compute mean strain for duplicates
        unique_stress, indices, _ = np.unique(stress_sorted, return_inverse=True, return_counts=True)
        unique_strain = np.zeros_like(unique_stress, dtype=float)

        for i in range(len(unique_stress)):
            unique_strain[i] = np.mean(strain_sorted[indices == i])  # Average strains for duplicate stresses

        stress_to_strain = interp1d(unique_stress, unique_strain, kind=interpolation_kind, bounds_error=False, fill_value=np.nan)

        force = self.q * self.l1  # Vectorized
        stress = force / self.A
        strain = stress_to_strain(stress)
        self.l0 = self.l1 / (1 + strain)
        self.l_scalar = np.min(self.l0/self.l1)
        return self.l0, self.l_scalar

    def optimize_vertices(self):
        """Find and update the optimal vertices of the network based on the error function.
        parameters:
        error_function: function
            Error function to minimize. If None, the standard error function is used.
        initial_guess: numpy array
            Initial guess for the vertices. If None, the current vertices are used.
        method: str
            Optimization method to use. Default is 'L-BFGS-B'.

        Returns
        -------
        Nothing. The vertices are updated in place.
        """
        self.vertices_optimized =  np.zeros_like(self.vertices)
        self.vertices_optimized[...,:2], self.l1, self.result = self.ShapeOptimizer.optimal_vertices()
        self.l_scalar = np.min(self.l0/self.l1)

    def scale_vertices(self, reference_point, l_scalar = None, account_for_leafs = False):
        """Scale the network based on a reference point and a scale factor.
        parameters:
        reference_point: list or numpy array
            Coordinates of the reference point
        scale_factor: float
            Scale factor
        account_for_leafs: bool
            If True, the leafs of the network are updated to ensure they have the correct length after scaling.
        """
        if l_scalar is None:
            l_scalar = self.l_scalar
        reference_point = np.array(reference_point)
        self.vertices_scaled = reference_point + l_scalar * (self.vertices_optimized - reference_point)
        self.l1_scaled = l_scalar * self.l1

        if account_for_leafs:
            if not self.leafs:
                self._set_leafs()
            for leaf_edge, leaf_end_node in zip(self.leaf_edges, self.leafs):
                leaf_nodes = self.edges[leaf_edge]
                leaf_start_node = leaf_nodes[np.where(leaf_nodes != leaf_end_node)[0][0]]
                L_scaled = np.linalg.norm(self.vertices_scaled[leaf_start_node] - self.vertices_scaled[leaf_end_node])
                # L = np.linalg.norm(self.vertices_optimized[leaf_start_node] - self.vertices_optimized[leaf_end_node])
                dL = (self.l0[leaf_edge] - L_scaled)
                dir_vec = self.vertices_scaled[leaf_end_node] - self.vertices_scaled[leaf_start_node]
                dir_vec /= np.linalg.norm(dir_vec)
                if np.isnan(dir_vec).any():
                    dir_vec = np.array([0, 0, 0])
                self.vertices_scaled[leaf_end_node] += dL * dir_vec
                self.l1_scaled[leaf_edge] = np.linalg.norm(self.vertices_scaled[leaf_start_node] - self.vertices_scaled[leaf_end_node])
        return

    def arc_param(self, D = None, L = None, interpolation_points = 1e6):
        """
        Calculate the radius and angle of an arc based on the arc length and the chord length. Note that the solution is found using interpolation to save time.
        parameters:
        --------
        D: float
            Chord length (For us the scaled length l_scalar * l1)
        L: float
            Arc length (For us the initial length l0)
        interpolation_points: int
            Number of interpolation points to use for the interpolation. The higher the number, the more accurate the solution, but the longer the computation time.

        Returns
        -------
        R : array
            Radius of the arc
        th : array
            Angle of the arc
        """
        if D is None:
            D = self.l1_scaled
        if L is None:
            L = self.l0
        
        C   = np.linspace(0,np.pi/2,int(interpolation_points + 1))
        C   = np.delete(C,0)
        DL  = np.sin(C)/C
        Cvec = np.interp(np.divide(D,L),np.flip(DL),np.flip(C))
        self.R       = np.divide(L, 2*Cvec)
        self.th      = np.divide(L, self.R)
        return self.R, self.th

    def arc_points(self, R = None, th = None, dir = None, n = 10):
        """Calculate n points on an arc defined by a starting point, a radius, an angle, and a direction.
        input:
            x0:         the starting coordinates of the arc
            R:          the radius of the arc
            angle:      the angle of the arc in radians
            n:      the number of points on the arc
            dir:        if -1, the arc is in the opposite direction
        output:
            X:          the points on the arc
        """
        if R is None:
            R = self.R
        if th is None:
            th = self.th
        if dir is None:
            dir = self.dir
        xyz = []
        for edge, R_i, th_i, dir_i in zip(self.edges, R, th, dir):
            p0 = self.vertices_scaled[edge[0]]
            p1 = self.vertices_scaled[edge[1]]
            arc_points = self._points_on_arc(p0[:2], p1[:2], R_i, dir_i, n)
            # arc_points = self._points_on_arc2(p0[:2], p1[:2], R_i, th_i, dir_i, n)
            p = np.stack((arc_points[:,0], arc_points[:,1], np.zeros(n)), axis = 1)
            xyz.append(p)
        self.xyz_vec = xyz
        
        # Update paths and intersections 
        self._update_intersections()
        self._update_paths()
        return xyz
    
    def _points_on_arc(self, p0, p1, R, dir=1, n=5, print_parameters = False):
        """Calculate n points on an arc defined by a starting point, ending point, a radius, and a direction.
        input:
            p0:         the starting coordinates of the arc
            p1:         the ending coordinates of the arc
            R:          the radius of the arc
            dir:        if -1, the arc is in the opposite direction
            n:          the number of points on the arc
        output:
            arc_points: the points on the arc
        """
        x0, y0 = p0
        x1, y1 = p1
        vec = np.array([x1 - x0, y1 - y0])
        chord_length = np.linalg.norm(vec)
        
        # Validate radius
        if chord_length > 2 * R:
            raise ValueError("Radius R is too small for the given points.")
        
        # Midpoint and perpendicular
        mid = np.array([(x0 + x1) / 2, (y0 + y1) / 2])
        perp_vec = dir * np.array([-vec[1], vec[0]]) / np.linalg.norm(vec)
        
        # Distance from midpoint to center
        h = np.sqrt(R**2 - (chord_length / 2)**2)
        center = mid + h * perp_vec
        
        # Angles
        theta0 = np.arctan2(y0 - center[1], x0 - center[0])
        theta1 = np.arctan2(y1 - center[1], x1 - center[0])
        
        # Normalize and interpolate angles
        theta0 = np.mod(theta0, 2 * np.pi)
        theta1 = np.mod(theta1, 2 * np.pi)
        if dir == 1 and theta1 < theta0:
            theta1 += 2 * np.pi
        elif dir == -1 and theta0 < theta1:
            theta0 += 2 * np.pi
        angles = np.linspace(theta0, theta1, n)
        if print_parameters:
            return center, theta0, theta1
        # Points on the arc
        arc_points = np.array([[center[0] + R * np.cos(angle), center[1] + R * np.sin(angle)] for angle in angles])
        return arc_points


    def flip_curve(self, edge_number, dir = None, n = None):
        """Flip the direction of the curve defined by the edge number.
        parameters:
        edge_number: int
            Index of the edge to flip
        dir: int
            Direction of the curve. If None, the direction will be flipped. dir = 1 or -1
        n: int
            Number of points to split the arc into. If None, the default value will be used.
        """
        if n is None:
            n = self.n_split
        if dir is None:
            self.dir[edge_number] = -self.dir[edge_number]
        else:
            self.dir[edge_number] = dir

        if isinstance(self.R, list) and not self.R:
            return

        p0 = self.vertices_scaled[self.edges[edge_number][0]]
        p1 = self.vertices_scaled[self.edges[edge_number][1]]
        arc_points = self._points_on_arc(p0[:2], p1[:2], self.R[edge_number], self.dir[edge_number], n)
        # arc_points = self._points_on_arc2(p0[:2], p1[:2], self.R[edge_number], self.th[edge_number], self.dir[edge_number], n)
        self.xyz_vec[edge_number] = np.stack((arc_points[:,0], arc_points[:,1], np.zeros(n)), axis = 1)

        # Update paths and intersections 
        self._update_intersections()
        self._update_paths()
        return
    
    def flip_curves(self, edges = None, dir = None, n = None):
        """Flip the direction of the curves defined by the edge numbers.
        parameters:
        edges: list of int
            Indices of the edges to flip
        dir: list of int or None
            Direction of the curves. If None, the direction will be flipped. dir is a list of 1 and -1s
        n: int or None
            Number of points to split the arcs into. If None, the default value will be used.
        """        
        if n is None:
            n = self.n_split
        
        if edges is None:
            edges = range(len(self.edges))
            
        if dir is None:
            dir_vec = [None] * len(edges)
            
        for edge, dir in zip(edges, dir_vec):
            self.flip_curve(edge, dir, n)
        return
    
    def auto_flip_curves(self, paths = None, n = None):
        """Automatically flip the direction of the curves. Loops paths and directions are set to 1, -1, 1, -1, ...
        parameters:
        paths: list of list of int or None
            List of paths to flip. If None, all paths are flipped.
        n: int or None
            Number of points to split the arcs into. If None, the default value will be used.
        """
        if n is None:
            n = self.n_split
        if paths is None:
            paths = self.path
        if  isinstance(paths[0], int):
            paths = [paths]
        for path in paths:
            for i, edge in enumerate(path):
                dir_goal = 1 if i % 2 == 0 else -1
                if self.dir[edge] != dir_goal:
                    self.flip_curve(edge, dir_goal, n)
        return
    
    def _update_paths(self):
        """Update lists of coordinates constituting printable paths paths based on the new points on the arcs. The first coordinate of each edge is removed to avoid overlap, except for the first edge."""
        self.paths_xyz = []
        for path in self.path:
            path_xyz = self.xyz_vec[path[0]].copy()
            for edge in path[1:]:
                path_xyz = np.concatenate((path_xyz, self.xyz_vec[edge][1:]), axis=0)
            self.paths_xyz.append(path_xyz)
        return
    
    def get_crossings(self):
        """Get the pairs of edges that cross each other in the path.
        returns:
        list of tuples
            Pairs of edges that cross each other
        """
        vertices_dict = {i: self.vertices_2d[i].tolist() + [0] for i in range(len(self.vertices_2d))}
        graph = Graph.from_nodes_and_edges(vertices_dict, self.edges.tolist())
        crossing_pairs = graph.find_crossings()

        # vertices_set_list = [set() for _ in range(len(self.path))]
        # for path, vertex_set  in zip(self.path, vertices_set_list):
        #     for edge in path:
        #         vertex_set.add(self.edges[edge][0])
        #         vertex_set.add(self.edges[edge][1])
        edge_set_list = [set(path) for path in self.path]
        for path, edge_set in zip(self.path, edge_set_list):
            for edge in path:
                edge_set.add(edge)

        crossing_encountered = [[False] * 2 for _ in range(len(crossing_pairs))]
        crossing_encountered = np.array(crossing_encountered)
        all_crossings = []
        for edge_set in edge_set_list:
            path_crossings = []
            for i, (vertex_pair_i, vertex_pair_j) in enumerate(crossing_pairs):
                edge_i = self.g_keys[vertex_pair_i]
                edge_j = self.g_keys[vertex_pair_j]
                if edge_i in edge_set:
                    crossing_encountered[i, 0] = True
                    if crossing_encountered[i, 1]:
                        self.crossings[edge_i] = edge_j
                        path_crossings.append([edge_i, edge_j])
                if edge_j in edge_set:
                    crossing_encountered[i, 1] = True
                    if crossing_encountered[i, 0]:
                        self.crossings[edge_j] = edge_i
                        path_crossings.append([edge_j, edge_i])
            all_crossings.append(path_crossings)
        self.all_crossings = all_crossings
        return

    def jump_at_crossings(self, crossing_width = 1, crossing_height = 1, interpolation_function = None):
        """
        Jump
        at the crossing of the paths. This loop checks how close a point in xyz_vec is to the crossing and changes the z-coordinate of the point based on the distance.
        """
        from compas.geometry import intersection_segment_polyline_xy
        if interpolation_function is None:
            interpolation_function = self.linear_interpolate
        
        self.crossing_width = crossing_width
        self.crossing_height = crossing_height
        self.interpolation_function = interpolation_function

        for crossing_pairs in self.all_crossings:
            for crossing_pair in crossing_pairs:
                xy_0 = self.xyz_vec[crossing_pair[0]][:,:2]
                xy_1 = self.xyz_vec[crossing_pair[1]][:,:2]
                for s0, s1 in zip(xy_0[:-1], xy_0[1:]) :
                    intersection = intersection_segment_polyline_xy([s0, s1], xy_1)
                    if intersection:
                        for k, p in enumerate(self.xyz_vec[crossing_pair[0]]):
                            d = distance_point_point(p, intersection)
                            self.xyz_vec[crossing_pair[0]][k,2] = interpolation_function(d, crossing_width, crossing_height)
                        break
        
        self._update_paths()
        return

    def _update_intersections(self):
        """determine which vertices are also in one or more previous paths.
        """
        vertices_set_list = [set() for _ in range(len(self.path))]
        for path, vertex_set  in zip(self.path, vertices_set_list):
            for edge in path:
                vertex_set.add(self.edges[edge][0])
                vertex_set.add(self.edges[edge][1])

        self.intersections = [[]]
        for i in range(1, len(vertices_set_list)):
            current_set = vertices_set_list[i]
            previous_union = set.union(*vertices_set_list[0:i])
            self.intersections.append(list(current_set.intersection(previous_union)))
        return
    
    def jump_at_intersection(self, intersection_width = 1, intersection_height = 1, interpolation_function = None):
        """
        Jump at the intersection of the paths. This loop checks how close a point in xyz_vec is to the intersection and changes the z-coordinate of the point based on the distance.
        """
        if interpolation_function is None:
            interpolation_function = self.linear_interpolate
        
        self.intersection_width = intersection_width
        self.intersection_height = intersection_height
        self.interpolation_function = interpolation_function

        for path, intersecting_nodes in zip(self.path[1:], self.intersections[1:]):
            for edge_n in path:
                vertex0, vertex1 = self.edges[edge_n]
                if vertex0 in intersecting_nodes and vertex1 in intersecting_nodes:
                    for k, p in enumerate(self.xyz_vec[edge_n]):
                        d = np.min([distance_point_point(p, self.vertices_scaled[vertex0]), distance_point_point(p, self.vertices_scaled[vertex1])])
                        self.xyz_vec[edge_n][k,2] = interpolation_function(d, intersection_width, intersection_height)
                elif vertex0 in intersecting_nodes:
                    for k, p in enumerate(self.xyz_vec[edge_n]):
                        d = distance_point_point(p, self.vertices_scaled[vertex0])
                        self.xyz_vec[edge_n][k,2] = interpolation_function(d, intersection_width, intersection_height)
                elif vertex1 in intersecting_nodes:
                    for k, p in enumerate(self.xyz_vec[edge_n]):
                        d = distance_point_point(p, self.vertices_scaled[vertex1])
                        self.xyz_vec[edge_n][k,2] = interpolation_function(d, intersection_width, intersection_height)
            
        # Closed paths always have an intersection at the start and end. We need to raise the end of these paths as well.
        for path in self.path:
            if self.is_closed(path):
                edge_n = path[-1]
                vertex0, vertex1 = self.edges[edge_n]
                for k, p in enumerate(self.xyz_vec[edge_n]):
                    d = distance_point_point(p, self.vertices_scaled[vertex1])
                    self.xyz_vec[edge_n][k,2] = np.max([interpolation_function(d, intersection_width, intersection_height), self.xyz_vec[edge_n][k,2]])

        self._update_paths()
        return
    
    def linear_interpolate(self, d, dl, z_height):
        """A
        Linear interpolation of the z-coordinate of the point on the network.
        param d: The distance from the point to intersection
        param dl: width of the intersection
        param z_height: the height of the intersection (jump)
        """
        d = np.asarray(d)  # Ensure 'd' is a numpy array
        return np.maximum((-z_height/dl) * d + z_height, 0)
    
    def rotate_vector(self, theta, vector):
        """Rotate a vector by an angle theta.
        parameters:
        theta: float
            Angle in radians
        vector: numpy array
            Vector to rotate
        """
        if len(vector)>2:
            R = np.array([[np.cos(theta), -np.sin(theta), 0],
                        [np.sin(theta),  np.cos(theta), 0],
                        [0, 0, 1]])
        else:
            R = np.array([[np.cos(theta), -np.sin(theta)],
                        [np.sin(theta),  np.cos(theta)]])
        return R.dot(vector)

    def find_loop_points(self, p0, angle0, L_loop = 10, alpha_loop = 20, n_points = 30):
        """ Find the points of a loop based on the starting point, starting angle, length, and angle of the loop.
        parameters:
        p0: list or numpy array
            Starting point
        angle0: float
            Starting angle in radians
        L_loop: float
            Length of the loop
        alpha_loop: float
            Angle of the loop in radians
        n_points: int
            Number of points on the loop
        """
        x0 = np.linspace(0,L_loop, n_points)
        x1 = np.flip(x0)[1:]
        x0 = x0[:-1]
        y0 = np.tan(alpha_loop/2)*x0
        y1 = np.tan(-alpha_loop/2)*x1
        R  = y0[-1]
        th = np.linspace(np.pi/2,-np.pi/2, n_points // 3)[1:]
        xr = R*np.cos(th) + x0[-1]
        yr = R*np.sin(th)

        d1 = np.linalg.norm(np.array([x1,y1]), axis = 0)
        dr = np.linalg.norm(np.array([xr,yr]), axis = 0)

        z0 = np.zeros_like(x0)
        z1 = self.interpolation_function(d1, self.intersection_width, self.intersection_height)
        zr = self.interpolation_function(dr, self.intersection_width, self.intersection_height)

        x_all = np.concatenate([x0, xr, x1])
        y_all = np.concatenate([y0, yr, y1])
        z_all = np.concatenate([z0, zr, z1])


        x_all, y_all = self.rotate_vector(angle0, np.array([x_all, y_all])) + np.array(p0[:2])[:, np.newaxis]
        loop_points = np.array([x_all, y_all, z_all]).T
        return loop_points

    def all_loop_to_path(self, start_loop_bools, end_loop_bools, L_loop, alpha_loop, n_points):
        """Convert all loop points to paths.
        parameters:
        loop_start_bools: list of bool
            List of bools indicating if a loop should be added at the start of the path
        loop_end_bools: list of bool
            List of bools indicating if a loop should be added at the end of the path
        L_loop: float
            Length of the loop
        alpha_loop: float
            Angle of the loop in radians
        n_points: int
            Number of points on the loop
        """
        # Add loop points at the start of the paths
        for i, (path_xyz, start, end) in enumerate(zip(self.paths_xyz, start_loop_bools, end_loop_bools)):
            if start:
                p0 = path_xyz[0]  # First point of the path
                p1 = path_xyz[1]  # Second point of the path to calculate angle
                angle0 = np.arctan2(p1[1] - p0[1], p1[0] - p0[0]) + np.pi  # Reversed direction for start
                loop_points = self.find_loop_points(p0, angle0, L_loop, alpha_loop, n_points)
                # Add loop_points at the start of the path
                self.paths_xyz[i] = np.vstack([loop_points, self.paths_xyz[i]])
            
            if end:
                p0 = path_xyz[-2]
                p1 = path_xyz[-1]
                angle0 = np.arctan2(p1[1] - p0[1], p1[0] - p0[0])
                loop_points = self.find_loop_points(p1, angle0, L_loop, alpha_loop, n_points)
                # Add loop_points at the end of the path
                self.paths_xyz[i] = np.vstack([self.paths_xyz[i], loop_points])
        return
    
    def add_running_start(self, start_bools, end_bools, L_running_start = 10, n_points = 10):
        """Add a running start to the paths.
        parameters:
        start_bools: list of bool
            List of bools indicating if a running start should be added at the start of the path
        end_bools: list of bool
            List of bools indicating if a running start should be added at the end of the path
        L_running_start: float
            Length of the running start
        """
        for i, (path_xyz, start, end) in enumerate(zip(self.paths_xyz, start_bools, end_bools)):
            if start:
                p0 = path_xyz[0]  # First point of the path
                p1 = path_xyz[1]  # Second point of the path to calculate angle
                direction = (p1 - p0)/np.linalg.norm(p1 - p0)
                end_point = p0 - direction * L_running_start
                running_start = np.linspace(end_point, p0, n_points)
                # Add running start at the start of the path
                self.paths_xyz[i] = np.vstack([running_start, self.paths_xyz[i]])

            if end:
                p0 = path_xyz[-2]
                p1 = path_xyz[-1]
                direction = (p1 - p0)/np.linalg.norm(p1 - p0)
                end_point = p1 + direction * L_running_start
                running_start = np.linspace(p1, end_point, n_points)
                # Add running start at the end of the path
                self.paths_xyz[i] = np.vstack([self.paths_xyz[i], running_start])
        return

    # Visualize the network
    def net_plot(self, color=False, elables=False, vlabels=False, plot_type='equilibrium'):
        """Plot the network using Plotly.
        parameters:
        color: bool
            If True, the edges will be colored based on the forces
        elables: bool
            If True, the edges will be labeled
        vlabels: bool
            If True, the vertices will be labeled
        plot_type: str
            Type of plot to display. Options are 'equilibrium', 'optimized', 'scaled', 'arcs', and 'projection'.
        """
        import plotly.graph_objects as go
        
        x, y, z, c = [], [], [], []
        text_positions = []
        elabels = []

        # Determine color scale if needed
        if color:
            minf, maxf = min(self.f), max(self.f)

        if plot_type == 'equilibrium':
            vertices = self.vertices
        elif plot_type == 'optimized':
            vertices = self.vertices_optimized
        elif plot_type == 'projection':
            vertices = np.hstack([self.vertices_2d, np.zeros((self.vertices_2d.shape[0], 1))])
        else:
            vertices = self.vertices_scaled
        
        if plot_type == 'arcs':
            color = False
            for i, (u, v) in enumerate(self.edges):
                x.extend(list(self.xyz_vec[i][:,0])+[None])
                y.extend(list(self.xyz_vec[i][:,1])+[None])
                z.extend(list(self.xyz_vec[i][:,2])+[None])

                if elables:
                    xyz_u = vertices[u]
                    xyz_v = vertices[v]
                    mid_point = [(xyz_u[j] + xyz_v[j]) / 2 for j in range(3)]
                    text_positions.append(mid_point)
                    elabels.append(str(i))

        else:
            for i, (u, v) in enumerate(self.edges):
                xyz_u = vertices[u]
                xyz_v = vertices[v]
                x.extend([xyz_u[0], xyz_v[0], None])
                y.extend([xyz_u[1], xyz_v[1], None])
                z.extend([xyz_u[2], xyz_v[2], None])
                if color:
                    c.extend([self.f[i], self.f[i], 0])
                if elables:
                    mid_point = [(xyz_u[j] + xyz_v[j]) / 2 for j in range(3)]
                    text_positions.append(mid_point)
                    elabels.append(str(i))

        if color:
            line_marker = dict(width=3, color=c, colorscale='Viridis', colorbar=dict(title='Force'), showscale=True)
        else:
            line_marker = dict(width=3, color='black')
        # Create lines
        lines = [go.Scatter3d(x=x, y=y, z=z, mode='lines', line=line_marker)]

        # Create text labels for edges if needed
        if elables:
            etexts = [go.Scatter3d(x=[pos[0]], y=[pos[1]], z=[pos[2]], text=[label], mode='text') for pos, label in zip(text_positions, elabels)]
            lines.extend(etexts)

        if vlabels:
            # if plot_type == 'equilibrium':
            vtexts = [go.Scatter3d(x=[pos[0]], y=[pos[1]], z=[pos[2]], text=[label], mode='text') for label, pos in enumerate(vertices.tolist())]
            # else:
            #     vtexts = [go.Scatter3d(x=[pos[0]], y=[pos[1]], z=[pos[2]], text=[str(i)], mode='text') for i, pos in enumerate(vertices.tolist())]
            lines.extend(vtexts)

        # Setup layout
        layout = go.Layout(title='Plotly Plot', scene=dict(aspectmode='data'))

        # Combine data and layout into a figure
        fig = go.Figure(data=lines, layout=layout)
        # Set axis limits
        # fig.update_layout(scene=dict(
        #     xaxis=dict(nticks=10, range=[-200, 200]),
        #     yaxis=dict(nticks=10, range=[-200, 200]),
        #     zaxis=dict(nticks=10, range=[-200, 200])
        # ))
        # Display the figure
        fig.show()

    def net_plot_mat(self, ax, elables=False, vlabels=False, plot_type='equilibrium', fp = False):
        """Plot the network using Matplotlib.
        parameters:
        color: bool
            If True, the edges will be colored based on the forces
        elables: bool
            If True, the edges will be labeled
        vlabels: bool
            If True, the vertices will be labeled
        plot_type: str
            Type of plot to display. Options are 'equilibrium', 'scaled' and arcs.
        """
        import matplotlib.pyplot as plt
        import matplotlib.patheffects as path_effects

        x, y = [], []
        text_positions = []
        elabels = []

        if plot_type == 'equilibrium':
            vertices = self.vertices
        else:
            vertices = self.vertices_scaled
        
        if plot_type == 'arcs':
            for i, (u, v) in enumerate(self.edges):
                x.extend(list(self.xyz_vec[i][:,0])+[None])
                y.extend(list(self.xyz_vec[i][:,1])+[None])

                if elables:
                    xyz_u = vertices[u]
                    xyz_v = vertices[v]
                    mid_point = [(xyz_u[j] + xyz_v[j]) / 2 for j in range(3)]
                    text_positions.append(mid_point)
                    elabels.append(str(i))

        else:
            for i, (u, v) in enumerate(self.edges):
                xyz_u = vertices[u]
                xyz_v = vertices[v]
                x.extend([xyz_u[0], xyz_v[0], None])
                y.extend([xyz_u[1], xyz_v[1], None])
                if elables:
                    mid_point = [(xyz_u[j] + xyz_v[j]) / 2 for j in range(3)]
                    text_positions.append(mid_point)
                    elabels.append(str(i))
                    ax.text(mid_point[0], mid_point[1], str(i), color='white', fontsize=10, fontweight='bold')
                    ax.text(mid_point[0], mid_point[1], str(i), color='black', fontsize=10)
        # Create lines
        ax.plot(x, y, 'g--', lw=1, alpha=0.6, label = 'Theoretical network')
        if fp:
            for fixed_point in self.fixed:
                coorx, coory = vertices[fixed_point][:2]
                ax.plot(coorx, coory, 'ro', markersize=10, label = 'Fixed point', alpha=0.6)
            
        if vlabels:
            for i, pos in enumerate(vertices.tolist()):
                ax.text(pos[0], pos[1], str(i), color='white', fontsize=10, fontweight='bold')
                ax.text(pos[0], pos[1], str(i), color='black', fontsize=10)
        return ax
    
    def save_network(self, path = None):
        """Save the network to a file.
        parameters:
        path: str
            Path to the file. If None, the network will be saved to the current directory.
        """
        import pickle
        if path is None:
            path = 'network.json'
        with open(path, 'wb') as f:
            pickle.dump(self, f)
        return
    
    def load_stress_strain_curve(self, file_path, convert_to_MPa = True):
        """Load the stress-strain curve from a file.
        parameters:
        """
        import pandas as pd
        df = pd.read_csv(file_path)
        # strain_columns = [col for col in df.columns if 'Strain' in col]
        # stress_columns = [col for col in df.columns if 'Stress' in col]

        # strain_data = np.array(df[strain_columns])[0]
        # stress_data = np.array(df[stress_columns])[0]
        strain_data = np.array(df['Strain'])
        stress_data = np.array(df['Stress (Pa)'])
        # stress measurement is in Pa. Rest of data is in mm and N, so we need to convert stress to N/m^2
        if convert_to_MPa:
            stress_data = stress_data*1e-6

        stress_mask = stress_data >= 0
        return stress_data[stress_mask], strain_data[stress_mask]


    @staticmethod
    def load_network(path = None):
        """Load the network from a file.
        parameters:
        path: str
            Path to the file. If None, the network will be loaded from the current directory.
        """
        import pickle
        if path is None:
            path = 'network.json'
        with open(path, 'rb') as f:
            self = pickle.load(f)
        return self


def replace_brackets(line, temperature_settings):
    """"Replace the brackets and the brackets contents with the corresponding variable from the dictionary"""
    for key, value in temperature_settings.items():
        key_b = '[' + key + ']'
        if key_b in line:
            line = line.replace(key_b, str(value))
    return line