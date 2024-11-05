from compas_fd.solvers import fd_constrained_numpy
from compas.geometry import distance_point_point
import numpy as np
import warnings as Warning

class Network_custom(object):
    def __init__(self):
        self.vertices = []              # Initial vertices of the network in equilibrium
        self.vertices_scaled = []       # Scaled vertices of the network
        self.edges = []                 # Edges of the network
        self.q = []                     # Force densities of the edges
        self.f = []                     # Forces in the edges
        self.fixed = []                 # Fixed vertices
        self.loads = []                 # Loads on the vertices
        self.l0 = []                    # Initial lengths of the edges
        self.l1 = []                    # Stressed lengths of the edges
        self.ls1 = []                   # Lengths of the edges after scaling
        self.R  = []                    # Radius of the arcs
        self.th = []                    # Angle of the arcs
        self.xyz_vec= []                # Points on the arcs
        self.path = []                  # Printable paths of connected edges
        self.paths_xyz = []             # Points on the paths
        self.constraints = []           # vertices that are constraint to a line (not tested yet)
        self.dir = []                   # Direction of the arcs either 1 or -1
        self.n_split = 5                # Number of points to split the arc into
        self.intersections = []         # list of list of vertices that are also in one or more previous paths

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
        self.f = np.array(result.forces).reshape(-1)
        self.l1 = np.array(result.lengths).reshape(-1)

    def is_closed(self, path):
        """Check if path is closed."""
        if self.edges[path[0]][0] == self.edges[path[-1]][1]:
            return True
        else:
            return False
    


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
        net.f = np.array(result.forces).reshape(-1)
        net.l1 = np.array(result.lengths).reshape(-1)
        
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
    

    def scale_vertices(self, reference_point, l_scalar = None):
        """Scale the network based on a reference point and a scale factor.
        parameters:
        reference_point: list or numpy array
            Coordinates of the reference point
        scale_factor: float
            Scale factor
        """
        if l_scalar is None:
            l_scalar = self.l_scalar
        reference_point = np.array(reference_point)
        self.vertices_scaled = reference_point + l_scalar * (self.vertices - reference_point)
        return self.vertices_scaled

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
            D = self.l_scalar * self.l1
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
        for edge, R_i, dir_i in zip(self.edges, R, dir):
            p0 = self.vertices_scaled[edge[0]]
            p1 = self.vertices_scaled[edge[1]]
            arc_points = self._points_on_arc(p0[:2], p1[:2], R_i, dir_i, n)
            p = np.stack((arc_points[:,0], arc_points[:,1], np.zeros(n)), axis = 1)
            xyz.append(p)
        self.xyz_vec = xyz
        
        # Update paths and intersections 
        self._update_intersections()
        self._update_paths()
        return xyz
            
    def _points_on_arc(self, p0, p1, R, dir = 1, n = 5):
        # import matplotlib.pyplot as plt
        """Calculate n points on an arc defined by a starting point, ending point, a radius, and a direction.
        input:
            p0:         the starting coordinates of the arc
            p1:         the ending coordinates of the arc
            R:          the radius of the arc
            n:          the number of points on the arc
            dir:        if -1, the arc is in the opposite direction
        output:
            arc_points: the points on the arc
        """
        # Step 1: Midpoint
        x0, y0 = p0
        x1, y1 = p1
        mid = np.array([(x0 + x1) / 2, (y0 + y1) / 2])
        
        # Step 2: Perpendicular vector
        vec = np.array([x1 - x0, y1 - y0])
        perp_vec = dir * np.array([-vec[1], vec[0]])
        
        # Step 3: Normalize and find center point
        perp_unit = perp_vec / np.linalg.norm(perp_vec)
        center = mid + np.sqrt(R**2 - np.linalg.norm(mid - np.array([x0, y0]))**2) * perp_unit
        
        # Step 4: Start and end angles
        theta0 = np.arctan2(y0 - center[1], x0 - center[0])
        theta1 = np.arctan2(y1 - center[1], x1 - center[0])

        # Step 5: Generate n angles between theta0 and theta1. If the angle difference is larger than pi, the arc is in the opposite direction
        angle_diff = theta1 - theta0
        if np.abs(theta0-theta1) > np.pi:
            angle_diff = 2*np.pi - np.abs(angle_diff)
            angles = np.linspace(theta0, theta0 - angle_diff, n)
        else:
            angles = np.linspace(theta0-theta1, 0, n) + theta1
        
        # Step 6: Calculate points
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

        p0 = self.vertices_scaled[self.edges[edge_number][0]]
        p1 = self.vertices_scaled[self.edges[edge_number][1]]
        arc_points = self._points_on_arc(p0[:2], p1[:2], self.R[edge_number], self.dir[edge_number], n)
        self.xyz_vec[edge_number] = np.stack((arc_points[:,0], arc_points[:,1], np.zeros(n)), axis = 1)

        # Update paths and intersections 
        self._update_intersections()
        self._update_paths()
        return
    
    def flip_curves(self, edges, dir = None, n = None):
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

        if dir is None:
            dir_vec = [None] * len(edges)
            
        for edge, dir in zip(edges, dir_vec):
            self.flip_curve(edge, dir, n)
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
                loop_points = self.find_loop_points(p0, angle0, L_loop, alpha_loop, n_points)
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
            Type of plot to display. Options are 'equilibrium', 'scaled' and arcs.
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
            if plot_type == 'equilibrium':
                vtexts = [go.Scatter3d(x=[pos[0]], y=[pos[1]], z=[pos[2]], text=[label], mode='text') for label, pos in enumerate(self.vertices.tolist())]
            else:
                vtexts = [go.Scatter3d(x=[pos[0]], y=[pos[1]], z=[pos[2]], text=[str(i)], mode='text') for i, pos in enumerate(self.vertices_scaled.tolist())]
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

def replace_brackets(line, temperature_settings):
    """"Replace the brackets and the brackets contents with the corresponding variable from the dictionary"""
    for key, value in temperature_settings.items():
        key_b = '[' + key + ']'
        if key_b in line:
            line = line.replace(key_b, str(value))
    return line