import numpy as np
import warnings as Warning

from scipy.optimize import minimize

class ShapeOptimizer(object):
    """

    Class for optimizing the shape of a network of vertices and edges.

    """

    def __init__(self, vertices, edges, l0):
        """
        Parameters
        ----------
        vertices : ndarray
            Array of vertex coordinates.
        edges : ndarray
            Array of edge indices.
        l0 : ndarray
            Array of edge lengths.
        """
        self.error_function_type = "standard"  # Default error function
        self.error_params = {"a": 1, "b": 1}  # Default parameters for sigmoid
        self.vertices = vertices
        self.edges = edges
        self.l0 = l0
        self.result = None
        self.method = 'L-BFGS-B'
        self.history = []
        self.history_error = []
        self.options = {'maxiter': 1000}

    def set_error_function(self, function_type, **params):
        """
        Set the active error function and its parameters.

        Parameters
        ----------
        function_type : str
            Type of error function ('standard' or 'sigmoid').
        params : dict
            Parameters for the error function (e.g., a, b for sigmoid).
        """
        if function_type not in ["standard", "sigmoid", "no optimization"]:
            raise ValueError("Invalid error function type. Choose 'standard' or 'sigmoid'.")
        self.error_function_type = function_type
        self.error_params.update(params)
    
    def set_optimization_method(self, method):
        """
        Set the optimization method.

        Parameters
        ----------
        method : str
            Optimization method to use.
        """
        self.method = method

    def set_options(self, options):
        """
        Set the maximum number of iterations for the optimization.

        Parameters
        ----------
        setting : dict
            Optimization settings.
        """
        self.options = options

    def _compute_error(self, flat_coords):
        """
        Compute the error using the active error function.

        Parameters
        ----------
        flat_coords : ndarray
            Flattened array of vertex coordinates.

        Returns
        -------
        error : float
            Computed error.
        """
        if self.error_function_type == "standard":
            return self.standard_error_function(flat_coords)
        elif self.error_function_type == "sigmoid":
            return self.sigmoid_error_function(flat_coords, **self.error_params)
        elif self.error_function_type == "no optimization":
            return 0
        
    def optimal_vertices(self, initial_guess = None):
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
        if self.error_function_type == "no optimization":
            l1 = np.array([np.linalg.norm(self.vertices[i] - self.vertices[j]) for i, j in self.edges])
            return self.vertices[...,:2], l1, None
    
        if initial_guess is None:
            initial_guess = self.vertices[...,:2].flatten()
        elif initial_guess.shape == (len(self.vertices), 2):
            initial_guess = initial_guess.flatten()
        elif initial_guess.shape == (len(self.vertices), 3):
            initial_guess = initial_guess[...,:2].flatten()
            Warning.warn("Initial guess has 3 dimensions. Only the first two dimensions were used.")
        elif len(initial_guess) != 2 * len(self.vertices):
            raise ValueError("Initial guess does not match the number of vertices.")

        # Minimize the error

        if self.method == 'Gauss-Seidel':
            coord = self.gauss_seidel_form_finding(initial_guess, self.options)
            l1 = np.array([np.linalg.norm(coord[i] - coord[j]) for i, j in self.edges])
            return coord.reshape(-1, 2), l1, None
        else:
            result = minimize(self._compute_error, initial_guess, method=self.method, callback=self.optimization_callback, options=self.options)
        if not result.success:
            raise RuntimeError(f"Optimization failed: {result.message}")
        vertices = result.x.reshape(-1, 2)
        self.history = np.array(self.history)
        l1 = np.array([np.linalg.norm(vertices[i] - vertices[j]) for i, j in self.edges])
        return vertices, l1, result


    def standard_error_function(self, flat_coords):
        coords = flat_coords.reshape(-1, 2)
        error = 0
        for (i, j), L_ij in zip(self.edges, self.l0):
            dist = np.linalg.norm(coords[i] - coords[j])
            # error += (dist - L_ij) ** 2
            error += ((dist - L_ij) / L_ij) ** 2
        return error

    def sigmoid_error_function(self, flat_coords, a = 1, b = 1):
        """
        Error function with sigmoid weighting to penalize longer/shorter edges more than shorter/longer ones.

        Parameters
        ----------
        flat_coords : ndarray
            Flattened array of vertex coordinates.
        a : float, optional
            Scaling factor for the weight. Default is 1.
        b : float, optional
            Controls the steepness of the sigmoid function. Default is 1. if b is positive, shorter edges are penalized more. If b is negative, longer edges are penalized more.

        Returns
        -------
        error : float
            Weighted error based on the sigmoid function.
        """
        coords = flat_coords.reshape(-1, 2)
        error = 0
        for (i, j), L_ij in zip(self.edges, self.l0):
            dist = np.linalg.norm(coords[i] - coords[j])
            diff = dist - L_ij
            weight = a / (1 + np.exp(-b * diff))
            error += weight * (diff ** 2)
        return error
    
    def gauss_seidel_form_finding(self, initial_guess, options):
        max_iter = options.get('maxiter', 1000)
        damping = options.get('damping', 1.0)
        correction_scalar = options.get('correction_scalar', 1.0)
        tol = options.get('tol', 1e-6)
        coords = initial_guess.reshape(-1, 2)
        prev_error = float('inf')

        for _ in range(max_iter):
            random_order = np.random.permutation(len(self.edges))
            for (i, j), L_ij in zip(self.edges[random_order], self.l0[random_order]):
                diff = coords[i] - coords[j]
                dist = np.linalg.norm(diff)

                if dist > tol/100:
                    # correction = (L_ij - dist) * (diff / dist) / 2  # Move both nodes halfway
                    correction_factor = (L_ij - dist) / L_ij  # Relative error
                    if dist > L_ij:
                        correction_factor *= correction_scalar  # Damping for short edges
                    correction = correction_factor * diff/2
                    coords[i] += correction*damping
                    coords[j] -= correction*damping
            total_error = sum((np.linalg.norm(coords[i] - coords[j]) - L_ij) ** 2 for (i, j), L_ij in zip(self.edges, self.l0))
            # Convergence check
            self.history.append(coords.copy().reshape(-1, 2))
            self.history_error.append(total_error)
            if _ % 100 == 0:
                print(f"Iteration {_}: Current error = {total_error}")
            if abs(prev_error - total_error) < tol:  # Check if error stops decreasing
                print(f"Converged after {_} iterations.")
                break
            prev_error = total_error
        print(f"Final error: {total_error}")
        self.history = np.array(self.history)
        return coords


    def optimization_callback(self, xk):
        """
        Callback function for optimization. Prints the current error and increments the iteration counter. Saves shape history.
        """
        if not hasattr(self, "_iteration"):
            self._iteration = 0  # Initialize counter on first call
        print(f"Iteration {self._iteration}: Current error = {self._compute_error(xk)}")
        self._iteration += 1  # Increment counter
        self.history.append(xk.reshape(-1, 2))
        self.history_error.append(self._compute_error(xk))
    
    def animate_optimization(self, frame_range=None):
        """
        Animate the optimization process.
        """
        import matplotlib.pyplot as plt
        import matplotlib.animation as animation
        fig, ax = plt.subplots()
        ax.set_aspect("equal")
        ax.set_xlim(np.min(self.history[:, :, 0]) - 5, np.max(self.history[:, :, 0]) +5)
        ax.set_ylim(np.min(self.history[:, :, 1]) - 15, np.max(self.history[:, :, 1]) + 5)
        ax.axis('off')

        def new_network_arrays(frame):
            vertex_array_x = []
            vertex_array_y = []
            
            for i, j in self.edges:
                vertices = self.history[frame]
                vertex_array_x.extend([vertices[i, 0], vertices[j, 0], None])
                vertex_array_y.extend([vertices[i, 1], vertices[j, 1], None])
            
            return np.array(vertex_array_x), np.array(vertex_array_y)

        vertex_array_x, vertex_array_y = new_network_arrays(0)
        network_plot, = ax.plot(vertex_array_x, vertex_array_y, 'k-')
        text_plot = ax.text(0.02, 0.01, f"it #: {0}, error: {self.history_error[0]:.4f}", transform=ax.transAxes)

        def update(frame):
            vertex_array_x, vertex_array_y = new_network_arrays(frame)
            network_plot.set_xdata(vertex_array_x)
            network_plot.set_ydata(vertex_array_y)
            text_plot.set_text(f"it #: {frame}, error: {self.history_error[frame]:.4f}")
            return network_plot, text_plot

        if frame_range is None:
            frame_range = range(len(self.history))

        ani = animation.FuncAnimation(fig, update, frames=frame_range, interval=200)
        plt.show()
        return ani