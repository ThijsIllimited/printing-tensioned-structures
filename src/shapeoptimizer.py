import numpy as np
import warnings as Warning

from scipy.optimize import minimize

class ShapeOptimizer(object):
    """
    
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
        result = minimize(self._compute_error, initial_guess, method=self.method, callback=self.optimization_callback, options=self.options)
        if not result.success:
            raise RuntimeError(f"Optimization failed: {result.message}")
        vertices = result.x.reshape(-1, 2)
        l1 = np.array([np.linalg.norm(vertices[i] - vertices[j]) for i, j in self.edges])
        return vertices, l1, result

    def optimization_callback(self, xk):
        if not hasattr(self, "_iteration"):
            self._iteration = 0  # Initialize counter on first call
        print(f"Iteration {self._iteration}: Current error = {self._compute_error(xk)}")
        self._iteration += 1  # Increment counter

    def standard_error_function(self, flat_coords):
        coords = flat_coords.reshape(-1, 2)
        error = 0
        for (i, j), L_ij in zip(self.edges, self.l0):
            dist = np.linalg.norm(coords[i] - coords[j])
            error += (dist - L_ij) ** 2
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