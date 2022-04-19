# Linear Planning Model Class

import numpy as np
from scipy.io import loadmat

class LPM:
    """Linear Planning Model class.

    This class represents the Linear Planning Model, which converts trajectory parameters 
    to trajectories. It contains attributes for the trajectory parameterization.

    Attributes
    ----------
    t_peak : float
        Time in trajectory corresponding to peak velocity
    t_total : float
        Total time duration of trajectory
    t_sample : float
        Sample time of trajectory (discretization)
    time : np.array
        Time vector for trajectory
    p_mat : np.array
    v_mat : np.array
    a_mat : np.array
    
    Methods
    -------
    compute_trajectory(k) 
        Compute nominal trajectory from a given trajectory parameter
    solve_trajectory(v_0, a_0, p_goal)

    """
    # Class variables
    t_peak = 0
    t_total = 0
    t_sample = 0
    time = []
    p_mat = []
    v_mat = []
    a_mat = []


    def __init__(self, mat_file):
        """Construct LPM object from .mat file.

        Parameters
        ----------
            mat_file : .mat
                .mat file containing all the LPM parameters

        """
        # Load the .mat file
        lpm = loadmat(mat_file)
        lpm = lpm['LPM']

        # Extract variables, convert arrays to numpy arrays
        self.t_peak = lpm['t_peak'][0,0][0][0]
        self.t_total = lpm['t_total'][0,0][0][0]
        self.t_sample = lpm['t_sample'][0,0][0][0]
        self.time = np.array(lpm['time'][0,0])[0]
        self.p_mat = np.array(lpm['position'][0,0])
        self.v_mat = np.array(lpm['velocity'][0,0])
        self.a_mat = np.array(lpm['acceleration'][0,0])


    def compute_trajectory(self, k):
        """Compute nominal trajectory from a given trajectory parameter.

        Parameters
        ----------
        k : np.array
            trajectory parameter k = (v_0, a_0, v_peak), n x 3 where n is workspace dimension
        
        Returns
        -------
        Tuple 
            Tuple of the form (p,v,a) where each of p,v,a are n x N where n is the workspace dimension.
        
        """
        p = np.dot(k, self.p_mat)
        v = np.dot(k, self.v_mat)
        a = np.dot(k, self.a_mat)
        return p,v,a


    def solve_trajectory(self, v_0, a_0, p_goal):
        """Solve for the peak velocity which reaches a desired goal position.

        Note that this solution does not account for max velocity constraints

        Parameters
        ----------
        v_0 : np.array
            Initial velocity (1 x n row vector)
        a_0 : np.array
            Initial acceleration (1 x n row vector)
        p_goal : np.array
            Goal position (1 x n row vector)

        Returns
        -------
        np.array 
            Peak velocity (1 x n row vector).
            
        """
        # change to column vectors
        v_0 = np.reshape(v_0, (3,1))
        a_0 = np.reshape(a_0, (3,1))
        # position component due to v_0 and a_0
        p_from_ic = np.dot(np.hstack((v_0, a_0)), self.p_mat[0:2,-1])
        # solve for v_peak
        v_peak = (p_goal - p_from_ic) / self.p_mat[2,-1]
        return v_peak