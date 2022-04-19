import numpy as np
from LPM import LPM
import params

class plan:
    """Plan class

    Planned trajectory
    
    Attributes
    ----------
    n_plan : int
    t : np.array (1 x N)
    p : np.array (1 x N)
    v : np.array (1 x N)
    a : np.array (1 x N)

    """

    def __init__(self):
        pass


class rtd_planner:
    """RTD Planner class

    Trajectory planner which recomputes collision-free trajectories to follow 
    in a receding-horizon fashion.

    Attributes
    ----------

    Methods
    -------

    """

    def __init__(self, lpm_file):
        # Initialize LPM object
        self.lpm = LPM(lpm_file)
        self.N_T_PLAN = len(self.lpm.time)  # planned trajectory length
        self.DT = self.lpm.t_sample  # trajectory discretization time interval

        # Initial conditions [m],[m/s],[m/s^2]
        self.p_0 = np.zeros((params.N_DIM,1))
        self.v_0 = np.zeros((params.N_DIM,1))
        self.a_0 = np.zeros((params.N_DIM,1))

        # Current plan ([p,v,a] x N)
        self.plan = np.zeros((1+3*params.N_DIM, self.N_T_PLAN))

        # Initialize with stationary plan
        self.plan[0,:] = self.lpm.time
        self.plan[1:4,:] = np.tile(self.p_0, (1,self.N_T_PLAN))

        # Goal position [m]
        self.p_goal = np.zeros((2,1))

        # Obstacles
        self.obstacles = []


    def traj_opt(self):
        """Trajectory Optimization

        Attempt to find a collision-free plan (v_peak) which brings the agent 
        closest to its goal.
        
        """

        # Generate potential v_peak samples
        V_peak = rand_in_bounds(params.V_BOUNDS, params.N_PLAN_MAX)
        # Eliminate samples that exceed the max velocity and max delta from initial velocity
        V_peak = prune_vel_samples(V_peak, self.v_0, params.V_MAX_NORM, params.DELTA_V_PEAK_MAX)
        


    def replan(self):
        """Replan
        
        """




def check_obs_collision(plan, obs, r_collision):
    """Check a plan against a single obstacle for collision.

    Obstacles are cylinders represented as (center, radius), and 
    assumed to have infinite height for now.

    Parameters
    ----------
    plan : np.array
    obs : tuple

    Returns
    -------
    bool
        True if the plan is safe, False is there is a collision

    """
    c_obs,r_obs = obs
    traj = plan[1:4,:]
    d_vec = np.linalg.norm(traj[:2,:] - c_obs[:,None], 2, 0)
    if any(d_vec <= r_collision + r_obs):
        return False
    else:
        return True


def rand_in_bounds(bounds, n):
    """Generate random samples within specified bounds

    Parameters
    ----------
    bounds : list
        List of min and max values for each dimension.
    n : int
        Number of points to generate.

    Returns
    -------
    np.array 
        Random samples

    """
    x_pts = np.random.uniform(bounds[0], bounds[1], n)
    y_pts = np.random.uniform(bounds[2], bounds[3], n)
    # 2D 
    if len(bounds) == 4:
        return np.vstack((x_pts, y_pts))
    # 3D
    elif len(bounds) == 6:
        z_pts = np.random.uniform(bounds[4], bounds[5], n)
        return np.vstack((x_pts, y_pts, z_pts))
    else:
        raise ValueError('Please pass in bounds as either [xmin xmax ymin ymax] '
                            'or [xmin xmax ymin ymax zmin zmax] ')


def prune_vel_samples(V, v_0, max_norm, max_delta):
    """Prune Velocity Samples
    
    """
    V_mag = np.linalg.norm(V, ord=2, axis=0)
    delta_V = np.linalg.norm(V - v_0, ord=2, axis=0)
    keep_idx = np.logical_and(V_mag < max_norm, delta_V < max_delta)
    return V[:,keep_idx]