import numpy as np
from scipy.io import loadmat

# load LPM from .mat file
lpm = loadmat('quadrotor_linear_planning_model.mat')
lpm = lpm['LPM']

# extract variables, convert arrays to numpy arrays
t_peak = lpm['t_peak'][0,0][0][0]
t_total = lpm['t_total'][0,0][0][0]
t_sample = lpm['t_sample'][0,0][0][0]
time = np.array(lpm['time'][0,0])
position = np.array(lpm['position'][0,0])
velocity = np.array(lpm['velocity'][0,0])
acceleration = np.array(lpm['acceleration'][0,0])

# specify some arbitrary trajectory parameters
v_x_0 = 1.0; v_y_0 = 1.0
a_x_0 = -2.0; a_y_0 = -2.0
v_x_peak = 2.3; v_y_peak = 2.3

k = np.array([[v_x_0, a_x_0, v_x_peak],
              [v_y_0, a_y_0, v_y_peak]])

# compute nominal trajectory
traj = np.dot(k, position)
print(traj) 
