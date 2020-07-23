# LQG-RTD
Reachability-based Trajectory Design (RTD) for a Linear-Quadratic-Gaussian (LQG) agent 

## agents
- `linear_agent.m` - 2D agent with continuous linear Gaussian dynamics (LQR feedback)
- `LQG_agent.m` - 2D agent with continuous LQG dynamics and measurement model (Kalman filter and LQR feedback)
- `discrete_agent.m` - 2D agent with discrete LQG dynamics and measurement model (Kalman filter and LQR feedback)

### low_level_controllers
- `LQR_LLC.m` - 

## FRS
- `generate_linear_dynamics.m` - create matlab linear dynamics function to be used by CORA for reachability computation
- `linear_FRS.m` - "trajectory" FRS computation for linear dynamics (no error included)

## other
- `run_linear_agent.m` - run linear agent
- `run_LQG_agent.m` - run LQG agent (comment in/out for discrete or continuous version)
