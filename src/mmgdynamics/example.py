import matplotlib
import numpy as np
import mmgdynamics.calibrated_vessels as cvs

from mmgdynamics.maneuvers import *
from mmgdynamics.structs import Vessel

matplotlib.rcParams['font.family'] = ['sans-serif']
font = {'weight': 'normal','size': 18}
matplotlib.rc('font', **font)


"""In here you find some basic tests to validate the 
    MMG model.
    Uncomment the desired test and run the script.
    
    Tests:
        - 35° turning tests for either starboard or port direction
        - ZigZag test 
    
    Note:
        See the docs for each function for additional information
        or use help(somefunction)
    
"""

# Initial values for the full scale kvlcc2
ivs = InitialValues(
    u     = 3.85, # Longitudinal vessel speed [m/s]
    v     = 0.0, # Lateral vessel speed [m/s]
    r     = 0.0, # Yaw rate acceleration [rad/s]
    delta = 0.0, # Rudder angle [rad]
    nps   = 1.05 # Propeller revs [s⁻¹]
)
# # Initial values for the 7-meter model
# ivs = InitialValues(
#     u     = 1.128, # Longitudinal vessel speed [m/s]
#     v     = 0.0, # Lateral vessel speed [m/s]
#     r     = 0.0, # Yaw rate acceleration [rad/s]
#     delta = 0.0, # Rudder angle [rad]
#     nps   = 13.4 # Propeller revs [s⁻¹]
# )

# Initial values for 1/5 kvlcc2
ivs = InitialValues(
    u     = 4.0, # Longitudinal vessel speed [m/s]
    v     = 0.0, # Lateral vessel speed [m/s]
    r     = 0.0, # Yaw rate acceleration [rad/s]
    delta = 0.0, # Rudder angle [rad]
    nps   = 3.0 # Propeller revs [s⁻¹]
)

# Use a pre-calibrated vessel
vessel = Vessel(**cvs.kvlcc2)

iters = 3000

p = turning_maneuver(ivs, vessel, iters, "port", water_depth=None)
q = turning_maneuver(ivs, vessel, iters, "port", water_depth=8)
plot_trajecory([p,q], vessel)
plot_r([p,q])
z, l = zigzag_maneuver(ivs, vessel, dir=1, dps=4.76, max_deg=10,wd=[1.2*vessel.d,None])
plot_zigzag(z, l, vessel, ivs)