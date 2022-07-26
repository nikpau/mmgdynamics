import json
import matplotlib
import numpy as np
import mmgdynamics.calibrated_vessels as cvs

from mmgdynamics.maneuvers import *
from mmgdynamics.structs import Vessel
from mmgdynamics.dynamics import calibrate

matplotlib.rcParams['font.family'] = ["Tahoma", 'sans-serif']

font = {'weight': 'normal',
        'size': 18}

matplotlib.rc('font', **font)


"""In here you find some basic tests to validate the 
    MMG model.
    
    Tests:
        - N° turning tests for either starboard or port direction
        - ZigZag test 
    
    Note:
        See the docs for each function for additional information
        or use help(somefunction)
    
"""

# Some initial values
ivs = InitialValues(
    u     = 7.97, # Longitudinal vessel speed [m/s]
    v     = 0.0, # Lateral vessel speed [m/s]
    r     = 0.0, # Yaw rate acceleration [rad/s]
    delta = 0.0, # Rudder angle [rad]
    nps   = 2 # Propeller revs [s⁻¹]
)
# # Some initial values
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

# Uncomment to calibrate a vessel from the minimal dict above
#vessel = calibrate(FUb,rho = 1000)

# Use a pre-calibrated vessel
vessel = Vessel(**cvs.kvlcc2_full)
vessel = Vessel(**cvs.kvlcc2)

ZIGZAG: bool = False

iters = 1500
#s = turning_maneuver(ivs, vessel, iters, "starboard",maxdeg=35)
#p = turning_maneuver(ivs, vessel, iters, "starboard",water_depth=20)
#q = turning_maneuver(ivs, vessel, iters, "starboard",water_depth=10)
#plot_trajecory([s], vessel)
#plot_r(s)

#free_flow_test(vessel, ivs)
angles = np.arange(180)
angles = angles/180*np.pi
current_test(vessel, ivs, 1000, 45/180*math.pi)
#static_current_test(vessel,angles)


if ZIGZAG:
    z, l = zigzag_maneuver(ivs, vessel, dir=1, dps=4.76, max_deg=10,wd=[1.2*vessel.d,None])
    plot_zigzag(z, l, vessel, ivs)

# Print the vessel dict to output
print(json.dumps(vessel.__dict__,sort_keys=True, indent=4))