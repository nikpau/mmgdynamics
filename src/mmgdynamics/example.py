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
ivs = np.array([7.97, # Longitudinal vessel speed [m/s]
                0.0, # Lateral vessel speed [m/s]
                0.0, # Yaw rate acceleration [rad/s]
                0.0, # Rudder angle [rad]
                2] # Propeller revs [s⁻¹]
               )
# Some initial values
ivs = np.array([1.128, # Longitudinal vessel speed [m/s]
                0.0, # Lateral vessel speed [m/s]
                0.0, # Yaw rate acceleration [rad/s]
                0.0, # Rudder angle [rad]
                13.4] # Propeller revs [s⁻¹]
               )

# Initial values for 1/5 kvlcc2
ivs = np.array([4.0, # Longitudinal vessel speed [m/s]
                0.0, # Lateral vessel speed [m/s]
                0.0, # Yaw rate acceleration [rad/s]
                0.0, # Rudder angle [rad]
                3.0] # Propeller revs [s⁻¹]
               )

# Uncomment to calibrate a vessel from the minimal dict above
#vessel = calibrate(FUb,rho = 1000)

# Use a pre-calibrated vessel
vessel = Vessel(**cvs.kvlcc2)

ZIGZAG: bool = True

iters = 800
#s = turning_maneuver(ivs, vessel, iters, "starboard",maxdeg=35)
#p = turning_maneuver(ivs, vessel, iters, "starboard",water_depth=20)
#q = turning_maneuver(ivs, vessel, iters, "starboard",water_depth=10)
#plot_trajecory([s], vessel)
#plot_r(s)

#free_flow_test(vessel, ivs)
angles = np.arange(180)
angles = angles/180*np.pi
#current_test(vessel, ivs, 400, 0/180*math.pi)
#static_current_test(vessel,angles)


if ZIGZAG:
    z, l = zigzag_maneuver(ivs, vessel, dir=1, dps=4.76, max_deg=10,wd=[1.2*vessel.d,None])
    plot_zigzag(z, l, vessel, ivs)

# Print the vessel dict to output
print(json.dumps(vessel.__dict__,sort_keys=True, indent=4))