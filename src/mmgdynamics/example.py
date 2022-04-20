import json
import matplotlib
import numpy as np
import calibrated_vessels as cvs

from maneuvers import *
from mmgdynamics.structs import Vessel
from mmgdynamics.dynamics import calibrate

matplotlib.rcParams['font.family'] = ["Noto Sans CJK TC", 'sans-serif']

font = {'weight': 'bold',
        'size': 14}

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

# Example minimal dict for a European fishing vessel
FUb = {
    "m":        642.4, # Displacement
    "d":        4.5, # Draft
    "A_R":      4.36, # Rudder Area
    "B":        10.0, # Width
    "Lpp":      26.16, # Length
    "C_b":      0.57, # Block coefficient
    "D_p":      2.5, # Propeller diameter
    "eta":      0.838 # Ratio of propeller diameter to rudder span
}

# Some initial values
ivs = np.array([7.0, # Longitudinal vessel speed [m/s]
                0.0, # Lateral vessel speed [m/s]
                0.0, # Yaw rate acceleration [rad/s]
                0.0, # Rudder angle [rad]
                6.0] # Propeller revs [s⁻¹]
               )

# Uncomment to calibrate a vessel from the minimal dict above
#vessel = calibrate(FUb,rho = 1000)

# Use a pre-calibrated vessel
vessel = Vessel(**cvs.kvlcc2)

ZIGZAG: bool = True

iters = 1000
#s = turning_maneuver(ivs, vessel, iters, "starboard")
#p = turning_maneuver(ivs, vessel, iters, "starboard",water_depth=20)
#q = turning_maneuver(ivs, vessel, iters, "starboard",water_depth=10)
#plot_trajecory([s, p, q], vessel)
#plot_r(p)

free_flow_test(vessel)

if ZIGZAG:
    z, l = zigzag_maneuver(ivs, vessel, rise_time=10, max_deg=10,wd=5)
    plot_zigzag(z, l)

# Print the vessel dict to output
print(json.dumps(vessel.__dict__,sort_keys=True, indent=4))