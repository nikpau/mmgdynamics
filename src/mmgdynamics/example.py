import json
import numpy as np

from maneuvers import *
import calibrated_vessels as cvs

"""In here you find some basic tests to validate the 
    MMG model.
    
    Tests:
        - N° turning tests for either starboard or port direction
        - ZigZag test 
    
    Note:
        See the docs for each function for additional information
        or use help(somefunction)
    
"""

# Example minimal dict for a German Großmotorgüterschiff (inland transport vessel)
fully_loaded_GMS = {
    "m":        3614.89, # Displacement
    "d":        3.1891, # Draft
    "A_R":      5.29, # Rudder Area
    "B":        11.45, # Width
    "Lpp":      110, # Length
    "C_b":      0.9, # Block coefficient
    "t_P":      0.2, # Thrust deduction factor
    "D_p":      1.751, # Propeller diameter
    "eta":      0.960, # Ratio of propeller diameter to rudder span
    "w_P0":     0.4, # Wake fraction coefficient
    "f_alpha":  2.45 # Rudder lift gradient coefficient
}

# Some initial values
ivs = np.array([5.0, # Longitudinal vessel speed [m/s]
                0.0, # Lateral vessel speed [m/s]
                0.0, # Yaw rate acceleration [rad/s]
                0.0, # Rudder angle [rad]
                5.0] # Propeller revs [s⁻¹]
               )

# Uncomment to calibrate a vessel from the minimal dict above
#vs = calibrate(fully_loaded_GMS,rho = 1000)

# Use a pre-calibrated vessel
vessel = cvs.GMS1

ZIGZAG: bool = False

iters = 600
s = turning_maneuver(ivs, vessel, iters, "starboard")
p = turning_maneuver(ivs, vessel, iters, "starboard", water_depth=150.)
plot_trajecory([s, p], vessel)
plot_r(p)

if ZIGZAG:
    z, l = zigzag_maneuver(ivs, vessel, rise_time=30, max_deg=20)
    plot_zigzag(z, l)

# Print the vessel dict to output
print(json.dumps(vessel,sort_keys=True, indent=4))