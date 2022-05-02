import json
import matplotlib
import numpy as np
import mmgdynamics.calibrated_vessels as cvs

from mmgdynamics.maneuvers import *
from mmgdynamics.structs import InlandVessel

matplotlib.rcParams['font.family'] = ["Helevtica", 'sans-serif']

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
europaschiff = {
    "m":        1545.6, # Displacement
    "d":        1.781, # Draft
    "A_R":      6.0, # Rudder Area
    "B":        9.5, # Width
    "Lpp":      105, # Length
    "w_P0":     0.4,
    "t_P":      0.2, # Thrust deduction factor
    "C_b":      0.7, # Block coefficient
    "D_p":      1.5, # Propeller diameter
    "f_alpha":  2.4,
    "eta":      0.938 # Ratio of propeller diameter to rudder span
}

# Some initial values
ivs = np.array([4.0, # Longitudinal vessel speed [m/s]
                0.0, # Lateral vessel speed [m/s]
                0.0, # Yaw rate acceleration [rad/s]
                0.0, # Rudder angle [rad]
                5.0] # Propeller revs [s⁻¹]
               )

# Uncomment to calibrate a vessel from the minimal dict above
#vessel = calibrate(MinimalVessel(**europaschiff),rho = 1000)

# Use a pre-calibrated vessel
vessel = InlandVessel(**cvs.SPTRR1)

ZIGZAG: bool = False

iters = 200
s = turning_maneuver(ivs, vessel, iters, "starboard",maxdeg=15)
#p = turning_maneuver(ivs, vessel, iters, "starboard",water_depth=10,maxdeg=5)
#q = turning_maneuver(ivs, vessel, iters, "starboard",water_depth=10)
plot_trajecory([s], vessel)
#plot_r(p)

#free_flow_test(vessel)
angles = np.arange(180)
angles = angles/180*np.pi
#current_test(vessel, 400, 45/180*math.pi)
#static_current_test(vessel,angles)


if ZIGZAG:
    z, l = zigzag_maneuver(ivs, vessel, rise_time=5, max_deg=10,wd=[5,None])
    plot_zigzag(z, l)

# Print the vessel dict to output
print(json.dumps(vessel.__dict__,sort_keys=True, indent=4))