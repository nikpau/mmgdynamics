import numpy as np
import json
from dynamics import turning_maneuver, plot_r, plot_trajecory, plot_zigzag, zigzag_maneuver
import calibrated_vessels as cvs

"""In here you find some basic tests to validate the 
    MMG model.
    
    Tests:
        - N° turning tests for either starboard or port direction
        - ZigZag test 
    
    Note:
        See the docs for each function for addidtional information
    
"""

# Example river dict
some_river = {
    "wd_avg": 14.,
    "B_K":    200.,
}

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

# Vessel parameter dictionary calculated from minimal dict
vs = cvs.get_coef_dict(fully_loaded_GMS,1000)

# Uncomment next line to use a pre-calibrated vessel
#vs = cvs.GMS1

# Print the vessel dict to output
print(json.dumps(vs,sort_keys=True, indent=4))


iters = 1000
s = turning_maneuver(ivs, vs, iters, "starboard")
p = turning_maneuver(ivs, vs, iters, "port")
z, l = zigzag_maneuver(ivs, vs, rise_time=10, max_deg=20)

plot_trajecory([s, p], vs)
plot_r(s)
plot_zigzag(z, l)