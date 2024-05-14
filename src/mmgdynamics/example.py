import matplotlib
import mmgdynamics.calibrated_vessels as cvs
from dataclasses import dataclass

from mmgdynamics.maneuvers import *
from mmgdynamics.structs import Vessel

matplotlib.rcParams['font.family'] = ['Liberation Sans', 'sans-serif']
font = {'weight': 'normal','size': 14}
matplotlib.rc('font', **font)

"""
In here you find some basic tests to validate the 
    MMG model.
    Uncomment the desired test and run the script.
    
    Tests:
        - 35° turning tests for either starboard or port direction
        - ZigZag test 
    
    Note:
        See the docs for each function for additional information
        or use help(somefunction)
    
"""
@dataclass
class KVLCC2Inits:
    full_scale = InitialValues(
        u     = 3.85, # Longitudinal vessel speed [m/s]
        v     = 0.0, # Lateral vessel speed [m/s]
        r     = 0.0, # Yaw rate acceleration [rad/s]
        delta = 0.0, # Rudder angle [rad]
        nps   = 1.05 # Propeller revs [s⁻¹]
    )
    
    l_64 = InitialValues(
        u     = 4.0, # Longitudinal vessel speed [m/s]
        v     = 0.0, # Lateral vessel speed [m/s]
        r     = 0.0, # Yaw rate acceleration [rad/s]
        delta = 0.0, # Rudder angle [rad]
        nps   = 3.0 # Propeller revs [s⁻¹]
    )
    
    l_7 = InitialValues(
        u     = 1.128, # Longitudinal vessel speed [m/s]
        v     = 0.0, # Lateral vessel speed [m/s]
        r     = 0.0, # Yaw rate acceleration [rad/s]
        delta = 0.0, # Rudder angle [rad]
        nps   = 13.4 # Propeller revs [s⁻¹]
    )


# Use a pre-calibrated vessel
vessel = Vessel(**cvs.kvlcc2_l64)

iters = 3000
def test_turning_maneuver():
    s = turning_maneuver(KVLCC2Inits.l_64, vessel, iters, "port", water_depth=None)
    p = turning_maneuver(KVLCC2Inits.l_64, vessel, iters, "port", water_depth=25)
    plot_trajecory([s,p], vessel)
    plot_r([s,p])


def test_zig_zag():
    zigzag_maneuver(
        KVLCC2Inits.l_64, 
        vessel, 
        dir=1, 
        dps=5, 
        max_deg=20,
        wd=[1.2*vessel.d,None]
    )
# Print the vessel dict to output
print(json.dumps(vessel.__dict__,sort_keys=True, indent=4))