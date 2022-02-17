# MMG Standard model with extensions for currents


In here you find an implementation of the [MMG standard model](https://doi.org/10.1007/s00773-014-0293-y) by Yasukawa, H., Yoshimura, Y. (2015).

## Dynamics


#### Note: Shallow water adaption is currently not working.

The dynamics of the system can be found in the `dynamics.py` file. The model can be used straight out of the box with any vessel found in the `calibrated_vessels.py` file. If you want to embed this model into your own framework you just need the `mmg_step()` function.

### Calibrate custom vessel
In order to calibrate a vessel that is not present in the `calibrated_vessels.py` file, you can define a minimal dict with basic information about the vessel and input it into the `get_coef_dict()` function. Several empirical formulas will then be used to estimate the relevant hydrodynamic derivatives for your vessel and return back a dict which can be used as an input to the `mmg_step()` function. 

#### Form of the dict:
The empirical estimations need at least the following information:
```
my_vessel = {
  "m":        0.0, # Vessel displacement [m³]
  "B":        0.0, # Vessel Breadth (width)
  "Lpp":      0.0, # Length between perpendiculars
  "C_b":      0.0, # Block coefficient
  "t_P":      0.0, # Thrust deduction factor
  "D_p":      0.0, # Propeller Diameter
  "eta":      0.0, # Ratio of propeller diameter to rudder span
  "f_alpha" : 0.0 # Rudder lift gradient coefficient (If not given you will be asked for the rudder aspect ratio)
}
```

### Extension for currents

This model currently only implements currents in x-direction seen from the global coordinate system. Positive current velocity flow in postive x direction, negative ones in negative x direction. 

##### Note: If you do not use currents, please set `fl_vel = None` in the `mmg_step()` function as zero-valued currents can not be handled at the moment.


