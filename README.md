# MMG Standard model with extensions for currents and shallow water

In here you find an implementation of the [MMG standard model](https://doi.org/10.1007/s00773-014-0293-y) by Yasukawa, H., Yoshimura, Y. (2015).

## Installation

Install the package via pip:

```bash
pip install git+https://github.com/nikpau/mmgdynamics
```

## Dynamics

The dynamics of the system can be found in the `dynamics.py` file. The model can be used straight out of the box with any vessel found in the `calibrated_vessels.py` file. If you want to embed this model into your own framework you just need the `step()` function.

### Calibrate custom vessel

In order to calibrate a vessel that is not present in the `calibrated_vessels.py` file, you can define a minimal dict with basic information about the vessel and use `calibrate()` to make it usable in the `step()` function. Several empirical formulas will be used to estimate the relevant hydrodynamic derivatives for your vessel and return back a dict which can be used as an input to the `step()` function.

#### Form of the dict:

The empirical estimations need at least the following information:

```python
my_vessel = {
  "m":        0.0, # Vessel displacement [m³]
  "B":        0.0, # Vessel Breadth (width)
  "Lpp":      0.0, # Length between perpendiculars
  "C_b":      0.0, # Block coefficient (< 0.8)
  "t_P":      0.0, # Thrust deduction factor
  "D_p":      0.0, # Propeller Diameter
  "eta":      0.0, # Ratio of propeller diameter to rudder span
  "f_alpha":  0.0  # Rudder lift gradient coefficient 
                   # (If not given you will be asked for the rudder aspect ratio)
}
```

### Extension for currents

Current forces are calculated according to [Budak and Beji, 2020](https://doi.org/10.1016/j.oceaneng.2020.108126)
The angle of attack for currents is set as an angle from the global reference frame. 0° current are parallel to the x-axis.

## Examples

You can find common test cases for vessel maeuvering such as the ZigZag or turning maneuver test in the `example.py` file. In there you also find the definition and use of the minimal dict from above.
