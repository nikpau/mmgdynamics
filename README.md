# MMG Standard model with extensions for winds, currents and shallow water

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
> Disclaimer: The quality of the empirical estimations for hydrodynamic derivatives vary greatly for different ships. Please consider comprehensive testing before using a custom vessel.

#### Calibration process:

Under `src/structs.py` you find the dataclasses responsible for modeling the vessel objects. For using a minmal dict as a vessel, you must define it as seen below and then pass it into the `calibrate()` function which returns a full vessel object

The empirical estimations need at least the following information:

```python
from mmgdynamics.structs import MinimalVessel
from mmgdynamics.dynamics import calibrate

my_vessel = {
  "m":        0.0, # Vessel displacement [m³]
  "B":        0.0, # Vessel Breadth (width)
  "Lpp":      0.0, # Length between perpendiculars
  "C_b":      0.0, # Block coefficient (< 0.8)
  "D_p":      0.0, # Propeller Diameter
  "eta":      0.0, # Ratio of propeller diameter to rudder span
  "f_alpha":  0.0  # Rudder lift gradient coefficient 
                   # (If not given you will be asked for the rudder aspect ratio)
}

# To create a complete vessel object, you must pass
# the minimal dict and the water density of your environment 
# into the calibrate function as a minimal Vessel:
full_vessel = calibrate(MinimalVessel(**my_vessel),rho = 1000)
```

### Extension for winds and currents

Current and wind forces are calculated according to [Fossen, 2011](https://doi.org/10.1002/9781119994138). 
The angle of attack for currents is set as an angle from the global reference frame. 0° current are parallel to the x-axis. Angles rotate clockwise, directions are modeled as `coming from`. (Wind direction of 90° means wind flows from east to west.)

### Shallow water adaption

The effects of shallow water are incorporated using various semi-empirical formulas summarized in [Taimuri et. al. (2020)](https://doi.org/10.1016/j.oceaneng.2020.108103)

## Examples

You can find common test cases for vessel maeuvering such as the ZigZag or turning maneuver test in the `example.py` file. In there you also find the definition and use of the minimal dict from above.


## Citation

If you use this code in one of your projects or papers, please cite it as follows.

```bibtex
@misc{mmgdynamics,
  author = {Niklas Paulig},
  title = {MMG standard model for ship maneuvering},
  year = {2022},
  publisher = {GitHub},
  journal = {GitHub Repository},
  howpublished = {\url{https://github.com/nikpau/mmgdynamics}}
}
```
