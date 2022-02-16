# MMG Standard model with extensions for currents

======

In here you find an implementation of the [MMG standard model](https://doi.org/10.1007/s00773-014-0293-y) by Yasukawa, H., Yoshimura, Y. (2015).

## Dynamics

======

#### Please use the stable version on the [stable branch](https://github.com/nikpau/mmgdynamics/tree/stable)

The dynamics of the system can be found in the `dynamics.py` file. The model can be used straight out of the box with any model found in the `calibrated_vessels.py` file. If you want to embed this model into your own framework you just need the `mmg_step()` function.
