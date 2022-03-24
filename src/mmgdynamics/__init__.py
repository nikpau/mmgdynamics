import copy
import numpy as np

from warnings import warn
from typing import Optional, Any
from scipy.integrate import solve_ivp

from dynamics import (
    LogicError, shallow_water_hdm,
    mmg_dynamics
)


__version__ = "1.0.3"
__author__ = "Niklas Paulig <niklas.paulig@tu-dresden.de>"
__all__ = ["step"]


def step(*, X: np.ndarray,params: dict, sps: float, nps_old: float, 
         delta_old: float, fl_psi: Optional[float] = None, 
         fl_vel: Optional[float] = None, water_depth: Optional[float] = None, 
         debug: bool = False, atol: float = 1e-5, rtol: float = 1e-5,
         **sol_options # Additional options (see help(solve_ivp))
         )->Any:
    """Solve the MMG system for a given vessel for an arbitrarily long timestep
    
    Args:
        X (np.ndarray): Initial values
        params (dict): Vessel dict
        sps (float): Seconds to simulate per timestep
        nps_old (float): propeller rotations [s⁻¹] last timestep
        delta_old (float): rudder angle last timestep
        fl_psi (float): Attack angle of current relative to heading [rad]
        fl_vel (Optional[float]): Fluid velocity (Current velocity)
        water_depth( Optional[float]): Water depth if vessel is simulated in shallow water
        debug (bool): Print all intermediate time steps of the solver

    Raises:
        LogicError: - If the waterdepth is less than the 
                      ships draft. Ship would run aground. It is recommended 
                      to have at least (1.2*draft) meters of water under the vessel.
                    - If a current velocity has been set but no 
                      attack angle for it. 

    Returns:
        Any: Solver output of two modes:
            Normal: Return the result of the last solver timestep
            Debug : Return also all intermediate solver values.
    """

    if fl_vel is not None and fl_psi is None:
        raise LogicError("No current direction specified. "
                         "Use the `fl_psi` keyword to set it.")
    if fl_vel is None and fl_psi is not None:
        warn("Current attack angle is given but current velocity is "
             "turned off. Attack angle is ignored.")
    
    # Correct for shallow water if a water depth is given. 
    # If none is given, open water with infinite depth is assumed
    if water_depth is not None:
        if water_depth < params["d"]:
            raise LogicError(
                "Water depth cannot be less than ship draft.\n"
                "Water depth:{} | Ship draft: {}".format(water_depth, params["d"])
                )
        sh_params = copy.deepcopy(params) # params is a pointer to the underlying array. But we need a new one
        shallow_water_hdm(sh_params, water_depth)

    solution = solve_ivp(fun=mmg_dynamics,
                         t_span=(float(0), float(sps)), # Calculate the system dynamics for this time span
                         y0=X,
                         t_eval=np.array([float(sps)]), # Evaluate the result at the final time
                         args=(params if water_depth is None else sh_params, # Order is important! Do not change
                               fl_psi,fl_vel, 
                               nps_old, delta_old),
                         method="RK45",
                         rtol = rtol,
                         atol = atol,
                         **sol_options
                        )

    # Return the complete solver output, with all steps
    if debug:
        return solution

    # Just return the relevant derivatives
    return solution.y