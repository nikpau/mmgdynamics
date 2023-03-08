import copy
import numpy as np

from warnings import warn
from typing import Optional, Any
from scipy.integrate import solve_ivp

from .structs import Vessel

from .dynamics import _shallow_water_hdm,mmg_dynamics

class LogicError(Exception):
    pass

__version__ = "1.0.7"
__author__ = "Niklas Paulig <niklas.paulig@tu-dresden.de>"
__all__ = ["step"]


def step(*, X: np.ndarray,vessel: Vessel, dT: float, nps: float, 
         delta: float, psi: float, fl_psi: Optional[float] = None,
         fl_vel: Optional[float] = None, w_vel: Optional[float] = None, 
         beta_w: Optional[float] = None, water_depth: Optional[float] = None, 
         debug: bool = False, atol: float = 1e-5, rtol: float = 1e-5, 
         **sol_options # Additional options (see help(solve_ivp))
        )->Any:
    """Solve the MMG system for a given vessel for an arbitrarily long timestep
    
    Args:
        X (np.ndarray): Initial values (u,v,r)
        params (dict): Vessel dict
        sps (float): Seconds to simulate per timestep
        nps (float): Propeller revolutions per second
        delta (float): Rudder angle in [rad]
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
        if water_depth < vessel.d:
            raise LogicError(
                "Water depth cannot be less than ship draft.\n"
                f"Water depth: {np.round(water_depth,3)} | Ship draft: {vessel.d}"
                )
        sh_vessel = copy.deepcopy(vessel)
        _shallow_water_hdm(sh_vessel, water_depth)

    solution = solve_ivp(
        fun=mmg_dynamics,
        t_span=(float(0), float(1)), # Calculate the system dynamics for this time span
        y0=X,
        t_eval=np.array([float(1)]), # Evaluate the result at the final time
        args=(vessel if water_depth is None else sh_vessel, # Order is important! Do not change
            psi,delta,nps,
            fl_psi,fl_vel,
            w_vel,beta_w,dT),
        method="RK23",
        rtol = rtol,
        atol = atol,
        **sol_options
    )

    # Return the complete solver output, with all steps
    if debug:
        return solution

    # Just return the relevant derivatives
    return solution.y