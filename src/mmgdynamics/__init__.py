import copy
import numpy as np

from warnings import warn
from typing import Optional

from .structs import Vessel
from .dynamics import _shallow_water_hdm,mmg_dynamics

class LogicError(Exception):
    pass

__version__ = "1.0.7"
__author__ = "Niklas Paulig <niklas.paulig@tu-dresden.de>"
__all__ = ["step", "pstep"]

# Vessel fixed coordinate system type aliases
Surge = float # Longitudinal velocity
Sway = float # Lateral velocity
YawRate = float # Yaw rate

Nu = tuple[Surge, Sway, YawRate]

# Earth fixed coordinate system type aliases
Northing = float # North coordinate
Easting = float # East coordinate
Psi = float # Heading angle [rad]

EarthFixedPos = tuple[Northing, Easting]
Eta = tuple[Northing, Easting, Psi]

def angle_to_two_pi(angle: float) -> float:
    """Converts an angle to the range [0, 2*pi]."""
    return angle % (2 * np.pi)

def rotpsi(psi: float) -> np.ndarray:        
    """Computes a rotation matrix for given heading (in rad)."""
    return np.array(
        [
            [np.cos(psi), -np.sin(psi), 0.],
            [np.sin(psi),  np.cos(psi), 0.],
            [0.,           0.,          1.]
        ]
    )

def step(*, X: np.ndarray, vessel: Vessel, dT: float, nps: float, 
         delta: float, psi: float, fl_psi: Optional[float] = None,
         fl_vel: Optional[float] = None, w_vel: Optional[float] = None, 
         beta_w: Optional[float] = None, water_depth: Optional[float] = None
         ) -> np.ndarray[Nu]:
    """Solve the MMG system for a given vessel for an arbitrarily long timestep
    
    Args:
        X (np.ndarray): Initial values: np.array([u,v,r])
        params (dict): Vessel dict
        nps (float): Propeller revolutions per second (not per minute!)
        delta (float): Rudder angle in [rad]
        fl_psi (float): Attack angle of current relative to heading [rad]
        fl_vel (Optional[float]): Fluid velocity (Current velocity)
        w_vel (Optional[float]): Wind velocity
        beta_w (Optional[float]): Wind attack angle
        water_depth( Optional[float]): Water depth if vessel is simulated in shallow water

    Raises:
        LogicError: - If the water depth is less than the 
                      ships draft. Ship would run aground. It is recommended 
                      to have at least (1.2*draft) meters of water under the vessel.
                    - If a current velocity has been set but no 
                      attack angle for it. 

    Returns:
        Derivatives of u,v and r in the vessel fixed coordinate system
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
        sh_vessel = _shallow_water_hdm(copy.deepcopy(vessel), water_depth)

    # Develop equation of motion
    uvr_dot = mmg_dynamics(
        X = X,
        params = sh_vessel if water_depth is not None else vessel,
        nps = nps,
        delta = delta,
        psi = psi,
        fl_psi = fl_psi,
        fl_vel = fl_vel,
        w_vel = w_vel,
        beta_w = beta_w
    )

    # Just return the relevant derivatives
    return uvr_dot * dT

def pstep(*, X: np.ndarray, pos: EarthFixedPos, vessel: Vessel, dT: float, nps: float, 
         delta: float, psi: float, fl_psi: Optional[float] = None,
         fl_vel: Optional[float] = None, w_vel: Optional[float] = None, 
         beta_w: Optional[float] = None, water_depth: Optional[float] = None
         ) -> tuple[Nu, Eta]:
    """
    Same as step but, the transformation of the velocities to the 
    earth fixed coordinate system is done here automatically.
    
    Returns:
        np.ndarray: New surge, sway and yaw rate of the vessel (Nu)
        np.ndarray: New position and heading of the vessel (Eta)
    """
    old_uvr = X.copy()
    old_eta_dot = np.dot(rotpsi(psi), old_uvr)
    
    uvr_dot = step(
        X = X,
        vessel = vessel,
        dT = dT,
        nps = nps,
        delta = delta,
        psi = psi,
        fl_psi = fl_psi,
        fl_vel = fl_vel,
        w_vel = w_vel,
        beta_w = beta_w,
        water_depth = water_depth
    )
    
    # Get new velocities in the vessel fixed coordinate system
    new_uvr = old_uvr + uvr_dot
    
    # Find new eta_dot via rotation
    eta_dot_new = np.dot(rotpsi(psi), new_uvr)
    
    # Update position in earth fixed coordinate system
    eta = np.hstack((pos,psi))
    new_eta = eta + 0.5 * (old_eta_dot + eta_dot_new) * dT
    
    # Correct for overshooting heading angles
    new_eta[2] = angle_to_two_pi(new_eta[2])
    
    # (array[u,v,r], array[N,E,psi])
    return new_uvr, new_eta
    