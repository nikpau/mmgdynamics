import math
from itertools import cycle
from typing import Optional, Sequence
from numpy.typing import NDArray

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

from mmgdynamics import step, pstep
from mmgdynamics.structs import InitialValues, Vessel

Trajectories = Sequence[NDArray[np.float64]]

def turning_maneuver(ivs: InitialValues, vessel: Vessel,
                     time: int, dir: str = "starboard",
                     water_depth: float = None) -> np.ndarray:
    """
    Perform a turning maneuver to reproduce
    https://doi.org/10.1007/s00773-014-0293-y .
    The rudder steering
    rate has to be 1.75°/s (1.76 in the paper)

    Args:
        ivs (np.ndarray): Set of initial values
        vessel (dict): vessel parameter dict
        time (int): Time in [s] to reach the maximum rudder angle
        dir (str, optional): Turning direction. Defaults to "starboard".

    Raises:
        Exception: Only if an invalid 'dir' was applied

    Returns:
        np.ndarray: Array of absolute x and y coordinates for the trajectory

    Note:
        There is no option to input curent speeds at the moment
        apart from hard-coding it into the function. Feel free to add the feature
    """

    # res[0,:] = x-coord
    # res[1,:] = y-coord
    res = np.zeros((3, time))

    # If maxdeg is set to 35° you have approx
    # 35°/20s=1.75°/s
    delta_list = np.concatenate(
        [np.array([0.]),
         np.linspace(0.0, 35, 15),
         np.full(time-15, 35)]
    )

    if dir == "starboard":
        delta_list = delta_list * np.pi/180.
    elif dir == "port":
        delta_list = delta_list * np.pi/180. * -1
    else:
        raise Exception("Please choose either 'port' or 'starboard' as dir")

    # Starting angle
    psi = 0.
    
    # Starting Pos:
    pos = np.array([0., 0.])
    
    # Initial values
    uvr = np.array([ivs.u,ivs.v,ivs.r])

    # Simulate for `time` seconds
    for s in range(time):

        # Solve
        uvr, eta = pstep(
            X=uvr,
            pos = pos,
            psi= psi,   
            vessel=vessel,
            dT=1,
            nps=ivs.nps,
            delta=delta_list[s],
            fl_vel=None,
            w_vel=0,
            beta_w=0,
            water_depth=water_depth
        )

        # Vel in x and y direction (m/s), angular turning rate (rad/s)
        u, v, r = uvr

        # Unpack position and heading
        x, y, psi = eta
        
        # Repack position
        pos = np.array([x,y])

        # Update the results
        res[0, s] = x
        res[1, s] = y

        # Calculate drift angle
        drift = -math.atan(-v/u)*180/math.pi
        res[2, s] = drift

    return res


def plot_r(trajectories: Trajectories) -> None:
    """
    Plots a list of yaw rates for a given trajectory.
    The trajectory is assumed to be a list of ndarrays
    with the following structure:
    list[n][0] = x-coord
    list[n][1] = y-coord
    list[n][2] = yaw rate
    """
    plt.figure(figsize=(16, 10))
    for t in trajectories:
        plt.plot(np.arange(len(t[2])), t[2], linewidth=2.5)
    plt.xlabel(r"$t(s)$", fontsize=14)
    plt.ylabel(r"$r(-)$", fontsize=14)

    plt.title(
        "Yaw rate acceleration $r(-)$ for $\pm 35^{\circ}$ turning maneuver")
    plt.grid(True)
    plt.show()


def plot_trajecory(trajectories: Trajectories, vessel: Vessel) -> None:
    """Plot trajecories 

    Args:
        trajectories (Trajectories) : List of trajectories
        vessel (dict): parameter dict for vessel
    """

    plt.figure(figsize=(16, 10))

    plt.plot(
        trajectories[0][1]/vessel.Lpp, 
        trajectories[0][0]/vessel.Lpp, 
        linewidth=2.5, color="orange"
    )
    
    plt.plot(
        trajectories[1][1]/vessel.Lpp, 
        trajectories[1][0]/vessel.Lpp, 
        linewidth=2.5
    )

    # This is just for the arrows depicting the flow direction
    # For now this is just added manually
    x, y = np.meshgrid(np.linspace(-5, 5, 20), np.linspace(-2, 4, 20))
    u, v = 0, 1
    plt.quiver(x, y, u, v, scale=100., width=0.001, color="grey")
    plt.xlabel("$y_0/L$", fontsize=14)
    plt.ylabel("$x_0/L$", fontsize=14)
    plt.axhline(y=0, color='black', linestyle=':')
    plt.axvline(x=0, color='black', linestyle=':')
    plt.axis("equal")
    plt.grid(True)
    plt.show()


def zigzag_maneuver(ivs: InitialValues, vessel: Vessel, 
                    max_deg: int, dps: float, dir=1, 
                    wd: Optional[Sequence[float]] = None) -> np.ndarray:
    """Perform ZigZag maneuver

    Args:
        ivs (np.ndarray): Set of inital values
        vessel (dict): vessel parameter dict
        max_deg (int): Maximum rudder angle in degrees (Passed to the build_delta_zigzag() function)
        wd (Optional[Sequence[float]], optional): Water depths. Defaults to None.

    Returns:
        np.ndarray: Array of yaw angles and rudder angles for each timestep
    """
    assert dir in [1,-1], "Invalid direction. Choose either 1 or -1."
    if wd is None:
        wd = [1.2*vessel.d, None]
        print("No water depth specified. Using 1.2*draft and infinite depth")
    
    def a2pm(angle: float) -> float:
        """[0,2pi] -> [-pi,pi]"""
        return angle if angle < np.pi else angle - 2*np.pi

    result = [[], []]
    delta_list = [[], []]

    reset_ivs = np.array([ivs.u,ivs.v,ivs.r])

    # Degrees of rudder change per second
    dps = dps*math.pi/180
    max_deg = max_deg*math.pi/180

    for idx, depth in enumerate(wd):

        uvr = reset_ivs
        pos = np.array([0., 0.])
        psi = 0.0
        delta = 0.0
        first,second,third = True, True, True
        while True:
            
            # Solve the ODE system for one second at a time
            uvr, eta = pstep(
                X=uvr,
                pos=pos,
                vessel=vessel,
                dT=1,
                psi=psi,
                nps=ivs.nps,
                delta=delta,
                fl_vel=None,
                fl_psi=None,
                water_depth=depth
            )

            x,y,psi = eta
            
            pos = np.array([x,y])
            
            psi = a2pm(psi)

            result[idx].append(float(psi))
            delta_list[idx].append(float(delta))

            if psi < max_deg*dir and first:
                delta = min(delta+dps*dir,max_deg*dir)
                continue
            else:
                first = False

            if psi > -max_deg*dir and second:
                delta = max(delta-dps*dir,-max_deg*dir)
                continue
            else:
                second = False

            if psi < -max_deg*dir and third:
                delta = min(delta+dps*dir,max_deg*dir)
                continue
            else:
                third = False
            
            if psi > 0.0:
                break
        
    for i, res in enumerate(result):
        result[i] = np.array(res)
        result[i] = result[i] * 180/np.pi

    for i, delta in enumerate(delta_list):
        delta_list[i] = np.array(delta)
        delta_list[i] = delta_list[i] * 180/np.pi
        
    return result, delta_list

def plot_zigzag(
    trajectories: Trajectories, 
    delta_list: np.ndarray, 
    vessel: Vessel, 
    ivs: InitialValues) -> None:
    """
    Plots 
    """
    
    ls = [(0, (3, 1, 1, 1, 1, 1)),"-"]
    labels = ["h/d = $\infty$","h/d = 1.2"]
    L = vessel.Lpp
    delta_max = max(delta_list[0])

    fig = plt.figure(figsize=(8, 4))

    for i,v in enumerate(trajectories):
        plt.plot(
            np.arange(len(v))*ivs.u/L, 
            v, 
            linewidth=2.5, 
            color="k",
            linestyle=ls[i],
            label = labels[i]
        )

    for i,v in enumerate(delta_list):
        plt.plot(
            np.arange(len(v))*ivs.u/L, 
            v, 
            linewidth=2.5, 
            color="k",
            linestyle=ls[i]
        )

    plt.xlabel("$t*U_0/L$", fontsize=18)
    plt.ylabel("$Angle [^{\circ}]$", fontsize=18)
    plt.title(f"${int(delta_max)}/-{int(delta_max)}Z$",loc="left")

    plt.ylim(-40,40)
    plt.xlim(0,12)

    plt.grid(True,"major")
    plt.grid(True,"minor",linestyle = ":")
    plt.minorticks_on()

    plt.legend(loc="upper right")

    plt.tight_layout()
    plt.show()