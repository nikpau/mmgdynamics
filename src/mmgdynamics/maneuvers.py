
import math
from matplotlib.patches import Rectangle
import numpy as np
import matplotlib.pyplot as plt
from calibrated_vessels import kvlcc2

from typing import Dict, List
from mmgdynamics import step


def turning_maneuver(ivs: np.ndarray, vessel: dict,
                     time: int, dir: str = "starboard",
                     maxdeg: int = 20, water_depth: float = None) -> np.ndarray:
    """
    Perform a turning maneuver for a given set 
    of inital values and a time horizon. 
    The trajectory returns absolute x and 
    y coordinates assuming the vessel started at the origin.

    To reproduce https://doi.org/10.1007/s00773-014-0293-y the rudder steering
    rate has to be 1.75°/s (1.76 in the paper)

    Args:
        ivs (np.ndarray): Set of initial values
        vessel (dict): vessel parameter dict
        time (int): Time in [s] to reach the maximum rudder angle
        dir (str, optional): Turning direction. Defaults to "starboard".
        maxdeg (int, optional): Maximum rudder angle in degrees. Defaults to 20.

    Raises:
        Exception: Only if an invalid 'dir' was applied

    Returns:
        np.ndarray: Array of absolute x and y coordinates for the trajectory

    Note:
        There is no option to input curent speeds at the moment
        apart from hard-coding it into the function. Feel free to add the feature
    """

    # res[0:] = x-coord
    # res[1:] = y-coord
    res = np.zeros((3, time))

    # If maxdeg is set to 35° you have approx
    # 35°/20s=1.75°/s
    delta_list = np.concatenate(
        [np.array([0.]),
         np.linspace(0.00, maxdeg, 20),
         np.full(time-20, maxdeg)]
    )

    if dir == "starboard":
        delta_list = delta_list * np.pi/180.
    elif dir == "port":
        delta_list = delta_list * np.pi/180. * -1
    else:
        raise Exception("Please choose either 'port' or 'starboard' as dir")

    # Starting angle
    psi = 0.

    # Simulate for `time` seconds
    for s in range(time):

        # Set r to the next value in our list
        ivs[3] = delta_list[s+1]

        # Solve the ODE system for one second at a time
        sol = step(X=ivs,
                   params=vessel,
                   sps=1,
                   nps_old=ivs[4],
                   delta_old=delta_list[s],
                   fl_vel=None,
                   water_depth=water_depth)

        # Vel in x and y direction (m/s), angular turning rate (rad/s)
        u, v, r, _, _ = sol

        # Transform to earth-fixed coordinate system

        psi += r
        print(s, round(float(r), 4), psi*180/math.pi)
        v_x = math.cos(psi) * u - math.sin(psi) * v
        v_y = math.sin(psi) * u + math.cos(psi) * v

        if s == 0:
            res[0, s] = v_x
            res[1, s] = v_y
        else:
            res[0, s] = res[0, s-1] + v_x
            res[1, s] = res[1, s-1] + v_y

        res[2, s] = r

        # Set current solution as next initial values
        ivs = sol.flatten()

    return res


def plot_trajecory(t: List[np.ndarray], vessel: dict) -> None:
    """Plot trajecories 

    Args:
        t (list[np.ndarray]): List of ndarrays with list[n][0] being the x-coord
        and list[n][1] the y coord of the trajectory
        vessel (dict): parameter dict for vessel
    """

    plt.figure(figsize=(16, 10))

    for tr in t:
        plt.plot(tr[1]/vessel["Lpp"], tr[0]/vessel["Lpp"], linewidth=2.5)

    # This is just for the arrows depicting the flow direction
    # For now this is just added manually
    x, y = np.meshgrid(np.linspace(-5, 5, 20), np.linspace(-2, 4, 20))
    u, v = 0, 1
    plt.quiver(x, y, u, v, scale=100., width=0.001, color="grey")
    plt.xlabel(r"$y_0/L$", fontsize=14)
    plt.ylabel(r"$x_0/L$", fontsize=14)
    plt.title(
        r"Seiun Maru ($L=105m,B=18m$)$\pm 35^{\circ}$ turning maneuver | $U_0=5.0m/s (\approx 9.72kn), rpm=90$")
    plt.axhline(y=0, color='r', linestyle='-')
    plt.axvline(x=0, color='r', linestyle='-')
    plt.grid(True)
    # plt.savefig(PLOTDIR+"/"+"turning.pdf")
    plt.show()


def build_delta_zigzag(rise_time: int = 20, delta_max: float = 20.) -> np.ndarray:
    """Build a zigzag trajectory for the zigzag test. 

    Args:
        rise_time (int, optional): Time [s] for the rudder to reach 'delta_max'. Defaults to 20.
        delta_max (float, optional): Maximum rudder angle in degrees. Defaults to 20.

    Returns:
        np.ndarray: Array of rudder angles per timestep, length of array
    """

    # Init increase from 0 to delta max °
    init_inc = np.linspace(0, delta_max, rise_time)
    # Decrease from delta max° to -delta max°
    decr = np.linspace(delta_max, -delta_max, rise_time)
    # Increase from -delta max° to delta max°
    incr = np.linspace(-delta_max, delta_max, rise_time)

    # Build rudder angle list
    delta_zigzag = np.concatenate([
        init_inc,
        np.full(rise_time*3, delta_max),
        decr,
        np.full(rise_time*10, -delta_max),
        incr,
        np.full(rise_time*10, delta_max)
    ]).flatten()
    delta_zigzag = delta_zigzag * np.pi/180

    return delta_zigzag, len(delta_zigzag)


def zigzag_maneuver(ivs: np.ndarray, vessel: dict, max_deg: int, rise_time: int) -> np.ndarray:
    """Perform ZigZag maneuver

    Args:
        ivs (np.ndarray): Set of inital values
        vessel (dict): vessel parameter dict
        max_deg (int): Maximum rudder angle in degrees (Passed to the build_delta_zigzag() function)
        rise_time (int): Maximum rudder angle in degrees (Passed to the build_delta_zigzag() function)

    Returns:
        np.ndarray: Array of yaw angles for each timestep
    """

    delta_list, len_list = build_delta_zigzag(
        rise_time=rise_time, delta_max=max_deg)
    delta_list = np.concatenate([np.array([0.]), delta_list])
    res = np.empty(len_list)
    psi = 0.

    for s in range(len_list):

        # Set r to the next value in our list
        ivs[3] = delta_list[s+1]

        # Solve the ODE system for one second at a time
        sol = step(X=ivs,
                   params=vessel,
                   sps=1,
                   nps_old=ivs[4],
                   delta_old=delta_list[s],
                   fl_vel=None,
                   water_depth=None)

        # Angular turning rate (rad/s)
        _, _, r, _, _ = sol

        psi += r

        res[s] = psi

        ivs = sol.flatten()

    return res, delta_list


def plot_zigzag(t: list, delta_list: np.ndarray) -> None:

    t = t * 180/np.pi
    delta_list = delta_list * 180/np.pi
    delta_list = np.delete(delta_list, 0)

    plt.figure(figsize=(16, 10))
    plt.plot(np.arange(len(t)), t, linewidth=2.5)
    plt.plot(np.arange(len(t)), delta_list, linewidth=2.5)
    plt.xlabel(r"$t(s)$", fontsize=14)
    plt.ylabel(r"$\psi (blue) \qquad \delta (orange)$", fontsize=14)
    plt.title(r"ZigZag Maneuver with $\pm 20^{\circ}$")
    plt.grid(True)
    # plt.savefig(PLOTDIR+"/"+"zigzag.pdf")
    plt.show()


def plot_r(t: List[float]):
    plt.figure(figsize=(16, 10))
    plt.plot(np.arange(len(t[2])), t[2], linewidth=2.5)
    plt.xlabel(r"$t(s)$", fontsize=14)
    plt.ylabel(r"$r(-)$", fontsize=14)
    plt.title(
        r"Yaw rate acceleration $r(-)$ for $\pm 35^{\circ}$ turning maneuver")
    plt.grid(True)
    # plt.savefig(PLOTDIR+"/"+"r.pdf")
    plt.show()


def free_flow_test(vessel: Dict[str, float]):

    fig = plt.figure()
    fig.patch.set_facecolor("#212529")
    ax: plt.Axes = fig.add_subplot(1, 1, 1)

    def deg(a): return a*180/math.pi

    timestep = 0
    twopi = 2*math.pi

    head = 0.0

    qx, qy = np.linspace(-500, 500, 11), np.linspace(-500, 500, 11)

    agx, agy = 0.0, 0.0

    # Length and Breadth of the simulated vessel
    L, B = vessel["Lpp"], vessel["B"]

    ivs = np.array(
        [
            5.0,  # Longitudinal vessel speed [m/s]
            0.0,  # Lateral vessel speed [m/s]
            0.0,  # Yaw rate acceleration [rad/s]
            deg(10),  # Rudder angle [rad]
            5.0  # Propeller revs [s⁻¹]
        ]
    )
    nps = ivs[-1]

    while True:
        fl_psi = math.pi

        sol = step(
            X=ivs,
            params=vessel,
            psi=head,
            nps_old=nps,
            delta_old=deg(10),
            fl_psi=fl_psi,
            fl_vel=10.0,
            water_depth=6.0,
            sps=1,
            atol=1e-6,rtol=1e-6
        )

        timestep += 1
        print(round(float(ang_diff(head, fl_psi)*180/math.pi), 3), timestep)

        # Vel in x and y direction (m/s), angular turning rate (rad/s)
        u, v, r, _, _ = sol

        # Transform to earth-fixed coordinate system

        head =  (head + r) % twopi
        vy = math.cos(head) * u - math.sin(head) * v
        vx = math.sin(head) * u + math.cos(head) * v

        agx += vx
        agy += vy

        anchor = agx - B/2, agy - L/2

        # Set current solution as next initial values
        ivs = np.hstack(sol)

        # Rectangle of the heading transformed vessel
        vessel_rect = Rectangle(anchor,
                                width=vessel["B"],
                                height=vessel["Lpp"],
                                rotation_point="center",
                                angle=((2*math.pi)-head)*180/math.pi,
                                color="black")

        ax.clear()
        ax.axhline(y=0, color='r', linestyle='-')
        ax.axvline(x=0, color='r', linestyle='-')
        ax.quiver(qx, qy,
                  np.full((11, 11), math.sin(fl_psi)),
                  np.full((11, 11), math.cos(fl_psi)),
                  scale=20, headwidth=2)
        ax.add_patch(vessel_rect)
        plt.pause(0.01)


def abs_ang_diff(a1: float, a2: float) -> float:
    """Relative angle difference for an angle range of [0,2*pi]

    Args:
        a1 (float): Angle in radians or degrees
        a2 (float): Angle in radians or degrees

    Returns:
        float: absolute diff in angles or degs
    """
    def deg(a): return a*180/math.pi

    a1, a2 = deg(a1), deg(a2)
    a = a1-a2
    a = (a+180) % 360-180

    return a/180*math.pi


def ang_diff(head, str_dir):

    twopi = 2*math.pi
    return (str_dir - head) % twopi


free_flow_test(kvlcc2)
