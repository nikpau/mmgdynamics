import math
from itertools import cycle
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgba
from matplotlib.patches import Patch, Rectangle

from mmgdynamics import step
from mmgdynamics.structs import InitialValues, Vessel


def turning_maneuver(ivs: InitialValues, vessel: Vessel,
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
         np.linspace(0.0, maxdeg, 15),
         np.full(time-15, maxdeg)]
    )

    if dir == "starboard":
        delta_list = delta_list * np.pi/180.
    elif dir == "port":
        delta_list = delta_list * np.pi/180. * -1
    else:
        raise Exception("Please choose either 'port' or 'starboard' as dir")

    # Starting angle
    psi = 0.
    
    # Initial values
    uvr = np.array([ivs.u,ivs.v,ivs.r])

    # Simulate for `time` seconds
    for s in range(time):

        # Solve the ODE system for one second at a time
        sol = step(X=uvr,
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
        u, v, r = sol

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

        U = math.sqrt(u**2+v**2)
        drift = -math.atan(-v/u)*180/math.pi
        res[2, s] = drift

        # Set current solution as next initial values
        uvr = np.hstack(sol)

    return res

def plot_r(t: Sequence[Sequence[float]]) -> None:
    plt.figure(figsize=(16, 10))
    for list in t:
        plt.plot(np.arange(len(list[2])), list[2], linewidth=2.5)
    plt.xlabel("$t(s)$", fontsize=14)
    plt.ylabel("$r(-)$", fontsize=14)
    plt.title(
        "Yaw rate acceleration $r(-)$ for $\pm 35^{\circ}$ turning maneuver")
    plt.grid(True)
    # plt.savefig(PLOTDIR+"/"+"r.pdf")
    plt.show()


def plot_trajecory(t: Sequence[np.ndarray], vessel: Vessel) -> None:
    """Plot trajecories 

    Args:
        t (list[np.ndarray]): List of ndarrays with list[n][0] being the x-coord
        and list[n][1] the y coord of the trajectory
        vessel (dict): parameter dict for vessel
    """

    plt.figure(figsize=(16, 10))

    # for tr in t:
    #     plt.plot(tr[1]/vessel.Lpp, tr[0]/vessel.Lpp, linewidth=2.5)
    plt.plot(t[0][1]/vessel.Lpp, t[0][0]/vessel.Lpp, linewidth=2.5, color="orange")
    plt.plot(t[1][1]/vessel.Lpp, t[1][0]/vessel.Lpp, linewidth=2.5)

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
    # plt.savefig(PLOTDIR+"/"+"turning.pdf")
    plt.show()


def zigzag_maneuver(ivs: InitialValues, vessel: Vessel, 
                    max_deg: int, dps: float, dir=1, 
                    wd: Optional[Sequence[float]] = None) -> np.ndarray:
    """Perform ZigZag maneuver

    Args:
        ivs (np.ndarray): Set of inital values
        vessel (dict): vessel parameter dict
        max_deg (int): Maximum rudder angle in degrees (Passed to the build_delta_zigzag() function)
        rise_time (int): Maximum rudder angle in degrees (Passed to the build_delta_zigzag() function)

    Returns:
        np.ndarray: Array of yaw angles for each timestep
    """
    assert dir in [1,-1], "Invalid dir"

    result = [[], []]
    delta_list = [[], []]

    reset_ivs = np.array([ivs.u,ivs.v,ivs.r])

    # Degrees of rudder change per second
    dps = dps*math.pi/180
    max_deg = max_deg*math.pi/180

    for idx, depth in enumerate(wd):

        uvr = reset_ivs
        psi = 0.0
        delta = 0.0
        first,second,third = True, True, True
        while True:
            
            # Solve the ODE system for one second at a time
            sol = step(X=uvr,
                    vessel=vessel,
                    dT=1,
                    psi=psi,
                    nps=ivs.nps,
                    delta=delta,
                    fl_vel=None,
                    fl_psi=None,
                    water_depth=depth,
                    atol=1e-6,
                    rtol=1e-6)

            # Angular turning rate (rad/s)
            _, _, r = sol

            psi += r

            result[idx].append(float(psi))
            delta_list[idx].append(float(delta))

            uvr = np.hstack(sol)

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
        

    return plot_zigzag(result, delta_list,vessel, ivs)


def plot_zigzag(
    t: list, delta_list: np.ndarray, vessel: Vessel, ivs: InitialValues) -> None:
    
    ls = [(0, (3, 1, 1, 1, 1, 1)),"-"]
    labels = ["h/d = $\infty$","h/d = 1.2"]
    L = vessel.Lpp
    delta_max = max(delta_list[0])

    fig = plt.figure(figsize=(8, 4))

    for i,v in enumerate(t):
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
    plt.savefig(f"zigzag{int(delta_max)}{int(delta_max)}.pdf")


def free_flow_test(vessel:Vessel, ivs: InitialValues):

    fig = plt.figure(figsize=(6,6))
    ax: plt.Axes = fig.add_subplot(1, 1, 1)
    plt.xlabel("$y_0/L$", fontsize=14)
    plt.ylabel("$x_0/L$", fontsize=14)

    def rad(a): return a/180*math.pi

    timestep = 0
    twopi = 2*math.pi

    # Length and Breadth of the simulated vessel
    L, B = vessel.Lpp/vessel.Lpp, vessel.B/vessel.Lpp

    ivs_init = np.array([ivs.u,ivs.v,ivs.r])
    
    # Define here lists of Rudder angles in degrees 
    # that are performed by the vessel.
    delta_list_deep = np.concatenate(
        [np.array([0.]),
         np.linspace(0.0, 35, 7),
         np.full(600-7, 35)]
    )    
    delta_list_shallow = np.concatenate(
        [np.array([0.]),
         np.linspace(0.0, 35, 7),
         np.full(1200-7, 35)]
    )
    
    deltas = [delta_list_shallow,delta_list_deep]
    
    # A color cycler to color different experiments
    colors = cycle([
        "#001219", "#005f73", "#0a9396", 
        "#94d2bd", "#e9d8a6", "#ee9b00", 
        "#ca6702", "#bb3e03", "#ae2012", "#9b2226"
        ])
    d = vessel.d
    
    depths = [d*1.2,None]

    # There must be one list per experiment
    assert len(depths)==len(deltas)
    for color,depth,dl in zip(colors,depths,deltas):
        
        agx, agy = 0.0, 0.0
        head = 0.0
        timestep = 0
        uvr = ivs_init
        plt.grid(ls=":")
        
        for d in range(len(dl)):

            sol = step(
                X=uvr,
                vessel=vessel,
                psi=head,
                nps=ivs.nps,
                delta=rad(dl[d]),
                fl_psi=None,
                fl_vel=None,
                water_depth=depth,
                dT=1,
                atol=1e-6,rtol=1e-3
            )

            timestep += 1

            # Vel in x and y direction (m/s), angular turning rate (rad/s)
            u, v, r = sol

            # Transform to earth-fixed coordinate system

            head =  (head + r) % twopi
            vy = math.cos(head) * u - math.sin(head) * v
            vx = math.sin(head) * u + math.cos(head) * v

            agx += vx/vessel.Lpp
            agy += vy/vessel.Lpp
 
            anchor = agx - B/2, agy - L/2

            # Set current solution as next initial values
            uvr = np.hstack(sol)

            # Rectangle of the heading transformed vessel
            vessel_rect = Rectangle(
                anchor,
                width=vessel.B/vessel.Lpp,
                height=vessel.Lpp/vessel.Lpp,
                rotation_point="center",
                angle=((2*math.pi)-head)*180/math.pi,
                edgecolor=to_rgba(color,1),
                facecolor = to_rgba(color,0.1)
            )

            # ADJUST HERE IF YOU WANT TO PLOT ONLY
            # EVERY NTH DATAPOINT 
            if timestep == 1 or timestep % 12 == 0:
                ax.add_patch(vessel_rect)
                
            print(timestep,depth,color)
    
    plt.xlim(-0.5,5.5)
    plt.ylim(-1,5)
    
    # Unfortunately we must generate our Legend
    # manually. 
    # So, for every experiment, add your legend
    # entry as a Patch() and add it to the legend
    # handler. PRs for improvements are welcome.
    handles, _ = ax.get_legend_handles_labels()
    shallow = Patch(color=colors[0], 
                 label="h/d = 1.2")
    deep = Patch(color=colors[1], 
                 label="h/d = $\infty$")
    
    
    handles.append(deep)
    handles.append(shallow)
    plt.legend(handles=handles)
    
    plt.savefig("freeflow.pdf")

def current_wind_test(vessel: Vessel, ivs: InitialValues, 
                      iters: int,fl_psi: float, fl_vel: float,
                      w_vel: float, beta_w:float) -> None:

    fig = plt.figure()
    ax: plt.Axes = fig.add_subplot(1, 1, 1)
    plt.xlabel(r"$x_0$", fontsize=14)
    plt.ylabel(r"$y_0$", fontsize=14)
    plt.grid()
    
    timestep = 0
    twopi = 2*math.pi

    agx, agy = 0.0, 0.0
    head = 0/180*math.pi

    qx, qy = np.linspace(-400, 400, 11), np.linspace(-400, 400, 11)

    # Length and Breadth of the simulated vessel
    L, B = vessel.Lpp, vessel.B
    
    nps = ivs.nps 
    uvr = np.array([ivs.u,ivs.v,ivs.r])
    
    colors = ["#9b2226","#001219"]
        
    for _ in range(iters):

        sol = step(
            X=uvr,
            vessel=vessel,
            psi=head,
            nps=nps,
            delta=0.0,
            fl_psi=fl_psi,
            fl_vel=fl_vel,
            w_vel=w_vel,
            beta_w=beta_w,
            water_depth=None,
            dT=1,
            atol=1e-5,rtol=1e-5
        )

        timestep += 1

        # Vel in x and y direction (m/s), angular turning rate (rad/s)
        u, v, r = sol

        # Transform to earth-fixed coordinate system
        head = float((head + r) % twopi)
        vy = math.cos(head) * u - math.sin(head) * v
        vx = math.sin(head) * u + math.cos(head) * v

        agx += vx
        agy += vy

        anchor = agx - B/2, agy - L/2
        
        # Set current solution as next initial values
        uvr = np.hstack(sol)
        
        speed = math.sqrt(u**2+v**2)
        print("Speed: {:.2f}".format(speed),f" | Timestep: {timestep}")

        # Rectangle of the heading transformed vessel
        vessel_rect = Rectangle(anchor,
                                width=vessel.B,
                                height=vessel.Lpp,
                                rotation_point="center",
                                angle=(twopi-head)*180/math.pi,
                                edgecolor=colors[0],
                                facecolor = to_rgba(colors[0],0.1))


        if timestep % 50 == 0:
            ax.add_patch(vessel_rect)

    ax.quiver(qx, qy,
                np.full((11, 11), -math.sin(beta_w)),
                np.full((11, 11), -math.cos(beta_w)),
                scale=20, headwidth=2)

        #print(timestep)
    plt.axis("equal")
    plt.show()

def _depr_turning_test(vessel:Vessel, ivs: InitialValues):
    """
    Currently not maintained.
    """

    fig = plt.figure(figsize=(6,6))
    #fig.patch.set_facecolor("#212529")
    ax: plt.Axes = fig.add_subplot(1, 1, 1)
    
    plt.xlabel(r"$y_0/L$", fontsize=14)
    #ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/vessel.Lpp)) 
    #ax.xaxis.set_major_formatter(ticks)
    
    plt.ylabel(r"$x_0/L$", fontsize=14)
    #ax.yaxis.set_major_formatter(ticks)

    
    def rad(a): return a/180*math.pi

    timestep = 0
    twopi = 2*math.pi

    qx, qy = np.linspace(-500, 500, 11), np.linspace(-500, 500, 11)

    # Length and Breadth of the simulated vessel
    L, B = vessel.Lpp/vessel.Lpp, vessel.B/vessel.Lpp

    ivs_init = np.array([ivs.u,ivs.v,ivs.r])
    
    delta_list_deep = np.concatenate(
        [np.array([0.]),
         np.linspace(0.0, 35, 7),
         np.full(600-7, 35)]
    )    
    delta_list_shallow = np.concatenate(
        [np.array([0.]),
         np.linspace(0.0, 35, 7),
         np.full(1200-7, 35)]
    )
    
    deltas = [delta_list_shallow,delta_list_deep]
    
    colors = ["#c1121f","#003049"]
    d = vessel.d
    depths = [d*1.2,None]

    for color,depth,dl in zip(colors,depths,deltas):
        
        agx, agy = 0.0, 0.0
        head = 0.0
        timestep = 0
        uvr = ivs_init
        plt.grid(ls=":")
        
        for d in range(len(dl)):

            sol = step(
                X=uvr,
                vessel=vessel,
                psi=head,
                nps=ivs.nps,
                delta=rad(dl[d]),
                fl_psi=None,
                fl_vel=None,
                water_depth=depth,
                dT=1,
                atol=1e-6,rtol=1e-3
            )

            timestep += 1

            # Vel in x and y direction (m/s), angular turning rate (rad/s)
            u, v, r = sol

            # Transform to earth-fixed coordinate system

            head =  (head + r) % twopi
            vy = math.cos(head) * u - math.sin(head) * v
            vx = math.sin(head) * u + math.cos(head) * v

            agx += vx/vessel.Lpp
            agy += vy/vessel.Lpp
 
            anchor = agx - B/2, agy - L/2

            # Set current solution as next initial values
            uvr = np.hstack(sol)

            # Rectangle of the heading transformed vessel
            vessel_rect = Rectangle(anchor,
                                    width=vessel.B/vessel.Lpp,
                                    height=vessel.Lpp/vessel.Lpp,
                                    rotation_point="center",
                                    angle=((2*math.pi)-head)*180/math.pi,
                                    edgecolor=to_rgba(color,1),
                                    facecolor = to_rgba(color,0.1))

            #ax.clear()
            # ax.quiver(qx, qy,
            #           np.full((11, 11), math.sin(fl_psi)),
            #           np.full((11, 11), math.cos(fl_psi)),
            #           scale=20, headwidth=2)
            if timestep == 1 or timestep % 12 == 0:
                ax.add_patch(vessel_rect)
                
            print(timestep,depth,color)
    
    plt.xlim(-0.5,5.5)
    plt.ylim(-1,5)
    
    handles, _ = ax.get_legend_handles_labels()
    shallow = Patch(color=colors[0], 
                 label=r"h/d = 1.2")
    deep = Patch(color=colors[1], 
                 label=r"h/d = $\infty$")
    
    
    handles.append(deep)
    handles.append(shallow)
    plt.legend(handles=handles)
    
    plt.savefig("turning_circles.pdf")
