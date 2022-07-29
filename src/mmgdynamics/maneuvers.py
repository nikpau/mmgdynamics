import csv
import math
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgba
from matplotlib.patches import Patch, Rectangle

import mmgdynamics
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
                   sps=1,
                   nps=ivs.nps,
                   delta=delta_list[s],
                   fl_vel=None,
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

def plot_r(t: list[Sequence[float]]):
    plt.figure(figsize=(16, 10))
    for list in t:
        plt.plot(np.arange(len(list[2])), list[2], linewidth=2.5)
    plt.xlabel(r"$t(s)$", fontsize=14)
    plt.ylabel(r"$r(-)$", fontsize=14)
    plt.title(
        r"Yaw rate acceleration $r(-)$ for $\pm 35^{\circ}$ turning maneuver")
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
    plt.xlabel(r"$y_0/L$", fontsize=14)
    plt.ylabel(r"$x_0/L$", fontsize=14)
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
                    sps=1,
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
        

    return result, delta_list


def plot_zigzag(t: list, delta_list: np.ndarray, vessel: Vessel, ivs: np.ndarray) -> None:
    
    colors = ["#bc6c25","#283618"]
    ls = ["-.","-"]
    L = vessel.Lpp

    fig = plt.figure(figsize=(10, 4))
    #fig.patch.set_facecolor("#212529")
    ax: plt.Axes = fig.add_subplot(1, 1, 1)

    for i,v in enumerate(t):
        plt.plot(
            np.arange(len(v))*ivs[0]/L, 
            v, 
            linewidth=3.5, 
            color=colors[i],
            linestyle=ls[1]
        )

    for i,v in enumerate(delta_list):
        plt.plot(
            np.arange(len(v))*ivs[0]/L, 
            v, 
            linewidth=3.5, 
            color=colors[i],
            linestyle=ls[0]
        )

    plt.xlabel(r"$t*U_0/L$", fontsize=24)
    plt.ylabel(r"$Angle [^{\circ}]$", fontsize=24)
    plt.title(r"$10/10Z$",loc="left")

    plt.ylim(-40,40)
    plt.xlim(-1,28)

    plt.grid(True,"major")
    plt.grid(True,"minor",linestyle = ":")
    plt.minorticks_on()

    handles, _ = ax.get_legend_handles_labels()
    shallow = Patch(color=colors[0], 
                 label=r"h/T = 1.2")
    deep = Patch(color=colors[1], 
                 label=r"h/T = $\infty$")
    
    handles.append(shallow)
    handles.append(deep)
    plt.legend(handles=handles, loc="upper right")

    plt.tight_layout()
    plt.savefig("zigzag1010.pdf")


def free_flow_test(vessel:Vessel, ivs: InitialValues):

    fig = plt.figure(figsize=(6,6))
    #fig.patch.set_facecolor("#212529")
    ax: plt.Axes = fig.add_subplot(1, 1, 1)
    
    plt.xlabel(r"$x_0/L$", fontsize=14)
    #ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/vessel.Lpp)) 
    #ax.xaxis.set_major_formatter(ticks)
    
    plt.ylabel(r"$y_0/L$", fontsize=14)
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
    
    colors = ["#bc6c25","#283618"]
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
                sps=1,
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
            if timestep == 0 or timestep % 8 == 0:
                ax.add_patch(vessel_rect)
                
            print(timestep,depth,color)
    
    plt.xlim(-0.5,7)
    plt.ylim(-1,6.5)
    
    handles, _ = ax.get_legend_handles_labels()
    shallow = Patch(color=colors[0], 
                 label=r"h/T = 1.2")
    deep = Patch(color=colors[1], 
                 label=r"h/T = $\infty$")
    
    handles.append(shallow)
    handles.append(deep)
    plt.legend(handles=handles)
    
    plt.savefig("turning_circles.pdf")

def current_wind_test(vessel: Vessel, ivs: InitialValues, 
                      iters: int,fl_psi: float, fl_vel: float,
                      w_vel: float, beta_w:float) -> None:

    fig = plt.figure()
    #fig.patch.set_facecolor("#212529")
    ax: plt.Axes = fig.add_subplot(1, 1, 1)
    
    plt.xlabel(r"$x_0$", fontsize=14)
    #ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/vessel.Lpp)) 
    #ax.xaxis.set_major_formatter(ticks)
    
    plt.ylabel(r"$y_0$", fontsize=14)
    #ax.yaxis.set_major_formatter(ticks)

    plt.grid()
    
    def rad(a): return a/180*math.pi

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
            sps=1,
            atol=1e-6,rtol=1e-3
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

def static_current_test(vessel: Vessel, ang_list: Sequence[float]) -> None:

    fig = plt.figure()
    ax: plt.Axes = fig.add_subplot(1, 1, 1)
    
    plt.xlabel(r"Angle of attack $\gamma_c$ [deg] for current (relative bow)", fontsize=14)

    plt.grid(linestyle="-.")

    colors = ["#9b2226","#43aa8b","#001219","#005f73","#ee9b00","#ca6702"]
    cx,cy,cn = [], [], []
        
    for angle in ang_list:
        cx.append(mmgdynamics.dynamics._C_X(angle))
        cy.append(mmgdynamics.dynamics._C_Y(angle))
        cn.append(mmgdynamics.dynamics._C_N(angle))

    ax.plot(np.arange(len(ang_list)),cx,c=colors[2],lw=3,label=r"$C_X$")
    ax.plot(np.arange(len(ang_list)),cy,c=colors[1],lw=3,label=r"$C_Y$")
    ax.plot(np.arange(len(ang_list)),cn,c=colors[0],lw=3,label=r"$C_N$")

    fossen_plot = []

    with open("wpd_datasets.csv") as csvfile:
        reader = csv.reader(csvfile,delimiter=",")
        for row in reader:
            fossen_plot.append(row)

    cxf = [float(x[1]) for x in fossen_plot]
    cyf = [float(x[3]) for x in fossen_plot]
    cnf = [float(x[5]) for x in fossen_plot]

    ax.plot(np.arange(0,181,10),cxf,c=colors[5],lw=3,marker="o",label=r"$C_X$ (experimental)")
    ax.plot(np.arange(0,181,10),cyf,c=colors[4],lw=3,marker="v",label=r"$C_Y$ (experimental)")
    ax.plot(np.arange(0,181,10),cnf,c=colors[3],lw=3,marker="s",label=r"$C_N$ (experimental)")

    plt.legend()

    plt.show()
