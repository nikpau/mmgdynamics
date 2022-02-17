import math
from typing import Callable, Optional
import copy
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.misc import derivative
import matplotlib
import matplotlib.pyplot as plt
import crossflow as cf

"""
Set up the System of ODEs for vessel maneuvering prediction after 
Yasukawa, H., Yoshimura, Y. (2015) Introduction of MMG standard method for ship maneuvering predictions.
J Mar Sci Technol 20, 37-52 https://doi.org/10.1007/s00773-014-0293-y


On each step, the longituidal, lateral and yaw-rate accelerations are calculated
from a set of initial conditions X for an arbitrary time horizon. 
"""

# This is the only dynamic function we need to export really
__all__ = ["mmg_step"]

# Plotting directory folder name
PLOTDIR="maneuver_plots"

# Gravitational constant
GRAVITY: float = 9.81

# System of ODEs after Yasukawa, H., Yoshimura, Y. (2015)
def mmg_dynamics(t: np.ndarray, # Timestep
                 X: np.ndarray, # Initial values
                 params: dict, # Vessel dict 
                 sps: float, # Seconds per timestep [-]
                 psi: float, # Current heading angle [rad] (global coordinate system)
                 fl_vel: float, # Velocity of current [m/s]
                 nps_old: float, # propeller rotation last timestep [s⁻¹]
                 delta_old: float, # rudder angle last timestep [rad]
                 ) -> np.ndarray:
    """System of ODEs after Yasukawa, H., Yoshimura, Y. (2015)
    for the MMG standard model

    Args:
        t (np.ndarray): time
        X (np.ndarray): Initial values
        params (dict):  Vessel dict 
        sps (float): Seconds per timestep [-]
        psi (float): Current heading angle [rad] (global coordinate system)
        fl_vel (float): Velocity of current [m/s]
        nps_old (float): propeller rotation last timestep [s⁻¹]
        delta_old (float): rudder angle last timestep [rad]
        water_depth (Optional[float]): Water depth at current timestep [m] .Defaults to None
        mode (str):  Mode to calculcate. Either "river" or "freeflow". Defaults to "freeflow".
        r_params (Optional[dict]): Optional River dict . Defaults to None.

    Raises:
        RuntimeError: If undefined mode was selected

    Returns:
        np.ndarray: Derivatives of the ODE for the current timestep (mostly for solver)
    """

    # Shorten the parameter dict to avoid clutter
    p = params

    # u: long vel, v:lat. vel, r:yaw rate (d_psi/dt),delta: rudder angle, nps: propeller revs per second
    u, v, r, delta, nps = X

    # v = v_m + x_G*r | lateral velocity midship + x coord center of gravity * yaw rate
    v_m = v - p["x_G"] * r
    U = math.sqrt(u**2 + v_m**2)  # Overall speed of the vessel

    if U == 0.0:  # No vessel movement. Velocity in all directions = 0
        beta = 0.0
        v_dash = 0.0
        r_dash = 0.0
    else:
        beta = np.arctan2(-v_m, u)  # Drift angle at midship position
        v_dash = v_m / U  # Non-dimensionalized lateral velocity | TODO shouldn't this be v_m?
        r_dash = r * p["Lpp"] / U  # Non-dimensionalized yaw rate

    # Redefine
    w_P = p["w_P0"] * np.exp(-4.0 * (beta - (p["x_P"]/p["Lpp"]) * r_dash)**2)

    if nps == 0.0:  # No propeller movement, no advance ratio
        J = 0.0
    else:
        J = (1 - w_P) * u / (nps * p["D_p"])  # Propeller advance ratio

    if all(k in p for k in ("k_0","k_1","k_2")):
        # Propeller thrust open water characteristic
        K_T = p["k_0"] + (p["k_1"] * J) + (p["k_2"] * J**2)
    else:
        # Inferred slope + intercept dependent on J (empirical)
        K_T = p["J_slo"] * J + p["J_int"]

    # Effective inflow angle to rudder in maneuvering motions
    beta_R = beta - p["l_R"] * r_dash

    # Flow straightening coefficient
    if "gamma_R" in p:
        gamma_R = p["gamma_R"]
    else:
        if beta_R < 0.0:
            gamma_R = p["gamma_R_minus"]
        else:
            gamma_R = p["gamma_R_plus"]

    # Lateral inflow velocity components to rudder
    v_R = U * gamma_R * beta_R

    # Longitudinal inflow velocity components to rudder
    if J == 0.0:
        u_R = math.sqrt(p["eta"] * (p["kappa"] * p["epsilon"] *
                      8.0 * p["k_0"] * nps ** 2 * p["D_p"]**4 / np.pi)**2)
    else:
        u_R = u * (1 - w_P) * p["epsilon"] * math.sqrt(
            p["eta"] * (1.0 + p["kappa"] * (
                math.sqrt(1.0 + 8.0 * K_T / (np.pi * J**2)) - 1))**2 + (1 - p["eta"])
        )
    # Rudder inflow velocity
    U_R = math.sqrt(u_R**2 + v_R**2)

    # Rudder inflow angle
    alpha_R = delta - np.arctan2(v_R, u_R)

    # Normal force on rudder
    if "A_R" in p:
        F_N = 0.5 * p["A_R"] * p["rho"] * \
            p["f_alpha"] * (U_R**2) * math.sin(alpha_R)
    else:
        F_N = 0.5 * p["A_R/Ld_em"] * (p["Lpp"]*p["d"]*p["rho"]) * \
            p["f_alpha"] * (U_R**2) * math.sin(alpha_R)


    # Longitudinal surge force around midship acting on ship hull
    X_H = (0.5 * p["rho"] * p["Lpp"] * p["d"] * (U**2) * (
        - p["R_0_dash"]
        + p["X_vv_dash"] * (v_dash**2)
        + p["X_vr_dash"] * v_dash * r_dash
        + p["X_rr_dash"] * (r_dash**2)
        + p["X_vvvv_dash"] * (v_dash**4)
    )
    )

    # Longitudinal surge force around midship by steering
    X_R = -(1 - p["t_R"]) * F_N * math.sin(delta)

    # Longitudinal Surge force due to propeller
    X_P = (1 - p["t_P"]) * p["rho"] * K_T * nps**2 * p["D_p"]**4

    # Longitudinal surge force around midship acting on ship hull
    Y_H = (0.5 * p["rho"] * p["Lpp"] * p["d"] * (U**2) * (
        p["Y_v_dash"] * v_dash
        + p["Y_r_dash"] * r_dash
        + p["Y_vvv_dash"] * (v_dash**3)
        + p["Y_vvr_dash"] * (v_dash**2) * r_dash
        + p["Y_vrr_dash"] * v_dash * (r_dash**2)
        + p["Y_rrr_dash"] * (r_dash**3)
    )
    )

    # Lateral surge force by steering
    Y_R = -(1 + p["a_H"]) * F_N * math.cos(delta)

    # Yaw moment around midship acting on ship hull
    N_H = (0.5 * p["rho"] * (p["Lpp"]**2) * p["d"] * (U**2) * (
        p["N_v_dash"] * v_dash
        + p["N_r_dash"] * r_dash
        + p["N_vvv_dash"] * (v_dash**3)
        + p["N_vvr_dash"] * (v_dash**2) * r_dash
        + p["N_vrr_dash"] * v_dash * (r_dash**2)
        + p["N_rrr_dash"] * (r_dash**3)
    )
    )

    # Redimensionalize x_H
    x_H = p["x_H_dash"] * p["Lpp"]

    # yaw moment around midship by steering
    N_R = -((-p["Lpp"]/2) + p["a_H"] * x_H) * F_N * math.cos(delta)

    # Forces related to currents:
    # In this simulation, since we simulate a river, 
    # the current only has a x component.
    # Flow speed in y0 direction is always 0
    if fl_vel is not None:
    
        # Longitudinal velocity of current dependent on ship heading
        u_c = math.cos(psi) * fl_vel

        # Lateral velocity of current dependent on ship heading
        v_c = math.sin(psi) * fl_vel

        U_c = math.sqrt((v_m - v_c)**2 + (u - u_c)**2)

        # Angle of attack for current
        aoa = np.arctan2(v_m - v_c, u - u_c)
        
        # Wetted surface area [m²]. Calculated as (L*B + 2*d(L+B))*
        S = (p["Lpp"] * p["B"] + 2*p["d"]*(p["Lpp"] + p["B"])) * p["C_b"]
        
        # Longitudinal current force
        X_C = 0.5 *p["rho"] * p["Lpp"] * p["d"] * U_c**2 * C_1c(aoa,S,fl_vel,p)

        # Lateral current force
        Y_C = 0.5 *p["rho"] * p["Lpp"] * p["d"] * U_c**2 * C_2c(aoa,p)
        
        # Current Moment
        N_C = 0.5 *p["rho"] * p["Lpp"]**2 * p["d"] * U_c**2 * C_6c(aoa,p)
    
    else: 
        X_C, Y_C, N_C = 0.0, 0.0, 0.0
    
    # Added masses and added moment of inertia
    m_x = p["m_x_dash"] * (0.5 * p["rho"] * (p["Lpp"]**2) * p["d"])
    m_y = p["m_y_dash"] * (0.5 * p["rho"] * (p["Lpp"]**2) * p["d"])
    J_z = p["J_z_dash"] * (0.5 * p["rho"] * (p["Lpp"]**4) * p["d"])
    m = p["m"]
    I_zG = p["I_zG"]

    # Longitudinal acceleration
    d_u = ((X_H + X_R + X_P + X_C) + (m + m_y) * v *
           r + p["x_G"] * m * (r**2)) / (m + m_x)

    # Lateral acceleration
    f = (I_zG + J_z + (p["x_G"]**2) * m)

    d_v = ((Y_H+Y_R+Y_C) - (m+m_x)*u*r - ((p["x_G"]*m*(N_H + N_R + N_C))/(f)) + ((p["x_G"]**2*m**2*u*r)/(f)))\
        / ((m+m_y)-((p["x_G"]**2*m**2)/(f)))

    # Yaw rate acceleration
    d_r = ((N_H + N_R + N_C) - (p["x_G"] * m * d_v + p["x_G"] * m * u * r)) / \
        (I_zG + J_z + (p["x_G"]**2) * m)

    time_list = np.arange(sps+1)
    # Change in rudder position
    delta_list = np.linspace(delta_old, delta, num=int(sps)+1)
    d_interpol = interp1d(time_list, delta_list, fill_value="extrapolate")
    d_delta = derivative(d_interpol, t)

    # Change in propeller revolutions
    nps_list = np.linspace(nps_old, nps, num=int(sps)+1)
    nps_interpol = interp1d(time_list, nps_list, fill_value="extrapolate")
    d_nps = derivative(nps_interpol, t)

    return np.array([d_u, d_v, d_r, d_delta, d_nps])

# Step function. Solve the mmg system of ODEs for one timestep
# given a vessel and initial values
# sps = seconds per step
def mmg_step(dynamics: Callable, # Dynamics to solve over time
              X: np.ndarray, # Initial values
              params: dict, # Vessel dict
              sps: float, # Seconds to simulate per timestep
              nps_old: float, # propeller rotations/s last timestep 
              delta_old: float, # rudder angle last timestep
              psi: float, # Yaw angle [rad]
              fl_vel: Optional[float] = None, # Fluid velocity (Current velocity)
              water_depth: Optional[float] = None, # Water depth if vessel is simulated in shallow water
              debug: bool = False, # Print all intermediate time steps of the solver
              **sol_options # Additional options (see help(solve_ivp))
              ):

    # Correct for shallow water if a water depth is given. 
    # If none is given, open water with infinite depth is assumed
    if water_depth is not None:
        sh_params = copy.deepcopy(params) # params is just a pointer to the underlying array. But we need a new one
        shallow_water_hdm(sh_params, water_depth)

    solution = solve_ivp(fun=dynamics,
                         t_span=(float(0), float(sps)), # Calculate the system dynamics for this time span
                         y0=X,
                         t_eval=np.array([float(sps)]), # Evaluate the result at the final time
                         args=(params if water_depth is None else sh_params, # Order is important! Do not change
                               sps, psi,fl_vel, 
                               nps_old, delta_old),
                         method="RK45",
                         rtol = 1e-5,
                         atol=1e-5,
                         **sol_options
                        )

    # Return the complete solver output, with all steps
    if debug:
        return solution

    # Just return the relevant derivatives
    return solution.y


def shallow_water_hdm(v: dict, water_depth: float) -> None:
    """Correct the hydrodynamic derivatives and
    hydrodynamic masses for shallow water conditions.
    
    Sources:
        Tang et. al (2020) https://doi.org/10.1016/j.oceaneng.2019.106679
        Furukawa et. al. (2016) https://dx.doi.org/10.18451/978-3-939230-38-0_33
        

    Args:
        v (dict): vessel parameter dict
        water_depth (float): water depth around vessel [m]
    """
    
    # Length, Beam (width), Draft, Block Coefficient
    L, B, d, Cb = v["Lpp"], v["B"], v["d"], v["C_b"]
    
    h = water_depth
    
    # Shallow water adaption for longitudinal hydrodynamic derivatives
    # Source: Furukawa et. al. (2016)
    v["X_vv_dash"]  *= -70.8*(d/h)**4 + 27.7*(d/h)**2 + 1
    v["X_vr_dash"]  *= 1.3*(d/h)**2 + 1
    
    # Hull frictional resistance coefficient
    v["R_0_dash"]   *= 0.388*(d/h)**2 + 1
    
    
    # Correction by Tang et al 2020
    k = (2*d)/L
    k_e = lambda q: k/(d/(2*h) + ((np.pi*d)/(2*h) * 1/math.tan((np.pi*d)/(2*h)))**q)
    
    v["Y_v_dash"] = -(0.5*np.pi*k_e(2.3)+1.4*Cb*B/L)
    v["Y_r_dash"] = 0.25*np.pi * k_e(0.7)
    v["N_v_dash"] = -k_e(1.7)
    v["N_r_dash"] = -(0.54*k_e(0.7) - k_e(0.7)**2)
    
    v["N_vvr_dash"] *= 1+6*(d/h)**2.5
    v["N_vrr_dash"] *= 1+6*(d/h)**2.5
    
    # Thrust deduction factor
    oneminustp = 1 - v["t_P"]
    oneminustp *= 1/(1-0.2*(d/h) + 0.7295*(d/h)**2)
    v["t_P"] = -oneminustp + 1
    
    # Wake fraction coefficient
    oneminuswp = 1 - v["w_P0"]
    oneminuswp *= math.cos(1.4*Cb*d/h)
    v["w_P0"] = -oneminuswp + 1
    
    # Flow-straightening coefficient
    flow_str_coef = 1 + 0.0161*(d/h) + 4.4222*(d/h)**2 - 4.9825*(d/h)**3
    if all(k in v for k in ["gamma_R_minus","gamma_R_plus"]):
        v["gamma_R_minus"] *= flow_str_coef
        v["gamma_R_plus"]  *= flow_str_coef
    else:
        v["gamma_R"] *= flow_str_coef
    
    
    # Correction term would not make sense otherwise
    if L > 60:
        L = 70.
    
    coef = (water_depth/d) - 1
    
    mxshallow = (coef**1.3+3.77+1.14*(B/d)-0.233*(L/d)-3.43*Cb)/(coef**1.3)
    myshallow = (coef**0.82+0.413+0.0320*(B/d)+0.129*(B/d)**2)/(coef**0.82)
    Jzshallow = (coef**0.82+0.413+0.0192*(B/d)+0.00554*(B/d)**2)/(coef**0.82)
    
    v["m_x_dash"] *= mxshallow 
    v["m_y_dash"] *= myshallow
    v["J_z_dash"] *= Jzshallow

    
# Get the hydrodynamic coefficients due to currents. Sources: 
def C_1c(psi: float, S:float, fl_vel: float, p: dict) -> float:
    """Longitudinal forces due to currents
    
    Sources:
        https://doi.org/10.1016/j.amc.2004.12.005
        https://doi.org/10.1016/S0141-1187(98)00002-9

    Args:
        psi (float): Angle of attack for current [rad]
        S (float): Wetted surface area
        fl_vel (float): fluid velocity (current velocity)
        p (dict): parameter dict for vessel

    Returns:
        float: Non-dimensionalized longitudinal current force coef
    """
    
    L, d, rho = p["Lpp"], p["d"], p["rho"]
    
    # Reynolds Number calculated for characteristic length of vessel. 
    # Dynamic viscosity of water at 20°C ~ 1.00e-6
    Re = abs(fl_vel)*L*rho/1.00e-6
    
    val = (0.9375*S/(((math.log(Re)-2)**2)*d*L)) * math.cos(psi) + \
        1/8 * np.pi*d/L*(math.cos(3*psi) - math.cos(psi))
        
    return val 

def C_2c(psi: float, p: dict) -> float:  
    """Lateral forces due to currents

    Args:
        psi (float): Angle of attack for current [rad]
        p (dict): parameter dict for vessel

    Returns:
        float: Non-dimensionalized lateral current force coef
    """
    
    L, B, d, Cb = p["Lpp"], p["B"], p["d"], p["C_b"]
    C_y = cf.interpolateY(B/(2*d),cf.cross_flow) # Fig 1 in https://doi.org/10.1016/S0141-1187(98)00002-9
    
    val1 = (C_y - (np.pi*d/(2*L))) * math.sin(psi)*abs(math.sin(psi)) + (np.pi*d/(2*L))*(math.sin(psi)**3)
    val2 = (np.pi*d/L)*(1+0.4*Cb*B/d)*math.sin(psi)*abs(math.cos(psi))
    
    return val1 + val2

def C_6c(psi: float, p: dict) -> float:
    """Additional moment force due to currents

    Args:
        psi (float): Angle of attack for current [rad]
        p (dict): parameter dict for vessel

    Returns:
        float: Non-dimensionalized moment current force coef
    """
    
    L, B, d, x_G = p["Lpp"], p["B"], p["d"], p["x_G"]
    C_y = cf.interpolateY(B/(2*d),cf.cross_flow)
    
    val1 = -x_G/L*(C_y - (np.pi*d/(2*L)))* math.sin(psi)*abs(math.sin(psi))
    val2 = np.pi*d/L*math.sin(psi) * math.cos(psi)
    val3 = ((1+abs(math.cos(psi)))/2)**2 * (np.pi*d/L) * (0.5-2.4*d/L)*math.sin(psi)*abs(math.cos(psi))
    
    return val1 - val2 - val3


def turning_maneuver(ivs: np.ndarray, vessel: dict, 
                     time: int, dir: str = "starboard", 
                     maxdeg: int=20, water_depth: float = None) -> np.ndarray:
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

        print(s)

        # Set r to the next value in our list
        ivs[3] = delta_list[s+1]

        # Solve the ODE system for one second at a time
        sol = mmg_step(mmg_dynamics, 
                        X=ivs, 
                        params=vessel, 
                        sps=1, 
                        nps_old=ivs[4], 
                        delta_old=delta_list[s],
                        fl_vel=None,
                        water_depth=water_depth,
                        r_params=None,
                        psi=psi)

        # Vel in x and y direction (m/s), angular turning rate (rad/s)
        u, v, r, _, _ = sol

        # Transform to earth-fixed coordinate system
        psi += r
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


def plot_trajecory(t: list[np.ndarray], vessel: dict) -> None:
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
    x,y=np.meshgrid(np.linspace(-5,5,20),np.linspace(-2,4,20))
    u, v = 0,1
    plt.quiver(x,y,u,v, scale=100., width = 0.001, color="grey")
    plt.xlabel(r"$y_0/L$", fontsize=14)
    plt.ylabel(r"$x_0/L$", fontsize=14)
    plt.title(r"Seiun Maru ($L=105m,B=18m$)$\pm 35^{\circ}$ turning maneuver | $U_0=5.0m/s (\approx 9.72kn), rpm=90$")
    plt.grid(True)
    #plt.savefig(PLOTDIR+"/"+"turning.pdf")
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
        sol = mmg_step(mmg_dynamics, 
                        X=ivs, 
                        params=vessel, 
                        sps=1, 
                        nps_old=ivs[4], 
                        delta_old=delta_list[s],
                        fl_vel=None,
                        water_depth=None,                        
                        r_params=None,
                        psi=psi)

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
    #plt.savefig(PLOTDIR+"/"+"zigzag.pdf")
    plt.show()


def plot_r(t: list[float]):
    plt.figure(figsize=(16, 10))
    plt.plot(np.arange(len(t[2])), t[2], linewidth=2.5)
    plt.xlabel(r"$t(s)$", fontsize=14)
    plt.ylabel(r"$r(-)$", fontsize=14)
    plt.title(r"Yaw rate acceleration $r(-)$ for $\pm 35^{\circ}$ turning maneuver")
    plt.grid(True)
    #plt.savefig(PLOTDIR+"/"+"r.pdf")
    plt.show()