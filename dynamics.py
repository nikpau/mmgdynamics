import math
from typing import Callable, Optional
import matplotlib
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.misc import derivative
import matplotlib
import matplotlib.pyplot as plt
import calibrated_vessels as cvs

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
                 water_depth: Optional[float] = None, # Water depth at current timestep [m]
                 mode: str = "river", # "river" or "freeflow"
                 r_params: Optional[dict] = None # River dict
                 ) -> np.ndarray:

    # Shorten the parameter dict to avoid clutter
    p = params
    if r_params is not None:
        rp = r_params

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

    # For river mode the longitudinal force on the ship hull is 
    # calculated according to the BAW FaRAO Model. TODO: Can this be disclosed?
    if mode == "river":
        assert rp is not None and water_depth is not None, "Please supply water depth and river property dict."

        # Hull resistance due to Pressure and surface friction due to fuselage flows
        cw = getcw(p["d"], water_depth)
        A_S = p["B"] * p["d"] * p["C_b"] # Main frame cross sectional area
        v_str = math.cos(psi) * fl_vel # Fluid velocity component in longitudinal direction 
        
        F_U = 0.5 * p["rho"] * cw * A_S * v_str**2
        
        # Hull resistance between ship hull and river ground (is this true?)
        O_S = (p["B"] + 2*p["d"]) * p["C_b"]
        lambda_s = 4*(1.89-1.62*(p["k_s"]/p["Lpp"]))**(-2.5)
        vr = getvr(rp,p,U,psi,fl_vel,water_depth)
        
        F_R = O_S * p["rho"] * (lambda_s/8) * (v_str - vr)**2
        
        X_H = F_R + F_U 
    
    elif mode == "freeflow":
        # Longitudinal surge force around midship acting on ship hull
        X_H = (0.5 * p["rho"] * p["Lpp"] * p["d"] * (U**2) * (
            - p["R_0_dash"]
            + p["X_vv_dash"] * (v_dash**2)
            + p["X_vr_dash"] * v_dash * r_dash
            + p["X_rr_dash"] * (r_dash**2)
            + p["X_vvvv_dash"] * (v_dash**4)
        )
        )
    else:
        raise RuntimeError("Mode string could not be matched. Avail. options: ['river','freeflow'].")

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
        # X_C = 0.5 *p["rho"] * p["Lpp"] * p["d"] * U_c**2 * C_1c(aoa,S,fl_vel,p)
        X_C = 0.

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
              mode: str = "river", # River or freeflow mode
              r_params: Optional[dict] = None, # Parameter dict for river 
              debug: bool = False, # Print all intermediate time steps of the solver
              **sol_options # Additional options (see help(solve_ivp))
              ):

    # Correct for shallow water if a water depth is given. 
    # If none is given, open water with infinite depth is assumed
    if water_depth is not None:
        params = shallow_water_hdm(params, water_depth)

    solution = solve_ivp(fun=dynamics,
                         t_span=(float(0), float(sps)), # Calculate the system dynamics for this time span
                         y0=X,
                         t_eval=np.array([float(sps)]), # Evaluate the result at the final time
                         args=(params, 
                               sps, 
                               psi,
                               fl_vel, 
                               nps_old, 
                               delta_old,
                               water_depth,
                               mode,
                               r_params),
                         method="RK45",
                         rtol = 1e-4,
                         atol=1e-4,
                         **sol_options
                        )

    # Return the complete solver output, with all steps
    if debug:
        return solution

    # Just return the relevant derivatives
    return solution.y

# Correct the hydrodynamic derivatives for shallow water conditions.
# Correction formulas come from https://doi.org/10.1016/j.oceaneng.2019.106679 (NOT WORKING)
def shallow_water_hdm(v: dict, water_depth: float) -> dict:
    
    # Length, Beam (width), Draft, Block Coefficient
    L, B, d, Cb = v["Lpp"], v["B"], v["d"], v["C_b"]
    
    # Correction term would not make sense otherwise
    if L > 80:
        L = 80.
    
    coef = (water_depth/d) - 1
    
    mxshallow = (coef**1.3+3.77+1.14*(B/d)-0.233*(L/d)-3.43*Cb)/(coef**1.3)
    myshallow = (coef**0.82+0.413+0.0320*(B/d)+0.129*(B/d)**2)/(coef**0.82)
    Jzshallow = (coef**0.82+0.413+0.0192*(B/d)+0.00554*(B/d)**2)/(coef**0.82)
    
    v["m_x_dash"] *= mxshallow
    v["m_y_dash"] *= myshallow
    v["J_z_dash"] *= Jzshallow
    
    return v

# Resistance coefficient
def getcw(d: float, h: float) -> float:

    cw0, cw1, cw2 = 0.1, 0.2, 2
    cw = cw0 + cw1*(d/h)**cw2
    
    return cw

def getvr(rp: dict, p: dict, U: float, psi: float, fl_vel: float, h: float) -> float:
    
    vsdw = U - math.cos(psi) * fl_vel
    
    ccr = 0.8
    A_S = p["B"] * p["d"]
    A_K = rp["B_K"] * rp["wd_avg"]
    
    vcr = 0.58 * ccr * math.sqrt(GRAVITY * h) * (A_K / A_S - 1)**(1/8)
    
    s = (vsdw/vcr)**3 * math.sqrt(p["d"]/h)

    A_V = A_S + s * p["B"] # TODO: Impl the effective width of the vessel here (Plane intersecting a cuboid)
    
    vr = vsdw * (A_V /(A_K - A_V)) # Check if A_K is really correct here

    return vr
    
    
# Get the hydrodynamic coefficients due to currents. Sources: 
# https://doi.org/10.1016/j.amc.2004.12.005
# https://doi.org/10.1016/S0141-1187(98)00002-9
# psi = Vessel yaw angle [rad]
# S = Wetted surface area [m²]. Calculated as L*B + 2*d(L+B)
# fl_vel = fluid velocity around vessel (current velocity)
def C_1c(psi: float, S:float, fl_vel: float, p: dict) -> float:
    
    L, d, rho = p["Lpp"], p["d"], p["rho"]
    
    # Reynolds Number calculated for characteristic length of vessel. 
    # Dynamic viscosity of water at 20°C ~ 1.00e-6
    Re = abs(fl_vel)*L*rho/1.00e-6
    
    val = (0.9375*S/(((math.log(Re)-2)**2)*d*L)) * math.cos(psi) + \
        1/8 * np.pi*d/L*(math.cos(3*psi) - math.cos(psi))
        
    return val 

def C_2c(psi: float, p: dict) -> float:
    
    L, B, d, Cb = p["Lpp"], p["B"], p["d"], p["C_b"]
    C_y = 0.6 # Fig 1 in https://doi.org/10.1016/S0141-1187(98)00002-9
    
    val1 = (C_y - (np.pi*d/(2*L))) * math.sin(psi)*abs(math.sin(psi)) + (np.pi*d/(2*L))*(math.sin(psi)**3)
    val2 = (np.pi*d/L)*(1+0.4*Cb*B/d)*math.sin(psi)*abs(math.cos(psi))
    
    return val1 + val2

def C_6c(psi: float, p: dict) -> float:
    
    L, B, d, x_G = p["Lpp"], p["B"], p["d"], p["x_G"]
    C_y = 0.6
    
    val1 = -x_G/L*(C_y - (np.pi*d/(2*L)))* math.sin(psi)*abs(math.sin(psi))
    val2 = np.pi*d/L*math.sin(psi) * math.cos(psi)
    val3 = ((1+abs(math.cos(psi)))/2)**2 * (np.pi*d/L) * (0.5-2.4*d/L)*math.sin(psi)*abs(math.cos(psi))
    
    return val1 - val2 - val3

# Calculate a trajectory for a given set of inital values and a time horizon (s)
# The trajectory returns absolute x and y coordinates assuming the vessel started
# at the origin
# To reproduce https://doi.org/10.1007/s00773-014-0293-y the rudder steering
# rate will be 1.75°/s (1.76 in the paper)
def turning_maneuver(ivs: np.ndarray, vessel: dict, time: int, dir: str = "starboard") -> np.ndarray:

    # res[0:] = x-coord
    # res[1:] = y-coord
    res = np.zeros((3, time))

    maxdeg = 20
    # 35/20=1.75
    delta_list = np.concatenate(
        [np.array([0.]), np.linspace(0.00, maxdeg, 20), np.full(time-20, maxdeg)])
    #delta_list = np.full((time+1,),35.)

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
                        water_depth=5.,
                        r_params=some_river,
                        mode="freeflow",
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


def build_delta_zigzag(rise_time: int = 3, delta_max: float = 20.) -> np.ndarray:

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

# ------------------------------------------------------------------------------
# TEST RANGE

# Example river dict
some_river = {
    "wd_avg": 14.,
    "B_K":    200.,
}

fully_loaded_GMS = {
    "m": 3848.,
    "d": 3.63,
    "A_R":5.29,
    "B": 11.45,
    "Lpp":110,
    "C_b": 0.9,
    "t_P": 0.2,
    "D_p": 1.6,
    "eta": 0.960,
    "w_P0": 0.4
}

# Some initial values
ivs = np.array([5.0, 0.0, 0.0, 0.0, 5.0])
vs = cvs.get_coef_dict(fully_loaded_GMS,1000)
#vs = cvs.GMS1
print(json.dumps(vs,sort_keys=True, indent=4))
iters = 1000
s = turning_maneuver(ivs, vs, iters, "starboard")
p = turning_maneuver(ivs, vs, iters, "port")
# z, l = zigzag_maneuver(ivs, vs, rise_time=10, max_deg=20)

plot_trajecory([s, p], vs)
plot_r(s)
#plot_zigzag(z, l)
