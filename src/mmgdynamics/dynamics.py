import math
import numpy as np

from . import crossflow as cf

from typing import Dict

"""
Set up the System of ODEs for vessel maneuvering prediction after 
Yasukawa, H., Yoshimura, Y. (2015) Introduction of MMG standard method for ship maneuvering predictions.
J Mar Sci Technol 20, 37-52 https://doi.org/10.1007/s00773-014-0293-y


On each step, the longituidal, lateral and yaw-rate accelerations are calculated
from a set of initial conditions X for an arbitrary time horizon. 
"""

# Export only selected functions if used as wildcard import
__all__ = [
    "turning_maneuver", "zigzag_maneuver",
    "plot_r", "plot_trajecory",
    "plot_zigzag", "calilbrate"
]

# System of ODEs after Yasukawa, H., Yoshimura, Y. (2015)


def mmg_dynamics(t: np.ndarray, X: np.ndarray, params: dict, psi:float,
                 fl_psi: float, fl_vel: float, nps_old: float,
                 delta_old: float,) -> np.ndarray:
    """System of ODEs after Yasukawa, H., Yoshimura, Y. (2015)
    for the MMG standard model

    Args:
        t (np.ndarray): time
        X (np.ndarray): Initial values
        params (dict):  Vessel dict 
        psi (float): vessel heading in the global frame
        fl_psi (float): Attack angle of current relative
                        longitudinal axis of motion [rad]
        fl_vel (float): Velocity of current [m/s]
        nps_old (float): propeller rotation last timestep [s⁻¹]
        delta_old (float): rudder angle last timestep [rad]

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
        beta = math.atan2(-v_m, u)  # Drift angle at midship position
        v_dash = v_m / U  # Non-dimensionalized lateral velocity
        r_dash = r * p["Lpp"] / U  # Non-dimensionalized yaw rate

    # Redefine
    w_P = p["w_P0"] * math.exp(-4.0 * (beta - (p["x_P"]/p["Lpp"]) * r_dash)**2)

    if nps == 0.0:  # No propeller movement, no advance ratio
        J = 0.0
    else:
        J = (1 - w_P) * u / (nps * p["D_p"])  # Propeller advance ratio

    if all(k in p for k in ("k_0", "k_1", "k_2")):
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
    alpha_R = delta - math.atan2(v_R, u_R)

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
    if fl_vel is not None and fl_vel != 0.:

        # Longitudinal velocity of current dependent on ship heading
        u_c = fl_vel * math.cos(fl_psi - psi)
        u_rc = (u - u_c)

        # Lateral velocity of current dependent on ship heading
        v_c = fl_vel * math.sin(fl_psi - psi)
        v_rc = (v_m - v_c)

        V_Rc = math.sqrt(u_rc**2 + v_rc**2)

        gamma_rc = -math.atan2(v_rc,u_rc)

        # Longitudinal current force
        A_Fc = p["B"] * p["d"] * p["C_b"]
        X_C = 0.5 * p["rho"] * A_Fc * C_X(gamma_rc) * V_Rc**2

        # Lateral current force
        A_Lc = p["Lpp"] * p["d"] * p["C_b"]
        Y_C = 0.5 * p["rho"] * A_Lc * C_Y(gamma_rc) * V_Rc**2

        # Current Moment
        N_C = 0.5 * p["rho"] * A_Lc * p["Lpp"] * C_N(gamma_rc) * V_Rc**2

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

    # Derivatives for delta and nps
    d_delta = (delta - delta_old) * t
    d_nps = (nps - nps_old) * t

    return np.array([d_u, d_v, d_r, d_delta, d_nps])


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
    v["X_vv_dash"] *= -70.8*(d/h)**4 + 27.7*(d/h)**2 + 1
    v["X_vr_dash"] *= 1.3*(d/h)**2 + 1

    # Hull frictional resistance coefficient
    v["R_0_dash"] *= 0.388*(d/h)**2 + 1

    # Correction by Tang et al 2020
    k = (2*d)/L
    def k_e(q): return k/(d/(2*h) + ((np.pi*d) /
                                     (2*h) * 1/math.tan((np.pi*d)/(2*h)))**q)

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
    if all(k in v for k in ["gamma_R_minus", "gamma_R_plus"]):
        v["gamma_R_minus"] *= flow_str_coef
        v["gamma_R_plus"] *= flow_str_coef
    else:
        v["gamma_R"] *= flow_str_coef

    coef = (water_depth/d) - 1

    mxshallow = (coef**1.3+3.77+1.14*(B/d)-0.233*(L/d)-3.43*Cb)/(coef**1.3)
    myshallow = (coef**0.82+0.413+0.0320*(B/d)+0.129*(B/d)**2)/(coef**0.82)
    Jzshallow = (coef**0.82+0.413+0.0192*(B/d)+0.00554*(B/d)**2)/(coef**0.82)

    v["m_x_dash"] *= mxshallow
    v["m_y_dash"] *= myshallow
    v["J_z_dash"] *= Jzshallow


def C_X(gamma_rc: float) -> float:

    return -0.0665*gamma_rc**5 + 0.5228*gamma_rc**4 - 1.4365*gamma_rc**3 \
    + 1.6024*gamma_rc**2 - 0.2967*gamma_rc - 0.4691

def C_Y(gamma_rc: float) -> float:

    return 0.1273*gamma_rc**4 - 0.802*gamma_rc**3 + 1.3216*gamma_rc**2\
        - 0.1799*gamma_rc

def C_N(gamma_rc: float) -> float:

    return -0.014*gamma_rc**5 + 0.1131*gamma_rc**4 - 0.2757*gamma_rc**3\
        + 0.1617*gamma_rc**2 + 0.0728*gamma_rc

# In here all the important hydrodynamic derivatives are calculated via empirical
# formulas from Suras and Sakir Bal (2019)


def calibrate(v: Dict[str, float], rho: float) -> Dict[str, float]:
    """Calculate relevant hydrodynamic derivatives based on a minimal
    dict and a given water density

    Sources:
        Suras and Sakir Bal (2019)

    Args:
        v (dict): Minimal dict of vessel parameters
        rho (float): water density

    Raises:
        KeyError: If your dict is incomplete, calculations will be aborted

    Returns:
        dict: Dict with all relevant information to be used as input in the step() function
    """

    # Length, Breadth, draft, Block coef, mass, Rudder area
    mininfo = ["B", "Lpp", "d", "C_b", "m", "A_R"]

    if not all(k in v for k in mininfo):
        raise KeyError(
            f"Mandatory keys not found. Need at least ['B','Lpp','d','C_b','A_R']")

    # Unpack initial dict
    L, B, d, Cb, m = v["Lpp"], v["B"], v["d"], v["C_b"], v["m"]

    # X-Coordinate of the center of gravity (m)
    if "x_G" not in v:
        v["x_G"] = float(0)

    # X-Coordinate of the propeller (assumed as -0.5*Lpp)
    v["x_P"] = -0.5*L

    # Add current water density to dict
    v["rho"] = rho

    # Masses and Moment of Intertia
    nondim_M = 0.5 * v["rho"] * L**2 * d
    nondim_N = 0.5 * v["rho"] * L**4 * d
    v["m"] = m * v["rho"]  # Displacement * water density
    v["I_zG"] = v["m"] * (0.25*L)**2
    v["m_x_dash"] = 0.05*v["m"] / nondim_M
    v["m_y_dash"] = (0.882-0.54*Cb*(1-1.6*d/B)-0.156*(1-0.673*Cb)*L/B +
                     0.826*((d*L)/(B**2))*(1-0.678*d/B)-0.638*Cb*((d*L)/(B**2))*(1-0.669*d/B))*v["m"] / nondim_M
    v["J_z_dash"] = ((0.01*(33-76.85*Cb*(1-0.784*Cb)+3.43 *
                     L/B*(1-0.63*Cb)))**2)*v["m"] / nondim_N

    # Hydrodynamic derivatives
    # Surge
    v["X_vv_dash"] = 1.15*(Cb/(L/B)) - 0.18  # Yoshimura and Masumoto (2012)
    v["X_vvvv_dash"] = -6.68*(Cb/(L/B)) + 1.1  # Yoshimura and Masumoto (2012)
    v["X_rr_dash"] = (-0.0027+0.0076*Cb*d/B)*L/d  # Lee et al. (1998)
    v["X_vr_dash"] = v["m_y_dash"] - 1.91 * \
        (Cb/(L/B)) + 0.08  # Yoshimura and Masumoto (2012)

    # Sway
    v["Y_v_dash"] = -(0.5*math.pi*(2*d/L)+1.4*(Cb/(L/B))
                      )  # Yoshimura and Masumoto (2012)
    v["Y_vvv_dash"] = -0.185*L/B + 0.48  # Yoshimura and Masumoto (2012)
    v["Y_r_dash"] = v["m_x_dash"] + 0.5*Cb*B/L  # Yoshimura and Masumoto (2012)
    v["Y_rrr_dash"] = -0.051
    v["Y_vrr_dash"] = -(0.26*(1-Cb)*L/B + 0.11)
    v["Y_vvr_dash"] = -0.75  # TODO: Change this

    # Yaw
    v["N_v_dash"] = -2*d/L  # Yoshimura and Masumoto (2012)
    v["N_vvv_dash"] = -(-0.69*Cb+0.66)  # Yoshimura and Masumoto (2012)
    v["N_r_dash"] = -0.54*(2*d/L) + (2*d/L)**2  # Yoshimura and Masumoto (2012)
    v["N_rrr_dash"] = ((0.25*Cb)/(L/B))-0.056  # Yoshimura and Masumoto (2012)
    v["N_vrr_dash"] = -0.075*(1-Cb)*L/B-0.098  # Yoshimura and Masumoto (2012)
    v["N_vvr_dash"] = ((1.55*Cb)/(L/B))-0.76  # Yoshimura and Masumoto (2012)

    # Wake coefficient
    if "w_P0" not in v:
        v["w_P0"] = 0.5*Cb - 0.05

    # Thrust deduction factor
    if "t_P" not in v:
        v["t_P"] = -0.27

    # Rudder force incease factor
    v["a_H"] = 0.627*Cb-0.153  # Quadvlieg (2013)

    # Longituinal acting point of longitudinal force
    v["x_H_dash"] = -0.45  # Aoki et. al. (2006)

    # Steering resistance deduction factor
    v["t_R"] = 0.39  # Yoshimura and Masumoto (2012)

    # Espilon [Lee and Shin (1998)]
    v["epsilon"] = -2.3281+8.697*Cb-3.78 * \
        ((2*d)/L)+1.19*Cb**2+292*((2*d)/L)**2 - \
        81.51*Cb*((2*d)/L)  # Lee and Shin (1998)
    # v["epsilon"] = -156.5*(Cb*B/L)**2 + 41.6*(Cb*B/L) - 1.76 # Kijima (1990)
    # v["epsilon"] = 2.26*1.82*(1-v["w_P0"]) # Yoshimura and Masumoto (2012)

    # Kappa
    v["kappa"] = 0.55-0.8*Cb*B/L

    # Correction of flow straightening factor to yaw-rate
    v["l_R"] = -0.9  # Yoshimura and Masumoto (2012)

    # Flow straightening coefficient
    v["gamma_R"] = 2.06*Cb*B/L+0.14

    # Additional assumptions
    v["J_int"] = 0.4
    v["k_0"] = 0.4
    v["J_slo"] = -0.5

    # Rudder aspect ratio
    if "f_alpha" not in v:
        while True:
            asp = input(
                "No rudder aspect ratio found. Please enter manually: ")
            try:
                asp = float(asp)
                break
            except ValueError as e:
                print("Wrong input type. Only float or int allowed. Err:", e)

        asp = float(asp)
        falpha = (6.13*asp)/(asp+2.25)
        v["f_alpha"] = falpha

        print(
            f"Aspect ration saved. Rudder lift coef calculated to be {falpha}")

    # Frictional resistance coefficent [ARBITRARY, no calculations so far]
    v["R_0_dash"] = 0.025

    return v


class LogicError(Exception):
    pass
