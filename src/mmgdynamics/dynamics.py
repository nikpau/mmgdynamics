import math
import numpy as np

from .structs import InlandVessel

"""
Set up the System of ODEs for vessel maneuvering prediction after 
Yasukawa, H., Yoshimura, Y. (2015) Introduction of MMG standard method for ship maneuvering predictions.
J Mar Sci Technol 20, 37-52 https://doi.org/10.1007/s00773-014-0293-y


On each step, the longituidal, lateral and yaw-rate accelerations are calculated
from a set of initial conditions X for an arbitrary time horizon. 
"""

# Export only selected functions if used as wildcard import
__all__ = ["calilbrate"]

# System of ODEs after Yasukawa, H., Yoshimura, Y. (2015)


def mmg_dynamics(t: np.ndarray, X: np.ndarray, params: InlandVessel, psi:float,
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
    v_m = v - p.x_G * r
    U = math.sqrt(u**2 + v_m**2)  # Overall speed of the vessel

    if U == 0.0:  # No vessel movement. Velocity in all directions = 0
        beta = 0.0
        r_dash = 0.0
    else:
        beta = math.atan2(-v_m, u)  # Drift angle at midship position
        r_dash = r * p.Lpp / U  # Non-dimensionalized yaw 
    # Redefine
    beta_P = beta - (p.x_P/p.Lpp) * r_dash
    w_P = p.w_P0 * math.exp(-4.0 * (beta_P)**2)

    if nps == 0.0:  # No propeller movement, no advance ratio
        J = 0.0
    else:
        J = (1 - w_P) * u / (nps * p.D_p)  # Propeller advance ratio

    # Propeller thrust coefficient (starboard/port)
    K_Ts = p.k_0s + (p.k_1s * J) + (p.k_2s * J**2)
    K_Tp = p.k_0p + (p.k_1p * J) + (p.k_2p * J**2)

    # Effective inflow angle to rudder in maneuvering motions
    beta_R = beta - p.l_R * r_dash

    # Lateral inflow velocity components to rudder
    v_Rs = U * p.gamma_R * beta_R
    v_Rp = U * p.gamma_R * beta_R

    # Wake fraction at the position of the rudder during maneuvering motions
    w_R = p.w_R0*w_P/p.w_P0

    # Longitudinal inflow velocity components to rudder
    if J == 0.0:
        u_R = math.sqrt(p.eta * (p.kappa  *
                                    8.0 * p.k_0 * nps ** 2 * p.D_p**4 / np.pi)**2)
    else:
        u_Rs = (u - p.y_rs*r) * (1 - w_R)  * math.sqrt(
            p.eta * (1.0 + p.kappa * (
                math.sqrt(1.0 + 8.0 * K_Ts / (np.pi * J**2)) - 1))**2 + (1 - p.eta)
        )
        u_Rp = (u - p.y_rp*r) * (1 - w_R)  * math.sqrt(
            p.eta * (1.0 + p.kappa * (
                math.sqrt(1.0 + 8.0 * K_Tp / (np.pi * J**2)) - 1))**2 + (1 - p.eta)
        )
        
    # Rudder inflow velocity
    U_Rs = math.sqrt(u_Rs**2 + v_Rs**2)
    U_Rp = math.sqrt(u_Rp**2 + v_Rp**2)

    # Rudder inflow angle
    alpha_Rs = delta - math.atan2(v_Rs, u_Rs)
    alpha_Rp = delta - math.atan2(v_Rp, u_Rp)
    delta_hs = math.atan2(v_Rs, u_Rs)
    delta_hp = math.atan2(v_Rp, u_Rp)

    # Coefficients for rudder lift an drag coefficients
    C_Ls = 6.175*math.sin(alpha_Rs) * p.asp/(p.asp+2.25)
    C_Lp = 6.175*math.sin(alpha_Rp) * p.asp/(p.asp+2.25)
    C_Ds = 0.032*math.sin(alpha_Rs) * p.asp/(p.asp+2.25)
    C_Dp = 0.032*math.sin(alpha_Rp) * p.asp/(p.asp+2.25)

    # Rudder lift and drag forces (*2 for two rudders)
    F_Ls = 0.5 * p.rho * p.A_R * C_Ls * U_Rs**2
    F_Lp = 0.5 * p.rho * p.A_R * C_Lp * U_Rp**2
    F_Ds = 0.5 * p.rho * p.A_R * C_Ds * U_Rs**2
    F_Dp = 0.5 * p.rho * p.A_R * C_Dp * U_Rp**2

    # Longituinal and lateral rudder force
    F_RXs = F_Ls * math.sin(delta_hs) + F_Ds * math.cos(delta_hs)
    F_RXp = F_Lp * math.sin(delta_hp) + F_Dp * math.cos(delta_hp)
    F_RYs = F_Ls * math.cos(delta_hs) - F_Ds * math.sin(delta_hs)
    F_RYp = F_Lp * math.cos(delta_hp) - F_Dp * math.sin(delta_hp)

    # Added masses and added moment of inertia
    m_x = p.m_x_dash * (0.5 * p.rho * (p.Lpp**2) * p.d)
    m_y = p.m_y_dash * (0.5 * p.rho * (p.Lpp**2) * p.d)
    J_z = p.J_z_dash * (0.5 * p.rho * (p.Lpp**4) * p.d)
    m = p.m
    I_zG = p.I_zG

    # Longitudinal surge force around midship acting on ship hull
    X_H = ((0.5 * p.rho * p.Lpp * p.d * (U**2)) * (
        p.X_0_dash + 
        p.X_bb_dash * beta**2 +
        (p.X_br_dash - p.m_y_dash)*beta*r_dash +
        (p.X_rr_dash + p.x_G/p.Lpp*p.m_y_dash)*r_dash**2 + 
        p.X_bbbb_dash * beta**4
    ))

    # Total rudder forces (Longitudinal, lateral, yaw moment)
    X_Rs = -(1 - p.t_R) * F_RXs
    X_Rp = -(1 - p.t_R) * F_RXp
    Y_Rs = -(1 + p.a_H) * F_RYs
    Y_Rp = -(1 + p.a_H) * F_RYp
    N_Rs = -(p.x_R*p.Lpp + p.a_H*p.x_H_dash*p.Lpp)*F_RYs + p.y_rs*(1 - p.t_R) * F_RXs
    N_Rp = -(p.x_R*p.Lpp + p.a_H*p.x_H_dash*p.Lpp)*F_RYp + p.y_rp*(1 - p.t_R) * F_RXp
    
    X_R = X_Rs + X_Rp
    Y_R = Y_Rs + Y_Rp
    N_R = N_Rs + N_Rp

    # Longitudinal Surge force due to propeller
    X_Ps = (1 - p.t_P) * p.rho * K_Ts * nps**2 * p.D_p**4
    X_Pp = (1 - p.t_P) * p.rho * K_Tp * nps**2 * p.D_p**4
    X_P = (X_Ps + X_Pp)

    N_Ps = -p.y_ps * (1 - p.t_P) * p.rho * K_Ts * nps**2 * p.D_p**4
    N_Pp = -p.y_pp * (1 - p.t_P) * p.rho * K_Tp * nps**2 * p.D_p**4
    N_P = (N_Ps + N_Pp)

    # Longitudinal surge force around midship acting on ship hull
    Y_H = (0.5 * p.rho * p.Lpp * p.d * (U**2)) * (
        p.Y_b_dash*beta+
        p.Y_r_dash*r_dash+
        p.Y_bb_dash*beta*abs(beta)+
        p.Y_rr_dash*r_dash*abs(r_dash)+
        (p.Y_bbr_dash*beta + p.Y_brr_dash*r_dash)*beta*r_dash
    )


    # Yaw moment around midship acting on ship hull
    N_H = (0.5 * p.rho * (p.Lpp**2) * p.d * (U**2)) * (
        p.N_b_dash*beta+
        p.N_r_dash*r_dash+
        p.N_bb_dash*beta*abs(beta)+
        p.N_rr_dash*r_dash*abs(r_dash)+
        (p.N_bbr_dash*beta + p.N_brr_dash*r_dash)*beta*r_dash
    )

    # Forces related to currents:
    if fl_vel is not None and fl_vel != 0.:
        
        # Longitudinal velocity of current dependent on ship heading
        u_c = -fl_vel * math.cos(fl_psi - psi)
        u_rc = u - u_c

        # Lateral velocity of current dependent on ship heading
        v_c = fl_vel * math.sin(fl_psi - psi)
        v_rc = v_m - v_c

        g_rc = abs(-math.atan2(v_rc,u_rc))

        # Longitudinal current force
        A_Fc = p.B * p.d * p.C_b
        X_C = 0.5 * p.rho * A_Fc * _C_X(g_rc) * abs(u_rc) * u_rc

        # Lateral current force
        A_Lc = p.Lpp * p.d * p.C_b
        Y_C = 0.5 * p.rho * A_Lc * _C_Y(g_rc) * abs(v_rc) * v_rc

        # Current Moment
        N_C = 0.5 * p.rho * A_Lc * p.Lpp * _C_N(g_rc) * abs(v_rc) * v_rc

    else:
        X_C, Y_C, N_C = 0.0, 0.0, 0.0

    # Longitudinal acceleration
    d_u = ((X_H + X_R + X_P + X_C) + (m + m_y) * v *
           r + p.x_G * m * (r**2)) / (m + m_x)

    # Lateral acceleration
    f = (I_zG + J_z + (p.x_G**2) * m)

    d_v = ((Y_H+Y_R+Y_C) - (m+m_x)*u*r - ((p.x_G*m*(N_H + N_R + N_P + N_C))/(f)) + 
           ((p.x_G**2*m**2*u*r)/(f)))/((m+m_y)-((p.x_G**2*m**2)/(f)))

    # Yaw rate acceleration
    d_r = ((N_H + N_R + N_P + N_C) - (p.x_G * m * d_v + p.x_G * m * u * r)) / \
        (I_zG + J_z + (p.x_G**2) * m)

    # Derivatives for delta and nps
    d_delta = (delta - delta_old) * t
    d_nps = (nps - nps_old) * t

    return np.array([d_u, d_v, d_r, d_delta, d_nps])

def _shallow_water_hdm(v: InlandVessel, water_depth: float) -> None:
    """Correct the hydrodynamic derivatives and
    hydrodynamic masses for shallow water conditions.

    Sources:
        Maimun et. al. (2011) https://doi.org/10.1016/j.oceaneng.2011.05.011
        Furukawa et. al. (2016) https://dx.doi.org/10.18451/978-3-939230-38-0_33
        Tang et. al (2020) https://doi.org/10.1016/j.oceaneng.2019.106679

    Args:
        v (InlandVessel): vessel parameter class
        water_depth (float): water depth below keel [m]
    """

    # Length, Beam (width), Draft, Block Coefficient
    L, B, d, Cb = v.Lpp, v.B, v.d, v.C_b

    h = d/water_depth
    k = 2*d/L

    # Correction by Furukawa et. al. (2016)
    # Hull frictional resistance coefficient
    v.X_0_dash *= 0.388*(h)**2 + 1

    # Correction by Maimun et. al. (2011)
    # First Correction
    def f(n: float) -> float:
        return pow(1-h,-n) - h

    v.Y_b_dash *= f(0.4*Cb*B/d)
    v.Y_bb_dash *= f(-0.26*Cb*B/d+1.74)
    v.Y_brr_dash *= f(-2.13*Cb/B+1.8)
    v.N_b_dash *= f(0.425*Cb*B/d)
    v.N_r_dash *= f(-7.14*k+1.5)

    # WARNING: f(.) now changes
    def f(a1: float,a2: float,a3: float) -> float:
        return 1+a1*h+a2*h**2+a3*h**3

    v.Y_rr_dash *= f(
        a1=-5.5*(Cb*B/d)**2+26*Cb*(B/d)-31.5,
        a2= 37*(Cb*(B/d))**2-185*Cb*(B/d)+230,
        a3= -38*(Cb*(B/d))**2+197*Cb*(B/d)-250
    )

    v.Y_rr_dash *= f(
        a1= -0.15*105*(1-Cb)**5,
        a2=1.16*105*(1-Cb)**5,
        a3=-1.28*105*(1-Cb)**5
    )

    # v.Y_bbr_dash *= f(
    #     a1=2.15*10**4*((d*(1-Cb)/B)**2)-0.48*10**4*d*(1-Cb)/B+220,
    #     a2=-4.08*10**4*((d*(1-Cb)/B)**2)-0.75*10**4*d*(1-Cb)/B-274,
    #     a3=-9.08*10**4*((d*(1-Cb)/B)**2)+2.55*10**4*d*(1-Cb)/B-1400
    # )

    # v.N_bb_dash *= f(
    #     a1=-0.24*10**3*(1-Cb)+57,
    #     a2=1.77*10**3*(1-Cb)-413,
    #     a3=-1.98*10**5*(1-Cb)+467
    # )

    # v.N_rr_dash *= f(
    #     a1=-0.196*10**4*((d*(1-Cb)/B)**2)+448*d*(1-Cb)/B-25,
    #     a2=1.222*10**4*((d*(1-Cb)/B)**2)-2720*d*(1-Cb)/B+446,
    #     a3=-1.216*10**4*((d*(1-Cb)/B)**2)+2650*d*(1-Cb)/B-137
    # )

    # v.N_bbr_dash *= f(
    #     a1=0.91*10**2*Cb*(d/B)-25,
    #     a2=-5.15*10**2*Cb*(d/B)+144,
    #     a3=5.08*10**2*Cb*(d/B)-143
    # )

    # v.N_brr_dash *= f(
    #     a1=0.4*10**3*Cb*(B/d)-88,
    #     a2=-2.95*10**3*Cb*(B/d)+645,
    #     a3=3.12*10**3*Cb*(B/d)-678
    # )

    # Correction by Tang et al 2020
    # Thrust deduction factor
    oneminustp = 1 - v.t_P
    oneminustp *= 1/(1-0.2*(h) + 0.7295*(h)**2)
    v.t_P = -oneminustp + 1

    # Wake fraction coefficient
    oneminuswp = 1 - v.w_P0
    oneminuswp *= math.cos(1.4*Cb*h)
    v.w_P0 = -oneminuswp + 1

    # Flow-straightening coefficient
    flow_str_coef = 1 + 0.0161*(h) + 4.4222*(h)**2 - 4.9825*(h)**3
    v.gamma_R *= flow_str_coef


def _C_X(g_rc: float) -> float:

    return (
       -0.0665*g_rc**5 + 
        0.5228*g_rc**4 - 
        1.4365*g_rc**3 + 
        1.6024*g_rc**2 - 
        0.2967*g_rc - 
        0.4691)


def _C_Y(g_rc: float) -> float:

    # return (
    #     0.1273*g_rc**4 - 
    #     0.8020*g_rc**3 + 
    #     1.3216*g_rc**2 - 
    #     0.1799*g_rc)
    return (
      0.05930686*g_rc**4 -
      0.37522028*g_rc**3 +
      0.46812233*g_rc**2 +
      0.39114522*g_rc -
      0.00273578
    )

def _C_N(g_rc: float) -> float:

    return (
       -0.0140*g_rc**5 + 
        0.1131*g_rc**4 -
        0.2757*g_rc**3 + 
        0.1617*g_rc**2 + 
        0.0728*g_rc)

def pow(base: float,exponent: float) -> float:
    return base**exponent