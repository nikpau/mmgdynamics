import math
import numpy as np

from .structs import MinimalVessel, Vessel

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


def mmg_dynamics(t: np.ndarray, X: np.ndarray, params: Vessel, 
                psi:float,delta: float, nps: float, fl_psi: float, 
                fl_vel: float,w_vel:float, beta_w: float, dT: float) -> np.ndarray:
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
    u, v_m, r  = X
    
    if fl_vel is not None:
        u_c = fl_vel * math.cos(fl_psi - psi - math.pi)
        v_c = fl_vel * math.sin(fl_psi - psi - math.pi)
        
        # Relative velocities
        u_r = u - u_c
        v_r = v_m - v_c

        # Hydrodynamic forces with relative velocities
        u_pure = u
        vm_pure = v_m
        u = u_r
        v_m = v_r

    U = math.sqrt(u**2 + v_m**2)  # Overall speed of the vessel

    if U == 0.0:  # No vessel movement. Velocity in all directions = 0
        beta = 0.0
        v_dash = 0.0
        r_dash = 0.0
    else:
        beta = math.atan(-v_m/u)  # Drift angle at midship position
        v_dash = v_m / U  # Non-dimensionalized lateral velocity
        r_dash = r * p.Lpp / U  # Non-dimensionalized yaw rate

    # Redefine
    beta_P = beta - (p.x_P/p.Lpp) * r_dash
    if all(getattr(p,k) is not None for k in ("C_1","C_2_plus","C_2_minus")):
        C_2 = p.C_2_plus if beta_P >= 0 else p.C_2_minus
        tmp = 1-math.exp(-p.C_1*abs(beta_P))*(C_2-1)
        w_P = 1-(1-p.w_P0)*(1+tmp)
    else:
        w_P = p.w_P0 * math.exp(-4.0 * beta_P**2)

    if nps == 0.0:  # No propeller movement, no advance ratio
        J = 0.0
    else:
        J = (1 - w_P) * u / (nps * p.D_p)  # Propeller advance ratio

    if all(getattr(p,k) is not None for k in ("k_0", "k_1", "k_2")):
        # Propeller thrust open water characteristic
        K_T = p.k_0 + (p.k_1 * J) + (p.k_2 * J**2)
    else:
        # Inferred slope + intercept dependent on J (empirical)
        K_T = p.J_slo * J + p.J_int

    # Effective inflow angle to rudder in maneuvering motions
    beta_R = beta - p.l_R * r_dash

    # Flow straightening coefficient
    if p.gamma_R is not None:
        gamma_R = p.gamma_R
    else:
        if beta_R < 0.0:
            gamma_R = p.gamma_R_minus
        else:
            gamma_R = p.gamma_R_plus

    # Lateral inflow velocity components to rudder
    v_R = U * gamma_R * beta_R

    # Longitudinal inflow velocity components to rudder
    if J == 0.0:
        u_R = 0.0
    else:
        u_R = u * (1 - w_P) * p.epsilon * math.sqrt(
            p.eta * (1.0 + p.kappa * (
                math.sqrt(1.0 + 8.0 * K_T / (math.pi * J**2)) - 1))**2 + (1 - p.eta)
        )
    # Rudder inflow velocity
    U_R = math.sqrt(u_R**2 + v_R**2)

    # Rudder inflow angle
    if p.delta_prop is not None:
        alpha_R = delta - p.delta_prop - v_R/u_R
    else:
        alpha_R = delta - math.atan2(v_R, u_R)

    # Normal force on rudder
    if p.A_R is not None:
        F_N = 0.5 * p.A_R * p.rho * \
            p.f_alpha * (U_R**2) * math.sin(alpha_R)
    else:
        F_N = 0.5 * p.A_R_Ld_em * (p.Lpp * p.d * p.rho) * \
            p.f_alpha * (U_R**2) * math.sin(alpha_R)

    # Longitudinal surge force around midship acting on ship hull
    X_H = (0.5 * p.rho * p.Lpp * p.d * (U**2) * (
        - p.R_0_dash
        + p.X_vv_dash * (v_dash**2)
        + p.X_vr_dash * v_dash * r_dash
        + p.X_rr_dash * (r_dash**2)
        + p.X_vvvv_dash * (v_dash**4)
    )
    )

    # Longitudinal surge force around midship by steering
    X_R = -(1 - p.t_R) * F_N * math.sin(delta)

    # Longitudinal Surge force due to propeller
    X_P = (1 - p.t_P) * p.rho * K_T * nps**2 * p.D_p**4

    # Longitudinal surge force around midship acting on ship hull
    Y_H = (0.5 * p.rho * p.Lpp * p.d * (U**2) * (
        p.Y_v_dash * v_dash
        + p.Y_r_dash * r_dash
        + p.Y_vvv_dash * (v_dash**3)
        + p.Y_vvr_dash * (v_dash**2) * r_dash
        + p.Y_vrr_dash * v_dash * (r_dash**2)
        + p.Y_rrr_dash * (r_dash**3)
    )
    )

    # Lateral surge force by steering
    Y_R = -(1 + p.a_H) * F_N * math.cos(delta)

    # Yaw moment around midship acting on ship hull
    N_H = (0.5 * p.rho * (p.Lpp**2) * p.d * (U**2) * (
        p.N_v_dash * v_dash
        + p.N_r_dash * r_dash
        + p.N_vvv_dash * (v_dash**3)
        + p.N_vvr_dash * (v_dash**2) * r_dash
        + p.N_vrr_dash * v_dash * (r_dash**2)
        + p.N_rrr_dash * (r_dash**3)
    )
    )

    # Redimensionalize x_H
    x_H = p.x_H_dash * p.Lpp

    # yaw moment around midship by steering
    N_R = -(-0.5*p.Lpp + p.a_H * x_H) * F_N * math.cos(delta)

        #-------------------------------- Wind ---------------------------------
    if w_vel != 0.0 and w_vel is not None:
        
        # relative velocity computations
        u_w = w_vel * math.cos(beta_w - psi)
        v_w = w_vel * math.sin(beta_w - psi)

        u_rw = u - u_w
        v_rw = v_m - v_w

        V_rw_sq = u_rw**2 + v_rw**2
        g_rw = -math.atan2(v_rw, u_rw)

        # forces               
        X_W = 1/2 * p.rho_air * V_rw_sq * _C_X_wind(g_rw) * p.A_Fw
        Y_W = 1/2 * p.rho_air * V_rw_sq * _C_Y_wind(g_rw) * p.A_Lw
        N_W = 1/2 * p.rho_air * V_rw_sq * _C_N_wind(g_rw) * p.A_Lw * p.Lpp

    else:
        X_W, Y_W, N_W = 0.0, 0.0, 0.0


    # Added masses and added moment of inertia
    m_x = p.m_x_dash * (0.5 * p.rho * (p.Lpp**2) * p.d)
    m_y = p.m_y_dash * (0.5 * p.rho * (p.Lpp**2) * p.d)
    J_z = p.J_z_dash * (0.5 * p.rho * (p.Lpp**4) * p.d)
    m = p.displ*p.rho
    I_zG = m*(0.25*p.Lpp)**2
    
    # Mass matrices
    M_RB = np.array([[m, 0.0, 0.0],
                    [0.0, m, m * p.x_G],
                    [0.0, m * p.x_G, I_zG]])
    M_A = np.array([[m_x, 0.0, 0.0],
                    [0.0, m_y, 0.0],
                    [0.0, 0.0, J_z + (p.x_G**2) * m]])
    M_inv = np.linalg.inv(M_RB + M_A)

    FX = X_H + X_R + X_P + X_W
    FY = Y_H + Y_R + Y_W
    FN = N_H + N_R + N_W
    
    F = np.array([FX,FY,FN])
    
    if fl_vel != 0.0 and fl_vel is not None:
        nu_c_dot = np.array([v_c*r, -u_c*r, 0.0])
        return np.dot(
            M_inv, F - np.dot(_C_RB(m=m,x_G=p.x_G,r=r), np.array([u_pure, vm_pure, r]))     
            - np.dot(_C_A(vm=v_m,u=u,m_x=m_x,m_y=m_y), np.array([u_r, v_r, r])  
            + np.dot(M_A, nu_c_dot))
        )
    else:
        return np.dot(
            M_inv, F - np.dot((_C_RB(m=m,x_G=p.x_G,r=r) + _C_A(vm=v_m,u=u,m_x=m_x,m_y=m_y)), 
                              np.array([u, v_m, r]))
            )



def _shallow_water_hdm(v: Vessel, water_depth: float) -> None:
    """Correct the hydrodynamic derivatives and
    hydrodynamic masses for shallow water conditions.

    Sources:
        Taimuri et. al. (2020) https://doi.org/10.1016/j.oceaneng.2020.108103


    Args:
        v (dict): vessel parameter dict
        water_depth (float): water depth around vessel [m]
    """
    frac = lambda x,y: x/y

    # Length, Beam (width), Draft, Block Coefficient
    L, B, T, Cb = v.Lpp, v.B, v.d, v.C_b

    H = water_depth

    HT = H/T -1
    TH = T/H
    K0 = 1+frac(0.0775,HT**2)-frac(0.011,HT**3)+frac(0.000068,HT**5)
    K1 = -frac(0.0643,HT)+frac(0.0724,HT**2)-frac(0.0113,HT**3)+frac(0.0000767,HT**5)
    K2 = frac(0.0342,HT) if B/T <= 4 else frac(0.137*B,HT*T)
    B1 = Cb*B*(1+frac(B,L))**2
    B2 = 0.83*frac(B1,Cb)

    A1Yr   = -5.5*(frac(Cb*B,T))**2+26*frac(Cb*B,T)-31.5
    A1Yrr  = -15600*(1-Cb)**5
    A1Yvrr = 21500*(((1-Cb)*frac(T,B))**2)-4800*((1-Cb)*frac(T,B))+220
    A1Nvv  = -240*(1-Cb)+57
    A1Nrr  = -1960*((1-Cb)*frac(T,B))**2+448*(1-Cb)*frac(T,B)-25
    A1Nvvr = 91*Cb*frac(T,B)-25
    A1Nvrr = 40*Cb*frac(B,T)-88
    A2Yr   = 37*frac(Cb*B,T)**2-185*frac(Cb*B,T)+230
    A2Yrr  = 116000*(1-Cb)**5
    A2Yvrr = -40800*((1-Cb)*frac(T,B))**2+7500*(1-Cb)*frac(T,B)-274
    A2Nvv  = 1770*(1-Cb)-413
    A2Nrr  = 12220*((1-Cb)*frac(T,B))**2-2720*(1-Cb)*frac(T,B)+146
    A2Nvvr = -515*Cb*frac(T,B)+144
    A2Nvrr = -295*Cb*frac(B,T)+645
    A3Yr   = -38*frac(Cb*B,T)**2+197*frac(Cb*B,T)-250
    A3Yrr  = -128000*(1-Cb)**5
    A3Yvvr = -90800*((1-Cb)*frac(T,B))**2+25500*(1-Cb)*frac(T,B)-1400
    A3Nvv  = -1980*(1-Cb)+467
    A3Nrr  = -12160*((1-Cb)*frac(T,B))**2+2650*(1-Cb)*frac(T,B)-137
    A3Nvvr = 508*Cb*frac(T,B)-143
    A3Nvrr = 312*Cb*frac(B,T)-678

    gv = K0+frac(2,3)*K1*frac(B1,T)+frac(8,15)*K2*(frac(B1,T))**2
    gnr = K0+frac(8,15)*K1*frac(B1,T)+frac(40,105)*K2*(frac(B1,T))**2
    fyr = K0+frac(2,5)*K1*frac(B1,T)+frac(24,105)*K2*(frac(B1,T))**2
    fnr = K0+frac(1,2)*K1*frac(B1,T)+frac(1,3)*K2*(frac(B1,T))**2
    fyv = 1.5*fnr-0.5
    fnv = K0+K1*frac(B1,T)+K2*frac(B1,T)**2

    # Corrections 

    v.X_vv_dash *= fyv
    v.X_vvvv_dash *= fyv
    v.X_rr_dash *= fnr
    v.X_vr_dash *= fyr
    v.Y_vvv_dash *= fyv
    v.N_v_dash *= fnv
    v.N_vvv_dash *= fyv
    v.Y_rrr_dash *= gnr
    v.N_rrr_dash *= gnr
    v.Y_vvr_dash *= fyv
    v.Y_vrr_dash *= fyv

    v.Y_v_dash *= (-TH+frac(1,(1-TH)**(frac(0.4*Cb*B,T))))
    v.Y_r_dash *= (1+A1Yr*TH+A2Yr*TH**2+A3Yr*TH**3)
    v.N_r_dash *= (-TH+frac(1,(1-TH)**(frac(-14.28*T,L)+1.5)))

    v.N_vvr_dash *= (1+A1Nvvr*TH+A2Nvvr*TH**2+A3Nvvr*TH**3)
    v.N_vrr_dash *= (1+A1Nvrr*TH+A2Nvrr*TH**2+A3Nvrr*TH**3)

    # Corrections for wake fraction, thrust deduction,
    # and flow-straighening coefficients
    v.w_P0 *= (1+(-4.932+0.6425*frac(Cb*L,T)-0.0165*(frac(Cb*L,T)**2))*TH**1.655)
    ctp = 1+((29.495-14.089*frac(Cb*L,B)+1.6486*frac(Cb*L,B)**2)*(frac(1,250)-frac(7*TH,200)-frac(13*TH**2,125)))
    v.t_P = 1-ctp*(1-v.t_P)
    cgr1 = 1+((frac(-5129,500)+178.207*frac(Cb*B,L)-frac(2745,4)*frac(Cb*B,L)**2)*(frac(-1927,500)+frac(2733*TH,200)-frac(2617*TH**2,250)))
    cgr2 = 1+(frac(-541,4)+2432.95*frac(Cb*B,L)-10137.7*frac(Cb*B,L)**2)*TH**4.81
    if TH <= (-0.332*frac(T,B)+0.581):
        v.gamma_R_minus *= cgr2
        v.gamma_R_plus *= cgr2
    else:
        v.gamma_R_minus *= cgr1
        v.gamma_R_plus *= cgr1


def _C_X_wind(g_w, cx=0.9):
    return -cx*math.cos(g_w)

def _C_Y_wind(g_w, cy=0.95):
    return cy*math.sin(g_w)

def _C_N_wind(g_w, cn=0.2):
    return cn*math.sin(2*g_w)    

def _C_RB(*,m,x_G, r):
    return np.array(
        [[0.0, -m * r, -m * x_G * r],
        [m * r, 0.0, 0.0],
        [m * x_G * r, 0.0, 0.0]])

def _C_A(*,m_x,m_y,u, vm):
    return np.array(
        [[0.0, 0.0, -m_y * vm],
        [0.0, 0.0, m_x * u],
        [0.0, 0.0, 0.0]])

def calibrate(v: MinimalVessel, rho: float) -> Vessel:
    """Calculate relevant hydrodynamic derivatives based on a minimal
    dict and a given water density

    Sources:
        Suras and Sakir Bal (2019)

    Args:
        v (MinimalVessel): Minimal vessel parameters
        rho (float): water density

    Raises:
        KeyError: If your dict is incomplete, calculations will be aborted

    Returns:
        dict: Dict with all relevant information to be used as input in the step() function
    """

    # Unpack initial dict
    L, B, d, Cb, m = v.Lpp, v.B, v.d, v.C_b, v.m

    # X-Coordinate of the center of gravity (m)
    if v.x_G is None:
        v.x_G = float(0)

    # X-Coordinate of the propeller (assumed as -0.5*Lpp)
    v.x_P = -0.5*L

    # Add current water density to dict
    v.rho = rho

    # Masses and Moment of Inertia
    nondim_M = 0.5 * v.rho * L**2 * d
    nondim_N = 0.5 * v.rho * L**4 * d
    v.m = m * v.rho  # Displacement * water density
    v.I_zG = v.m * (0.25*L)**2
    v.m_x_dash = 0.05*v.m / nondim_M
    v.m_y_dash = (0.882-0.54*Cb*(1-1.6*d/B)-0.156*(1-0.673*Cb)*L/B +
                     0.826*((d*L)/(B**2))*(1-0.678*d/B)-0.638*Cb*((d*L)/(B**2))*(1-0.669*d/B))*v.m / nondim_M
    v.J_z_dash = (L/100*(33-76.85*Cb*(1-0.784*Cb)+3.43*L/B*(1-0.63*Cb)))**2*m/nondim_N

    # Hydrodynamic derivatives
    # Surge
    #v.X_vv_dash = 1.15*(Cb/(L/B)) - 0.18  # Yoshimura and Masumoto (2012)
    v.X_vv_dash = 0.0014-0.1975*d*(1-Cb)/B*L/d # Lee et. al. (1998)
    v.X_vvvv_dash = -6.68*(Cb/(L/B)) + 1.1  # Yoshimura and Masumoto (2012)
    #v.X_rr_dash = (-0.0027+0.0076*Cb*d/B)*L/d  # Lee et al. (1998)
    v.X_rr_dash = (-0.0027+0.0076*Cb*d/B)*L/d # Lee et al. (1998)
    #v.X_vr_dash = v.m_y_dash - 1.91 * (Cb/(L/B)) + 0.08  # Yoshimura and Masumoto (2012)
    v.X_vr_dash = (m+0.1176*v.m_y_dash*(0.5+Cb))*L/d

    # Sway
    #v.Y_v_dash = -(0.5*math.pi*(2*d/L)+1.4*(Cb/(L/B)))  # Yoshimura and Masumoto (2012)
    v.Y_v_dash = -(0.5*math.pi*2*d/L+1.4*Cb*B/L) # Kijima et. al. (1990)
    #v.Y_vvv_dash = -0.185*L/B + 0.48  # Yoshimura and Masumoto (2012)
    v.Y_vvv_dash = (-0.6469*(1-Cb)*d/B+0.0027)*L/d # Lee et al. (1998)
    #v.Y_r_dash = v.m_x_dash + 0.5*Cb*B/L  # Yoshimura and Masumoto (2012)
    v.Y_r_dash = (-0.115*Cb*B/L+0.0024)*L/d # Lee et al. (1998)
    #v.Y_rrr_dash = -0.051
    v.Y_rrr_dash = (-0.233*Cb*d/B+0.0063)*L/d # Lee et al. (1998)
    #v.Y_vrr_dash = -(0.26*(1-Cb)*L/B + 0.11)
    v.Y_vrr_dash = -(5.95*d*(1-Cb)/d) # Kijima et. al. (1990)
    #v.Y_vvr_dash = -0.75  # TODO: Change this
    v.Y_vvr_dash = 1.5*d*Cb/B-0.65 # Kijima et. al. (1990)

    # Yaw
    v.N_v_dash = -2*d/L  # Yoshimura and Masumoto (2012)
    #v.N_vvv_dash = -(-0.69*Cb+0.66)  # Yoshimura and Masumoto (2012)
    v.N_vvv_dash = (0.0348-0.5283*(1-Cb)*d/B)*L/d # Lee et al. (1998)
    #v.N_r_dash = -0.54*(2*d/L) + (2*d/L)**2  # Yoshimura and Masumoto (2012)
    v.N_r_dash = -0.54*2*d/L+(2*d/L)**2 # Kijima et. al. (1990)
    #v.N_rrr_dash = ((0.25*Cb)/(L/B))-0.056  # Yoshimura and Masumoto (2012)
    v.N_rrr_dash = (-0.0572+0.03*Cb*d/L)*L/d # Lee et al. (1998)
    #v.N_vrr_dash = -0.075*(1-Cb)*L/B-0.098  # Yoshimura and Masumoto (2012)
    v.N_vrr_dash = (0.5*d*Cb/B)-0.05 # Kijima et. al. (1990)
    #v.N_vvr_dash = ((1.55*Cb)/(L/B))-0.76  # Yoshimura and Masumoto (2012)
    v.N_vvr_dash = -(57.5*(Cb*B/L)**2-18.4*(Cb*B/L)+1.6 ) # Kijima et. al. (1990)

    # Wake coefficient
    if v.w_P0 is None:
        v.w_P0 = 0.5*Cb - 0.05

    # Thrust deduction factor
    if v.t_P is None:
        v.t_P = -0.27

    # Rudder force increase factor
    v.a_H = 0.627*Cb-0.153  # Quadvlieg (2013)

    # Longitudinal acting point of longitudinal force
    v.x_H_dash = -0.37  # Khanﬁr et al. (2011)

    # Steering resistance deduction factor
    v.t_R = 0.39  # Yoshimura and Masumoto (2012)

    # Espilon [Lee and Shin (1998)]
    v.epsilon = -2.3281+8.697*Cb-3.78 * \
        ((2*d)/L)+1.19*Cb**2+292*((2*d)/L)**2 - \
        81.51*Cb*((2*d)/L)  # Lee and Shin (1998)
    # v["epsilon"] = -156.5*(Cb*B/L)**2 + 41.6*(Cb*B/L) - 1.76 # Kijima (1990)
    # v["epsilon"] = 2.26*1.82*(1-v.w_P0) # Yoshimura and Masumoto (2012)

    # Kappa
    v.kappa = 0.55-0.8*Cb*B/L

    # Correction of flow straightening factor to yaw-rate
    v.l_R = -0.9  # Yoshimura and Masumoto (2012)

    # Flow straightening coefficient
    v.gamma_R = 2.06*Cb*B/L+0.14

    # Additional assumptions
    v.J_int = 0.4
    v.k_0 = 0.4
    v.J_slo = -0.5

    # Rudder aspect ratio
    if v.f_alpha is None:
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
        v.f_alpha = falpha

        print(
            f"Aspect ratio saved. Rudder lift coef calculated to be {falpha}")

    # Frictional resistance coefficent [ARBITRARY, no calculations so far]
    v.R_0_dash = 0.025

    return Vessel(**v.__dict__)
