from dataclasses import dataclass
from typing import Optional

@dataclass
class InlandVessel:
    
    # Water density
    rho: float # Water density
    
    # Vessel particulars
    C_b: float # Block Coefficient
    Lpp: float # Length over pependiculars (m)/
    B: float # Overall width
    d: float # Ship draft (Tiefgang)
    w_P0: float # Wake fraction coefficient
    x_G: float # X-Coordinate of the center of gravity (m)
    D_p: float # Diameter of propeller (m)
    l_R: float # correction of flow straightening factor to yaw-rate
    eta: float # Ratio of propeller diameter to rudder span
    kappa: float # An experimental constant for expressing "u_R"
    A_R: float # Moveable rudder area
    t_R: float # Steering resistance deduction factor
    t_P: float # Thrust deduction factor
    x_H_dash: float # Longitudinal coordinate of acting point of the additional lateral force
    x_R: float # 
    x_P: float # 
    a_H: float # Rudder force increase factor
    w_R0: float
    y_ps: float
    y_pp: float
    y_rs: float
    y_rp: float
    # MMG hydrodynamic derivatives
    X_0_dash: float # frictional resistance coefficient
    
    X_bb_dash:   float
    X_br_dash:   float
    X_rr_dash:   float
    X_bbbb_dash: float
    Y_b_dash:    float
    Y_r_dash:    float
    Y_bb_dash:   float
    Y_rr_dash:   float
    Y_brr_dash:  float
    Y_bbr_dash:  float
    N_b_dash:    float
    N_r_dash:    float
    N_bb_dash:   float
    N_rr_dash:   float
    N_brr_dash:  float
    N_bbr_dash:  float
    
    # Masses and added masses
    m: float # Mass of ship as calculated by ▽*rho (displacement * water density)
    m_x_dash: float # Non dimensionalized added masses coefficient in x direction
    m_y_dash: float # Non dimensionalized added masses coefficient in y direction
    
    # Moment of inertia and added moment of inertia
    I_zG: float # Moment of inertia of ship around center of gravity (m*(0.25*Lpp)**2) (Point mass Inertia)
    J_z_dash: float # Added moment of inertia coefficient
    
    # Wake change coefficients and propeller advance ratio polynomial
    k_0s: float
    k_1s: float
    k_2s: float
    k_0p: float
    k_1p: float
    k_2p: float

    
    # Optional parameters depending on specification
    gamma_R: float # Flow straightening coefficient
    asp: float