from dataclasses import dataclass
from typing import Optional, TypeVar

import numpy as np


@dataclass
class Vessel:
    
    rho: float # Water density
    rho_air: float # Air density
    
    # Vessel particulars
    C_b: float # Block Coefficient
    Lpp: float # Length over pependiculars (m)/
    B: float # Overall width
    d: float # Ship draft (Tiefgang)
    w_P0: float # Wake fraction coefficient
    x_G: float # X-Coordinate of the center of gravity (m)
    x_P: float # X-Coordinate of the propeller (-0.5*Lpp)
    D_p: float # Diameter of propeller (m)
    l_R: float # correction of flow straightening factor to yaw-rate
    eta: float # Ratio of propeller diameter to rudder span
    kappa: float # An experimental constant for expressing "u_R"
    A_R: float # Moveable rudder area
    epsilon: float # Ratio of wake fraction at propeller and rudder positions ((1 - w_R) / (1 - w_P))
    t_R: float # Steering resistance deduction factor
    t_P: float # Thrust deduction factor
    x_H_dash: float # Longitudinal coordinate of acting point of the additional lateral force
    a_H: float # Rudder force increase factor

    # Experimental wind projected areas+
    A_Fw: float # Frontal projected area of the wind
    A_Lw: float # Lateral projected area of the wind
    
    # MMG hydrodynamic derivatives
    R_0_dash: float # frictional resistance coefficient
    
    # Hull derivatives for longitudinal forces
    X_vv_dash: float # Hull derivatives
    X_vr_dash: float # Hull derivatives
    X_rr_dash: float # Hull derivatives
    X_vvvv_dash: float # Hull derivatives
    
    # Hull derivatives for lateral forces
    Y_v_dash: float # Hull derivatives
    Y_r_dash: float # Hull derivatives
    Y_vvv_dash: float # Hull derivatives
    Y_vvr_dash: float # Hull derivatives
    Y_vrr_dash: float # Hull derivatives
    Y_rrr_dash: float # Hull derivatives
    
    # Hull derivatives for yaw moment
    N_v_dash: float # Hull derivatives
    N_r_dash: float # Hull derivatives
    N_vvv_dash: float # Hull derivatives
    N_vvr_dash: float # Hull derivatives
    N_vrr_dash: float # Hull derivatives
    N_rrr_dash: float # Hull derivatives
    
    # Masses and added masses
    displ: float # Displacement in [m³]
    m_x_dash: float # Non dimensionalized added masses coefficient in x direction
    m_y_dash: float # Non dimensionalized added masses coefficient in y direction
    
    # Moment of inertia and added moment of inertia
    J_z_dash: float # Added moment of inertia coefficient
    
    # Wake change coefficients and propeller advance ratio polynomial
    k_0: Optional[float] = None
    k_1: Optional[float] = None
    k_2: Optional[float] = None
    C_1: Optional[float] = None
    C_2_plus: Optional[float] = None
    C_2_minus: Optional[float] = None
    J_slo: Optional[float] = None
    J_int: Optional[float] = None
    
    # Optional parameters depending on specification
    gamma_R_plus: Optional[float] = None # Flow straightening coefficient for positive rudder angles
    gamma_R_minus: Optional[float] = None # Flow straightening coefficient for negative rudder angles
    gamma_R: Optional[float] = None # Flow straightening coefficient
    A_R_Ld_em: Optional[float] = None # Fraction of moveable Rudder area to length*draft
    f_alpha: Optional[float] = None # Rudder lift gradient coefficient

@dataclass
class MinimalVessel:
        
    # Vessel particulars
    m: float # Displacent of ship
    C_b: float # Block Coefficient
    Lpp: float # Length over pependiculars (m)/
    B: float # Overall width
    d: float # Ship draft (Tiefgang)
    eta: float # Ratio of propeller diameter to rudder span
    A_R: float # Rudder Area
    D_p: float # Propeller diameter
    f_alpha: Optional[float] = None # Rudder lift gradient coefficient
    x_G: Optional[float] = None # X-Coordinate of the center of gravity (m)
    w_P0: Optional[float] = None # Wake fraction coefficient
    t_P: Optional[float] = None # Thrust deduction factor


Surge = TypeVar("Surge",bound=float)
Sway = TypeVar("Sway",bound=float)
YawRate = TypeVar("YawRate",bound=float)
RudderAngle = TypeVar("RudderAngle",bound=float)
RevPerSecond = TypeVar("RevPerSecond",bound=float)

@dataclass
class InitialValues:
    u: Surge
    v: Sway
    r: YawRate
    delta: RudderAngle
    nps: RevPerSecond