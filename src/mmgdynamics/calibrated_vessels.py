# Dicts of tested and measured vessels. 
# Data from Yoshimura, Yasuo, and Yumiko Masumoto. "Hydrodynamic database 
# and manoeuvring prediction method with medium high-speed merchant ships and fishing vessels." 
# International MARSIM Conference. 2012.
#
# The SEIUN MARU is a Japanese trainig vessel of 105m length and ~18m width.
# Notes:
#   - Mass is calculated from the displacement value of a 1:24.5 scale model of the vessel:
#       - ▽ = 0.3877 m³ (Displacement)
#       - Scale 1/24.5
#       - Mass = ▽ * 24.5³ * rho (water density) = 5_701_000kg
#
#   - L/B Ratio: 5.866
#   - Block coefficient: 0.508
#
#   - Max Propeller Revs: 176 rpm (2.93 rps) (https://doi.org/10.2534/jjasnaoe1968.1991.170_111)
#
seiunmaru = {
    "C_b":          0.508, # Block Coeffiient
    "Lpp":          105.0, # Length over pependiculars (m)
    "B":            17.91, # Overall width
    "m":            5701.565*1000, # Mass of ship as calculated by ▽*rho (displacement * water density)
    "w_P0":         0.204, # Wake fraction coefficient Kijima et. al. (1990)
    "J_int":        0.4, # Intercept for the calculation of K_T (https://doi.org/10.1615/ICHMT.2012.ProcSevIntSympTurbHeatTransfPal.500)
    "J_slo":       -0.4, # Slope for the calculation of K_T
    "x_G":         -2.249, # X-Coordinate of the center of gravity (m)
    "x_P":         -52.5, # X-Coordinate of the propeller (-0.5*Lpp)
    "D_p":          4.7, # Diameter of propeller (m)
    "k_0":          0.4, # Same value as "J_int" | Propeller open water coefficients. 
    "l_R":         -0.947, # correction of flow straightening factor to yaw-rate
    "gamma_R_plus": 0.394, # Flow straightening coefficient for positive rudder angles
    "gamma_R_minus":0.312, # Flow straightening coefficient for negative rudder angles
    "eta":          0.960, # Ratio of propeller diameter to rudder span
    "kappa":        0.676, # An experimental constant for expressing "u_R"
    "A_R":          13.405, # Moveable rudder area
    "epsilon":      1.033, # Ratio of wake fraction at propeller and rudder positions ((1 - w_R) / (1 - w_P))
    "A_R_Ld_em":    1/46.8, # Fraction of moveable Rudder area to length*draft
    "f_alpha":      2.88, # Rudder lift gradient coefficient (assumed rudder aspect ratio = 2)
    "rho":          1000, # Water density of freshwater
    "t_R":          0.341, # Steering resistance deduction factor (Yoshimura and Masumoto (2012))
    "t_P":          0.21, # Thrust deduction factor. Kulzyk (1995)
    "x_H_dash":    -0.4, # Longitudinal coordinate of acting point of the additional lateral force (Yoshimura and Masumoto (2012))
    "d":            5.966, # Ship draft (Tiefgang)
    #"m_x_dash":     0.0087, # Added masses in x direction (m*0.05)(Clarke et. al.(1983))
    "m_x_dash":     0.0167, # Added masses in x direction (m*0.05)(Clarke et. al.(1983))
    "m_y_dash":     0.1576, # Added mass in y direction (m*0.9089) (Zhou et. al.(1983))
    "J_z_dash":     0.009,# Added moment of inertia (m* 0.0537) (Zhou et. al.(1983))
    "R_0_dash":     0.016, # frictional resistance coefficient TODO Estimate this via Schoenherr's formula
    "X_vv_dash":   -0.0804, # Hull derivatives (Yoshimura and Masumoto (2012))
    "X_vr_dash":    0.072, # Hull derivatives (Yoshimura and Masumoto (2012))
    #"X_rr_dash":    0.355, # Hull derivatives (Yoshimura and Masumoto (2012))
    "X_rr_dash":    0.022, # Hull derivatives (Yoshimura and Masumoto (2012))
    "X_vvvv_dash":  0.439, # Hull derivatives (Yoshimura and Masumoto (2012))
    "Y_v_dash":    -0.317, # Hull derivatives (Yoshimura and Masumoto (2012))
    "Y_r_dash":     0.052, # Hull derivatives (Yoshimura and Masumoto (2012))
    "Y_vvv_dash":  -1.604, # Hull derivatives (Yoshimura and Masumoto (2012))
    "Y_vvr_dash":  -0.75, # Hull derivatives (Yoshimura and Masumoto (2012))
    "Y_vrr_dash":  -0.8599, # Hull derivatives (Yoshimura and Masumoto (2012))
    "Y_rrr_dash":  -0.051, # Hull derivatives (Yoshimura and Masumoto (2012))
    "N_v_dash":    -0.1136, # Hull derivatives (Yoshimura and Masumoto (2012))
    "N_r_dash":    -0.0795, # Hull derivatives (Yoshimura and Masumoto (2012))
    "N_vvv_dash":  -0.3095, # Hull derivatives (Yoshimura and Masumoto (2012))
    "N_vvr_dash":  -0.6257, # Hull derivatives (Yoshimura and Masumoto (2012))
    "N_vrr_dash":  -0.3143, # Hull derivatives (Yoshimura and Masumoto (2012))
    "N_rrr_dash":  -0.0502, # Hull derivatives (Yoshimura and Masumoto (2012))
    "I_zG":         3_928_734_632, # Moment of inertia of ship around center of gravity (m*(0.25*Lpp)**2) (Point mass Inertia)
    "a_H":          0.236 # Rudder force increase factor (Quadvlieg (2013))
}

# Großmotorschiff nach BAW Flottenstruktur Niederrhein (Leer)
GMS1 = {
    "C_b":          0.9, # Block Coeffiient
    "Lpp":          110.0, # Length over pependiculars (m)/
    "B":            11.45, # Overall width
    "m":            1848.82*1000, # Mass of ship as calculated by ▽*rho (displacement * water density)
    "w_P0":         0.40, # Wake fraction coefficient BAW(2020)
    "J_int":        0.4, # Intercept for the calculation of K_T (https://doi.org/10.1615/ICHMT.2012.ProcSevIntSympTurbHeatTransfPal.500)
    "J_slo":       -0.4, # Slope for the calculation of K_T
    "x_G":          0.0, # X-Coordinate of the center of gravity (m)
    "x_P":         -55.0, # X-Coordinate of the propeller (-0.5*Lpp)
    "D_p":          1.7, # Diameter of propeller (m)
    "k_0":          0.4, # Same value as "J_int" | Propeller open water coefficients. 
    "l_R":         -0.9, # correction of flow straightening factor to yaw-rate
    "gamma_R_plus": 0.394, # Flow straightening coefficient for positive rudder angles
    "gamma_R_minus":0.312, # Flow straightening coefficient for negative rudder angles
    "gamma_R":      0.333,
    "eta":          0.960, # Ratio of propeller diameter to rudder span
    "kappa":        0.475, # An experimental constant for expressing "u_R"
    "A_R":          4.402, # Moveable rudder area
    "epsilon":      2.467, # Ratio of wake fraction at propeller and rudder positions ((1 - w_R) / (1 - w_P))
    "A_R_Ld_em":    1/46.8, # Fraction of moveable Rudder area to length*draft
    "f_alpha":      2.45, # Rudder lift gradient coefficient (assumed rudder aspect ratio = 2)
    "rho":          1000, # Water density of freshwater
    "t_R":          0.39, # Steering resistance deduction factor (Yoshimura and Masumoto (2012))
    "t_P":          0.2, # Thrust deduction factor. BAW(2020)
    "x_H_dash":    -0.4, # Longitudinal coordinate of acting point of the additional lateral force (Yoshimura and Masumoto (2012))
    "d":            1.631, # Ship draft (Tiefgang)
    "m_x_dash":     0.009, # Added masses in x direction (m*0.05)(Clarke et. al.(1983))
    "m_y_dash":     0.0423, # Added mass in y direction (m*0.2261) (Zhou et. al.(1983))
    "J_z_dash":     1.121e-6,# Added moment of inertia (m* 0.0724) (Zhou et. al.(1983))
    "R_0_dash":     0.026, # frictional resistance coefficient TODO Estimate this via Schoenherr's formula
    "X_vv_dash":   -0.0723, # Hull derivatives (Yoshimura and Masumoto (2012))
    "X_vr_dash":   -0.0566, # Hull derivatives (Yoshimura and Masumoto (2012))
    "X_rr_dash":   -0.116, # Hull derivatives (Yoshimura and Masumoto (2012))
    "X_vvvv_dash":  0.474, # Hull derivatives (Yoshimura and Masumoto (2012))
    "Y_v_dash":    -0.177, # Hull derivatives (Yoshimura and Masumoto (2012))
    "Y_r_dash":     0.056, # Hull derivatives (Yoshimura and Masumoto (2012))
    "Y_vvv_dash":  -1.297, # Hull derivatives (Yoshimura and Masumoto (2012))
    "Y_vvr_dash":  -0.75, # Hull derivatives (Yoshimura and Masumoto (2012))
    "Y_vrr_dash":  -0.3598, # Hull derivatives (Yoshimura and Masumoto (2012))
    "Y_rrr_dash":  -0.051, # Hull derivatives (Yoshimura and Masumoto (2012))
    "N_v_dash":    -0.029, # Hull derivatives (Yoshimura and Masumoto (2012))
    "N_r_dash":    -0.0151, # Hull derivatives (Yoshimura and Masumoto (2012))
    "N_vvv_dash":  -0.039, # Hull derivatives (Yoshimura and Masumoto (2012))
    "N_vvr_dash":  -0.614, # Hull derivatives (Yoshimura and Masumoto (2012))
    "N_vrr_dash":   0.014, # Hull derivatives (Yoshimura and Masumoto (2012))
    "N_rrr_dash":  -0.0326, # Hull derivatives (Yoshimura and Masumoto (2012))
    "I_zG":         1_398_170_125.0,# 5_592_680_500, # Moment of inertia of ship around center of gravity (m*(0.25*Lpp)**2) (Point mass Inertia)
    "a_H":          0.4113 # Rudder force increase factor (Quadvlieg (2013))
    
}

# L7 model of the KVLCC2 Tanker ship. Originally 320m length, 58m width.
# Scale of model 1/5
kvlcc2 = {
    "C_b":          0.810, # Block Coeffiient
    "Lpp":          64, # Length over pependiculars (m)
    "B":            11.6, # Overall width
    "m":            2500.8*1000, # Mass of ship as calculated by ▽*rho (displacement * water density)
    "w_P0":         0.35, # Assumed wake fraction coefficient
    "J_int":        0.4, # Intercept for the calculation of K_T (https://doi.org/10.1615/ICHMT.2012.ProcSevIntSympTurbHeatTransfPal.500)
    "J_slo":       -0.5, # Slope for the calculation of K_T
    "x_G":          2.24, # X-Coordinate of the center of gravity (m)
    "x_P":         -32.0, # X-Coordinate of the propeller (-0.5*Lpp)
    "D_p":          1.972, # Diameter of propeller (m)
    "k_0":          0.2931, # Same value as "J_int" | Propeller open water coefficients. 
    "k_1":         -0.2753,
    "k_2":         -0.1359,
    "C_1":          2.0,
    "C_2_plus":     1.6,
    "C_2_minus":    1.1,
    "l_R":         -0.710, # correction of flow straightening factor to yaw-rate
    "gamma_R_plus": 0.640, # Flow straightening coefficient for positive rudder angles
    "gamma_R_minus":0.395, # Flow straightening coefficient for negative rudder angles
    "eta":          0.626, # Ratio of propeller diameter to rudder span
    "kappa":        0.50, # An experimental constant for expressing "u_R"
    "A_R":          4.5, # Moveable rudder area
    "epsilon":      1.09, # Ratio of wake fraction at propeller and rudder positions ((1 - w_R) / (1 - w_P))
    "A_R_Ld_em":    1/46.8, # Fraction of moveable Rudder area to length*draft
    "f_alpha":      2.747, # Rudder lift gradient coefficient (assumed rudder aspect ratio = 2)
    "rho":          1000, # Water density of freshwater
    "t_R":          0.387, # Steering resistance deduction factor
    "t_P":          0.220, # Thrust deduction factor. TODO give this more than an arbitrary value
    "x_H_dash":    -0.464, # Longitudinal coordinate of acting point of the additional lateral force
    "d":            4.16, # Ship draft (Tiefgang)
    "m_x_dash":     0.022, # Non dimensionalized added masses coefficient in x direction
    "m_y_dash":     0.223, # Non dimensionalized added masses coefficient in y direction
    "R_0_dash":     0.022, # frictional resistance coefficient TODO Estimate this via Schoenherr's formula
    "X_vv_dash":   -0.040, # Hull derivatives
    "X_vr_dash":    0.002, # Hull derivatives
    "X_rr_dash":    0.011, # Hull derivatives
    "X_vvvv_dash":  0.771, # Hull derivatives
    "Y_v_dash":    -0.315, # Hull derivatives
    "Y_r_dash":     0.083, # Hull derivatives
    "Y_vvv_dash":  -1.607, # Hull derivatives
    "Y_vvr_dash":   0.379, # Hull derivatives
    "Y_vrr_dash":  -0.391, # Hull derivatives
    "Y_rrr_dash":   0.008, # Hull derivatives
    "N_v_dash":    -0.137, # Hull derivatives
    "N_r_dash":    -0.049, # Hull derivatives
    "N_vvv_dash":  -0.030, # Hull derivatives
    "N_vvr_dash":  -0.294, # Hull derivatives
    "N_vrr_dash":   0.055, # Hull derivatives
    "N_rrr_dash":  -0.013, # Hull derivatives
    "I_zG":         640204800, # Moment of inertia of ship around center of gravity (m*(0.25*Lpp)**2) (Point mass Inertia)
    "J_z_dash":     0.011, # Added moment of inertia coefficient
    "a_H":          0.312 # Rudder force increase factor
}

kvlcc2_full = {
    "C_b":          0.810, # Block Coeffiient
    "Lpp":          320.0, # Length over pependiculars (m)
    "B":            58., # Overall width
    "m":            312_600*1000, # Mass of ship as calculated by ▽*rho (displacement * water density)
    "w_P0":         0.35, # Assumed wake fraction coefficient
    "J_int":        0.4, # Intercept for the calculation of K_T (https://doi.org/10.1615/ICHMT.2012.ProcSevIntSympTurbHeatTransfPal.500)
    "J_slo":       -0.5, # Slope for the calculation of K_T
    "x_G":          11.2, # X-Coordinate of the center of gravity (m)
    "x_P":         -160.0, # X-Coordinate of the propeller (-0.5*Lpp)
    "D_p":          9.86, # Diameter of propeller (m)
    "k_0":          0.2931, # Same value as "J_int" | Propeller open water coefficients. 
    "k_1":         -0.2753,
    "k_2":         -0.1359,    
    "C_1":          2.0,
    "C_2_plus":     1.6,
    "C_2_minus":    1.1,
    "l_R":         -0.710, # correction of flow straightening factor to yaw-rate
    "gamma_R_plus": 0.640, # Flow straightening coefficient for positive rudder angles
    "gamma_R_minus":0.395, # Flow straightening coefficient for negative rudder angles
    "eta":          0.626, # Ratio of propeller diameter to rudder span
    "kappa":        0.50, # An experimental constant for expressing "u_R"
    "A_R":          112.5, # Moveable rudder area
    "epsilon":      1.09, # Ratio of wake fraction at propeller and rudder positions ((1 - w_R) / (1 - w_P))
    "A_R_Ld_em":    1/46.8, # Fraction of moveable Rudder area to length*draft
    "f_alpha":      2.747, # Rudder lift gradient coefficient (assumed rudder aspect ratio = 2)
    "rho":          1000, # Water density of freshwater
    "t_R":          0.387, # Steering resistance deduction factor
    "t_P":          0.220, # Thrust deduction factor. TODO give this more than an arbitrary value
    "x_H_dash":    -0.464, # Longitudinal coordinate of acting point of the additional lateral force
    "d":            20.8, # Ship draft (Tiefgang)
    "m_x_dash":     0.022, # Non dimensionalized added masses coefficient in x direction
    "m_y_dash":     0.223, # Non dimensionalized added masses coefficient in y direction
    "R_0_dash":     0.022, # frictional resistance coefficient TODO Estimate this via Schoenherr's formula
    "X_vv_dash":   -0.040, # Hull derivatives
    "X_vr_dash":    0.002, # Hull derivatives
    "X_rr_dash":    0.011, # Hull derivatives
    "X_vvvv_dash":  0.771, # Hull derivatives
    "Y_v_dash":    -0.315, # Hull derivatives
    "Y_r_dash":     0.083, # Hull derivatives
    "Y_vvv_dash":  -1.607, # Hull derivatives
    "Y_vvr_dash":   0.379, # Hull derivatives
    "Y_vrr_dash":  -0.391, # Hull derivatives
    "Y_rrr_dash":   0.008, # Hull derivatives
    "N_v_dash":    -0.137, # Hull derivatives
    "N_r_dash":    -0.049, # Hull derivatives
    "N_vvv_dash":  -0.030, # Hull derivatives
    "N_vvr_dash":  -0.294, # Hull derivatives
    "N_vrr_dash":   0.055, # Hull derivatives
    "N_rrr_dash":  -0.013, # Hull derivatives
    "I_zG":         2e12, # Moment of inertia of ship around center of gravity (m*(0.25*Lpp)**2) (Point mass Inertia)
    "J_z_dash":     0.011, # Added moment of inertia coefficient
    "a_H":          0.312 # Rudder force increase factor
}