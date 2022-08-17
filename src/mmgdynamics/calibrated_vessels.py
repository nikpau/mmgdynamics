# L7 model of the KVLCC2 Tanker ship. Originally 320m length, 58m width.
# Scale of model 1/5
kvlcc2 = {
    "C_b":          0.810, # Block Coeffiient
    "Lpp":          64, # Length over pependiculars (m)
    "B":            11.6, # Overall width
    "displ":        2500.8, # Displacement in [m³]
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
    "f_alpha":      2.747, # Rudder lift gradient coefficient (assumed rudder aspect ratio = 2)
    "rho":          1000, # Water density of freshwater
    "rho_air":      1.225, # Air density
    "A_Fw":         240,   # Frontal wind area [m²]
    "A_Lw":         720,   # Lateral wind area [m²]
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
    "J_z_dash":     0.011, # Added moment of inertia coefficient
    "a_H":          0.312 # Rudder force increase factor
}

kvlcc2_full = {
    "C_b":          0.810, # Block Coeffiient
    "Lpp":          320.0, # Length over pependiculars (m)
    "B":            58.0, # Overall width
    "displ":        312_600, # Displacement in [m³]
    "w_P0":         0.35, # Assumed wake fraction coefficient
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
    "f_alpha":      2.747, # Rudder lift gradient coefficient (assumed rudder aspect ratio = 2)
    "rho":          1030, # Water density of freshwater
    "rho_air":      1.225, # Air density
    "A_Fw":         1200,  # Frontal wind area [m²]
    "A_Lw":         3600,  # Lateral wind area [m²]
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
    "J_z_dash":     0.011, # Added moment of inertia coefficient
    "a_H":          0.312 # Rudder force increase factor
}

kvlcc2_l7 = {
    "C_b":          0.810, # Block Coeffiient
    "Lpp":          7.0, # Length over pependiculars (m)
    "B":            1.27, # Overall width
    "displ":        3.2724, # Displacement in [m³]
    "w_P0":         0.40, # Assumed wake fraction coefficient
    "J_int":        0.4, # Intercept for the calculation of K_T (https://doi.org/10.1615/ICHMT.2012.ProcSevIntSympTurbHeatTransfPal.500)
    "J_slo":       -0.5, # Slope for the calculation of K_T
    "x_G":          0.244, # X-Coordinate of the center of gravity (m)
    "x_P":         -3.36, # X-Coordinate of the propeller (-0.5*Lpp)
    "D_p":          0.204, # Diameter of propeller (m)
    "k_0":          0.2931, # Same value as "J_int" | Propeller open water coefficients. 
    "k_1":         -0.2753,
    "k_2":         -0.1385,    
    "C_1":          2.0,
    "C_2_plus":     1.6,
    "C_2_minus":    1.1,
    "l_R":         -0.710, # correction of flow straightening factor to yaw-rate
    "gamma_R_plus": 0.640, # Flow straightening coefficient for positive rudder angles
    "gamma_R_minus":0.395, # Flow straightening coefficient for negative rudder angles
    "eta":          0.626, # Ratio of propeller diameter to rudder span
    "kappa":        0.50, # An experimental constant for expressing "u_R"
    "A_R":          0.0654, # Moveable rudder area
    "epsilon":      1.09, # Ratio of wake fraction at propeller and rudder positions ((1 - w_R) / (1 - w_P))
    "A_R_Ld_em":    1/46.8, # Fraction of moveable Rudder area to length*draft
    "f_alpha":      2.747, # Rudder lift gradient coefficient (assumed rudder aspect ratio = 2)
    "rho":          1030, # Water density of freshwater
    "rho_air":      1.225, # Air density
    "t_R":          0.387, # Steering resistance deduction factor
    "t_P":          0.220, # Thrust deduction factor. TODO give this more than an arbitrary value
    "x_H_dash":    -0.464, # Longitudinal coordinate of acting point of the additional lateral force
    "d":            0.455, # Ship draft (Tiefgang)
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
    "J_z_dash":     0.011, # Added moment of inertia coefficient
    "a_H":          0.312 # Rudder force increase factor
}