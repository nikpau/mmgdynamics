# Fig 1 in https://doi.org/10.1016/S0141-1187(98)00002-9
# Points have been digitalized using an online graph digitizer

cross_flow = {
    "x": [0.03968253968254,0.138888888888889,
          0.21031746031746,0.253968253968254,
          0.321428571428571,0.253968253968254,
          0.416666666666667,0.373015873015873,
          0.448412698412698,0.452380952380952,
          0.452380952380952,0.46031746031746,
          0.464285714285714,0.464285714285714,
          0.464285714285714,0.468253968253968,
          0.468253968253968,0.472222222222222,
          0.456349206349206,0.48015873015873,
          0.51984126984127,0.563492063492063,
          0.607142857142857,0.638888888888889,
          0.734126984126984,0.869047619047619,
          0.932539682539682,1.36507936507937,
          1.83333333333333,2.29761904761905,
          2.67460317460317,3.65079365079365],
    "y": [1.96177370030581,1.96177370030581,
          1.95259938837921,1.93425076452599,
          1.90672782874618,1.93425076452599,
          1.82721712538226,1.86697247706422,
          1.79357798165138,1.74770642201835,
          1.65290519877676,1.58562691131498,
          1.53363914373089,1.47553516819572,
          1.42354740061162,1.35932721712538,
          1.30428134556575,1.25535168195719,
          1.70489296636086,1.20948012232416,
          1.1605504587156,1.10550458715596,
          1.05351681957187,1.00764525993884,
          0.958715596330275,0.879204892966361,
          0.839449541284404,0.753822629969419,
          0.655963302752294,0.628440366972477,
          0.603975535168196,0.570336391437309]
}

def interpolateY(x: float, vals: dict) -> float:
    """Interpolate y values from a discrete function table

    Args:
        x (float): input value for which an interpolation shall be returned
        vals (dict): Dictionary with discrete function values of form:
                    {"x": [x1,...,xn],"y": [y1,...,yn]}
    
    Raises:
        RuntimeError: Raised if input value is out of interpolation range

    Returns:
        float: Interpolated y value
    """
    
    # Get x and y vals from the function to be interpolated
    func_x, func_y = vals["x"], vals["y"]
    
    # Get index of closest value to input
    for idx, val in enumerate(func_x):
        if  val > x:
            func_ind = idx-1
            break
    
    if "func_ind" not in locals():
        raise RuntimeError(f"Value out of interpolation range. Choose a value between {min(func_x)} and {max(func_x)}")
    
    # Define the two points to interpolate in between
    xa, xb = func_x[func_ind+1], func_x[func_ind]
    ya, yb = func_y[func_ind+1], func_y[func_ind]
    
    # Linear equation
    m = (ya - yb) / (xa - xb)
    yc = (x - xb) * m + yb
    
    return yc