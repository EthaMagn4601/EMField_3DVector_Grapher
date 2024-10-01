import numpy as np
import sympy as sp

RHO = "\u03C1"
THETA = "\u03b8"
PHI = "\u03C6"

rhohat = RHO + u"\u0302" 
thetahat = THETA + u"\u0302"
phihat = PHI + u"\u0302"
xhat = "x" + u"\u0302"
yhat = "y" + u"\u0302"
zhat = "z" + u"\u0302"
rhat = "r" + u"\u0302"


def vector_cart_to_cart(x_input, y_input, z_input):

    x_output = str(x_input)
    y_output = str(y_input)
    z_output = str(z_input)

    return {"xhat": x_output, "yhat": y_output, "zhat": z_output}

def vector_cart_to_cyli(x_input, y_input, z_input):
    
    rho_output = f"{str(x_input)}*cos({PHI}) + {str(y_input)}*sin({PHI})"
    phi_output = f"{str(x_input)}*-sin({PHI}) + {str(y_input)}*cos({PHI})"
    z_output = f"{str(z_input)}"
    
    return {"rhohat": rho_output, "phihat": phi_output, "zhat": z_output}

def vector_cart_to_spher(x_input, y_input, z_input):

    r_output = f"{str(x_input)}*sin({THETA})*cos({PHI}) + {str(y_input)}*sin({THETA})*sin({PHI}) + {str(z_input)}*cos({PHI})"
    theta_output = f"{str(x_input)}*cos({THETA})*cos({PHI}) + {str(y_input)}*cos({THETA})*sin({PHI}) + {str(z_input)}*-sin({THETA})"
    phi_output = f"{str(x_input)}*-sin({PHI}) + {str(y_input)}*cos({PHI})"

    return {"rhat": r_output, "thetahat": theta_output, "phihat": phi_output}

def vector_cyli_to_cart(rho_input, phi_input, z_input):

    x_output = f"{str(rho_input)}*cos({PHI}) + {str(phi_input)}*-sin({PHI})"
    y_output = f"{str(rho_input)}*sin({PHI}) + {str(phi_input)}*cos({PHI})"
    z_output = f"{str(z_input)}"

    return {"xhat": x_output, "yhat": y_output, "zhat": z_output}

def vector_cyli_to_cyli(rho_input, phi_input, z_input):

    rho_output = str(rho_input)
    phi_output = str(phi_input)
    z_output = str(z_input)

    return {"rhohat": rho_output, "phihat": phi_output, "zhat": z_output}

def vector_cyli_to_spher(rho_input, phi_input, z_input):

    r_output = f"{str(rho_input)}*sin({THETA}) + {str(z_input)}*cos({THETA})"
    theta_output = f"{str(rho_input)}*cos({THETA}) + {str(z_input)}*-sin({THETA})"
    phi_output = f"{str(phi_input)}"

    return {"rhat": r_output, "thetahat": theta_output, "phihat": phi_output}

def vector_spher_to_cart(r_input, theta_input, phi_input):

    x_output = f"{str(r_input)}*sin({THETA})*cos({PHI}) + {str(theta_input)}*cos({THETA})*cos{PHI} + {str(phi_input)}*-sin({PHI})"
    y_output = f"{str(r_input)}*sin({THETA})*sin({PHI}) + {str(theta_input)}*cos({THETA})*sin{PHI} + {str(phi_input)}*cos({PHI})"
    z_output = f"{str(r_input)}*cos({THETA}) + {str(theta_input)}*-sin({THETA})"

    return {"xhat": x_output, "yhat": y_output, "zhat": z_output}

def vector_spher_to_cyli(r_input, theta_input, phi_input):

    rho_output = f"{str(r_input)}*sin({THETA}) + {str(theta_input)}*cos({THETA})"
    phi_output = str(phi_input)
    z_output = f"{str(r_input)}*cos({THETA}) + {str(theta_input)}*-sin({THETA})"

    return {"rhohat": rho_output, "phihat": phi_output, "zhat": z_output}

def vector_spher_to_spher(r_input, theta_input, phi_input):

    r_output = str(r_input)
    theta_output = str(theta_input)
    phi_output = str(phi_input)

    return {"rhat": r_output, "thetahat": theta_output, "phihat": phi_output}

# equationssolver() solves any equation defined in it, cart_equation_solver() exists to specifically solve for the vectors within the cartesian plane
def equationssolver(equation_to_solve):

    # Declaration of values and count 
    values = np.zeros(((10+1)**3, 7))
    count = 0

    # This loop gives all the possible values of x, y, and z for the values matrix
    for i in range(0, 11):
        for ii in range(0, 11):
            for iii in range (0, 11):
                values[count][0] = i
                values[count][1] = ii
                values[count][2] = iii
                if count < (10+1)**3:
                        count = count + 1

    # This loop gives all the possible values of rho, phi, theta, and r for the values matrix
    for i in range(0, 11**3):
        rho_val = 0
        phi_val = 0
        theta_val = 0
        r_val = 0
        x_val = values[i][0]
        y_val = values[i][1]
        z_val = values[i][2]
        rho_val = np.sqrt((x_val ** 2) + (y_val**2))
        if x_val == 0:
            phi_val = 90 # arcos(0) = 90
        else:
            phi_val = np.rad2deg(np.arctan(y_val/x_val))
        if z_val == 0:
            theta_val = 90 # arcos(0) = 90
        else:
            theta_val = np.rad2deg(np.arctan(rho_val/z_val))
        r_val = np.sqrt((rho_val**2) + (z_val**2))
        values[i][3] = rho_val
        values[i][4] = phi_val
        values[i][5] = theta_val
        values[i][6] = r_val
    
        
    # Variable declaration so sympy knows what to sub in values for
    x, y, z, rho_, phi_, theta_, r = sp.symbols(f'x y z {RHO} {PHI} {THETA} r')

    # Declaration of equation to evaluate
    equation1_input = equation_to_solve["rhohat"]
    equation1 = sp.sympify(equation1_input)
    equation2_input = equation_to_solve["phihat"]
    equation2 = sp.sympify(equation2_input)
    equation3_input = equation_to_solve["zhat"]
    equation3 = sp.sympify(equation3_input)

    result1 = np.zeros(11**3)
    result2 = np.zeros(11**3)
    result3 = np.zeros(11**3)

    for i in range (0, 11**3):
        values_dict = {x: values[i][0], y: values[i][1], z: values[i][2], rho_: values[i][3], phi_: values[i][4], theta_: values[i][5], r: values[i][6]}

        result1[i] = str(equation1.subs(values_dict))
        result2[i] = str(equation2.subs(values_dict))
        result3[i] = str(equation3.subs(values_dict))

    return {"r1": result1, "r2": result2, "r3": result3}

def cart_equation_solver(cart_equation, clutter_size):
    
    m_size = clutter_size

    # Declaration of values and count 
    values = np.zeros(((m_size)**3, 7))
    count = 0

    # This loop gives all the possible values of x, y, and z for the values matrix
    for i in range(0, m_size):
        for ii in range(0, m_size):
            for iii in range (0, m_size):
                values[count][0] = i
                values[count][1] = ii
                values[count][2] = iii
                if count < (m_size)**3:
                        count = count + 1
                
                # Debug code
                # if count == 100:
                #     for i in range (0,20):
                #         print(values[i])

    # This loop gives all the possible values of rho, phi, theta, and r for the values matrix
    for i in range(0, m_size**3):
        rho_val = 0
        phi_val = 0
        theta_val = 0
        r_val = 0
        x_val = values[i][0]
        y_val = values[i][1]
        z_val = values[i][2]
        rho_val = np.sqrt((x_val ** 2) + (y_val**2))
        r_val = np.sqrt((rho_val**2) + (z_val**2))
        if x_val == 0:
            if y_val == 0:
                phi_val = 0
            elif y_val > 0:
                phi_val = 90
            else:
                phi_val = 270
        elif y_val == 0:
            if x_val > 0:
                phi_val = 0
            else:
                phi_val = 180
        else:
            if y_val > 0:
                phi_val = np.rad2deg(np.arccos(x_val/rho_val))
            else:
                phi_val = 360 - np.rad2deg(np.arccos(x_val/rho_val))
        if z_val == 0:
            if r_val == 0:
                theta_val = 0
            else:
                theta_val = 90
        else:
                theta_val = np.rad2deg(np.arccos(z_val/r_val))
        values[i][3] = rho_val
        values[i][4] = phi_val
        values[i][5] = theta_val
        values[i][6] = r_val

    # Variable declaration so sympy knows what to sub in values for
    x, y, z, rho_, phi_, theta_, r = sp.symbols(f'x y z {RHO} {PHI} {THETA} r')

    # Declaration of equation to evaluate
    equation1_input = cart_equation["xhat"]
    equation1 = sp.sympify(equation1_input)
    equation2_input = cart_equation["yhat"]
    equation2 = sp.sympify(equation2_input)
    equation3_input = cart_equation["zhat"]
    equation3 = sp.sympify(equation3_input)

    u_dir = np.zeros((m_size,m_size,m_size))
    v_dir = np.zeros((m_size,m_size,m_size))
    w_dir = np.zeros((m_size,m_size,m_size))

    ii = 0
    iii = 0
    iv = 0

    for i in range (0, m_size**3):
        
        values_dict = {x: values[i][0], y: values[i][1], z: values[i][2], rho_: values[i][3], phi_: values[i][4], theta_: values[i][5], r: values[i][6]}

        if ii == m_size:
            ii = 0
            iii = iii + 1
        if iii == m_size:
            iii = 0
            iv = iv + 1
        
        u_dir[iv, iii, ii] = equation1.subs(values_dict)
        v_dir[iv, iii, ii] = equation2.subs(values_dict)
        w_dir[iv, iii, ii] = equation3.subs(values_dict)
        
        ii = ii + 1

    return {"xhat_out": u_dir, "yhat_out": v_dir, "zhat_out": w_dir}