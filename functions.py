import numpy as np
import sympy as sp
import tkinter as tk
from tkinter import ttk
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

RHO = "\u03C1"
THETA = "\u03b8"
PHI = "\u03C6"
# PI = "\u03C0"

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
    
    rho_output = f"({str(x_input)})*cos({PHI}) + ({str(y_input)})*sin({PHI})"
    phi_output = f"({str(x_input)})*-sin({PHI}) + ({str(y_input)})*cos({PHI})"
    z_output = f"({str(z_input)})"
    
    return {"rhohat": rho_output, "phihat": phi_output, "zhat": z_output}

def vector_cart_to_spher(x_input, y_input, z_input):

    r_output = f"({str(x_input)})*sin({THETA})*cos({PHI}) + ({str(y_input)})*sin({THETA})*sin({PHI}) + ({str(z_input)})*cos({THETA})"
    theta_output = f"({str(x_input)})*cos({THETA})*cos({PHI}) + ({str(y_input)})*cos({THETA})*sin({PHI}) + ({str(z_input)})*-sin({THETA})"
    phi_output = f"({str(x_input)})*-sin({PHI}) + ({str(y_input)})*cos({PHI})"

    return {"rhat": r_output, "thetahat": theta_output, "phihat": phi_output}

def vector_cyli_to_cart(rho_input, phi_input, z_input):

    x_output = f"({str(rho_input)})*cos({PHI}) + ({str(phi_input)})*-sin({PHI})"
    y_output = f"({str(rho_input)})*sin({PHI}) + ({str(phi_input)})*cos({PHI})"
    z_output = f"({str(z_input)})"

    return {"xhat": x_output, "yhat": y_output, "zhat": z_output}

def vector_cyli_to_cyli(rho_input, phi_input, z_input):

    rho_output = str(rho_input)
    phi_output = str(phi_input)
    z_output = str(z_input)

    return {"rhohat": rho_output, "phihat": phi_output, "zhat": z_output}

def vector_cyli_to_spher(rho_input, phi_input, z_input):

    r_output = f"({str(rho_input)})*sin({THETA}) + ({str(z_input)})*cos({THETA})"
    theta_output = f"({str(rho_input)})*cos({THETA}) + ({str(z_input)})*-sin({THETA})"
    phi_output = f"({str(phi_input)})"

    return {"rhat": r_output, "thetahat": theta_output, "phihat": phi_output}

def vector_spher_to_cart(r_input, theta_input, phi_input):

    x_output = f"({str(r_input)})*sin({THETA})*cos({PHI}) + ({str(theta_input)})*cos({THETA})*cos({PHI}) + ({str(phi_input)})*-sin({PHI})"
    y_output = f"({str(r_input)})*sin({THETA})*sin({PHI}) + ({str(theta_input)})*cos({THETA})*sin({PHI}) + ({str(phi_input)})*cos({PHI})"
    z_output = f"({str(r_input)})*cos({THETA}) + ({str(theta_input)})*-sin({THETA})"

    return {"xhat": x_output, "yhat": y_output, "zhat": z_output}

def vector_spher_to_cyli(r_input, theta_input, phi_input):

    rho_output = f"({str(r_input)})*sin({THETA}) + ({str(theta_input)})*cos({THETA})"
    phi_output = str(phi_input)
    z_output = f"({str(r_input)})*cos({THETA}) + ({str(theta_input)})*-sin({THETA})"

    return {"rhohat": rho_output, "phihat": phi_output, "zhat": z_output}

def vector_spher_to_spher(r_input, theta_input, phi_input):

    r_output = str(r_input)
    theta_output = str(theta_input)
    phi_output = str(phi_input)

    return {"rhat": r_output, "thetahat": theta_output, "phihat": phi_output}

def cart_equation_solver(cart_equation, x_inc_amt, y_inc_amt, z_inc_amt, grid_size):


    u = np.arange(-1*grid_size, grid_size + 1, x_inc_amt)
    v = np.arange(-1*grid_size, grid_size + 1, y_inc_amt)
    w = np.arange(-1*grid_size, grid_size + 1, z_inc_amt)
    u_len = len(u)
    v_len = len(v)
    w_len = len(w)

    # Declaration of values and count 
    values = np.zeros((u_len*v_len*w_len, 7))
    count = 0

    # This loop gives all the possible values of x, y, z, rho, phi, theta, and r for the values matrix
    for i in range(0, u_len):
        for ii in range(0, v_len):
            for iii in range(0, w_len):
                rho_val = 0
                phi_val = 0
                theta_val = 0
                r_val = 0
                x_val = -1*grid_size + i*x_inc_amt
                y_val = -1*grid_size + ii*y_inc_amt
                z_val = -1*grid_size + iii*z_inc_amt
                rho_val = np.sqrt((x_val ** 2) + (y_val**2))
                r_val = np.sqrt((rho_val**2) + (z_val**2))
                if x_val == 0:
                    if y_val == 0:
                        phi_val = np.deg2rad(0)
                    elif y_val > 0:
                        phi_val = np.deg2rad(90)
                    else:
                        phi_val = np.deg2rad(270)
                elif x_val < 0:
                    if y_val == 0:
                        phi_val = np.deg2rad(180)
                    elif y_val < 0:
                        phi_val = np.deg2rad(360 - np.rad2deg(np.arccos(x_val/rho_val)))
                    else:
                        phi_val = np.arccos(x_val/rho_val)
                else:
                    if y_val == 0:
                        phi_val = np.deg2rad(0)
                    elif y_val < 0:
                        phi_val = np.deg2rad(360 - np.rad2deg(np.arccos(x_val/rho_val)))
                    else:
                        phi_val = np.arccos(x_val/rho_val)
                if z_val == 0:
                    if r_val == 0:
                        theta_val = np.deg2rad(0)
                    else:
                        theta_val = np.deg2rad(90)
                else:
                        theta_val = np.arccos(z_val/r_val)
                values[count][0] = x_val
                values[count][1] = y_val
                values[count][2] = z_val
                values[count][3] = rho_val
                values[count][4] = phi_val
                values[count][5] = theta_val
                values[count][6] = r_val
                count = count + 1
                
    # Variable declaration so sympy knows what to sub in values for
    x, y, z, rho_, phi_, theta_, r_ = sp.symbols(f'x y z {RHO} {PHI} {THETA} r')

    # Declaration of equations to evaluate
    equation1_input = cart_equation["xhat"]
    equation1 = sp.sympify(equation1_input)
    equation2_input = cart_equation["yhat"]
    equation2 = sp.sympify(equation2_input)
    equation3_input = cart_equation["zhat"]
    equation3 = sp.sympify(equation3_input)

    u_dir = np.zeros((u_len,v_len,w_len))
    v_dir = np.zeros((u_len,v_len,w_len))
    w_dir = np.zeros((u_len,v_len,w_len))

    ii = 0
    iii = 0
    iv = 0

    for i in range (0, u_len*v_len*w_len):
        
        values_dict = {x: values[i][0], y: values[i][1], z: values[i][2], rho_: values[i][3], phi_: values[i][4], theta_: values[i][5], r_: values[i][6]}

        if ii == w_len:
            ii = 0
            iii = iii + 1
        if iii == v_len:
            iii = 0
            iv = iv + 1

        u_dir[iv, iii, ii] = equation1.subs(values_dict)
        v_dir[iv, iii, ii] = equation2.subs(values_dict)
        w_dir[iv, iii, ii] = equation3.subs(values_dict)
        
        ii = ii + 1

    return {"xhat_out": u_dir, "yhat_out": v_dir, "zhat_out": w_dir}

switch_var = 1

def main(spacerpady, font_size, text_font):
    root = tk.Tk()
    root.title("Cartesian 3D Vector Field Plotter")
    root.geometry("1500x900")

    fig = plt.Figure(figsize=(7,5), dpi=100)
    ax = fig.add_subplot(1,1,1,projection="3d")

    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().pack(side="right",fill="both", expand=True)

    entry_length = 6

    def on_closing():
        plt.close()
        root.destroy()

    def root_gui():

        def vector_plot(plot_type):
            
            def initialize_plot():
                x_text_field.delete(0, tk.END)
                y_text_field.delete(0, tk.END)
                z_text_field.delete(0, tk.END)

                rho_text_field.delete(0, tk.END)
                phi_c_text_field.delete(0, tk.END)
                zcyl_text_field.delete(0, tk.END)

                r_text_field.delete(0, tk.END)
                theta_text_field.delete(0, tk.END)
                phi_s_text_field.delete(0, tk.END)
                cart_equation = vector_cart_to_cart(x_input, y_input, z_input)
                results = cart_equation_solver(cart_equation, x_inc_amt, y_inc_amt, z_inc_amt, grid_size)
                return results
            
            ax.clear()
            try:
            # initial conditions for the plot to work & also default values in the event nothing is input into the entry fields
                x_input = 0
                y_input = 0
                z_input = 0
                rho_input = 0
                phicyl_input = 0
                zcyl_input = 0
                r_input = 0
                theta_input = 0
                phispher_input = 0
                view_size = 10
                clutter_size = "3, 3, 3"

                view_text = view_size_entry.get()
                clutter_text = clutter_size_entry.get()
                x_text = x_text_field.get()
                y_text = y_text_field.get()
                z_text = z_text_field.get()
                rho_text = rho_text_field.get()
                phicyl_text = phi_c_text_field.get()
                zcyl_text = zcyl_text_field.get()
                r_text = r_text_field.get()
                theta_text = theta_text_field.get()
                phispher_text = phi_s_text_field.get()

                if switch_var == 1:
                    if x_text:
                        x_input = x_text
                    if y_text:
                        y_input = y_text
                    if z_text:
                        z_input = z_text
                elif switch_var == 2:
                    if rho_text:
                        rho_input = rho_text
                    if phicyl_text:
                        phicyl_input = phicyl_text
                    if zcyl_text:
                        zcyl_input = zcyl_text
                elif switch_var == 3:
                    if r_text:
                        r_input = r_text
                    if theta_text:
                        theta_input = theta_text
                    if phispher_text:
                        phispher_input = phispher_text
                if view_text:
                    view_size = float(view_text)
                    view_size_entry.delete(0, tk.END)
                    view_size_entry.insert(tk.END, view_text)
                else:
                    view_size_entry.delete(0, tk.END)
                    view_size_entry.insert(tk.END, view_size)
                if clutter_text:
                    clutter_size = clutter_text
                    clutter_size_entry.delete(0, tk.END)
                    clutter_size_entry.insert(tk.END, clutter_text)
                else:
                    clutter_size_entry.delete(0, tk.END)
                    clutter_size_entry.insert(tk.END, clutter_size)

                grid_size = float(view_size)

                clutter_size = clutter_size.split(",")
                x_clutter = int(clutter_size[0])
                y_clutter = int(clutter_size[1])
                z_clutter = int(clutter_size[2])
                
                x_inc_amt = grid_size/x_clutter
                y_inc_amt = grid_size/y_clutter
                z_inc_amt = grid_size/z_clutter
                    
                x, y, z = np.meshgrid(np.arange(-1*grid_size, grid_size + 1, x_inc_amt),
                                np.arange(-1*grid_size, grid_size + 1, y_inc_amt),
                                np.arange(-1*grid_size, grid_size + 1, z_inc_amt),
                                indexing='ij')
                
                if switch_var == 1:
                    if (x_input != 0) or (y_input != 0) or (z_input != 0):
                        cart_equation = vector_cart_to_cart(x_input, y_input, z_input)
                        cyl_equation = vector_cart_to_cyli(x_input, y_input, z_input)
                        spher_equation = vector_cart_to_spher(x_input, y_input, z_input)
                        results = cart_equation_solver(cart_equation, x_inc_amt, y_inc_amt, z_inc_amt, grid_size)

                        rho_text_field.delete(0, tk.END)
                        rho_text_field.insert(tk.END, cyl_equation["rhohat"])
                        phi_c_text_field.delete(0, tk.END)
                        phi_c_text_field.insert(tk.END, cyl_equation["phihat"])
                        zcyl_text_field.delete(0, tk.END)
                        zcyl_text_field.insert(tk.END, cyl_equation["zhat"])

                        r_text_field.delete(0, tk.END)
                        r_text_field.insert(tk.END, spher_equation["rhat"])
                        theta_text_field.delete(0, tk.END)
                        theta_text_field.insert(tk.END, spher_equation["thetahat"])
                        phi_s_text_field.delete(0, tk.END)
                        phi_s_text_field.insert(tk.END, spher_equation["phihat"])
                    else:    
                        results = initialize_plot()
                
                if switch_var == 2:
                    if (rho_input != 0) or (phicyl_input != 0) or (zcyl_input != 0):
                        cart_equation = vector_cyli_to_cart(rho_input, phicyl_input, zcyl_input)
                        spher_equation = vector_cyli_to_spher(rho_input, phicyl_input, zcyl_input)
                        x_input = cart_equation["xhat"]
                        y_input = cart_equation["yhat"]
                        z_input = cart_equation["zhat"]
                        results = cart_equation_solver(cart_equation, x_inc_amt, y_inc_amt, z_inc_amt, grid_size)

                        x_text_field.delete(0, tk.END)
                        x_text_field.insert(tk.END, cart_equation["xhat"])
                        y_text_field.delete(0, tk.END)
                        y_text_field.insert(tk.END, cart_equation["yhat"])
                        z_text_field.delete(0, tk.END)
                        z_text_field.insert(tk.END, cart_equation["zhat"])

                        r_text_field.delete(0, tk.END)
                        r_text_field.insert(tk.END, spher_equation["rhat"])
                        theta_text_field.delete(0, tk.END)
                        theta_text_field.insert(tk.END, spher_equation["thetahat"])
                        phi_s_text_field.delete(0, tk.END)
                        phi_s_text_field.insert(tk.END, spher_equation["phihat"])
                    else:    
                        results = initialize_plot()
                if switch_var == 3:
                    if (r_input != 0) or (theta_input != 0) or (phispher_input != 0):
                        cart_equation = vector_spher_to_cart(r_input, theta_input, phispher_input)
                        cyl_equation = vector_spher_to_cyli(r_input, theta_input, phispher_input)
                        x_input = cart_equation["xhat"]
                        y_input = cart_equation["yhat"]
                        z_input = cart_equation["zhat"]
                        results = cart_equation_solver(cart_equation, x_inc_amt, y_inc_amt, z_inc_amt, grid_size)

                        x_text_field.delete(0, tk.END)
                        x_text_field.insert(tk.END, cart_equation["xhat"])
                        y_text_field.delete(0, tk.END)
                        y_text_field.insert(tk.END, cart_equation["yhat"])
                        z_text_field.delete(0, tk.END)
                        z_text_field.insert(tk.END, cart_equation["zhat"])

                        rho_text_field.delete(0, tk.END)
                        rho_text_field.insert(tk.END, cyl_equation["rhohat"])
                        phi_c_text_field.delete(0, tk.END)
                        phi_c_text_field.insert(tk.END, cyl_equation["phihat"])
                        zcyl_text_field.delete(0, tk.END)
                        zcyl_text_field.insert(tk.END, cyl_equation["zhat"])
                    else:    
                        results = initialize_plot()

                # Defines magnitude/direction vector arrays
                u = results["xhat_out"]
                v = results["yhat_out"]
                w = results["zhat_out"]

                if plot_type == "Norm":

                    # plot the plot
                    ax.quiver(x, y, z, u, v, w, length=view_size/10, pivot='middle', normalize=True)

                elif plot_type == "Rel_mag":
                    
                    # Sets default value for relative magnitude check button
                    rel_mag_length = 0.1
                    rel_mag_text = rel_mag_entry_field.get()

                    # if statement to check if there is anything actually in the entry field to compute (i.e. not empty)
                    if rel_mag_text:
                        rel_mag_length = float(rel_mag_text)
                    
                    # inputs 0.1 below relative magnitude check button if nothing is input, otherwise will replace 0.1 with user input
                    rel_mag_entry_field.delete(0, tk.END)
                    rel_mag_entry_field.insert(tk.END, rel_mag_length)

                    # plot the plot
                    ax.quiver(x, y, z, u, v, w, pivot='middle', length=rel_mag_length, normalize=False)

                # checks if an input is zero
                elif ((x_input == 0) and (y_input == 0) and (z_input == 0)) and plot_type == "Colors":
                    
                    # plot the plot
                    ax.quiver(x, y, z, u, v, w, length=view_size/10, pivot='middle', normalize=True)

                elif plot_type == "Colors":
                    magnitudes = np.sqrt((u**2) + (v**2) + (w**2))

                    #Error correction for divide by zero case
                    magnitudes = np.nan_to_num(magnitudes)

                    norm = mpl.colors.Normalize(vmin=0, vmax=magnitudes.max())
                    colours = mpl.cm.jet(norm(magnitudes))
                    ax.quiver(x, y, z, u, v, w, colors=colours.reshape(-1, 4), length=view_size/10, pivot="middle", normalize=True)
                    cax = ax.inset_axes([0.2, 0, 0.6, 0.05])
                    plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap="jet"), cax=cax, orientation="horizontal", label="Magnitudes")
                else:
                    tk.messagebox.showinfo("PlotTypeError", "Something went wrong with the plot_type variable declaration")
                
                ax.view_init(elev=30, azim=60)
                
                # ax.set_title()
                ax.set_xlabel("X Axis")
                ax.set_ylabel("Y Axis")
                ax.set_zlabel("Z Axis")
                
                canvas.draw()
            except ZeroDivisionError:
                tk.messagebox.showinfo("DivisionByZero", "Please make sure you've input a whole number value greater than zero into the Clutter Size Entry Field.")

            except IndexError:
                tk.messagebox.showinfo("IndexError", "Please make sure you've input the clutter size as a three integer string separated by commas (i.e. like 3, 4, 5)")

            except TypeError:
                tk.messagebox.showinfo("TypeError", "Please make sure you've input the appropriate varibles into respective Entry Fields (i.e. no f, g, etc.).")
            
            except ValueError:
                tk.messagebox.showinfo("ValueError", "Please make sure you've input a numerical value into the Grid Size & Clutter Size Entry Fields (i.e. no f, g, h, x, a, etc.).")

        def plot_type_check():
            if norm_var.get() == 1:
                vector_plot("Norm")
            elif rel_mag_var.get() == 1:
                vector_plot("Rel_mag")
            elif col_mag_var.get() == 1:
                vector_plot("Colors")
        
        def norm_button_toggle():
            if norm_var.get() == 1:
                vector_plot("Norm")
                rel_mag_entry_field.pack_forget()
                rel_mag_var.set(0)
                col_mag_var.set(0)
        
        def rel_mag_button_toggle():
            if rel_mag_var.get() == 1:
                vector_plot("Rel_mag")
                rel_mag_entry_field.pack(side="top")
                norm_var.set(0)
                col_mag_var.set(0)
        
        def col_mag_button_toggle():
            if col_mag_var.get() == 1:
                vector_plot("Colors")
                rel_mag_entry_field.pack_forget()
                norm_var.set(0)
                rel_mag_var.set(0)

        def zoom_in():
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            zlim = ax.get_zlim()
            
            ax.set_xlim(xlim[0]/2,xlim[1]/2)
            ax.set_ylim(ylim[0]/2,ylim[1]/2)
            ax.set_zlim(zlim[0]/2,zlim[1]/2)

            canvas.draw()
        
        def zoom_out():
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            zlim = ax.get_zlim()

            ax.set_xlim(xlim[0]*2,xlim[1]*2)
            ax.set_ylim(ylim[0]*2,ylim[1]*2)
            ax.set_zlim(zlim[0]*2,zlim[1]*2)

            canvas.draw()

        def reset_view():
            view_size = 10
            view_text = view_size_entry.get()

            if view_text:
                view_size = int(view_text)
            
            ax.set_xlim(-(view_size + (view_size/8)), view_size + (view_size/8))
            ax.set_ylim(-(view_size + (view_size/8)), view_size + (view_size/8))
            ax.set_zlim(-(view_size + (view_size/8)), view_size + (view_size/8))

            ax.view_init(elev=30, azim=60)

            canvas.draw()
        
        def gui_initial():
            global switch_var
            norm_var.set(1)
            switch_var = 1
            vector_plot("Norm")

        def lrf_pack_forget():
            lrf_top_space_fill.pack_forget()
            xyz_text.pack_forget()
            xyz_frame.pack_forget()
            xyz_frame_spacer.pack_forget()
            rhophiz_text.pack_forget()
            rhophiz_frame.pack_forget()
            cyl_input_buttons_frame.pack_forget()
            rhophiz_frame_spacer.pack_forget()
            rthetaphi_text.pack_forget()
            rthetaphi_frame.pack_forget()
            spher_input_buttons_frame.pack_forget()
            rthetaphi_frame_spacer.pack_forget()
            lrf_bottom_frame_text.pack_forget()
            lrf_bottom_entry_frame.pack_forget()
            lrf_bottom_space_fill.pack_forget()

        def cart_gui():
            global switch_var
            lrf_pack_forget()
            root.title("Cartesian 3D Vector Field Plotter")
            lrf_top_space_fill.pack(side="top",fill="x",expand=True)

            xyz_text.config(text="Enter vector x, y, z components:")
            xyz_text.pack(side="top",fill="x",padx=20, expand=False)
            xyz_frame.pack(side="top",padx=20, fill="x",expand=False)
            xyz_frame_spacer.pack(side="top",ipady=spacerpady)
            
            rhophiz_text.config(text=f"Vector {RHO}, {PHI}, and z components:")
            rhophiz_text.pack(side="top",fill="x", padx=20, expand=False)
            rhophiz_frame.pack(side="top",fill="x",padx=20,expand=False)
            rhophiz_frame_spacer.pack(side="top",ipady=spacerpady)

            rthetaphi_text.config(text=f"Vector r, {THETA}, and {PHI} components:")
            rthetaphi_text.pack(side="top",fill="x",padx=20,expand=False)
            rthetaphi_frame.pack(side="top",fill="x",padx=20,expand=False)
            rthetaphi_frame_spacer.pack(side="top",ipady=spacerpady)

            lrf_bottom_frame_text.pack(side="top", fill="x",padx=20, expand=False)
            lrf_bottom_entry_frame.pack(side="top", fill="x", padx=20, expand=False)
            lrf_bottom_space_fill.pack(side="top",fill="x",expand=True)
            switch_var = 1

        def cyl_gui():
            global switch_var
            lrf_pack_forget()
            root.title("Cylindrical 3D Vector Field Plotter")
            lrf_top_space_fill.pack(side="top",fill="x",expand=True)

            rhophiz_text.config(text=f"Enter vector {RHO}, {PHI}, and z components:")
            rhophiz_text.pack(side="top",fill="x", padx=20, expand=False)
            rhophiz_frame.pack(side="top",fill="x",padx=20, expand=False)
            cyl_input_buttons_frame.pack(side="top", fill="x",padx=20,expand=False)
            rhophiz_frame_spacer.pack(side="top",ipady=spacerpady)

            xyz_text.config(text="Vector x, y, z components:")
            xyz_text.pack(side="top",fill="x",padx=20, expand=False)
            xyz_frame.pack(side="top",padx=20, fill="x",expand=False)
            xyz_frame_spacer.pack(side="top",ipady=spacerpady)
            
            rthetaphi_text.config(text=f"Vector r, {THETA}, and {PHI} components:")
            rthetaphi_text.pack(side="top",fill="x",padx=20,expand=False)
            rthetaphi_frame.pack(side="top",fill="x",padx=20,expand=False)
            rthetaphi_frame_spacer.pack(side="top",ipady=spacerpady)

            lrf_bottom_frame_text.pack(side="top", fill="x",padx=20, expand=False)
            lrf_bottom_entry_frame.pack(side="top", fill="x", padx=20, expand=False)
            lrf_bottom_space_fill.pack(side="top",fill="x",expand=True)
            switch_var = 2

        def spher_gui():
            global switch_var
            lrf_pack_forget()
            root.title("Spherical 3D Vector Field Plotter")
            lrf_top_space_fill.pack(side="top",fill="x",expand=True)

            rthetaphi_text.config(text=f"Enter vector r, {THETA}, and {PHI} components:")
            rthetaphi_text.pack(side="top",fill="x",padx=20,expand=False)
            rthetaphi_frame.pack(side="top",fill="x",padx=20,expand=False)
            spher_input_buttons_frame.pack(side="top", fill="x",padx=20,expand=False)
            rthetaphi_frame_spacer.pack(side="top",ipady=spacerpady)

            xyz_text.config(text="Vector x, y, z components:")
            xyz_text.pack(side="top",fill="x",padx=20, expand=False)
            xyz_frame.pack(side="top",padx=20, fill="x",expand=False)
            xyz_frame_spacer.pack(side="top",ipady=spacerpady)
            
            rhophiz_text.config(text=f"Vector {RHO}, {PHI}, and z components:")
            rhophiz_text.pack(side="top",fill="x", padx=20, expand=False)
            rhophiz_frame.pack(side="top",fill="x",padx=20, expand=False)
            rhophiz_frame_spacer.pack(side="top",ipady=spacerpady)

            lrf_bottom_frame_text.pack(side="top", fill="x",padx=20, expand=False)
            lrf_bottom_entry_frame.pack(side="top", fill="x", padx=20, expand=False)
            lrf_bottom_space_fill.pack(side="top",fill="x",expand=True)
            switch_var = 3
        
        def rho_rho_insert():
            rho_text = rho_text_field.get()
            rho_text_field.delete(0, tk.END)
            if rho_text:
                rho_text_field.insert(tk.END, rho_text+f"*{RHO}")
            else:
                rho_text_field.insert(tk.END, RHO)
        
        def rho_phi_insert():
            rho_text = rho_text_field.get()
            rho_text_field.delete(0, tk.END)
            if rho_text:
                rho_text_field.insert(tk.END, rho_text+f"*{PHI}")
            else:
                rho_text_field.insert(tk.END, PHI)

        def rho_z_insert():
            rho_text = rho_text_field.get()
            rho_text_field.delete(0, tk.END)
            if rho_text:
                rho_text_field.insert(tk.END, rho_text+f"*z")
            else:
                rho_text_field.insert(tk.END, "z")

        def phi_c_rho_insert():
            phi_c_text = phi_c_text_field.get()
            phi_c_text_field.delete(0, tk.END)
            if phi_c_text:
                phi_c_text_field.insert(tk.END, phi_c_text+f"*{RHO}")
            else:
                phi_c_text_field.insert(tk.END, RHO)

        def phi_c_phi_insert():
            phi_c_text = phi_c_text_field.get()
            phi_c_text_field.delete(0, tk.END)
            if phi_c_text:
                phi_c_text_field.insert(tk.END, phi_c_text+f"*{PHI}")
            else:
                phi_c_text_field.insert(tk.END, PHI)

        def phi_c_z_insert():
            phi_c_text = phi_c_text_field.get()
            phi_c_text_field.delete(0, tk.END)
            if phi_c_text:
                phi_c_text_field.insert(tk.END, phi_c_text+f"*z")
            else:
                phi_c_text_field.insert(tk.END, "z")

        def zcyl_rho_insert():
            zcyl_text = zcyl_text_field.get()
            zcyl_text_field.delete(0, tk.END)
            if zcyl_text:
                zcyl_text_field.insert(tk.END, zcyl_text+f"*{RHO}")
            else:
                zcyl_text_field.insert(tk.END, RHO)

        def zcyl_phi_insert():
            zcyl_text = zcyl_text_field.get()
            zcyl_text_field.delete(0, tk.END)
            if zcyl_text:
                zcyl_text_field.insert(tk.END, zcyl_text+f"*{PHI}")
            else:
                zcyl_text_field.insert(tk.END, PHI)
        
        def zcyl_z_insert():
            zcyl_text = zcyl_text_field.get()
            zcyl_text_field.delete(0, tk.END)
            if zcyl_text:
                zcyl_text_field.insert(tk.END, zcyl_text+f"*z")
            else:
                zcyl_text_field.insert(tk.END, "z")
        
        def r_r_insert():
            r_text = r_text_field.get()
            r_text_field.delete(0, tk.END)
            if r_text:
                r_text_field.insert(tk.END, r_text+f"*r")
            else:
                r_text_field.insert(tk.END, "r")
        
        def r_phi_insert():
            r_text = r_text_field.get()
            r_text_field.delete(0, tk.END)
            if r_text:
                r_text_field.insert(tk.END, r_text+f"*{PHI}")
            else:
                r_text_field.insert(tk.END, PHI)
        
        def r_theta_insert():
            r_text = r_text_field.get()
            r_text_field.delete(0, tk.END)
            if r_text:
                r_text_field.insert(tk.END, r_text+f"*{THETA}")
            else:
                r_text_field.insert(tk.END, THETA)

        def theta_r_insert():
            theta_text = theta_text_field.get()
            theta_text_field.delete(0, tk.END)
            if theta_text:
                theta_text_field.insert(tk.END, theta_text+f"*r")
            else:
                theta_text_field.insert(tk.END, "r")
        
        def theta_phi_insert():
            theta_text = theta_text_field.get()
            theta_text_field.delete(0, tk.END)
            if theta_text:
                theta_text_field.insert(tk.END, theta_text+f"*{PHI}")
            else:
                theta_text_field.insert(tk.END, PHI)
        
        def theta_theta_insert():
            theta_text = theta_text_field.get()
            theta_text_field.delete(0, tk.END)
            if theta_text:
                theta_text_field.insert(tk.END, theta_text+f"*{THETA}")
            else:
                theta_text_field.insert(tk.END, THETA)
        
        def phi_s_r_insert():
            phi_s_text = phi_s_text_field.get()
            phi_s_text_field.delete(0, tk.END)
            if phi_s_text:
                phi_s_text_field.insert(tk.END, phi_s_text+f"*r")
            else:
                phi_s_text_field.insert(tk.END, "r")
        
        def phi_s_phi_insert():
            phi_s_text = phi_s_text_field.get()
            phi_s_text_field.delete(0, tk.END)
            if phi_s_text:
                phi_s_text_field.insert(tk.END, phi_s_text+f"*{PHI}")
            else:
                phi_s_text_field.insert(tk.END, PHI)
        
        def phi_s_theta_insert():
            phi_s_text = phi_s_text_field.get()
            phi_s_text_field.delete(0, tk.END)
            if phi_s_text:
                phi_s_text_field.insert(tk.END, phi_s_text+f"*{THETA}")
            else:
                phi_s_text_field.insert(tk.END, THETA)

        menubar = tk.Menu(root)

        rootmenu = tk.Menu(menubar, tearoff=0, bg="#ededed", relief=tk.RAISED)
        rootmenu.add_command(label="Cartesian", command=cart_gui)
        rootmenu.add_command(label="Cylindrical", command=cyl_gui)
        rootmenu.add_command(label="Spherical", command=spher_gui)
        rootmenu.add_separator()
        rootmenu.add_command(label="Exit", command=on_closing)

        menubar.add_cascade(label="File", menu=rootmenu)
        root.config(menu=menubar)

        # Left Frame of Main Root Frame
        left_root_frame = ttk.Frame(root, relief=tk.SUNKEN)
        left_root_frame.pack(side='left', fill='both',expand=True)
        
        # Top Fill Space of Left Frame
        lrf_top_space_fill = ttk.Label(left_root_frame)
        lrf_top_space_fill.pack(side="top",fill="x",expand=True)


        # "Enter vector x,y,z components"
        xyz_text = ttk.Label(left_root_frame, text="Enter vector x, y, z components:")
        xyz_text.config(font=(text_font, font_size, "bold"))
        xyz_text.pack(side="top",fill="x",padx=20, expand=False)

        # Frame for xhat, yhat, and zhat components
        xyz_frame = ttk.Frame(left_root_frame)
        xyz_frame.pack(side="top",padx=20, fill="x",expand=False)

        # Text field for xhat component
        x_text_field = ttk.Entry(xyz_frame)
        x_text_field.pack(side="left",ipadx=entry_length)

        xhat_label = ttk.Label(xyz_frame,text=xhat)
        xhat_label.config(width=2, font=(text_font,font_size,"bold"))
        xhat_label.pack(side="left", padx=2)

        # Text field for yhat component
        y_text_field = ttk.Entry(xyz_frame)
        y_text_field.pack(side="left",ipadx=entry_length)

        yhat_label = ttk.Label(xyz_frame,text=yhat)
        yhat_label.config(width=2, font=(text_font,font_size,"bold"))
        yhat_label.pack(side="left", padx=2)

        # Text field for zhat component
        z_text_field = ttk.Entry(xyz_frame)
        z_text_field.pack(side="left",ipadx=entry_length)

        zhat_label = ttk.Label(xyz_frame,text=zhat)
        zhat_label.config(width=2, font=(text_font,font_size,"bold"))
        zhat_label.pack(side="left", padx=2)

        # Filler space to separate xyz_frame from rhophitheta_frame 
        xyz_frame_spacer = ttk.Label(left_root_frame)
        xyz_frame_spacer.pack(side="top",ipady=spacerpady)


        # "Enter vector rho,phi,theta components"
        rhophiz_text = ttk.Label(left_root_frame, text=f"Vector {RHO}, {PHI}, z components:")
        rhophiz_text.config(font=(text_font, font_size, "bold"))
        rhophiz_text.pack(side="top",fill="x", padx=20, expand=False)

        # Frame for rhohat, phihat, and zhat components
        rhophiz_frame = ttk.Frame(left_root_frame)
        rhophiz_frame.pack(side="top",fill="x",padx=20,expand=False)

        # Text field for rhohat component
        rho_text_field = ttk.Entry(rhophiz_frame)
        rho_text_field.pack(side="left",ipadx=entry_length)

        rhohat_label = ttk.Label(rhophiz_frame,text=rhohat)
        rhohat_label.config(width=2, font=(text_font,font_size,"bold"))
        rhohat_label.pack(side="left",padx=2)

        # Text field for phihat component
        phi_c_text_field = ttk.Entry(rhophiz_frame)
        phi_c_text_field.pack(side="left",ipadx=entry_length)

        phihat_c_label = ttk.Label(rhophiz_frame,text=phihat)
        phihat_c_label.config(width=2, font=(text_font,font_size,"bold"))
        phihat_c_label.pack(side="left",padx=2)

        # Text field for zhat component
        zcyl_text_field = ttk.Entry(rhophiz_frame)
        zcyl_text_field.pack(side="left",ipadx=entry_length)

        zcyl_label = ttk.Label(rhophiz_frame,text=zhat)
        zcyl_label.config(width=2, font=(text_font,font_size,"bold"))
        zcyl_label.pack(side="left",padx=2)

        cyl_input_width = 6
        #
        cyl_input_buttons_frame = ttk.Frame(left_root_frame)
        rho_rho_button = ttk.Button(cyl_input_buttons_frame, text=RHO, command=rho_rho_insert, width=cyl_input_width)
        rho_rho_button.pack(side="left")
        rho_phi_button = ttk.Button(cyl_input_buttons_frame, text=PHI, command=rho_phi_insert, width=cyl_input_width)
        rho_phi_button.pack(side="left")
        rho_theta_button = ttk.Button(cyl_input_buttons_frame, text="z", command=rho_z_insert, width=cyl_input_width)
        rho_theta_button.pack(side="left")
        c_rho_entry_spacer = ttk.Label(cyl_input_buttons_frame, width=3)
        c_rho_entry_spacer.pack(side="left")

        phi_c_rho_button = ttk.Button(cyl_input_buttons_frame, text=RHO, command=phi_c_rho_insert, width=cyl_input_width)
        phi_c_rho_button.pack(side="left")
        phi_c_phi_button = ttk.Button(cyl_input_buttons_frame, text=PHI, command=phi_c_phi_insert, width=cyl_input_width)
        phi_c_phi_button.pack(side="left")
        phi_c_theta_button = ttk.Button(cyl_input_buttons_frame, text="z", command=phi_c_z_insert, width=cyl_input_width)
        phi_c_theta_button.pack(side="left")
        c_phi_entry_spacer = ttk.Label(cyl_input_buttons_frame, width=3) 
        c_phi_entry_spacer.pack(side="left")

        zcyl_rho_button = ttk.Button(cyl_input_buttons_frame, text=RHO, command=zcyl_rho_insert, width=cyl_input_width)
        zcyl_rho_button.pack(side="left")
        zcyl_phi_button = ttk.Button(cyl_input_buttons_frame, text=PHI, command=zcyl_phi_insert, width=cyl_input_width)
        zcyl_phi_button.pack(side="left")
        zcyl_theta_button = ttk.Button(cyl_input_buttons_frame, text="z", command=zcyl_z_insert, width=cyl_input_width) 
        zcyl_theta_button.pack(side="left")
        c_z_entry_spacer = ttk.Label(cyl_input_buttons_frame, width=3)
        c_z_entry_spacer.pack(side="left")

        # Filler space to separate rhophiz_frame from rthetaphi_frame 
        rhophiz_frame_spacer = ttk.Label(left_root_frame)
        rhophiz_frame_spacer.pack(side="top",ipady=spacerpady)


        # "Enter vector r,theta,phi components"
        rthetaphi_text = ttk.Label(left_root_frame, text=f"Vector r, {THETA}, {PHI} components:")
        rthetaphi_text.config(font=(text_font, font_size, "bold"))
        rthetaphi_text.pack(side="top",fill="x",padx=20,expand=False)

        # Frame for rhat, thetahat, and phihat components
        rthetaphi_frame = ttk.Frame(left_root_frame)
        rthetaphi_frame.pack(side="top",fill="x",padx=20,expand=False)

        # Text field for rhat component
        r_text_field = ttk.Entry(rthetaphi_frame)
        r_text_field.pack(side="left",ipadx=entry_length)

        rhat_label = ttk.Label(rthetaphi_frame,text=rhat)
        rhat_label.config(width=2, font=(text_font,font_size,"bold"))
        rhat_label.pack(side="left",padx=2)

        # Text field for thetahat component
        theta_text_field = ttk.Entry(rthetaphi_frame)
        theta_text_field.pack(side="left",ipadx=entry_length)

        thetahat_label = ttk.Label(rthetaphi_frame,text=thetahat)
        thetahat_label.config(width=2, font=(text_font,font_size,"bold"))
        thetahat_label.pack(side="left",padx=2)

        # Text field for phihat component
        phi_s_text_field = ttk.Entry(rthetaphi_frame)
        phi_s_text_field.pack(side="left",ipadx=entry_length)

        phihat_s_label = ttk.Label(rthetaphi_frame,text=phihat)
        phihat_s_label.config(width=2, font=(text_font,font_size,"bold"))
        phihat_s_label.pack(side="left",padx=2)

        spher_input_width = 6
        #
        spher_input_buttons_frame = ttk.Frame(left_root_frame)
        r_rho_button = ttk.Button(spher_input_buttons_frame, text="r", command=r_r_insert, width=spher_input_width)
        r_rho_button.pack(side="left")
        r_phi_button = ttk.Button(spher_input_buttons_frame, text=PHI, command=r_phi_insert, width=spher_input_width)
        r_phi_button.pack(side="left")
        r_theta_button = ttk.Button(spher_input_buttons_frame, text=THETA, command=r_theta_insert, width=spher_input_width)
        r_theta_button.pack(side="left")
        s_rho_entry_spacer = ttk.Label(spher_input_buttons_frame, width=3)
        s_rho_entry_spacer.pack(side="left")

        theta_rho_button = ttk.Button(spher_input_buttons_frame, text="r", command=theta_r_insert, width=spher_input_width)
        theta_rho_button.pack(side="left")
        theta_phi_button = ttk.Button(spher_input_buttons_frame, text=PHI, command=theta_phi_insert, width=spher_input_width)
        theta_phi_button.pack(side="left")
        theta_theta_button = ttk.Button(spher_input_buttons_frame, text=THETA, command=theta_theta_insert, width=spher_input_width) 
        theta_theta_button.pack(side="left")
        s_z_entry_spacer = ttk.Label(spher_input_buttons_frame, width=3)
        s_z_entry_spacer.pack(side="left")

        phi_s_rho_button = ttk.Button(spher_input_buttons_frame, text="r", command=phi_s_r_insert, width=spher_input_width)
        phi_s_rho_button.pack(side="left")
        phi_s_phi_button = ttk.Button(spher_input_buttons_frame, text=PHI, command=phi_s_phi_insert, width=spher_input_width)
        phi_s_phi_button.pack(side="left")
        phi_s_theta_button = ttk.Button(spher_input_buttons_frame, text=THETA, command=phi_s_theta_insert, width=spher_input_width)
        phi_s_theta_button.pack(side="left")
        s_phi_entry_spacer = ttk.Label(spher_input_buttons_frame, width=3) 
        s_phi_entry_spacer.pack(side="left")

        # Filler space to separate rthetaphi_frame from lrf_bottom_frame 
        rthetaphi_frame_spacer = ttk.Label(left_root_frame)
        rthetaphi_frame_spacer.pack(side="top",ipady=spacerpady)


        # 
        lrf_bottom_frame_text = ttk.Frame(left_root_frame)
        lrf_bottom_frame_text.pack(side="top", fill="x",padx=20, expand=False)

        #
        view_size_text = ttk.Label(lrf_bottom_frame_text, text=f"View size [n*n*n]:")
        view_size_text.config(font=(text_font, font_size, "bold"))
        view_size_text.pack(side="left")

        view_size_text_spacer = ttk.Label(lrf_bottom_frame_text)
        view_size_text_spacer.pack(side="left", padx=42)

        clutter_size_text = ttk.Label(lrf_bottom_frame_text, text=f"Clutter size [n*m*v]:")
        clutter_size_text.config(font=(text_font, font_size, "bold"))
        clutter_size_text.pack(side="left")

        #
        lrf_bottom_entry_frame = ttk.Frame(left_root_frame)
        lrf_bottom_entry_frame.pack(side="top", fill="x", padx=20, expand=False)

        view_size_entry = ttk.Entry(lrf_bottom_entry_frame)
        view_size_entry.pack(side="left", ipadx=entry_length)
        
        view_size_entry_spacer = ttk.Label(lrf_bottom_entry_frame)
        view_size_entry_spacer.pack(side="left", padx=23)

        clutter_size_entry = ttk.Entry(lrf_bottom_entry_frame)
        clutter_size_entry.pack(side="left", ipadx=entry_length)

        clutter_size_entry_spacer = ttk.Label(lrf_bottom_entry_frame)
        clutter_size_entry_spacer.pack(side="left",padx=20)

        plot_button = ttk.Button(lrf_bottom_entry_frame, text="Plot", command=plot_type_check)
        plot_button.pack(side="left")

        #
        lrf_bottom_space_fill = ttk.Label(left_root_frame)
        lrf_bottom_space_fill.pack(side="top",fill="x",expand=True)
        
        
        # Right Frame of Main Root Frame
        right_root_frame = ttk.Frame(root, relief=tk.RAISED)
        right_root_frame.pack(side='right',fill='both',expand=False)

        #
        rrf_top_frame = ttk.Frame(right_root_frame)
        rrf_top_frame.pack(side="top", fill="x", expand=True)
        
        #
        checkmark_box_frame = ttk.Frame(rrf_top_frame)
        checkmark_box_frame.pack(side="left")

        norm_var = tk.IntVar()
        normalize_check_button = ttk.Checkbutton(checkmark_box_frame, text="Normalize", variable=norm_var, onvalue=1, offvalue=0, width=20, command=norm_button_toggle)
        normalize_check_button.pack(side="top", pady=5)
        
        rel_mag_var = tk.IntVar()
        rel_mag_button = ttk.Checkbutton(checkmark_box_frame, text="Relative Magnitude", variable=rel_mag_var, onvalue=1, offvalue=0, width=20,command=rel_mag_button_toggle)
        rel_mag_button.pack(side="top", pady=5)
        rel_mag_entry_field = ttk.Entry(checkmark_box_frame)

        col_mag_var = tk.IntVar()
        col_mag_button = ttk.Checkbutton(checkmark_box_frame, text="Color Magnitude", variable=col_mag_var, onvalue=1, offvalue=0, width=20, command=col_mag_button_toggle)
        col_mag_button.pack(side="bottom", pady=5)

        #


        #
        plot_view_buttons_frame = ttk.Frame(rrf_top_frame)
        plot_view_buttons_frame.pack(side="left")

        #
        zoom_buttons_frame = ttk.Frame(plot_view_buttons_frame)
        zoom_buttons_frame.pack(side="top")
        
        zoom_in_button = ttk.Button(zoom_buttons_frame, text="+", width=5, command=zoom_in)
        zoom_in_button.pack(side="top", ipadx=2, ipady=2, expand=False)
        
        
        zoom_out_button = ttk.Button(zoom_buttons_frame, text="-", width=5, command=zoom_out)
        zoom_out_button.pack(side="top", ipadx=2, ipady=2, expand=False)
        
        reset_view_frame = ttk.Frame(plot_view_buttons_frame)
        reset_view_frame.pack(side="top")

        reset_view_button = ttk.Button(reset_view_frame, text="\u2302", width=5, command=reset_view)
        reset_view_button.pack(side="top", ipadx=2, ipady=2, expand=False)

        gui_initial()

    root_gui()

    root.protocol("WM_DELETE_WINDOW", on_closing)
    root.mainloop()
