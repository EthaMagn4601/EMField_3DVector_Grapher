import functions as func
import numpy as np
import matplotlib.pyplot as plt

# Uncomment np.set_print to view full matrix output
np.set_printoptions(threshold=np.inf)

clutter_size = 11
grid_size = 11

plot_quiver_count = 0

# cartesian equation example
x_input = "1"
y_input = "1"
z_input = "1"
cart_equation = func.vector_cart_to_cart(x_input, y_input, z_input)

# cylindrical equation example
# rho_input = func.PHI
# phi_input = func.PHI
# z_input = func.PHI
# cart_equation = func.vector_cyli_to_cart(rho_input, phi_input, z_input)

# spherical equation example
# r_input = "4"
# theta_input = "r"
# phi_input = func.PHI
# cart_equation = func.vector_spher_to_cart(r_input, theta_input, phi_input)

# print(cart_equation["xhat"])
# print(cart_equation["yhat"])
# print(cart_equation["zhat"])

results = func.cart_equation_solver(cart_equation, clutter_size)

# print(results["xhat_out"])
# print(results["yhat_out"])
# print(results["zhat_out"])

# print(f"Result1: "+ str(results["r1"]))
# print(f"Result2: "+ str(results["r2"]))
# print(f"Result3: "+ str(results["r3"]))

# ax = plt.figure().add_subplot(projection='3d')

# # Make the grid
# x, y, z = np.meshgrid(np.arange(0, 11, 1),
#                       np.arange(0, 11, 1),
#                       np.arange(0, 11, 1), 
#                       indexing='ij')

# # Make the direction data for the arrows
# u = results["xhat_out"]
# v = results["yhat_out"]
# w = results["zhat_out"]

# ax.quiver(x, y, z, u, v, w, length=0.5, normalize=True)

# plt.show()
               

