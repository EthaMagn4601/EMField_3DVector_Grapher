import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np
import matplotlib.pyplot as plt
import main

def plot_quiver():
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        if main.plot_quiver_count == 0:

            # Generate sample data
            x, y, z = np.meshgrid(np.arange(0, 11, 1),
                                np.arange(0, 11, 1),
                                np.arange(0, 11, 1), 
                                indexing='ij')
            
            u = main.results["xhat_out"]
            v = main.results["yhat_out"]
            w = main.results["zhat_out"]

            ax.quiver(x, y, z, u, v, w, length=0.5, normalize=True)

            # Embed the plot in Tkinter
            canvas = FigureCanvasTkAgg(fig, master=root)
            canvas.draw()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        else:
            root.destroy()
            root = tk.Tk()
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
            x, y, z = np.meshgrid(np.arange(0, 11, 1),
                                np.arange(0, 11, 1),
                                np.arange(0, 11, 1), 
                                indexing='ij')
            
            u = main.results["xhat_out"]
            v = main.results["yhat_out"]
            w = main.results["zhat_out"]
            ax.quiver(x, y, z, u, v, w, length=0.5, normalize=True)

            # Embed the plot in Tkinter
            canvas = FigureCanvasTkAgg(fig, master=root)
            canvas.draw()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        main.plot_quiver_count = main.plot_quiver_count + 1

def root_quit():
     root.destroy()

try:
    while(True):
            root = tk.Tk()
            root.wm_title("3D Quiver in Tkinter")

            button = tk.Button(master=root, text="Plot Quiver", command=plot_quiver)
            button.pack(side=tk.TOP)

            button = tk.Button(master=root, text="Quit", command=root_quit)
            button.pack(side=tk.BOTTOM)
            tk.mainloop()
except:
    SystemExit

