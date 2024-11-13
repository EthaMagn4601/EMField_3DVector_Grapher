# EMField_3DVector_Grapher
How to use:
- Input your vector field equation that you would like to plot into the corresponding component entry fields.
- Input the view size that you wish the grid to display.
  - I.e. This will make the view of the graph display n units along the ±x-axis, ±y-axis, and ±z-axis.
  - Default value is 10 units.
- Input the number of sample spaces you'd like to see on the plot
  - Default value is 3, 3, 3 corresponding to 3 x-axis samples, 3 y-axis
  - Minimum value along each axis is 1
- Select how you want to display the 3D field with the checkboxes to the left of the 3D plot
  - Select Normalize to make the length of every vector the same.
  - Select Relative Magnitude to give a constant (input value) multiplied by the relative magnitude of each vector.
  - Select Color Magnitude to show relative magnitude with colors.
- To switch between cartesian, cylindrical, or spherical vector field plotters, click on the upper left-hand corner of the application to file and select your desired field plotter. 

# Warnings:
- Recommended gridsize [x*y*z]: x < 6, y < 6, z < 6 (max should be 10 x 10 x 10, due to matplotlib being very laggy past this)
- Don't input random variables into the text entry field areas (i.e. f, g, etc.)
