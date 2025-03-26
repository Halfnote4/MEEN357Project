#Submittable code for the experiment visualization

import subfunctions as sf
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as spi
import define_experiment as de

experiment = de.experiment1()[0]
alpha_dist = experiment['alpha_dist']
alpha_deg = experiment['alpha_deg']

'''Task 2'''
#This function will plot the results of the experiment / experimental terrain
#using 100 points, evenly spaced between the minimum and maximum distance values
#create a single figure of terrain angle vs. position with the data contained in the define_experiment.py file


spaced_dist = np.linspace(min(alpha_dist), max(alpha_dist), 100) #create 100 evenly spaced points



#To calculate terrain angles, we need to first interpolate the terrain data

alpha_fun = spi.interp1d(alpha_dist, alpha_deg, kind = 'cubic', fill_value='extrapolate') #fit the cubic spline - given code

 #find the corresponding degrees for the 100 points

'''
#Create a figure, with evaluate the terrain angle using 100 points, evenly spaced between the minimum and
maximum distance (experiment[‘alpha_dist’] contains distance information) in the experiment file. The
script should create a single figure of terrain angle vs. position with the data contained in
define_experiment.py plotted as star symbols, and the 100 evaluated terrain angles plotted as a
line. Axes should be labeled appropriately.
'''



plt.figure()
plt.plot(alpha_dist, alpha_deg, 'r*', label = 'Experimental Data') #plot the experimental data
plt.plot(spaced_dist, alpha_fun(spaced_dist), label = 'Interpolated Data') #plot the interpolated data
plt.xlabel('Distance (m)')
plt.ylabel('Angle (degrees)')
plt.title('Terrain Angle vs. Position')
plt.legend()
plt.show()
#This function will plot the results of the experiment / experimental terrain
