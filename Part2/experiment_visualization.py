#Submittable code for the experiment visualization

import subfunctions as sf
import numpy as np
import matplotlib.pyplot as plt
import define_experiment as de
import scipy.interpolate as interp

#Eddy will do this


'''Task 2'''
#This function will plot the results of the experiment / experimental terrain


#To calculate terrain angles, we need to first interpolate the terrain data
alpha_fun = interp1d(alpha_dist, alpha_deg, kind = 'cubic', fill_value='extrapolate') #fit the cubic spline - given code