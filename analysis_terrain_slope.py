from subfunctions import *
import numpy


import scipy.optimize as so

#GOAL: Plot of v_max vs slope_array_deg


#Assumptions according to specifications
Crr = 0.2
slope_array_deg = numpy.linspace(-10,35,25)

#set vmax as empty array
vmax = np.zeros(len(slope_array_deg), dtype = float)

#To find vmax -> root of fmax = 0 at each given angle
for i in range(len(slope_array_deg)):
    vmax[i] = so.root(F_net(omega,float(slope_array_deg[i]),rover,planet,Crr))