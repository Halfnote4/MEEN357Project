from subfunctions import *
import numpy as np
import scipy.optimize as so
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#Assumptions according to specifications
Crr_array = np.linspace(0.01,0.4,25)
slope_array_deg = np.linspace(-10,35,25)

omega_nl = rover['wheel_assembly']['motor']['speed_noload']


CRR, SLOPE = np.meshgrid(Crr_array, slope_array_deg)


VMAX = np.zeros(np.shape(CRR), dtype = float)



N = np.shape(CRR)[0]
for i in range(N):
    for j in range(N):
        Crr_sample = float(CRR[i,j])
        slope_sample = float(SLOPE[i,j])
        fun = lambda omega: F_net(omega, slope_sample, rover, planet, Crr_sample)
        
        try:
            sol = root_scalar(fun, method='bisect', bracket=[0, omega_nl])
            VMAX[i,j] = (sol.root * rover['wheel_assembly']['wheel']['radius'])/ get_gear_ratio(speed_reducer) 
        except ValueError:
            VMAX[i,j] = np.nan



figure = plt.figure()
ax = Axes3D(figure, elev = 0, azim = 0) # where N1 and N2 will control the 3D view
ax.plot_surface(CRR, SLOPE, VMAX)

            
        