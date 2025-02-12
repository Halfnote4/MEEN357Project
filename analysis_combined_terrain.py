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
        except:
            VMAX[i,j] = np.NAN

figure = plt.figure()
ax = plt.axes(projection = '3d')

#CODE AXES/TITLE

ax.plot_surface(CRR, SLOPE, VMAX)
ax.view_init(azim = 30, elev = 30)
ax.set_box_aspect(None, zoom=0.8)
ax.set_title('Combined Terrain Analysis')
ax.set_xlabel('Crr')
ax.set_ylabel('Rover Speed [m/s]')
ax.set_zlabel('Terrain angle (degrees)')

#plt.show()

 