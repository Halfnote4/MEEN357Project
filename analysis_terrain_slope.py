from subfunctions import *
import numpy as np
import scipy.optimize as so
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt

#GOAL: Plot of v_max vs slope_array_deg


#Assumptions according to specifications
Crr = 0.2
slope_list_deg = np.linspace(-10,35,25)
omega_max = np.zeros(len(slope_list_deg), dtype=float)
omega_nl = rover['wheel_assembly']['motor']['speed_noload']



#To find vmax -> root of fmax = 0 at each given angle
for ii in range(len(slope_list_deg)):
    fun = lambda omega: F_net(omega, float(slope_list_deg[ii]), rover, planet, Crr)
    sol = root_scalar(fun, method='bisect', bracket=[0, omega_nl])
    omega_max[ii] = rover['wheel_assembly']['wheel']['radius']*(sol.root/get_gear_ratio(rover['wheel_assembly']['speed_reducer']))
    

#omegain to out with gear ratio -> vel

#plotting vmax vs slope_array_deg
plt.figure(1)
plt.plot(slope_list_deg, omega_max)  
plt.title('Terrain Slope Analysis')
plt.xlabel('Terrain Angle [degrees]')
plt.ylabel('Rover Speed [m/s]')
#plt.show()