from subfunctions import *
import numpy as np
import scipy.optimize as so
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt


Crr_array = np.linspace(0.01,0.4,25) #different Crr values
omega_nl = rover['wheel_assembly']['motor']['speed_noload']
omega_max = np.zeros(len(Crr_array), dtype=float)

#To find vmax -> root of fmax = 0 at each given angle
for ii in range(len(Crr_array)):
    fun = lambda omega: F_net(omega, 0, rover, planet, float(Crr_array[ii]))
    sol = root_scalar(fun, method='bisect', bracket=[0, omega_nl])
    omega_max[ii] = rover['wheel_assembly']['wheel']['radius']*(sol.root/get_gear_ratio(rover['wheel_assembly']['speed_reducer']))

#GOAL: Plot of v_max vs Crr

#plotting vmax vs Crr_array
plt.figure(1)
plt.plot(Crr_array, omega_max)  
plt.title('Rolling Resistance Analysis')
plt.xlabel('Crr')
plt.ylabel('Rover Speed [m/s]')
#plt.show()