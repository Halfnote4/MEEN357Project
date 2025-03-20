#Submittable code for the efficiency visualization

import subfunctions as sf
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as spi
import define_rovers as dr


rover = dr.define_rover_1()[0]


# initialize rover



# interpolation function (cubic spline)
effcy_fun = spi.interp1d(rover['wheel_assembly']['motor']['effcy_tau'], rover['wheel_assembly']['motor']['effcy'], kind='cubic')

# torque between the minimum and maximum values
new_torque_values = np.linspace(min(rover['wheel_assembly']['motor']['effcy_tau']), max(rover['wheel_assembly']['motor']['effcy_tau']), 100)

# Calculate the corresponding efficiency values for the 100 points
new_effcy_values = effcy_fun(new_torque_values)

# Plot the efficiency vs. torque
plt.figure()

# Plot the original data points (star symbols)
plt.scatter(rover['wheel_assembly']['motor']['effcy_tau'], rover['wheel_assembly']['motor']['effcy'], color='r', label='Data Points', marker='*', s=100)

# Plot the interpolated efficiency curve
plt.plot(new_torque_values, new_effcy_values, label='Interpolated Efficiency Curve', color='b', linewidth=2)


plt.xlabel('Motor Torque (Nm)')
plt.ylabel('Motor Efficiency')
plt.title('Motor Torque vs. Efficiency')

plt.legend()
plt.show()