#Submittable code for the rover experiment1
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as spi
import define_experiment as de
import subfunctions as sf
import end_of_mission_event as em


'''
From the worksheet given:

For this task you will create a Python script, rover_experiment1.py to simulate the trajectory of
the rover using the simulate_rover function.
• You need to load the experiment and end_event dictionaries (download
define_experiment.py from Canvas first)
• You should also set end_event fields to the following values:
o end_event[‘max_distance’] = 1000
o end_event[‘max_time’] = 10000
o end_event[‘min_velocity’] = 0.01
#done
The script should create a single figure with three subfigures (in a 3x1 arrangement):
1. Position vs. time
2. Velocity vs. time
3. Power vs. time



For this task you should also include the figure in your write-up and explain what you observe. The
explanation must be based on what you know about the particular terrain that the rover must traverse.
Please discuss: (a) if the position vs time graph is linear or non-linear and why that is the case (b) The
relationship between velocity and power that you observe and the basis of that relationship (c) how the velocity profile is related to the experimental terrain from the experiment_visualization.py plot, and (d) If
the velocity graph is smooth or not and why that is the case.
Finally, please include a table of the rover[‘telemetry’] data for the fields: completion_time,
distance_traveled, max_velocity, average_velocity, battery_energy, and batt_energy_per_distance.




Finally, please include a table of the rover[‘telemetry’] data for the fields: completion_time,
distance_traveled, max_velocity, average_velocity, battery_energy, and batt_energy_per_distance.
'''

#Eddy/Eimaan as needed

omega = 1 # rad/s (motor shaft speed)
angle = 5 # degrees (terrain angle)
Crr = 0.1

from define_rovers import define_rover_1
rover, planet = define_rover_1() # this is what is given in Appendix A and B


'''
Fd = F_drive(omega, rover)
Frr = F_rolling(omega, angle, rover, planet, Crr)
Fg = F_gravity(angle, rover, planet)
Fnet = F_net(omega, angle, rover, planet, Crr)
'''

#End Event dictionary definition as given, to be used in the simulate rover function
end_event = {
    'max_distance' : 1000,
    'max_time' : 10000,
    'min_velocity' : 0.01,
}


events = em.end_of_mission_event(end_event) #to call in simulate rover
experiment = de.experiment1()[0] #to call in simulate rover

rov = sf.simulate_rover(rover, planet, experiment, end_event)['telemetry'] #create var to call the updated telemetry data of the rover


#Subfigures, as requested:

fig, axs = plt.subplots(3, 1, figsize=(10, 10))

#1. Position vs. time
axs[0].plot(rov['Time'], rov['position'])
axs[0].set_title('Position vs. Time')
axs[0].set_xlabel('Time (s)')
axs[0].set_ylabel('Position (m)')
#axs[0].grid()


#2. Velocity vs. time
axs[1].plot(rov['Time'], rov['velocity'])
axs[1].set_title('Velocity vs. Time')
axs[1].set_xlabel('Time (s)')
axs[1].set_ylabel('Velocity (m/s)')
#axs[1].grid()


#3. Power vs. time
axs[2].plot(rov['Time'], rov['power'])
axs[2].set_title('Power vs. Time')
axs[2].set_xlabel('Time (s)')
axs[2].set_ylabel('Power (W)')
#axs[2].grid()


plt.tight_layout()
plt.show()