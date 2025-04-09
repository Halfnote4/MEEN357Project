#study_parachute_size python script
from define_edl_system import*
from define_rovers import*
from define_planet import*
from main_edl_simulation import*
from define_mission_events import*
from subfunctions_EDL import*
from redefine_edl_system import*
import numpy as np
import matplotlib.pyplot as plt

#overrides preivous values with prescribed values
# initial conditions
edl_system['altitude'] = 11000    # [m] initial altitude
edl_system['velocity'] = -590     # [m/s] initial velocity
edl_system['parachute']['deployed'] = True   # our parachute is open
edl_system['parachute']['ejected'] = False   # and still attached
edl_system['rover']['on_ground'] = False # the rover has not yet landed

tmax = 2000   # [s] maximum simulated time

parachute_diam = np.arange(14, 19, 0.5) # [m] parachute diameters to test

t_touch = np.zeros(len(parachute_diam)) # time of touch down
v_touch = np.zeros(len(parachute_diam)) # velocity at touch down
success = np.zeros(len(parachute_diam)) # success of landing (1=success; 0=failure)

Task6= True # True if you want to run Task 6, False if you want to run Task 5
''' New Drag Function in Task 6'''


Mach = [0.25, 0.5, 0.65, 0.7, 0.8, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 1.9, 2.0, 2.2, 2.5, 2.6]

MEF = [1.0, 1.0, 1.0, 0.97, 0.91, 0.72, 0.66, 0.75, 0.90, 0.96, 0.990, 0.999, 0.992, 0.98, 0.91, 0.85, 0.82, 0.75, 0.64, 0.62]

#MEF Model using the given array of values using cubic interpolation
MEF = np.interp1d(Mach, MEF, kind='cubic', fill_value="extrapolate")




''' If Commented out, this will produce task 5 plot OR if you change value to Task6 = True, it will produce task 6 plot'''
if Task6:
    Task6 = True

for i in range(len(parachute_diam)):
    edl_system = redefine_edl_system(edl_system) #need to reset everything
    edl_system['parachute']['diameter'] = parachute_diam[i]

# the simulation. changing last argument to false turns off message echo
    [t, Y, edl_system] = simulate_edl(edl_system, mars, mission_events, tmax, True)
    #put interesting results into arrays for plotting
    #t[-1] is the last time step in the simulation
    #Y[-1] is the last state vector in the simulation
    #edl_system['rover']['on_ground'] is 1 if rover has landed successfully, 0 if not
    t_touch[i] = t[-1] # time of touch down
    v_touch[i] = Y[0,-1] # velocity at touch down

    #Success if: v=<-1
    if v_touch[i] <= -1:
        success[i] = 0
    else:
        success[i] = 1




#Plots the results of the simulation in 3 subplots
plt.close()
fig0, axs = plt.subplots(3)
#plt.tight_layout()

#Subplot 1 - Simulated time (i.e., time at termination of simulation) vs. parachute diameter

axs[0].plot(parachute_diam, t_touch, '--')
axs[0].set_ylabel('Time (s)')
axs[0].set_title('Initial', fontsize=10)
axs[0].grid()

#Subplot 2 - Rover speed (relative to ground) at simulation termination vs. parachute diameter


axs[1].plot(parachute_diam, v_touch, '--')
axs[1].set_ylabel('Rover Speed (m/s)')
axs[1].grid()


#Subplot 3- Rover landing success (1=success; 0=failure) vs. parachute diameter 

axs[2].plot(parachute_diam, success, 'o--')
axs[2].set_xlabel('Parachute Diameter (m)')
axs[2].set_ylabel('Rover Landing Success')


#Show plots
#plt.show()
