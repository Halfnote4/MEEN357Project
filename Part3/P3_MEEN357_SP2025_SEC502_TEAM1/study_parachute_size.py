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
import subfunctions_EDL 

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

Task6 = subfunctions_EDL.Task6
#Reminder to change Task6 value to True or False depending on which you want to run

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



if Task6:
    axs[0].plot(parachute_diam, t_touch, '--', color='green')
    axs[0].set_title('MEF Consideration', fontsize=10)
else:
    axs[0].plot(parachute_diam, t_touch, '--', color='blue')
    axs[0].set_title('Initial', fontsize=10)

axs[0].set_ylabel('Time (s)')
axs[0].grid()


#Subplot 2 - Rover speed (relative to ground) at simulation termination vs. parachute diameter

if Task6:
    axs[1].plot(parachute_diam, v_touch, '--', color='green')
else:
    axs[1].plot(parachute_diam, v_touch, '--', color='blue')
axs[1].set_ylabel('Rover Speed (m/s)')
axs[1].grid()


#Subplot 3- Rover landing success (1=success; 0=failure) vs. parachute diameter 

if Task6:
    axs[2].plot(parachute_diam, success, 'o--', color='green')
else:
    axs[2].plot(parachute_diam, success, 'o--', color='blue')

axs[2].set_xlabel('Parachute Diameter (m)')
axs[2].set_ylabel('Rover Landing Success')


#Show plots
plt.show()
