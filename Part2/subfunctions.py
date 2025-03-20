"""###########################################################################
#   This file contains subfunctions for Phase 1 of the TAMU MEEN 357 project
#
#   Created by: MEEN 357 Instructional Team
#   Last Modified: 4 October 2023
###########################################################################"""

import math
import numpy as np
import scipy.interpolate as spi
import scipy.integrate as spitegrate


def get_mass(rover):
    """
    Inputs:  rover:  dict      Data structure containing rover parameters
    
    Outputs:     m:  scalar    Rover mass [kg].
    """
    
    # Check that the input is a dict
    if type(rover) != dict:
        raise Exception('Input must be a dict')
    
    # add up mass of chassis, power subsystem, science payload, 
    # and components from all six wheel assemblies
    m = rover['chassis']['mass'] \
        + rover['power_subsys']['mass'] \
        + rover['science_payload']['mass'] \
        + 6*rover['wheel_assembly']['motor']['mass'] \
        + 6*rover['wheel_assembly']['speed_reducer']['mass'] \
        + 6*rover['wheel_assembly']['wheel']['mass'] \
    
    return m


def get_gear_ratio(speed_reducer):
    """
    Inputs:  speed_reducer:  dict      Data dictionary specifying speed
                                        reducer parameters
    Outputs:            Ng:  scalar    Speed ratio from input pinion shaft
                                        to output gear shaft. Unitless.
    """
    
    # Check that the input is a dict
    if type(speed_reducer) != dict:
        raise Exception('Input must be a dict')
    
    # Check 'type' field (not case sensitive)
    if speed_reducer['type'].lower() != 'reverted':
        raise Exception('The speed reducer type is not recognized.')
    
    # Main code
    d1 = speed_reducer['diam_pinion']
    d2 = speed_reducer['diam_gear']
    
    Ng = (d2/d1)**2
    
    return Ng


def tau_dcmotor(omega, motor):
    """
    Inputs:  omega:  numpy array      Motor shaft speed [rad/s]
             motor:  dict             Data dictionary specifying motor parameters
    Outputs:   tau:  numpy array      Torque at motor shaft [Nm].  Return argument
                                      is same size as first input argument.
    """
    
    # Check that the first input is a scalar or a vector
    if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(omega, np.ndarray):
        omega = np.array([omega],dtype=float) # make the scalar a numpy array
    elif len(np.shape(omega)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')

    # Check that the second input is a dict
    if type(motor) != dict:
        raise Exception('Second input must be a dict')
        
    # Main code
    tau_s    = motor['torque_stall']
    tau_nl   = motor['torque_noload']
    omega_nl = motor['speed_noload']
    
    # initialize
    tau = np.zeros(len(omega),dtype = float)
    for ii in range(len(omega)):
        if omega[ii] >= 0 and omega[ii] <= omega_nl:
            tau[ii] = tau_s - (tau_s-tau_nl)/omega_nl *omega[ii]
        elif omega[ii] < 0:
            tau[ii] = tau_s
        elif omega[ii] > omega_nl:
            tau[ii] = 0
        
    return tau
    
    


def F_rolling(omega, terrain_angle, rover, planet, Crr):
    """
    Inputs:           omega:  numpy array     Motor shaft speed [rad/s]
              terrain_angle:  numpy array     Array of terrain angles [deg]
                      rover:  dict            Data structure specifying rover 
                                              parameters
                    planet:  dict            Data dictionary specifying planetary 
                                              parameters
                        Crr:  scalar          Value of rolling resistance coefficient
                                              [-]
    
    Outputs:           Frr:  numpy array     Array of forces [N]
    """
    
    # Check that the first input is a scalar or a vector
    if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(omega, np.ndarray):
        omega = np.array([omega],dtype=float) # make the scalar a numpy array
    elif len(np.shape(omega)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the second input is a scalar or a vector
    if (type(terrain_angle) != int) and (type(terrain_angle) != float) and (not isinstance(terrain_angle, np.ndarray)):
        raise Exception('Second input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(terrain_angle, np.ndarray):
        terrain_angle = np.array([terrain_angle],dtype=float) # make the scalar a numpy array
    elif len(np.shape(terrain_angle)) != 1:
        raise Exception('Second input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the first two inputs are of the same size
    if len(omega) != len(terrain_angle):
        raise Exception('First two inputs must be the same size')
    
    # Check that values of the second input are within the feasible range  
    if max([abs(x) for x in terrain_angle]) > 75:    
        raise Exception('All elements of the second input must be between -75 degrees and +75 degrees')
        
    # Check that the third input is a dict
    if type(rover) != dict:
        raise Exception('Third input must be a dict')
        
    # Check that the fourth input is a dict
    if type(planet) != dict:
        raise Exception('Fourth input must be a dict')
        
    # Check that the fifth input is a scalar and positive
    if (type(Crr) != int) and (type(Crr) != float):
        raise Exception('Fifth input must be a scalar')
    if Crr <= 0:
        raise Exception('Fifth input must be a positive number')
        
    # Main Code
    m = get_mass(rover)
    g = planet['g']
    r = rover['wheel_assembly']['wheel']['radius']
    Ng = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    
    v_rover = r*omega/Ng
    
    Fn = np.array([m*g*math.cos(math.radians(x)) for x in terrain_angle],dtype=float) # normal force
    Frr_simple = -Crr*Fn # simple rolling resistance
    
    Frr = np.array([math.erf(40*v_rover[ii]) * Frr_simple[ii] for ii in range(len(v_rover))], dtype = float)
    
    return Frr


def F_gravity(terrain_angle, rover, planet):
    """
    Inputs:  terrain_angle:  numpy array   Array of terrain angles [deg]
                     rover:  dict          Data structure specifying rover 
                                            parameters
                    planet:  dict          Data dictionary specifying planetary 
                                            parameters
    
    Outputs:           Fgt:  numpy array   Array of forces [N]
    """
    
    # Check that the first input is a scalar or a vector
    if (type(terrain_angle) != int) and (type(terrain_angle) != float) and (not isinstance(terrain_angle, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(terrain_angle, np.ndarray):
        terrain_angle = np.array([terrain_angle],dtype=float) # make the scalar a numpy array
    elif len(np.shape(terrain_angle)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that values of the first input are within the feasible range  
    if max([abs(x) for x in terrain_angle]) > 75:    
        raise Exception('All elements of the first input must be between -75 degrees and +75 degrees')

    # Check that the second input is a dict
    if type(rover) != dict:
        raise Exception('Second input must be a dict')
    
    # Check that the third input is a dict
    if type(planet) != dict:
        raise Exception('Third input must be a dict')
        
    # Main Code
    m = get_mass(rover)
    g = planet['g']
    
    Fgt = np.array([-m*g*math.sin(math.radians(x)) for x in terrain_angle], dtype = float)
        
    return Fgt


def F_drive(omega, rover):
    """
    Inputs:  omega:  numpy array   Array of motor shaft speeds [rad/s]
             rover:  dict          Data dictionary specifying rover parameters
    
    Outputs:    Fd:  numpy array   Array of drive forces [N]
    """
    
    # Check that the first input is a scalar or a vector
    if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(omega, np.ndarray):
        omega = np.array([omega],dtype=float) # make the scalar a numpy array
    elif len(np.shape(omega)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')

    # Check that the second input is a dict
    if type(rover) != dict:
        raise Exception('Second input must be a dict')
    
    # Main code
    Ng = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    
    tau = tau_dcmotor(omega, rover['wheel_assembly']['motor'])
    tau_out = tau*Ng
    
    r = rover['wheel_assembly']['wheel']['radius']
    
    # Drive force for one wheel
    Fd_wheel = tau_out/r 
    
    # Drive force for all six wheels
    Fd = 6*Fd_wheel
    
    return Fd


def F_net(omega, terrain_angle, rover, planet, Crr):
    """
    Inputs:           omega:  list     Motor shaft speed [rad/s]
              terrain_angle:  list     Array of terrain angles [deg]
                      rover:  dict     Data structure specifying rover 
                                      parameters
                     planet:  dict     Data dictionary specifying planetary 
                                      parameters
                        Crr:  scalar   Value of rolling resistance coefficient
                                      [-]
    
    Outputs:           Fnet:  list     Array of forces [N]
    """
    
    # Check that the first input is a scalar or a vector
    if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
    # if (not isinstance(omega, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(omega, np.ndarray):
        omega = np.array([omega],dtype=float) # make the scalar a numpy array
    elif len(np.shape(omega)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the second input is a scalar or a vector
    if (type(terrain_angle) != int) and (type(terrain_angle) != float) and (not isinstance(terrain_angle, np.ndarray)):
        raise Exception('Second input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(terrain_angle, np.ndarray):
        terrain_angle = np.array([terrain_angle],dtype=float) # make the scalar a numpy array
    elif len(np.shape(terrain_angle)) != 1:
        raise Exception('Second input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the first two inputs are of the same size
    if len(omega) != len(terrain_angle):
        raise Exception('First two inputs must be the same size')
    
    # Check that values of the second input are within the feasible range  
    if max([abs(x) for x in terrain_angle]) > 75:    
        raise Exception('All elements of the second input must be between -75 degrees and +75 degrees')
        
    # Check that the third input is a dict
    if type(rover) != dict:
        raise Exception('Third input must be a dict')
        
    # Check that the fourth input is a dict
    if type(planet) != dict:
        raise Exception('Fourth input must be a dict')
        
    # Check that the fifth input is a scalar and positive
    if (type(Crr) != int) and (type(Crr) != float):
        raise Exception('Fifth input must be a scalar')
    if Crr <= 0:
        raise Exception('Fifth input must be a positive number')
    
    # Main Code
    Fd = F_drive(omega, rover)
    Frr = F_rolling(omega, terrain_angle, rover, planet, Crr)
    Fg = F_gravity(terrain_angle, rover, planet)
    
    Fnet = Fd + Frr + Fg # signs are handled in individual functions
    
    return Fnet

'''
#%% Hint for students
omega = 1 # rad/s (motor shaft speed)
angle = 5 # degrees (terrain angle)
Crr = 0.1

from define_rovers import define_rover_1
rover, planet = define_rover_1() # this is what is given in Appendix A and B

Fd = F_drive(omega, rover)
Frr = F_rolling(omega, angle, rover, planet, Crr)
Fg = F_gravity(angle, rover, planet)
Fnet = F_net(omega, angle, rover, planet, Crr)
'''



'''
PART 2

'''




def motorW(v, rover):

    """
    motorW

    Compute the rotational speed of the motor shaft [rad/s] given the translational velocity of the rover and the rover dictionary.
    Inputs:
        v : 1D numpy array [m/s]
        rover: dict, dictionary containing rover parameters

    Return:
        w : 1D numpy array [rad/s]

    Notes:
    -Should call get_gear_ratio() to get the gear ratio
    -Should validate (a) the first input is a scalar or vector. If a vector, it should be defined as a numpy array and matrices are not allowed
    -Should (b) the second input is a dict
    >Should any of these fail, need to raise an exception

    """

    #Input verification

    #a
    if not isinstance(v, (np.ndarray, float, int)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')    
    
    try:
        if len(v.shape) != 1:
            raise Exception('First input must be a 1D numpy array')
    except:
        pass    

    
    #b
    if type(rover) != dict:
        raise Exception('Second input must be a dict')
    

    #Main code
    r = rover['wheel_assembly']['wheel']['radius'] # get wheel radius
    Ng = get_gear_ratio(rover['wheel_assembly']['speed_reducer']) # get gear ratio
    
    w = v * Ng / r # motor shaft speed equation
    
    return w


def rover_dynamics(t, y, rover, planet, experiment):


    """
    rover_dynamics - Eddy

    Computes the derivative of the state vector ([velocity, position]) for the rover given its current state. 
    It requires rover and experiment dictionary input parameters. It is intended to be passed to an ODE
    solver.
    Inputs:
        t : scalar, time sample [s]
        y : 1D numpy array, state vector [m/s, m]
        rover : dict, dictionary containing rover parameters
        planet : dict, dictionary containing planet parameters
        experiment : dict, dictionary containing experiment parameters

    Return:
        dydt : 1D numpy array, acceleration [m/s^2] and rover velocity [m/s]


    Notes:
    -Should validate all meaningful inputs

    -Should call F_net() to get the net force
    -Should call MotorW() to get the motor speed
    -Should call get_mass() to get the rover mass
    -Create an interpolation function to calculate terrain angle and evaluate it at the current rover position (terrain_angle = alpha_fun(y[1]))

    -Assumes y[1] has been defined to be the rover position


    """

    #Input verification

    #y must be a 1D array
    if not isinstance(y, (np.ndarray, int, float)):
        raise Exception('y must be a numpy array or scalar')
    if isinstance(y, np.ndarray):
        if len(y.shape) != 1:
            raise Exception('y must be a 1D numpy array')
        
    #t must be scalar
    if not isinstance(t, (float, int)):
        raise Exception('t must be a scalar')
    
    #rover must be dict
    if type(rover) != dict:
        raise Exception('rover must be a dict')
    
    #planet must be dict
    if type(planet) != dict:
        raise Exception('planet input must be a dict')
    
    #experiment must be dict
    if type(experiment) != dict:
        raise Exception('experiment input must be a dict')
    
    
    #Main code

    #interpolation function
    alpha_fun = spi.interp1d(experiment['alpha_dist'], experiment['alpha_deg'], kind='cubic', fill_value='extrapolate')

    #Get the terrain angle - given definition in notes
    terrain_angle = float(alpha_fun(float(y[1])))

    #Get the motor speed
    omega = float(motorW(float(y[0]), rover))

    #Get the net force
    Fnet = F_net(omega, terrain_angle, rover, planet, experiment['Crr'])

    #Get the rover mass
    m = get_mass(rover)

    #Calculate acceleration
    a = float(Fnet / m)

    dydt = np.array([a, float(y[0])])

    return dydt


"""
mechpower - Eimaan

"""
def mechpower(v, rover):
    """
    mechpower - 
    
    Calculates the power output of a single DC motor for each point in a velocity profile.
    
    Inputs:
        v      : 1D numpy array or scalar [m/s] - Rover velocity data from the simulation.
        rover  : dict - Dictionary with rover details, including motor specs.
    
    Outputs:
        P      : 1D numpy array or scalar [W] - The instantaneous power from the motor.
    
    Notes:
    - First input should be a scalar or a 1D numpy array, not a matrix
    - The second input should be a dictionary with motor parameters.
    - Use `tau_dcmotor` for torque and `motorW` for angular velocity
    - Check units
    
    """
    

    # Validate inputs
    if not isinstance(v, (np.ndarray, float, int)):
        raise Exception("First input must be a scalar or a 1D numpy array.")
    
    if isinstance(v, np.ndarray):
        if len(v.shape) != 1:
            raise Exception("First input must be a 1D numpy array, not a matrix.")
    
    if not isinstance(rover, dict):
        raise Exception("Second input must be a dictionary containing rover definition.")
    

    
    # Compute angular velocity using motorW
    omega = motorW(v, rover)
    tau = tau_dcmotor(omega, rover['wheel_assembly']['motor'])
    
    # Compute power as P = torque * angular velocity
    P = tau * omega
  
    return P
 



"""
battenergy - Eimaan

"""


def battenergy(t, v, rover):
    
    
    """
    battenergy - 
    
    Computes how much electrical energy the rover battery uses over time, based on velocity data. 
    This function also takes into account the inefficiencies of converting electrical energy to mechanical energy in the motor.

    Calling Syntax:
    E = battenergy(t, v, rover)

    Inputs:
        t      : 1D numpy array [s] - Time samples from the rover simulation.
        v      : 1D numpy array [m/s] - Rover velocity data from the simulation.
        rover  : dict - Dictionary containing rover details.

    Outputs:
        E      : scalar [J] - Total electrical energy used by the rover's battery pack.

    Notes:
        - The inputs should be equal-length arrays for time and velocity, and the rover input should be a dictionary.
        - The function calls `mechpower` and `tau_dcmotor` to calculate the necessary power and torque.
        - Use interpolation for efficiency based on torque, so make sure the dictionary has that data.
        - Check units
     """
    #validate inputs
    if not (isinstance(t, np.ndarray) and isinstance(v, np.ndarray)):
        raise Exception("First and second inputs must be 1D numpy arrays.")
    if len(t) != len(v):
        raise Exception("Time and velocity arrays must be of equal length.")
    if len(v.shape) != 1:
        raise Exception("Inputs must be 1D arrays")
        
    if not isinstance(rover, dict):
        raise Exception("Third input must be a dictionary containing rover definition.")
    
    
    # mechanical power for a single motor
    P_mech = mechpower(v, rover)
    
    # torque
    omega = motorW(v, rover)
    tau = tau_dcmotor(omega, rover['wheel_assembly']['motor'])
    
    # interpolate  
    efficiency_fun = spi.interp1d(rover['wheel_assembly']['motor']['effcy_tau'], 
                              rover['wheel_assembly']['motor']['effcy'],
                              kind = "cubic")
    efficiency = efficiency_fun(tau)
    # electrical power (efficiency and 6 motors)
    P_elec = P_mech / efficiency
    
    #total energy
    Energy = spitegrate.trapezoid(P_elec, x=t)
    E = Energy * 6
    
    return E









#import end_of_mission_event as eome

def simulate_rover(rover, planet, experiment, end_event):



    """
    simulate_rover - Eddy

    Integrates the trajectory of the rover using input arguments.

    Inputs:
        rover : dict, dictionary containing rover parameters
        planet : dict, dictionary containing planet parameters
        experiment : dict, dictionary containing experiment parameters
        end_event : dict, dictionary containing termination parameters

    Return:
        rover: dict, dictionary containing rover parameters with updated state vector/trajectory

    Notes: *as written*
    This function integrates the trajectory of a rover according to the terrain and initial conditions contained in
    experiment1 (defined in the define_experiment.py file). It uses end_event to define the necessary and sufficient
    conditions to terminate the simulation. Note that you can use the Python ODE solvers. Please determine the best
    solver for the particular problem (Hint: Think about whether it is a stiff system and needs an implicit solver). Note
    that you need to provide the conditions to stop the simulation using the options structure to be fed into the ode
    solver. You will need to download end_of_mission_event function from Canvas and add it to your
    subfunctions.py file. This function requires the end_event structure. The return argument will be used as one of the
    inputs to your ODE solver.
    The function must check validity of inputs. The function should populate the rover[‘telemetry’] dictionary.

    """

    #Input verification

    #a
    if type(rover) != dict:
        raise Exception('First input must be a dict')
    #b
    if type(planet) != dict:
        raise Exception('Second input must be a dict')
    #c
    if type(experiment) != dict:
        raise Exception('Third input must be a dict')
    #d
    if type(end_event) != dict:
        raise Exception('Fourth input must be a dict')
    
    #Main code



    #Define the lambda fun
    fun = lambda t,y: rover_dynamics(t, y, rover, planet, experiment)

    #sol = (diff eq, time span, initial conditions, end condition/event, method-> for stiff system)
    sol = spitegrate.solve_ivp(fun, experiment['time_range'], experiment['initial_conditions'], method = 'BDF')
        #get y = ....
    
    #time
    time = experiment['time_span']

    #completion time
    if max(time) < end_event['max_time']:
        completion_time = max(time)
    else:
        completion_time = end_event['max_time']

    #velocity info
    velocity = sol.y[0]
    max_velocity = max(velocity)
    average_velocity = np.mean(velocity)
    
    #position info
    position = sol.y[1]
    distance_traveled = position[-1]


    #update rover telemetry
    rover['telemetry'] = {
        'time' : time,
        'completion_time' : completion_time,
        'velocity' : velocity,
        'position' : position,
        'distance_traveled' : distance_traveled,
        'max_velocity' : max_velocity,
        'average_velocity' : average_velocity,
        'power' : 1, #Eimaan
        'battery_energy' : 1, #Eimaan
        'energy_per_distance' : 1, #Eimaan
        }

    return rover
