

'''
all necessary functions for Part 1 of Project for MEEN 357 - Spring 2025
'''

# Constants for Mars Rover Calculations
# Nested dictionary of the rover
speed_reducer= {
    'type': 'reverted',
    'diam_pinion': 0.04,  # m
    'diam_gear': 0.07,    # m
    'mass': 1.5           # kg
}

wheel= {
    'radius': 0.30,  # m
    'mass': 1.0      # kg
}

motor= {
    'torque_stall': 170,    # N·m
    'torque_noload': 0,     # N·m
    'speed_noload': 3.80,   # rad/s
    'mass': 5.0             # kg
}

chassis= {
    'mass': 659  # kg
}

wheel_assembly = {
    'wheel': wheel,
    'speed_reducer': speed_reducer,
    'motor': motor
}


rover = {
    'wheel_assembly': wheel_assembly
    ,
    'chassis': chassis
    ,
    'science_payload': {
        'mass': 75  # kg
    },
    'power_subsys': {
        'mass': 90  # kg (RTG - Radioisotope Thermoelectric Generator)
    }
}



# Planet constants
planet = {
    'g': 3.72  # m/s^2 (acceleration due to gravity)
}




import numpy as np
from scipy import erf


def tau_dcmotor(omega, motor):
    torque_stall = rover['wheel_assembly']['motor']['torque_stall']
    tauNoLoad = rover['wheel_assembly']['motor']['torque_noload']
    omegaNoLoad = rover['wheel_assembly']['motor']['speed_noload']
    #checking for errors in input
    if (type(omega) is not int) and (type(omega) is not isinstance(motor, np.ndarray)):
        raise Exception('<omega is not vector or scalar>')
    if type(motor) is not dict:
        raise Exception('<Motor input is not a dict>')
    #verifying edge/corner cases
    if omega > omegaNoLoad:
        tau = 0
    if omega <= 0:
        tau = torque_stall
    else:
        #calculating tau from 2.1 formula
        tau = torque_stall - ((torque_stall - tauNoLoad)/omegaNoLoad) * omega
    return tau



def get_gear_ratio(speed_reducer):
    if type(speed_reducer) is not dict:
        raise Exception('<Speed reducer input is not a dict>')
    
    pinion = rover['wheel_assembly']['speed_reducer']['diam_pinion']
    gear = rover['wheel_assembly']['speed_reducer']['diam_gear']
    
    typeo = rover['wheel_assembly']['speed_reducer']['type']
    
    if typeo != 'reverted':
        raise Exception('<Invalid type called')
    
    Ng = (gear/pinion)**2

    return Ng

def get_mass(rover):

    '''
This function computes rover mass in kilograms. It accounts for the chassis, power subsystem, science payload,
and six wheel assemblies, which itself is comprised of a motor, speed reducer, and the wheel itself.
    '''
    
    if type(rover) is not dict:
        raise Exception('<Rover input is not a dict>')
    
    m_chassey = rover['chassis']['mass']
    m_payload = rover['science_payload']['mass']
    m_power = rover['power_subsys']['mass']
    
    m_wheel = rover['wheel_assembly']['wheel']['mass']  
    m_reducer = rover['wheel_assembly']['speed_reducer']['mass'] 
    m_motor = rover['wheel_assembly']['motor']['mass']     
    
    m_wheelass = m_wheel + m_reducer + m_motor
    
    m = 6 * m_wheelass + m_chassey + m_power + m_payload

    return m

def F_rolling(omega, terrain_angle, rover, planet, Crr):
    
    #a/c - verify scalars or vectors of same size and dicts
    
    if (type(omega) is not int) and (type(omega) is not isinstance(omega, np.ndarray)):
        raise Exception('<omega is not vector or scalar>')      
    if (type(terrain_angle) is not int) and (type(terrain_angle) is not isinstance(terrain_angle, np.ndarray)):
        raise Exception('<terrain_angle is not vector or scalar>')  
    if type(rover) is not dict:
        raise Exception('<Rover input is not a dict>')
    if type(planet) is not dict:
        raise Exception('<Planet input is not a dict>')        
    if (type(Crr) is not int) and Crr < 0 :
        raise Exception('<Crr is not a positive scalar>')
    
    #b - verify bounds
    
    if terrain_angle < -75 and terrain_angle > 75:
        raise Exception('<Terrain_angle not in bounds')    
    
    Fn = get_mass(rover)*planet['g']*np.cos(np.degrees(terrain_angle))
    
    Frr_simple = Crr*Fn
    
    omega_out = omega/get_gear_ratio(speed_reducer)
    
    v_rover = rover['wheel_assembly']['wheel']['radius']*omega_out
    
    Frr = erf(40*v_rover)*Frr_simple
    
    return Frr


def F_gravity(terrain_angle, rover, planet):
    # validate first input is scalar or vector
    if not isinstance(terrain_angle, (int, float, np.ndarray)):
        raise Exception('First input must be a scalar or vector.')

    # all elements of the first argument are between -75 and +75 degrees
    terrain_angle_array = np.asarray(terrain_angle)
    if np.any((terrain_angle_array < -75) | (terrain_angle_array > 75)):
        raise Exception('<Terrain_angle not in bounds>')

    # validate last two inputs are dictionaries
    if not isinstance(rover, dict):
        raise Exception('Second input must be a dictionary.')
    if not isinstance(planet, dict):
        raise Exception('Third input must be a dictionary.')

    #calculate Fg
    mass = get_mass(rover)  # Call get mass
    angle_rad = np.deg2rad(terrain_angle)  # angle to radians
    Fgt = -mass * planet['g'] * np.sin(angle_rad)  # Apply formula

    return Fgt



def F_drive(omega, rover):
    #validate first two inputs are scalars or vectors
    if (type(omega) is not int) and (not isinstance(omega, np.ndarray)):
        raise Exception('<omega is not vector or scalar>')
    
    # validate second input is a dict
    if not isinstance(rover, dict): 
        raise Exception('Second input must be a dictionary.')

    # vectorization
    omega = np.atleast_1d(omega)

    # Compute the drive force,
    Fd = 6 * get_gear_ratio(rover) * tau_dcmotor(omega, motor) / rover['wheel_radius']

    return Fd


def F_net(omega, terrain_angle, rover, planet, Crr):
    # validate thatfirst two inputs are scalars or vectors
    if (type(omega) is not int) and (not isinstance(omega, np.ndarray)):
        raise Exception('<omega is not vector or scalar>')
    if (type(terrain_angle) is not int) and (not isinstance(terrain_angle, np.ndarray)):
        raise Exception('<terrain angle is not vector or scalar>')

    # validate that both are vectors and have the same size
    if isinstance(omega, np.ndarray) and isinstance(terrain_angle, np.ndarray):
        if omega.shape != terrain_angle.shape:
            raise Exception('Omega and terrain angle must be vectors of the same size.')

    # validate all elements of terrain_angle are between -75 and +75 degrees
    if np.any((np.asarray(terrain_angle) < -75) | (np.asarray(terrain_angle) > 75)):
        raise Exception('<Terrain angle not in bounds>')

    # validate rover and planet are dictionaries
    if type(rover) is not dict:
        raise Exception('<Rover input is not a dict>')
    if type(planet) is not dict:
        raise Exception('<Planet input is not a dict>')

    #Validate Crr is a positive scalar
    if (type(Crr) is not int) and Crr < 0:
        raise Exception('<Crr is not a positive scalar>')

    # calc forces
    Fd = F_drive(omega, rover)                 
    Fg = F_gravity(terrain_angle, rover, planet) 
    Fr = F_rolling(terrain_angle, rover, planet, Crr)  

    #calc net force
    Fnet = Fd + Fg - Fr  

    return Fnet

