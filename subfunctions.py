

'''
all necessary functions for Part 1 of Project for MEEN 357 - Spring 2025
'''

# Constants for Mars Rover Calculations
# Nested dictionary of the rover
rover = {
    'wheel_assembly': {
        'wheel': {
            'radius': 0.30,  # m
            'mass': 1.0      # kg
        },
        'speed_reducer': {
            'type': 'reverted',
            'diam_pinion': 0.04,  # m
            'diam_gear': 0.07,    # m
            'mass': 1.5           # kg
        },
        'motor': {
            'torque_stall': 170,    # N·m
            'torque_noload': 0,     # N·m
            'speed_noload': 3.80,   # rad/s
            'mass': 5.0             # kg
        }
    },
    'chassis': {
        'mass': 659  # kg
    },
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



# def get_gear_ratio(speed_reducer):

#     return Ng

# def get_mass(rover):

#     return m

# def F_rolling(omega, terrain_angle, rover, planet, Crr):
#     return

# def F_gravity(terrain_angle, rover, planet):
#     return Fgt

# def F_drive(omega, rover):
#     return Fd

# def F_net(omega, terrain_angle, rover, planet, Crr):
#     return Fnet
