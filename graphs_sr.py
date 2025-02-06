import matplotlib.pyplot as plt
import numpy as np
from subfunctions import tau_dcmotor 
from subfunctions import get_gear_ratio  


motor= {
    'torque_stall': 170,    # N·m
    'torque_noload': 0,     # N·m
    'speed_noload': 3.80,   # rad/s
    'mass': 5.0             # kg
}

speed_reducer= {
    'type': 'reverted',
    'diam_pinion': 0.04,  # m
    'diam_gear': 0.07,    # m
    'mass': 1.5           # kg
}

#gear ratio
Ng = get_gear_ratio(speed_reducer)

#define torque range for motor
motor_torque = np.linspace(0, motor['torque_stall'], 100)  # [Nm]

#calculate motor shaft speed using tau_dcmotor function
motor_speed = tau_dcmotor(motor_torque, motor)  # [rad/s]

#calculate motor power
motor_power = motor_torque * motor_speed  # [W]

#calc speed reducer output values
sr_torque = motor_torque * Ng  # Output torque
sr_speed = motor_speed / Ng  # Output speed
sr_power = sr_torque * sr_speed  # Output power

#plot
plt.figure(figsize=(8, 12))

# speed reducer output shaft speed vs. torque
plt.subplot(3, 1, 1)
plt.plot(sr_torque, sr_speed)
plt.xlabel('Speed Reducer Output Shaft Torque [Nm]')
plt.ylabel('Speed Reducer Output Shaft Speed [rad/s]')
plt.grid(True)

#speed reducer output power vs. torque
plt.subplot(3, 1, 2)
plt.plot(sr_torque, sr_power)
plt.xlabel('Speed Reducer Output Shaft Torque [Nm]')
plt.ylabel('Speed Reducer Output Power [W]')
plt.grid(True)

# speed reducer output power vs. speed
plt.subplot(3, 1, 3)
plt.plot(sr_speed, sr_power)
plt.xlabel('Speed Reducer Output Shaft Speed [rad/s]')
plt.ylabel('Speed Reducer Output Power [W]')
plt.grid(True)


plt.tight_layout()
plt.show()


