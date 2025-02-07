import matplotlib.pyplot as plt
import numpy as np
from subfunctions import tau_dcmotor, get_gear_ratio

# def motor and speed reducer
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

# gear ratio
Ng = get_gear_ratio(speed_reducer)

# motor shaft speed range
motor_speed = np.linspace(0, motor['speed_noload'], 100)  # [rad/s]

# motor torque
motor_torque = tau_dcmotor(motor_speed, motor)  # [Nm]

# calc speed reducer output shaft torque and speed
sr_torque = motor_torque * Ng
sr_speed = motor_speed / Ng

# calc speed reducer power
sr_power = sr_torque * sr_speed  # [W]

# plot
plt.figure(figsize=(8, 12))

# speed reducer output shaft speed vs. torque
plt.subplot(3, 1, 1)
plt.plot(sr_torque, sr_speed)
plt.xlabel('Speed Reducer Output Shaft Torque [Nm]')
plt.ylabel('Speed Reducer Output Shaft Speed [rad/s]')
plt.grid(True)

# speed reducer output power vs. torque
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

#plt.show()

