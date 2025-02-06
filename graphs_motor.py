
import matplotlib.pyplot as plt
import numpy as np
from subfunctions import tau_dcmotor 


motor= {
    'torque_stall': 170,    # N·m
    'torque_noload': 0,     # N·m
    'speed_noload': 3.80,   # rad/s
    'mass': 5.0             # kg
}

#torque range for motor
torque = np.linspace(0, motor['torque_stall'], 100)  # [Nm]

#calc motor shaft speed with tau_dcmotor 
speed = tau_dcmotor(torque, motor)  # [rad/s]

#calculate motor power
power = torque * speed  # [W]
# Plotting
plt.figure(figsize=(8, 12))

#motor shaft speed vs. motor shaft torque
plt.subplot(3, 1, 1)
plt.plot(torque, speed)
plt.xlabel('Motor Shaft Torque [Nm]')
plt.ylabel('Motor Shaft Speed [rad/s]')
plt.grid(True)

# motor power vs. motor shaft torque
plt.subplot(3, 1, 2)
plt.plot(torque, power)
plt.xlabel('Motor Shaft Torque [Nm]')
plt.ylabel('Motor Power [W]')
plt.grid(True)

# motor power vs. motor shaft speed
plt.subplot(3, 1, 3)
plt.plot(speed, power)
plt.xlabel('Motor Shaft Speed [rad/s]')
plt.ylabel('Motor Power [W]')
plt.grid(True)


plt.tight_layout()

