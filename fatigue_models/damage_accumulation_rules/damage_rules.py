'''
Created on 25.04.2018

@author: abaktheer
'''

import math
from scipy.optimize import newton, bisect
import matplotlib.pyplot as plt
import numpy as np

# P-M rule
x_0 = np.array([0, 1])
y_0 = np.array([0, 1])

# Shah
x_1 = np.linspace(0, 1, 100)
y_1 = 1.14 * x_1**3 - 2.4 * x_1**2 + 2.26 * x_1


# Meyer
x_2 = np.linspace(0, 1, 100)
y_21 = x_2**(1.6 * 0.65)
y_22 = x_2**(1.6 * 0.7)
y_23 = x_2**(1.6 * 0.75)
y_24 = x_2**(1.6 * 0.8)
y_25 = x_2**(1.6 * 0.85)
y_26 = x_2**(1.6 * 0.9)
y_27 = x_2**(1.6 * 0.95)


# DLC
N1 = 179
N2 = 1336

x_30 = np.linspace(0, 1, 100)
x_knee = 1 - 0.65 * (0.5)**0.25
y_knee = 0.35 * (0.5)**0.25

x_3 = np.array([0, x_knee, 1])
y_3 = np.array([0, y_knee, 1])


ax1 = plt.subplot(111)
ax1.plot(x_0, y_0, '--k', linewidth=2, label='P-M')
ax1.plot(x_1, y_1, 'r', label='Shah')
ax1.plot(x_2, y_21, 'g', label='Meyer, S=0.65')
ax1.plot(x_2, y_22, 'g', label='Meyer, S=0.7')
ax1.plot(x_2, y_23, 'g', label='Meyer, S=0.75')
ax1.plot(x_2, y_24, 'g', label='Meyer, S=0.8')
ax1.plot(x_2, y_25, 'g', label='Meyer, S=0.85')
ax1.plot(x_2, y_26, 'g', label='Meyer, S=0.9')
ax1.plot(x_2, y_27, 'g', label='Meyer, S=0.95')
ax1.plot(x_3, y_3, 'b', linewidth=2, label='LCD')

plt.xlabel('cycle ratio (N/Nf)')
plt.ylabel('damage')
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.legend()

plt.show()
