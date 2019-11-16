'''
Created on 14.11.2019

@author: abaktheer
'''


'''
Created on 06.07.2018

@author: abaktheer
'''
import matplotlib.pyplot as plt
import numpy as np


ax1 = plt.subplot(111)
# Smin =0.2

#Exp-IMB

# H-L
n_1 = np.array([0.32 , 1.41, 2.09, 0.48, 0.51, 0.37, 1.3])
s_1 = np.array([0.24,  0.24, 0.24, 0.24, 0.47, 0.47, 0.47 ])
ax1.plot(n_1, s_1, 'ro', markersize=4, color='r')

# L-H
n_1 = np.array([0.09, 0.09, 0.09, 0.26, 0.26])
s_1 = np.array([2.29, 1.2, 3.51, 1.11, 2.33])
ax1.plot(n_1, s_1, 'ro', markersize=4, color='r')

# average
n_1 = np.array([1 , 0])
s_1 = np.array([0,  1])
ax1.plot(n_1, s_1,  color='r',)




#ax1.set_xlim(1, 7)
#ax1.set_ylim(0.6, 0.9)



plt.show()