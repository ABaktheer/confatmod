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
# Points
n_1 = np.array([6322, 13713, 3610, 10949, 818, 647, ])
s_1 = np.array([0.85,  0.85, 0.85, 0.85, 0.85, 0.85 ])
ax1.plot(np.log10(n_1), s_1, 'ro', markersize=4, color='r')

n_1 = np.array([74190, 416585, 253449, 142547, 538303, 2012546 ])
s_1 = np.array([0.75,  0.75, 0.75, 0.75, 0.75, 0.75 ])
ax1.plot(np.log10(n_1), s_1, 'ro', markersize=4, color='r')

# average
n_1 = np.array([6010 , 572937])
s_1 = np.array([0.85,  0.75])
ax1.plot(np.log10(n_1), s_1,  color='r',)


# FIB model code 2010
n_11 = np.array([32064472260,    2853787864,    22605643,    2011937,    179066,    15937,    1418,    126,    11,    1,])
s_11 = np.array([0.5,    0.55,    0.65,    0.7,    0.75,    0.8,    0.85,    0.9,    0.95,    1])
ax1.plot(np.log10(n_11), s_11, color='k')


ax1.set_xlim(1, 7)
ax1.set_ylim(0.6, 0.9)



plt.show()
