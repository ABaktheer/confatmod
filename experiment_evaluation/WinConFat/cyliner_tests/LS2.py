'''
Created on 27.03.2019

@author: abaktheer
'''


import matplotlib.pyplot as plt
import numpy as np

#===================================
"LS2"
#===================================

"#C40"
#F-U
u_1 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-53_WA-Mittelwert.npy')
F_1 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-53_Kraft.npy')


plt.subplot(331)
plt.plot((abs(u_1))/300, abs(F_1 *1000)/(np.pi * 150**2 /4), "k")

plt.xlim(0.0, 0.008)
plt.ylim(0.0, 80)
plt.xlabel('eps')
plt.ylabel('sigma')


eps = np.array([0,   0.0031,    0.0045,    0.0054,    0.0063])
eps_p = np.array([0,    0.0015,    0.0025,    0.0035,    0.0043])
E = np.array([31947.80546,    30202.79572,    14579.20851,    11622.73358,    8495.066223])
omega = np.array([0,    0.054620645,    0.5436554,    0.636196183,    0.734095469])

plt.subplot(334)
plt.plot(eps, eps_p, "ko--",  markersize=3)
plt.xlim(0.0, 0.008)
plt.ylim(0.0, 0.008)


plt.subplot(337)
plt.plot(eps, omega , "ro--",  markersize=3)
plt.xlim(0.0, 0.008)
plt.ylim(0, 1)

#===========================================================================
u_2 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-54_WA-Mittelwert.npy')
F_2 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-54_Kraft.npy')


plt.subplot(332)
plt.plot((abs(u_2))/300, abs(F_2 *1000)/(np.pi * 150**2 /4), "r")

plt.xlim(0.0, 0.008)
plt.ylim(0.0, 80)
plt.xlabel('eps')
plt.ylabel('sigma')

eps = np.array([0,    0.003074627,    0.003820895,    0.004726368,    0.005671642])
eps_p = np.array([0,    0.001293532,    0.002059701,    0.003124378,    0.004099502])
E = np.array([30534.54621,    29865.7682,   29690.01138,    21112.90841,    17448.3067])
omega = np.array([0,    0.02190234,    0.027658339,    0.308556667,    0.428571606])

plt.subplot(335)
plt.plot(eps, eps_p, "ko--",  markersize=3)
plt.xlim(0.0, 0.008)
plt.ylim(0.0, 0.008)


plt.subplot(338)
plt.plot(eps, omega , "ro--",  markersize=3)
plt.xlim(0.0, 0.008)
plt.ylim(0, 1)


#=============================================================================
u_3 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-55_WA-Mittelwert.npy')
F_3 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-55_Kraft.npy')


plt.subplot(333)
plt.plot((abs(u_3))/300, abs(F_3 *1000)/(np.pi * 150**2 /4), "g")

plt.xlim(0.0, 0.008)
plt.ylim(0.0, 80)
plt.xlabel('eps')
plt.ylabel('sigma')


eps = np.array([0,    0.003111753,    0.003921002,    0.004816955,    0.005732177])
eps_p = np.array([0,    0.00123314,    0.002138728,    0.003236994,    0.003689788])
E = np.array([29164.76307,    29083.41124,    28991.31572,    20749.62542,    12848.80997])
omega = np.array([0,    0.047524367,    0.050540476,    0.320454109,    0.579204162])

plt.subplot(336)
plt.plot(eps, eps_p, "ko--",  markersize=3)
plt.xlim(0.0, 0.008)
plt.ylim(0.0, 0.008)


plt.subplot(339)
plt.plot(eps, omega , "ro--",  markersize=3)
plt.xlim(0.0, 0.008)
plt.ylim(0, 1)



plt.show()