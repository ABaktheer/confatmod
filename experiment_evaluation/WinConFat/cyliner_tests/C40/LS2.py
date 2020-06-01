

import matplotlib.pyplot as plt
import numpy as np


u_1 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-40_00-100Zykl_001_avg_WA2 mm_WA3 mm.npy')
F_1 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-40_00-100Zykl_001_Kraft kN.npy')




plt.subplot(221)
plt.plot((abs(u_1)- 0.038)/300, abs(F_1 *1000)/(np.pi * 152**2 /4), "k")


plt.xlim(0.0, 0.0025)
plt.xlabel('eps')
plt.ylabel('sigma')




N_1 =  len(np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-40_00-100Zykl_001_avg_WA2 mm_WA3 mm_max.npy'))
N_max_1 = np.arange(1, N_1 +1, 1)
eps_max_1 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-40_00-100Zykl_001_avg_WA2 mm_WA3 mm_max.npy')
eps_min_1 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-40_00-100Zykl_001_avg_WA2 mm_WA3 mm_min.npy')

plt.subplot(222)
plt.plot(N_max_1[1:] *100/N_1, (abs(eps_max_1[1:])-0.038)/300, "k")
plt.plot(N_max_1[1:]*100/N_1, (abs(eps_min_1[1:])-0.038)/300, "k")

plt.ylim(0.0, 0.003)
plt.title('Fatigue creep curve normalized (H-L)')
plt.xlabel('N/Nf')
plt.ylabel('Displacement [mm]')


plt.show()