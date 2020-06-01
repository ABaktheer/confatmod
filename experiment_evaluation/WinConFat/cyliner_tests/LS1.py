'''
Created on 27.03.2019

@author: abaktheer
'''


import matplotlib.pyplot as plt
import numpy as np

#===================================
"LS1"
#===================================

"#C40"

u_1 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-50_WA-Mittelwert.npy')
F_1 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-50_Kraft.npy')


u_2 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-51_WA-Mittelwert.npy')
F_2 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-51_Kraft.npy')


u_3 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-52_WA-Mittelwert.npy')
F_3 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-52_Kraft.npy')

plt.subplot(231)
plt.plot((abs(u_1))/300, abs(F_1 *1000)/(np.pi * 150**2 /4), "k")
plt.plot((abs(u_2))/300, abs(F_2 *1000)/(np.pi * 150**2 /4), "k")
plt.plot((abs(u_3))/300, abs(F_3 *1000)/(np.pi * 150**2 /4), "k")

plt.xlim(0.0, 0.008)
plt.ylim(0.0, 130)
plt.xlabel('eps')
plt.ylabel('sigma')

#========================================================================

"#C80"

u_11 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\NPY\CT80-33_avg_WA_1_WA_2_WA_3.npy')
F_11 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\NPY\CT80-33_Kraft.npy')

u_22 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\NPY\CT80-34_avg_WA_1_WA_2_WA_3.npy')
F_22 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\NPY\CT80-34_Kraft.npy')

u_33 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\NPY\CT80-35_avg_WA_1_WA_2_WA_3.npy')
F_33 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\NPY\CT80-35_Kraft.npy')


u_44 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\NPY\CT80-36_avg_WA_1_WA_2_WA_3.npy')
F_44 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\NPY\CT80-36_Kraft.npy')

u_55 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\NPY\CT80-37_avg_WA_1_WA_2_WA_3.npy')
F_55 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\NPY\CT80-37_Kraft.npy')


u_66 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\NPY\CT80-38_avg_WA_1_WA_2_WA_3.npy')
F_66 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\NPY\CT80-38_Kraft.npy')


plt.subplot(232)
plt.plot((abs(u_11))/300, abs(F_11 *1000)/(np.pi * 100**2 /4), "k")
plt.plot((abs(u_22))/300, abs(F_22 *1000)/(np.pi * 100**2 /4), "k")
plt.plot((abs(u_33))/300, abs(F_33 *1000)/(np.pi * 100**2 /4), "k")
plt.plot((abs(u_44))/300, abs(F_44 *1000)/(np.pi * 100**2 /4), "k")
plt.plot((abs(u_55))/300, abs(F_55 *1000)/(np.pi * 100**2 /4), "k")
plt.plot((abs(u_66))/300, abs(F_66 *1000)/(np.pi * 100**2 /4), "k")

plt.xlim(0.0, 0.008)
plt.ylim(0.0, 130)
plt.xlabel('eps')
plt.ylabel('sigma')







#===========================================================
"#C120"
u_111 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\CT_120_1_13_LS1_avg_WA_1_WA_2_WA_3.npy')
F_111 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\CT_120_1_13_LS1_Kraft.npy')

u_222 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\CT_120_1_14_LS1_avg_WA_1_WA_2_WA_3.npy')
F_222 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\CT_120_1_14_LS1_Kraft.npy')


u_333 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\CT_120_1_15_LS1_avg_WA_1_WA_2_WA_3.npy')
F_333 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\CT_120_1_15_LS1_Kraft.npy')

plt.subplot(233)
plt.plot((abs(u_111))/300, abs(F_111 *1000)/(np.pi * 100**2 /4), "k")
plt.plot((abs(u_222))/300, abs(F_222 *1000)/(np.pi * 100**2 /4), "k")
plt.plot((abs(u_333))/300, abs(F_333 *1000)/(np.pi * 100**2 /4), "k")

plt.xlim(0.0, 0.008)
plt.ylim(0.0, 130)
plt.xlabel('eps')
plt.ylabel('sigma')










plt.show()