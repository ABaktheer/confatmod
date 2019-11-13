'''
Created on 11.03.2019

@author: abaktheer
'''

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from io import StringIO
    import os
    
    
'''
L-H-1
'''

#=============================================
''' L-H (S= 0.85 -0.75) (CT_80-54) '''
#=============================================
N_1 =  len(np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\Low-High\NPY\CT80-54_63388_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_1 = np.arange(1, N_1 +1, 1)
eps_max_1 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\Low-High\NPY\CT80-54_63388_Zykl_avg_WA_1_WA_2_WA_3_max.npy')


#=============================================
''' L-H  (S= 0.85 -0.75) (CT_80-55) '''
#=============================================
N_2 =  len(np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\Low-High\NPY\CT80-55_56805_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_2 = np.arange(1, N_2 +1, 1)
eps_max_2 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\Low-High\NPY\CT80-55_56805_Zykl_avg_WA_1_WA_2_WA_3_max.npy')


#=============================================
''' L-H  (S= 0.85 -0.75) (CT_80-56) '''
#=============================================
N_3 =  len(np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\Low-High\NPY\CT80-56_70691_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_3 = np.arange(1, N_3 +1, 1)
eps_max_3 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\Low-High\NPY\CT80-56_70691_Zykl_avg_WA_1_WA_2_WA_3_max.npy')




'''
L-H-2
'''
#=============================================
''' L-H  (S= 0.85 -0.75) (CT_80-57) '''
#=============================================
N_4 =  len(np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\Low-High\NPY\CT80-57_155538_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_4 = np.arange(1, N_4 +1, 1)
eps_max_4 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\Low-High\NPY\CT80-57_155538_Zykl_avg_WA_1_WA_2_WA_3_max.npy')

#=============================================
''' L-H  (S= 0.85 -0.75) (CT_80-58) '''
#=============================================
N_5 =  len(np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\Low-High\NPY\CT80-58_162843_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_5 = np.arange(1, N_5 +1, 1)
eps_max_5 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\Low-High\NPY\CT80-58_162843_Zykl_avg_WA_1_WA_2_WA_3_max.npy')




n_01 = 49615
n_02 = 49615
n_03 = 49615

n_04 = 148845
n_05 = 148845


n_hf = 6010.0
n_lf = 572937.0

# normalized accumulative damage

N_01 = np.hstack(
    (N_max_1[0:n_01 + 1] / n_lf, n_01 / n_lf + (N_max_1[n_01 + 1:] - N_max_1[n_01]) / n_hf))
N_02 = np.hstack(
    (N_max_2[0:n_02 + 1] / n_lf, n_02 / n_lf + (N_max_2[n_02 + 1:] - N_max_2[n_02]) / n_hf))
N_03 = np.hstack(
    (N_max_3[0:n_03 + 1] / n_lf, n_03 / n_lf + (N_max_3[n_03 + 1:] - N_max_3[n_03]) / n_hf))

N_04 = np.hstack(
    (N_max_4[0:n_04 + 1] / n_lf, n_04 / n_lf + (N_max_4[n_04 + 1:] - N_max_4[n_04]) / n_hf))
N_05 = np.hstack(
    (N_max_5[0:n_05 + 1] / n_lf, n_05 / n_lf + (N_max_5[n_05 + 1:] - N_max_5[n_05]) / n_hf))


#========================================
# Plotting
#========================================

# plt.subplot(221)
# plt.plot(N_max_1, abs(eps_max_1), 'k')
# plt.plot(N_max_2, abs(eps_max_2), 'k')
# plt.plot(N_max_3, abs(eps_max_3), 'k')
# plt.plot(N_max_4, abs(eps_max_4), 'k')
# plt.plot(N_max_5, abs(eps_max_5), 'k')
# 
# 
# plt.ylim(0.5, 1.2)
# plt.title('Fatigue creep curve (H-L)')
# plt.xlabel('N')
# plt.ylabel('Displacement [mm]')


plt.subplot(111)
plt.plot(N_01[1:], abs(eps_max_1[1:]), "k")
plt.plot(N_02[1:], abs(eps_max_2[1:]), "r")
plt.plot(N_03[1:], abs(eps_max_3[1:]), "b")
plt.plot(N_04[1:], abs(eps_max_4[1:]), "g")
plt.plot(N_05[1:], abs(eps_max_5[1:]), "y")



plt.ylim(0.5, 1.2)
plt.title('Fatigue creep curve normalized (H-L)')
plt.xlabel('N/Nf')
plt.ylabel('Displacement [mm]')


plt.show()
