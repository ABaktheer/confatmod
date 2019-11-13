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
H-L-1
'''

#=============================================
''' H-L (S= 0.85 -0.75) (CT_80-49) '''
#=============================================
N_1 =  len(np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\High-Low\NPY\CT80-49_185821_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_1 = np.arange(1, N_1 +1, 1)
eps_max_1 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\High-Low\NPY\CT80-49_185821_Zykl_avg_WA_1_WA_2_WA_3_max.npy')


#=============================================
''' H-L (S= 0.85 -0.75) (CT_80-50) '''
#=============================================
N_2 =  len(np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\High-Low\NPY\CT80-50_810229_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_2 = np.arange(1, N_2 +1, 1)
eps_max_2 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\High-Low\NPY\CT80-50_810229_Zykl_avg_WA_1_WA_2_WA_3_max.npy')


#=============================================
''' H-L (S= 0.85 -0.75) (CT_80-51) '''
#=============================================
N_3 =  len(np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\High-Low\NPY\CT80-51_1200277_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_3 = np.arange(1, N_3 +1, 1)
eps_max_3 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\High-Low\NPY\CT80-51_1200277_Zykl_avg_WA_1_WA_2_WA_3_max.npy')

#=============================================
''' H-L (S= 0.85 -0.75) (CT_80-41) '''
#=============================================
N_4 =  len(np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\High-Low\NPY\CT80-41_274256_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_4 = np.arange(1, N_4 +1, 1)
eps_max_4 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\High-Low\NPY\CT80-41_274256_Zykl_avg_WA_1_WA_2_WA_3_max.npy')


'''
H-L-2
'''
#=============================================
''' H-L (S= 0.85 -0.75) (CT_80-53) '''
#=============================================
N_5 =  len(np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\High-Low\NPY\CT80-53_292317_Zykl2_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_5 = np.arange(1, N_5 +1, 1)
eps_max_5 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\High-Low\NPY\CT80-53_292317_Zykl2_avg_WA_1_WA_2_WA_3_max.npy')

#=============================================
''' H-L (S= 0.85 -0.75) (CT_80-62) '''
#=============================================
N_6 =  len(np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\High-Low\NPY\CT80-62_292317_Zykl2_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_6 = np.arange(1, N_6 +1, 1)
eps_max_6 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\High-Low\NPY\CT80-62_292317_Zykl2_avg_WA_1_WA_2_WA_3_max.npy')

#=============================================
''' H-L (S= 0.85 -0.75) (CT_80-63) '''
#=============================================
N_7 =  len(np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\High-Low\NPY\CT80-63_747618_Zykl2_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_7 = np.arange(1, N_7 +1, 1)
eps_max_7 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\High-Low\NPY\CT80-63_747618_Zykl2_avg_WA_1_WA_2_WA_3_max.npy')



n_01 = 1416
n_02 = 1416
n_03 = 1416
n_04 = 1416
n_05 = 2833
n_06 = 2833
n_07 = 2833

n_hf = 6010.0
n_lf = 572937.0

# normalized accumulative damage
N_01 = np.hstack(
    (N_max_1[0:n_01 + 1] / n_hf, n_01 / n_hf + (N_max_1[n_01 + 1:] - N_max_1[n_01]) / n_lf))
N_02 = np.hstack(
    (N_max_3[0:n_02 + 1] / n_hf, n_02 / n_hf + (N_max_2[n_02 + 1:] - N_max_2[n_02]) / n_lf))
N_03 = np.hstack(
    (N_max_3[0:n_03 + 1] / n_hf, n_03 / n_hf + (N_max_3[n_03 + 1:] - N_max_3[n_03]) / n_lf))
N_04 = np.hstack(
    (N_max_4[0:n_04 + 1] / n_hf, n_04 / n_hf + (N_max_4[n_04 + 1:] - N_max_4[n_04]) / n_lf))

N_05 = np.hstack(
    (N_max_5[0:n_05 + 1] / n_hf, n_05 / n_hf + (N_max_5[n_05 + 1:] - N_max_5[n_05]) / n_lf))
N_06 = np.hstack(
    (N_max_6[0:n_06 + 1] / n_hf, n_06 / n_hf + (N_max_6[n_06 + 1:] - N_max_6[n_06]) / n_lf))
N_07 = np.hstack(
    (N_max_7[0:n_07 + 1] / n_hf, n_07 / n_hf + (N_max_7[n_07 + 1:] - N_max_7[n_07]) / n_lf))

#========================================
# Plotting
#========================================

# plt.subplot(221)
# plt.plot(N_max_1, abs(eps_max_1), 'k')
# plt.plot(N_max_2, abs(eps_max_2), 'k')
# plt.plot(N_max_3, abs(eps_max_3), 'k')
# plt.plot(N_max_4, abs(eps_max_4), 'k')
# plt.plot(N_max_5, abs(eps_max_5), 'k')
# plt.plot(N_max_6, abs(eps_max_6), 'k')
# plt.plot(N_max_7, abs(eps_max_7), 'k')
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
plt.plot(N_06[1:], abs(eps_max_6[1:]), "m")
plt.plot(N_07[1:], abs(eps_max_7[1:]), "-k")


plt.ylim(0.5, 1.2)
plt.title('Fatigue creep curve normalized (H-L)')
plt.xlabel('N/Nf')
plt.ylabel('Displacement [mm]')


plt.show()
