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



eps_max_high_1 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\CT40-64__000-527112_001_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')
eps_max_low_1 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\CT40-64_527112-1050000_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')
N_1 =  len(eps_max_high_1) + len(eps_max_low_1 ) 
N_max_1 = np.arange(1, N_1 +1, 1)
eps_max_1 = np.hstack((eps_max_high_1 ,  eps_max_low_1))


eps_min_high_11 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\CT40-64__000-527112_001_avg_WA1 mm_WA2 mm_WA3 mm_min.npy')
eps_min_low_11 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\CT40-64_527112-1050000_avg_WA1 mm_WA2 mm_WA3 mm_min.npy')
N_11 =  len(eps_min_high_11) + len(eps_min_low_11 ) 
N_min_1 = np.arange(1, N_11 +1, 1)
eps_min_1 = np.hstack((eps_min_high_11 ,  eps_min_low_11))

plt.subplot(131)
plt.plot(N_max_1[1:], (abs(eps_max_1[1:]))/300, "k")
plt.plot(N_min_1[1:], (abs(eps_min_1[1:]))/300, "k")

plt.ylim(0.0, 0.0045)
plt.xlim(-100, 750000)
plt.title('Fatigue creep curve normalized (H-L)')
plt.xlabel('N/Nf')
plt.ylabel('Displacement [mm]')

#========================================================================
"#C80"

# u_11 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 1\NPY\CT_80-5-86Zyk_avg_WA_1_WA_2_WA_3.npy')
# F_11 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 1\NPY\CT_80-5-86Zyk_Kraft.npy')
# 
# 
# plt.subplot(334)
# plt.plot((abs(u_11))/300, abs(F_11 *1000)/(np.pi * 100**2 /4), "k")
# 
# 
# plt.xlim(0.0, 0.004)
# plt.ylim(0.0, 120)
# plt.xlabel('eps')
# plt.ylabel('sigma')


#creep-fatigue
N_11 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 1\NPY\CT_80-7-469110Zyk_g_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_11 = np.arange(1, N_11 +1, 1)
eps_max_11 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 1\NPY\CT_80-7-469110Zyk_g_avg_WA_1_WA_2_WA_3_max.npy')
eps_min_11 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 1\NPY\CT_80-7-469110Zyk_g_avg_WA_1_WA_2_WA_3_min.npy')


plt.subplot(132)
plt.plot(N_max_11[1:], (abs(eps_max_11[1:])-0.21)/300, "k")
plt.plot(N_max_11[1:-1], (abs(eps_min_11[1:-1])-0.21)/300, "k")


plt.ylim(0.0, 0.0045)
plt.xlim(-100, 750000)
plt.title('Fatigue creep curve normalized (H-L)')
plt.xlabel('N/Nf')
plt.ylabel('Displacement [mm]')




#==========================================================================================
"#C120"
# #F-U
# u_111 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\CT_120_1_19_LS3_avg_WA_1_WA_2_WA_3.npy')
# F_111 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\CT_120_1_19_LS3_Kraft.npy')
# 
# 
# plt.subplot(331)
# #plt.plot((abs(u_111))/300, abs(F_111 *1000)/(np.pi * 100**2 /4), "k")
# 
# 
# plt.xlim(0.0, 0.004)
# plt.ylim(0.0, 130)
# plt.xlabel('eps')
# plt.ylabel('sigma')


#creep-fatigue
N_111 =  len(np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\CT_120_1_19_LS3_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_111 = np.arange(1, N_111 +1, 1)
eps_max_111 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\CT_120_1_19_LS3_avg_WA_1_WA_2_WA_3_max.npy')
eps_min_111 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\CT_120_1_19_LS3_avg_WA_1_WA_2_WA_3_min.npy')



plt.subplot(133)
plt.plot(N_max_111[1:] , (abs(eps_max_111[1:]))/300, "k")
plt.plot(N_max_111[1:], (abs(eps_min_111[1:]))/300, "k")


plt.ylim(0.0, 0.0045)
plt.xlim(-100, 750000)
plt.title('Fatigue creep curve normalized (H-L)')
plt.xlabel('N/Nf')
plt.ylabel('Displacement [mm]')








plt.show()