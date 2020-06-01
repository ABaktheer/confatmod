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
u_1 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-40_00-100Zykl_001_avg_WA2 mm_WA3 mm.npy')
F_1 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-40_00-100Zykl_001_Kraft kN.npy')

u_2 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\CT40-41_00-100Zykl_001_001_avg_WA1 mm_WA2 mm_WA3 mm.npy')
F_2 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\CT40-41_00-100Zykl_001_001_Kraft kN.npy')

plt.subplot(231)
plt.plot((abs(u_1)-0.042)/300, abs(F_1 *1000)/(np.pi * 150**2 /4), "k")
plt.plot((abs(u_2)-0.06)/300, abs(F_2 *1000)/(np.pi * 150**2 /4), "r")

plt.xlim(0.0, 0.004)
plt.ylim(0.0, 60)
plt.xlabel('eps')
plt.ylabel('sigma')


#creep-fatigue
N_1 =  len(np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-40_00-100Zykl_001_avg_WA2 mm_WA3 mm_max.npy'))
N_max_1 = np.arange(1, N_1 +1, 1)
eps_max_1 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-40_00-100Zykl_001_avg_WA2 mm_WA3 mm_max.npy')
eps_min_1 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-40_00-100Zykl_001_avg_WA2 mm_WA3 mm_min.npy')


#creep-fatigue
N_2 =  len(np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\CT40-41_00-100Zykl_001_001_avg_WA1 mm_WA2 mm_WA3 mm_max.npy'))
N_max_2 = np.arange(1, N_2 +1, 1)
eps_max_2 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\CT40-41_00-100Zykl_001_001_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')
eps_min_2 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\CT40-41_00-100Zykl_001_001_avg_WA1 mm_WA2 mm_WA3 mm_min.npy')


plt.subplot(234)
plt.plot(N_max_1[1:] *100/N_1, (abs(eps_max_1[1:])-0.042)/300, "k")
plt.plot(N_max_1[1:]*100/N_1, (abs(eps_min_1[1:])-0.042)/300, "k")

plt.plot(N_max_2[1:] *100/N_2, (abs(eps_max_2[1:])-0.06)/300, "r")
plt.plot(N_max_2[1:-1]*100/N_2, (abs(eps_min_2[1:])-0.06)/300, "r")

plt.ylim(0.0, 0.0033)
plt.title('Fatigue creep curve normalized (H-L)')
plt.xlabel('N/Nf')
plt.ylabel('Displacement [mm]')


#========================================================================
"#C80"

u_11 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 1\NPY\CT_80-5-86Zyk_avg_WA_1_WA_2_WA_3.npy')
F_11 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 1\NPY\CT_80-5-86Zyk_Kraft.npy')

u_22 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 1\NPY\CT_80-6-86Zyk_avg_WA_1_WA_2_WA_3.npy')
F_22 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 1\NPY\CT_80-6-86Zyk_Kraft.npy')


u_33 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 1\NPY\CT_80-8-91Zyk_avg_WA_1_WA_2_WA_3.npy')
F_33 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 1\NPY\CT_80-8-91Zyk_Kraft.npy')


plt.subplot(232)
plt.plot((abs(u_11))/300, abs(F_11 *1000)/(np.pi * 100**2 /4), "k")
plt.plot((abs(u_22))/300, abs(F_22 *1000)/(np.pi * 100**2 /4), "r")
plt.plot((abs(u_33))/300, abs(F_33 *1000)/(np.pi * 100**2 /4), "g")

plt.xlim(0.0, 0.004)
plt.ylim(0.0, 120)
plt.xlabel('eps')
plt.ylabel('sigma')


#creep-fatigue
N_11 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 1\NPY\CT_80-5-86Zyk_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_11 = np.arange(1, N_11 +1, 1)
eps_max_11 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 1\NPY\CT_80-5-86Zyk_avg_WA_1_WA_2_WA_3_max.npy')
eps_min_11 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 1\NPY\CT_80-5-86Zyk_avg_WA_1_WA_2_WA_3_min.npy')

N_22 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 1\NPY\CT_80-6-86Zyk_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_22= np.arange(1, N_22 +1, 1)
eps_max_22 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 1\NPY\CT_80-6-86Zyk_avg_WA_1_WA_2_WA_3_max.npy')
eps_min_22 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 1\NPY\CT_80-6-86Zyk_avg_WA_1_WA_2_WA_3_min.npy')


N_33 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 1\NPY\CT_80-8-91Zyk_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_33= np.arange(1, N_33 +1, 1)
eps_max_33 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 1\NPY\CT_80-8-91Zyk_avg_WA_1_WA_2_WA_3_max.npy')
eps_min_33 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 1\NPY\CT_80-8-91Zyk_avg_WA_1_WA_2_WA_3_min.npy')


plt.subplot(235)
plt.plot(N_max_11[1:] *100/N_11, (abs(eps_max_11[1:]))/300, "k")
plt.plot(N_max_11[1:-1]*100/N_11, (abs(eps_min_11[1:-1]))/300, "k")

plt.plot(N_max_22[1:] *100/N_22, (abs(eps_max_22[1:]))/300, "k")
plt.plot(N_max_22[1:]*100/N_22, (abs(eps_min_22[1:]))/300, "k")

plt.plot(N_max_33[1:] *100/N_33, (abs(eps_max_33[1:]))/300, "k")
plt.plot(N_max_33[1:-1]*100/N_33, (abs(eps_min_33[1:]))/300, "k")

plt.ylim(0.0, 0.004)
plt.title('Fatigue creep curve normalized (H-L)')
plt.xlabel('N/Nf')
plt.ylabel('Displacement [mm]')




#==========================================================================================
"#C120"
#F-U
u_111 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\CT_120_1_16_LS2_avg_WA_1_WA_2_WA_3.npy')
F_111 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\CT_120_1_16_LS2_Kraft.npy')

u_222 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\CT_120_1_17_LS2_avg_WA_1_WA_2_WA_3.npy')
F_222 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\CT_120_1_17_LS2_Kraft.npy')

plt.subplot(233)
plt.plot((abs(u_111))/300, abs(F_111 *1000)/(np.pi * 100**2 /4), "k")
plt.plot((abs(u_222))/300, abs(F_222 *1000)/(np.pi * 100**2 /4), "r")

plt.xlim(0.0, 0.004)
plt.ylim(0.0, 130)
plt.xlabel('eps')
plt.ylabel('sigma')


#creep-fatigue
N_111 =  len(np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\CT_120_1_16_LS2_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_111 = np.arange(1, N_111 +1, 1)
eps_max_111 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\CT_120_1_16_LS2_avg_WA_1_WA_2_WA_3_max.npy')
eps_min_111 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\CT_120_1_16_LS2_avg_WA_1_WA_2_WA_3_min.npy')

N_222 =  len(np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\CT_120_1_17_LS2_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_222 = np.arange(1, N_222 +1, 1)
eps_max_222 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\CT_120_1_17_LS2_avg_WA_1_WA_2_WA_3_max.npy')
eps_min_222 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\CT_120_1_17_LS2_avg_WA_1_WA_2_WA_3_min.npy')

plt.subplot(236)
plt.plot(N_max_111[1:] *100/N_111, (abs(eps_max_111[1:]))/300, "k")
plt.plot(N_max_111[1:-1]*100/N_111, (abs(eps_min_111[1:]))/300, "k")

plt.plot(N_max_222[1:] *100/N_222, (abs(eps_max_222[1:]))/300, "k")
plt.plot(N_max_222[1:]*100/N_222, (abs(eps_min_222[1:]))/300, "k")

plt.ylim(0.0, 0.004)
plt.title('Fatigue creep curve normalized (H-L)')
plt.xlabel('N/Nf')
plt.ylabel('Displacement [mm]')








plt.show()