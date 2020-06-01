'''
Created on 27.03.2019

@author: abaktheer
'''


import matplotlib.pyplot as plt
import numpy as np

#===================================
"LS5"
#===================================

"#C40"
#-------------------------
#creep-fatigue (L)
#C40_11
N_40_1 =  len(np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-11-181336Zyk_avg_WA1_WA2_WA3_max.npy'))
N_40_max_1 = np.arange(1, N_40_1 +1, 1)
eps_40_max_1 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-11-181336Zyk_avg_WA1_WA2_WA3_max.npy')
eps_40_min_1 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-11-181336Zyk_avg_WA1_WA2_WA3_min.npy')

#C40_12
N_40_2 =  len(np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-12-61761Zyk_avg_WA1_WA2_WA3_max.npy'))
N_40_max_2 = np.arange(1, N_40_2 +1, 1)
eps_40_max_2 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-12-61761Zyk_avg_WA1_WA2_WA3_max.npy')
eps_40_min_2 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-12-61761Zyk_avg_WA1_WA2_WA3_min.npy')

#C40_13
N_40_3 =  len(np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-13-122674Zyk_avg_WA1_WA2_WA3_max.npy'))
N_40_max_3 = np.arange(1, N_40_3 +1, 1)
eps_40_max_3 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-13-122674Zyk_avg_WA1_WA2_WA3_max.npy')
eps_40_min_3 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-13-122674Zyk_avg_WA1_WA2_WA3_min.npy')

#C40_37
N_40_4 =  len(np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-37-854634_avg_WA1 mm_WA2 mm_WA3 mm_max.npy'))
N_40_max_4 = np.arange(1, N_40_4 +1, 1)
eps_40_max_4 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-37-854634_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')
eps_40_min_4 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-37-854634_avg_WA1 mm_WA2 mm_WA3 mm_min.npy')

#C40_38
N_40_5 =  len(np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-38-518973_avg_WA1 mm_WA2 mm_WA3 mm_max.npy'))
N_40_max_5 = np.arange(1, N_40_5 +1, 1)
eps_40_max_5 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-38-518973_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')
eps_40_min_5 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-38-518973_avg_WA1 mm_WA2 mm_WA3 mm_min.npy')


#C40_39
N_40_6 =  len(np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-39-2847129_avg_WA1 mm_WA2 mm_WA3 mm_max.npy'))
N_40_max_6 = np.arange(1, N_40_6 +1, 1)
eps_40_max_6 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-39-2847129_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')
eps_40_min_6 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-39-2847129_avg_WA1 mm_WA2 mm_WA3 mm_min.npy')

plt.subplot(331)
plt.plot(N_40_max_1[1:] /N_40_1, (abs(eps_40_max_1[1:]))/300, "k")
plt.plot(N_40_max_2[1:] /N_40_2, (abs(eps_40_max_2[1:]))/300, "k")
plt.plot(N_40_max_3[1:] /N_40_3, (abs(eps_40_max_3[1:]))/300, "k")
plt.plot(N_40_max_4[1:] /N_40_4, (abs(eps_40_max_4[1:])-0.15)/300, "r")
plt.plot(N_40_max_5[1:] /N_40_5, (abs(eps_40_max_5[1:]))/300, "g")
plt.plot(N_40_max_6[1:] /N_40_6, (abs(eps_40_max_6[1:])-0.15)/300, "b")


plt.ylim(0.001, 0.0045)
# plt.title('Fatigue creep curve normalized (H-L)')
# plt.xlabel('N/Nf')
# plt.ylabel('Displacement [mm]')

#--------------------------
#creep-fatigue (H)
#C40_8
N_40_11 =  len(np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-8-6094Zyk_avg_WA1_WA2_WA3_max.npy'))
N_40_max_11 = np.arange(1, N_40_11 +1, 1)
eps_40_max_11 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-8-6094Zyk_avg_WA1_WA2_WA3_max.npy')
eps_40_min_11 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-8-6094Zyk_avg_WA1_WA2_WA3_min.npy')

#C40_9
N_40_22 =  len(np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-9-29009Zyk_avg_WA1_WA2_WA3_max.npy'))
N_40_max_22 = np.arange(1, N_40_22 +1, 1)
eps_40_max_22 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-9-29009Zyk_avg_WA1_WA2_WA3_max.npy')
eps_40_min_22 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-9-29009Zyk_avg_WA1_WA2_WA3_min.npy')

#C40_10
N_40_33 =  len(np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-10-21018Zyk_avg_WA1_WA2_WA3_max.npy'))
N_40_max_33 = np.arange(1, N_40_33 +1, 1)
eps_40_max_33 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-10-21018Zyk_avg_WA1_WA2_WA3_max.npy')
eps_40_min_33 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-10-21018Zyk_avg_WA1_WA2_WA3_min.npy')


#C40_14
N_40_44 =  len(np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-14-13356Zyk_avg_WA1_WA2_WA3_max.npy'))
N_40_max_44 = np.arange(1, N_40_44 +1, 1)
eps_40_max_44 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-14-13356Zyk_avg_WA1_WA2_WA3_max.npy')
eps_40_min_44 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-14-13356Zyk_avg_WA1_WA2_WA3_min.npy')

#C40_15
N_40_55 =  len(np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-15-16960Zyk_avg_WA1_WA2_WA3_max.npy'))
N_40_max_55 = np.arange(1, N_40_55 +1, 1)
eps_40_max_55 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-15-16960Zyk_avg_WA1_WA2_WA3_max.npy')
eps_40_min_55 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-15-16960Zyk_avg_WA1_WA2_WA3_min.npy')

#C40_16
N_40_66 =  len(np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-16-10606Zyk_avg_WA1_WA2_WA3_max.npy'))
N_40_max_66 = np.arange(1, N_40_66 +1, 1)
eps_40_max_66 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-16-10606Zyk_avg_WA1_WA2_WA3_max.npy')
eps_40_min_66 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-16-10606Zyk_avg_WA1_WA2_WA3_min.npy')


plt.subplot(332)
plt.plot(N_40_max_11[1:] /N_40_11, (abs(eps_40_max_11[1:]))/300, "k")
plt.plot(N_40_max_22[1:] /N_40_22, (abs(eps_40_max_22[1:]))/300, "k")
plt.plot(N_40_max_33[1:] /N_40_33, (abs(eps_40_max_33[1:]))/300, "k")
plt.plot(N_40_max_44[1:] /N_40_44, (abs(eps_40_max_44[1:]))/300, "k")
plt.plot(N_40_max_55[1:] /N_40_55, (abs(eps_40_max_55[1:]))/300, "k")
plt.plot(N_40_max_66[1:] /N_40_66, (abs(eps_40_max_66[1:]))/300, "k")

#plt.plot(N_max_1[1:]/N_1, (abs(eps_min_1[1:])-0.038)/300, "k")

plt.ylim(0.001, 0.004)

#------------------------
#Wohler curve


#Smax=0.75
N_40_1 = np.array([6094, 29009, 21018, 13356, 16960, 10606, 55140, 47806, 37864 ])
S_40_1 = np.array([0.75,  0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75  ])

#Smax=0.65
N_40_2 = np.array([181366, 61761, 122674, 18043, 173249, 133530, 854634, 518973, 2847129])
S_40_2 = np.array([0.65,  0.65, 0.65, 0.65,  0.65, 0.65, 0.65,  0.65, 0.65   ])

# average
N_40_3 = np.array([np.average(N_40_1) , np.average(N_40_2)])
S_40_3 = np.array([0.75,  0.65])


# FIB model code 2010
N_fib = np.array([89251862,    14307230,    367648,    58935,    9447,    1514,    243,    39,    6,    1])
S_fib = np.array([0.5,    0.55,    0.65,    0.7,    0.75,    0.8,    0.85,    0.9,    0.95,    1])


plt.subplot(333)
plt.plot(np.log10(N_40_1), S_40_1, 'ko', markersize=4, color='k')
plt.plot(np.log10(N_40_2), S_40_2, 'ko', markersize=4, color='k')
plt.plot(np.log10(N_40_3), S_40_3,  color='r',)
plt.plot(np.log10(N_fib), S_fib,  color='k',)

plt.ylim(0.55, 0.95)
plt.xlim(2, 8)




#==========================================================================================================

"#C80"

#creep-fatigue (L)
#C80_44
N_80_1 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low\NPY\CT80-44_74190_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_80_max_1 = np.arange(1, N_80_1 +1, 1)
eps_80_max_1 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low\NPY\CT80-44_74190_Zykl_avg_WA_1_WA_2_WA_3_max.npy')
eps_80_min_1 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low\NPY\CT80-44_74190_Zykl_avg_WA_1_WA_2_WA_3_min.npy')

#C80_45
N_80_2 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low\NPY\CT80-45_32437_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_80_max_2 = np.arange(1, N_80_2 +1, 1)
eps_80_max_2 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low\NPY\CT80-45_32437_Zykl_avg_WA_1_WA_2_WA_3_max.npy')
eps_80_min_2 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low\NPY\CT80-45_32437_Zykl_avg_WA_1_WA_2_WA_3_min.npy')

#C80_46
N_80_3 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low\NPY\CT80-46_416585_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_80_max_3 = np.arange(1, N_80_3 +1, 1)
eps_80_max_3 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low\NPY\CT80-46_416585_Zykl_avg_WA_1_WA_2_WA_3_max.npy')
eps_80_min_3 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low\NPY\CT80-46_416585_Zykl_avg_WA_1_WA_2_WA_3_min.npy')

#C80_47
N_80_4 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low\NPY\CT80-47_253449_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_80_max_4 = np.arange(1, N_80_4 +1, 1)
eps_80_max_4 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low\NPY\CT80-47_253449_Zykl_avg_WA_1_WA_2_WA_3_max.npy')
eps_80_min_4 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low\NPY\CT80-47_253449_Zykl_avg_WA_1_WA_2_WA_3_min.npy')

#C80_59
N_80_5 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low\NPY\CT80-59_142547_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_80_max_5 = np.arange(1, N_80_5 +1, 1)
eps_80_max_5 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low\NPY\CT80-59_142547_Zykl_avg_WA_1_WA_2_WA_3_max.npy')
eps_80_min_5 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low\NPY\CT80-59_142547_Zykl_avg_WA_1_WA_2_WA_3_min.npy')

#C80_60
N_80_6 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low\NPY\CT80-60_538303_Zykl2_avg_WA_1_WA_2_WA_3_max.npy'))
N_80_max_6 = np.arange(1, N_80_6 +1, 1)
eps_80_max_6 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low\NPY\CT80-60_538303_Zykl2_avg_WA_1_WA_2_WA_3_max.npy')
eps_80_min_6 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low\NPY\CT80-60_538303_Zykl2_avg_WA_1_WA_2_WA_3_min.npy')

plt.subplot(334)
plt.plot(N_80_max_1 [1:] /N_80_1, (abs(eps_80_max_1[1:]))/300, "r")
#plt.plot(N_120_max_1 [0:] /N_120_1, (abs(eps_120_min_1[1:]))/300, "g")

plt.plot(N_80_max_2 [1:] /N_80_2, (abs(eps_80_max_2[1:]))/300, "r")
#plt.plot(N_120_max_1 [0:] /N_120_1, (abs(eps_120_min_1[1:]))/300, "g")

plt.plot(N_80_max_3 [1:] /N_80_3, (abs(eps_80_max_3[1:]))/300, "r")
#plt.plot(N_120_max_1 [0:] /N_120_1, (abs(eps_120_min_1[1:]))/300, "g")

plt.plot(N_80_max_4 [1:] /N_80_4, (abs(eps_80_max_4[1:]))/300, "r")
#plt.plot(N_120_max_1 [0:] /N_120_1, (abs(eps_120_min_1[1:]))/300, "g")

plt.plot(N_80_max_5 [1:] /N_80_5, (abs(eps_80_min_5[1:-1]))/300, "r")
#plt.plot(N_120_max_1 [0:] /N_120_1, (abs(eps_120_min_1[1:]))/300, "g")

plt.plot(N_80_max_6 [1:] /N_80_6, (abs(eps_80_max_6[1:]))/300, "r")
#plt.plot(N_120_max_1 [0:] /N_120_1, (abs(eps_120_min_1[1:]))/300, "g")

plt.ylim(0.0010, 0.004)


#=============================================================
#creep-fatigue (H)
#C80_39
N_80_11 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High\NPY\CT80-39_6322_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_80_max_11 = np.arange(1, N_80_11 +1, 1)
eps_80_max_11 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High\NPY\CT80-39_6322_Zykl_avg_WA_1_WA_2_WA_3_max.npy')
eps_80_min_11 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High\NPY\CT80-39_6322_Zykl_avg_WA_1_WA_2_WA_3_min.npy')

#C80_40
N_80_22 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High\NPY\CT80-40_13713_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_80_max_22 = np.arange(1, N_80_22 +1, 1)
eps_80_max_22 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High\NPY\CT80-40_13713_Zykl_avg_WA_1_WA_2_WA_3_max.npy')
eps_80_min_22 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High\NPY\CT80-40_13713_Zykl_avg_WA_1_WA_2_WA_3_min.npy')

#C80_42
N_80_33 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High\NPY\CT80-42_3610_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_80_max_33 = np.arange(1, N_80_33 +1, 1)
eps_80_max_33 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High\NPY\CT80-42_3610_Zykl_avg_WA_1_WA_2_WA_3_max.npy')
eps_80_min_33 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High\NPY\CT80-42_3610_Zykl_avg_WA_1_WA_2_WA_3_min.npy')

#C80_43
N_80_44 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High\NPY\CT80-43_10949_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_80_max_44 = np.arange(1, N_80_44 +1, 1)
eps_80_max_44 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High\NPY\CT80-43_10949_Zykl_avg_WA_1_WA_2_WA_3_max.npy')
eps_80_min_44 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High\NPY\CT80-43_10949_Zykl_avg_WA_1_WA_2_WA_3_max.npy')

#C80_48
N_80_55 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High\NPY\CT80-48_818_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_80_max_55 = np.arange(1, N_80_55 +1, 1)
eps_80_max_55 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High\NPY\CT80-48_818_Zykl_avg_WA_1_WA_2_WA_3_max.npy')
eps_80_min_55 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High\NPY\CT80-48_818_Zykl_avg_WA_1_WA_2_WA_3_min.npy')

#C80_52
N_80_66 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High\NPY\CT80-52_647_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_80_max_66 = np.arange(1, N_80_66 +1, 1)
eps_80_max_66 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High\NPY\CT80-52_647_Zykl_avg_WA_1_WA_2_WA_3_max.npy')
eps_80_min_66 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High\NPY\CT80-52_647_Zykl_avg_WA_1_WA_2_WA_3_min.npy')


plt.subplot(335)
plt.plot(N_80_max_11 [1:] /N_80_11, (abs(eps_80_max_11[1:]))/300, "r")
#plt.plot(N_120_max_1 [0:] /N_120_1, (abs(eps_120_min_1[1:]))/300, "g")

plt.plot(N_80_max_22 [1:] /N_80_22, (abs(eps_80_max_22[1:]))/300, "r")
#plt.plot(N_120_max_1 [0:] /N_120_1, (abs(eps_120_min_1[1:]))/300, "g")

plt.plot(N_80_max_33 [1:-1] /N_80_33, (abs(eps_80_min_33[1:]))/300, "r")
#plt.plot(N_120_max_1 [0:] /N_120_1, (abs(eps_120_min_1[1:]))/300, "g")

plt.plot(N_80_max_44 [1:] /N_80_44, (abs(eps_80_max_44[1:]))/300, "r")
#plt.plot(N_120_max_1 [0:] /N_120_1, (abs(eps_120_min_1[1:]))/300, "g")

plt.plot(N_80_max_55 [1:] /N_80_55, (abs(eps_80_max_55[1:]))/300, "r")
#plt.plot(N_120_max_1 [0:] /N_120_1, (abs(eps_120_min_1[1:]))/300, "g")

plt.plot(N_80_max_66 [1:] /N_80_66, (abs(eps_80_max_66[1:]))/300, "r")
#plt.plot(N_120_max_1 [0:] /N_120_1, (abs(eps_120_min_1[1:]))/300, "g")


plt.ylim(0.0010, 0.004)

#=============================================================
#Wohler

#Smax=0.85
N_80_1 = np.array([6322, 13713, 3610, 10949, 818, 647, 1969, 916])
S_80_1 = np.array([0.85,  0.85, 0.85, 0.85, 0.85, 0.85,  0.85, 0.85])

#Smax=0.80
N_80_2 = np.array([7835,    17932,    32173,    11339 ])
S_80_2 = np.array([0.80,  0.80, 0.80, 0.80 ])

#Smax=0.75
N_80_3 = np.array([74190, 32437, 416585, 253449, 142547, 538303, 2012546, 163931,    84205,    234316])
S_80_3 = np.array([0.75,  0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75])

#Smax=0.65
N_80_4 = np.array([2160000,    2065037,    3074229])
S_80_4 = np.array([0.65,  0.65,  0.65 ])

# average
N_80_5 = np.array([np.average(N_80_1) , np.average(N_80_2) , np.average(N_80_3) , np.average(N_80_4)])
S_80_5 = np.array([0.85, 0.80, 0.75,  0.65])

# FIB model code 2010
N_fib = np.array([32064472260,    2853787864,    22605643,    2011937,    179066,    15937,    1418,    126,    11,    1,])
S_fib = np.array([0.5,    0.55,    0.65,    0.7,    0.75,    0.8,    0.85,    0.9,    0.95,    1])


plt.subplot(336)
plt.plot(np.log10(N_80_1), S_80_1, 'ro', markersize=4, color='r')
plt.plot(np.log10(N_80_2), S_80_2, 'ro', markersize=4, color='r')
plt.plot(np.log10(N_80_3), S_80_3, 'ro', markersize=4, color='r')
plt.plot(np.log10(N_80_4), S_80_4, 'ro', markersize=4, color='r')
plt.plot(np.log10(N_80_5), S_80_5,  color='r',)
plt.plot(np.log10(N_fib), S_fib,  color='k',)

plt.ylim(0.55, 0.95)
plt.xlim(2, 8)









#==========================================================================================================

"#C120"


#creep-fatigue (L)
#C120_20
N_120_1 =  len(np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\CT_120_1_20_LS4L_6750Zyk_avg_WA_1_WA_2_WA_3_max.npy'))
N_120_max_1 = np.arange(1, N_120_1 +1, 1)
eps_120_max_1 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\CT_120_1_20_LS4L_6750Zyk_avg_WA_1_WA_2_WA_3_max.npy')
eps_120_min_1 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\CT_120_1_20_LS4L_6750Zyk_avg_WA_1_WA_2_WA_3_min.npy')

#C120_22
N_120_2 =  len(np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\WinConFat_CT_120_1_22_Zyklen_avg_WA_1_WA_2_WA_3_max.npy'))
N_120_max_2 = np.arange(1, N_120_2 +1, 1)
eps_120_max_2 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\WinConFat_CT_120_1_22_Zyklen_avg_WA_1_WA_2_WA_3_max.npy')
eps_120_min_2 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\WinConFat_CT_120_1_22_Zyklen_avg_WA_1_WA_2_WA_3_min.npy')

#C120_23
N_120_3 =  len(np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\WinConFat_CT_120_1_23_0-5626000Zyklen_avg_WA_1_WA_2_WA_3_max.npy'))
N_120_max_3 = np.arange(1, N_120_3 +1, 1)
eps_120_max_3 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\WinConFat_CT_120_1_23_0-5626000Zyklen_avg_WA_1_WA_2_WA_3_max.npy')
#eps_120_min_2 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\WinConFat_CT_120_1_22_Zyklen_avg_WA_1_WA_2_WA_3_min.npy')

plt.subplot(337)
plt.plot(N_120_max_1 [1:] /N_120_1, (abs(eps_120_max_1[1:])-0.2)/300, "g")
#plt.plot(N_120_max_1 [0:] /N_120_1, (abs(eps_120_min_1[1:]))/300, "g")

plt.plot(N_120_max_2 [1:] /N_120_2, (abs(eps_120_max_2[1:]))/300, "g")

plt.plot(N_120_max_3 [1:] /N_120_3, (abs(eps_120_max_3[1:]))/300, "r")
#plt.plot(N_120_max_2 [1:-1] /N_120_2, (abs(eps_120_min_2[1:]))/300, "g")



plt.ylim(0.001, 0.004)



#=====================================
#creep-fatigue (H)
#C120_24
N_120_11 =  len(np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\WinConFat_CT_120_1_24_13504Zyklen_avg_WA_1_WA_2_WA_3_max.npy'))
N_120_max_11 = np.arange(1, N_120_11 +1, 1)
eps_120_max_11 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\WinConFat_CT_120_1_24_13504Zyklen_avg_WA_1_WA_2_WA_3_max.npy')
eps_120_min_11 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\WinConFat_CT_120_1_24_13504Zyklen_avg_WA_1_WA_2_WA_3_min.npy')

#C120_25
N_120_22 =  len(np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\WinConFat_CT_120_1_25_23256Zyklen_avg_WA_1_WA_2_WA_3_max.npy'))
N_120_max_22 = np.arange(1, N_120_22 +1, 1)
eps_120_max_22 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\WinConFat_CT_120_1_25_23256Zyklen_avg_WA_1_WA_2_WA_3_max.npy')
eps_120_min_22 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\WinConFat_CT_120_1_25_23256Zyklen_avg_WA_1_WA_2_WA_3_min.npy')

#C120_26
N_120_33 =  len(np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\WinConFat_CT_120_1_26_ca23000Zyklen_avg_WA_1_WA_2_WA_3_max.npy'))
N_120_max_33 = np.arange(1, N_120_33 +1, 1)
eps_120_max_33 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\WinConFat_CT_120_1_26_ca23000Zyklen_avg_WA_1_WA_2_WA_3_max.npy')
eps_120_min_33 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\WinConFat_CT_120_1_26_ca23000Zyklen_avg_WA_1_WA_2_WA_3_min.npy')

#C120_29
N_120_44 =  len(np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\WinConFat_CT_120_1_29_5755Zyklen_avg_WA_1_WA_2_WA_3_max.npy'))
N_120_max_44 = np.arange(1, N_120_44 +1, 1)
eps_120_max_44 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\WinConFat_CT_120_1_29_5755Zyklen_avg_WA_1_WA_2_WA_3_max.npy')
eps_120_min_44 = np.load( r'D:\Heimarbeit\Rohedata\C120\NPY\WinConFat_CT_120_1_29_5755Zyklen_avg_WA_1_WA_2_WA_3_min.npy')

plt.subplot(338)
plt.plot(N_120_max_11 [2:] /N_120_11, (abs(eps_120_max_11[2:])-0.22)/300, "g")
#plt.plot(N_120_max_1 [1:] /N_120_1, (abs(eps_120_min_1[1:]))/300, "g")

plt.plot(N_120_max_22 [2:] /N_120_22, (abs(eps_120_max_22[2:])-0.5)/300, "g")
#plt.plot(N_120_max_2 [1:-1] /N_120_2, (abs(eps_120_min_2[1:]))/300, "--g")

plt.plot(N_120_max_33 [2:] /N_120_33, (abs(eps_120_max_33[2:])-0.18)/300, "g")
#plt.plot(N_120_max_2 [1:-1] /N_120_2, (abs(eps_120_min_2[1:]))/300, "g")

plt.plot(N_120_max_44 [2:] /N_120_44, (abs(eps_120_min_44[2:-1])-0.6)/300, "g")
#plt.plot(N_120_max_2 [1:-1] /N_120_2, (abs(eps_120_min_2[1:]))/300, "g")

plt.ylim(0.0010, 0.004)











#Smax=0.85
N_120_1 = np.array([13504, 23256, 24881 , 5755 ])
S_120_1 = np.array([0.85,  0.85, 0.85 , 0.85])

#Smax=0.75
N_120_2 = np.array([6750, 193528, 5626000, 5597950])
S_120_2 = np.array([0.75,  0.75, 0.75 , 0.75 ])

# average
N_120_3 = np.array([np.average(N_120_1), np.average(N_120_2)]) 
S_120_3 = np.array([0.85,  0.75])


# FIB model code 2010
N_fib = np.array([32064472260,    2853787864,    22605643,    2011937,    179066,    15937,    1418,    126,    11,    1,])
S_fib = np.array([0.5,    0.55,    0.65,    0.7,    0.75,    0.8,    0.85,    0.9,    0.95,    1])


plt.subplot(339)
plt.plot(np.log10(N_120_1), S_120_1, 'go', markersize=4, color='r')
plt.plot(np.log10(N_120_2), S_120_2, 'go', markersize=4, color='r')
plt.plot(np.log10(N_120_3), S_120_3,  color='r',)
plt.plot(np.log10(N_fib), S_fib,  color='k',)

plt.ylim(0.55, 0.95)
plt.xlim(2, 8)







plt.show()