'''
Created on 27.03.2019

@author: abaktheer
'''


import matplotlib.pyplot as plt
import numpy as np



#===================================
"LS6"
#===================================

#============================================
"#C40"

#=============================================
''' H-L (S= 0.75 -0.65) (CT_40_23) '''
#=============================================
eps_max_high = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-23-High_000-4043_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')
eps_max_low = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-23-low4043-123977_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')
N_1 =  len(eps_max_high) + len(eps_max_low ) 
N_max_1 = np.arange(1, N_1 +1, 1)
eps_max_1 = np.hstack((eps_max_high ,  eps_max_low)) - 0.10

#=============================================
''' H-L (S= 0.75 -0.65) (CT_40_25) '''
#=============================================
eps_max_high_2 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-25_high_000_4043_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')
eps_max_low_2 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-25_low_4043_154014_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')
N_2 =  len(eps_max_high_2) + len(eps_max_low_2 ) 
N_max_2 = np.arange(1, N_2 +1, 1)
eps_max_2 = np.hstack((eps_max_high_2 ,  eps_max_low_2)) -0.065

#=============================================
''' H-L (S= 0.75 -0.65) (CT_40_26) '''
#=============================================
eps_max_high_3 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-26_High_000-4043_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')
eps_max_low_3 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-26_low_4043_251272_avg_WA1 mm_WA2 mm_WA3 mm_min.npy')
N_3 =  len(eps_max_high_3) + len(eps_max_low_3 ) 
N_max_3 = np.arange(1, N_3 +1, 1)
eps_max_3 = np.hstack((eps_max_high_3 ,  eps_max_low_3)) -0.166


#=============================================
''' H-L (S= 0.75 -0.65) (CT_40_27) '''
#=============================================
eps_max_high_4 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-27_high_000-4043_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')
eps_max_low_4 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-27_low_4043_595562_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')
N_4 =  len(eps_max_high_4) + len(eps_max_low_4 ) 
N_max_4 = np.arange(1, N_4 +1, 1)
eps_max_4 = np.hstack((eps_max_high_4 ,  eps_max_low_4)) -0.061



n_01 = 4043
n_02 = 4043
n_03 = 4043
n_04 = 4043


n_hf = 26428.0
n_lf = 545707.0

# normalized accumulative damage
N_01 = np.hstack(
    (N_max_1[0:n_01 + 1] / n_hf, n_01 / n_hf + (N_max_1[n_01 + 1:] - N_max_1[n_01]) / n_lf))

N_02 = np.hstack(
    (N_max_2[0:n_02 + 1] / n_hf, n_02 / n_hf + (N_max_2[n_02 + 1:] - N_max_2[n_02]) / n_lf))

N_03 = np.hstack(
    (N_max_3[0:n_03 + 1] / n_hf, n_03 / n_hf + (N_max_3[n_03 + 1:] - N_max_3[n_03]) / n_lf))

N_04 = np.hstack(
    (N_max_4[0:n_04 + 1] / n_hf, n_04 / n_hf + (N_max_4[n_04 + 1:] - N_max_4[n_04]) / n_lf))





#=============================================
''' H-L (S= 0.75 -0.65) (CT_40_43) '''
#=============================================
eps_max_high_5 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-43-318066_avg_WA1 mm_WA2 mm_WA3 mm_max.npy') - 0.092
N_5 =  len(eps_max_high_5) 
N_max_5 = np.arange(1, N_5 +1, 1)
eps_max_5 = np.hstack((eps_max_high_5))

#=============================================
''' H-L (S= 0.75 -0.65) (CT_40_44) '''
#=============================================
eps_max_high_6 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-44_00-1463573Zykl_avg_WA1 mm_WA2 mm_WA3 mm_max.npy') 
N_6 =  len(eps_max_high_6) 
N_max_6 = np.arange(1, N_6 +1, 1)
eps_max_6 = np.hstack((eps_max_high_6))




#=============================================
''' H-L (S= 0.75 -0.65) (CT_40_46) '''
#=============================================
eps_max_high_7 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-46_high_000-1099604_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')
eps_max_low_7 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-46_high_1099604-1542589_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')
N_7 =  len(eps_max_high_7) + len(eps_max_low_7 ) 
N_max_7 = np.arange(1, N_7 +1, 1)
eps_max_7 = np.hstack((eps_max_high_7 ,  eps_max_low_7))

#=============================================
''' H-L (S= 0.75 -0.65) (CT_40_47) '''
#=============================================
eps_max_high_8 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-47_low&high_0-2072915_avg_WA1 mm_WA2 mm_WA3 mm_max.npy') -0.065
N_8 =  len(eps_max_high_8) 
N_max_8 = np.arange(1, N_8 +1, 1)
eps_max_8 = np.hstack((eps_max_high_8))



#=============================================
''' H-L (S= 0.75 -0.65) (CT_40_48) '''
#=============================================
eps_max_high_9 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-48_high_0-23468_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')
eps_max_low_9 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-48_high_23468-1571981_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')
N_9 =  len(eps_max_high_9) + len(eps_max_low_9) 
N_max_9 = np.arange(1, N_9 +1, 1)
eps_max_9 = np.hstack((eps_max_high_9 ,  eps_max_low_9)) -0.09






n_05 = 8087
n_06 = 8087
n_07 = 23468
n_08 = 23468
n_09 = 23468

n_hf = 49671.0
n_lf = 2217025.0

N_05 = np.hstack(
    (N_max_5[0:n_05 + 1] / n_hf, n_05 / n_hf + (N_max_5[n_05 + 1:] - N_max_5[n_05]) / n_lf))
N_06 = np.hstack(
    (N_max_6[0:n_06 + 1] / n_hf, n_06 / n_hf + (N_max_6[n_06 + 1:] - N_max_6[n_06]) / n_lf))

N_07 = np.hstack(
    (N_max_7[0:n_07 + 1] / n_hf, n_07 / n_hf + (N_max_7[n_07 + 1:] - N_max_7[n_07]) / n_lf))

N_08 = np.hstack(
    (N_max_8[0:n_08 + 1] / n_hf, n_08 / n_hf + (N_max_8[n_08 + 1:] - N_max_8[n_08]) / n_lf))

N_09 = np.hstack(
    (N_max_9[0:n_09 + 1] / n_hf, n_09 / n_hf + (N_max_9[n_09 + 1:] - N_max_9[n_09]) / n_lf))



#========================================
# Plotting
#========================================


plt.subplot(331)
plt.plot(N_01[2:], abs(eps_max_1[2:])/300, "k")
plt.plot(N_02[2:], abs(eps_max_2[2:])/300, "k")
plt.plot(N_03[2:], abs(eps_max_3[2:])/300, "k")
plt.plot(N_04[2:], abs(eps_max_4[2:])/300, "k")
plt.plot(N_05[2:], abs(eps_max_5[2:])/300, "b")
plt.plot(N_06[2:], abs(eps_max_6[2:])/300, "b")
plt.plot(N_07[2:], abs(eps_max_7[2:])/300, "r")
plt.plot(N_08[2:], abs(eps_max_8[2:])/300, "r")
plt.plot(N_09[2:], abs(eps_max_9[2:])/300, "r")

plt.ylim(0.001, 0.0035)
plt.xlim(-0.1 , 2.0 )



#===================================================================================
'''
L-H-1
'''

#=============================================
''' L-H (S= 0.75 -0.65) (CT_40_28) '''
#=============================================
eps_max_high_1 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-28_low_000_28776_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')
eps_max_low_1 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-28_high_28776_46174_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')
N_1 =  len(eps_max_high_1) + len(eps_max_low_1 ) 
N_max_1 = np.arange(1, N_1 +1, 1)
eps_max_1 = np.hstack((eps_max_high_1 ,  eps_max_low_1)) - 0.183
print('N_1=', N_1)

#=============================================
''' L-H (S= 0.75 -0.65) (CT_40_29) '''
#=============================================
eps_max_high_2 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-29_low_28776_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')
eps_max_low_2 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-29_high_28776_55659_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')
N_2 =  len(eps_max_high_2) + len(eps_max_low_2 ) 
N_max_2 = np.arange(1, N_2 +1, 1)
eps_max_2 = np.hstack((eps_max_high_2 ,  eps_max_low_2)) - 0.06
print('N_2=', N_2)


#=============================================
''' L-H  (S= 0.85 -0.75) (CT_80-30) '''
#=============================================
eps_max_3 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-30-42053_avg_WA1 mm_WA2 mm_WA3 mm_max.npy') - 0.061
N_3 =  len(eps_max_3) 
N_max_3 = np.arange(1, N_3 +1, 1)
print('N_3=', N_3)

#=============================================
''' L-H  (S= 0.85 -0.75) (CT_80-31) '''
#=============================================
eps_max_4 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-31-61337_avg_WA1 mm_WA2 mm_WA3 mm_min.npy') - 0.11
N_4 =  len(eps_max_4) 
N_max_4 = np.arange(1, N_4 +1, 1)
print('N_4=', N_4)

#=============================================
''' L-H  (S= 0.85 -0.75) (CT_80-32) '''
#=============================================
eps_max_5 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-32-72867_avg_WA1 mm_WA2 mm_WA3 mm_max.npy') - 0.253
N_5 =  len(eps_max_5) 
N_max_5 = np.arange(1, N_5 +1, 1)
print('N_5=', N_5)

#=============================================
''' L-H  (S= 0.85 -0.75) (CT_80-33) '''
#=============================================
eps_max_6 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-33-91796_avg_WA1 mm_WA2 mm_WA3 mm_min.npy') - 0.05
N_6 =  len(eps_max_6) 
N_max_6 = np.arange(1, N_6 +1, 1)
print('N_6=', N_6)

n_01 = 28776
n_02 = 28776
n_03 = 28776
n_04 = 28776
n_05 = 28776
n_06 = 28776

n_hf = 26428.0
n_lf = 545707.0

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
N_06 = np.hstack(
    (N_max_6[0:n_06 + 1] / n_lf, n_06 / n_lf + (N_max_6[n_06 + 1:] - N_max_6[n_06]) / n_hf))




#=============================================
''' L-H (S= 0.75 -0.65) (CT_40_56) '''
#=============================================
eps_max_high_7 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-56_low_281382_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')
eps_max_low_7 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-56_high_281382_356224_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')
N_7 =  len(eps_max_high_7) + len(eps_max_low_7 ) 
N_max_7 = np.arange(1, N_7 +1, 1)
eps_max_7 = np.hstack((eps_max_high_7 ,  eps_max_low_7)) -0.083


#=============================================
''' L-H (S= 0.75 -0.65) (CT_40_57) '''
#=============================================
eps_max_8 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-57_low&high_350247_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')-0.25
N_8 =  len(eps_max_8) 
N_max_8 = np.arange(1, N_8 +1, 1)



#=============================================
''' L-H (S= 0.75 -0.65) (CT_40_58) '''
#=============================================
eps_max_9 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-58_lowhigh 454099_avg_WA1 mm_WA2 mm_WA3 mm_max.npy') - 0.15
N_9 =  len(eps_max_9) 
N_max_9 = np.arange(1, N_9 +1, 1)


#=============================================
''' L-H (S= 0.75 -0.65) (CT_40_59) '''
#=============================================
eps_max_10 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-59_lowhigh_000-450805_avg_WA1 mm_WA2 mm_WA3 mm_max.npy')
N_10 =  len(eps_max_10) 
N_max_10 = np.arange(1, N_10 +1, 1)


#=============================================
''' L-H (S= 0.75 -0.65) (CT_40_60) '''
#=============================================
eps_max_11 = np.load( r'D:\Heimarbeit\Rohedata\C40\NPY\C40-60_lowhigh_000-426683_avg_WA1 mm_WA2 mm_WA3 mm_max.npy') - 0.717
N_11 =  len(eps_max_11) 
N_max_11 = np.arange(1, N_11 +1, 1)



n_07 = 281382
n_08 = 281382
n_09 = 409280
n_10 = 409280
n_11 = 409280


n_hf = 49671.0
n_lf = 2217025.0

# normalized accumulative damage

N_07 = np.hstack(
    (N_max_7[0:n_07 + 1] / n_lf, n_07 / n_lf + (N_max_7[n_07 + 1:] - N_max_7[n_07]) / n_hf))
N_08 = np.hstack(
    (N_max_8[0:n_08 + 1] / n_lf, n_08 / n_lf + (N_max_8[n_08 + 1:] - N_max_8[n_08]) / n_hf))
N_09 = np.hstack(
    (N_max_9[0:n_09 + 1] / n_lf, n_09 / n_lf + (N_max_9[n_09 + 1:] - N_max_9[n_09]) / n_hf))
N_10 = np.hstack(
    (N_max_10[0:n_10 + 1] / n_lf, n_10 / n_lf + (N_max_10[n_10 + 1:] - N_max_10[n_10]) / n_hf))
N_11 = np.hstack(
    (N_max_11[0:n_11 + 1] / n_lf, n_11 / n_lf + (N_max_11[n_11 + 1:] - N_max_11[n_11]) / n_hf))










#========================================
# Plotting
#========================================
plt.subplot(332)
plt.plot(N_01[2:], abs(eps_max_1[2:])/300, "k")
plt.plot(N_02[2:], abs(eps_max_2[2:])/300, "k")
plt.plot(N_03[2:], abs(eps_max_3[2:])/300, "k")
plt.plot(N_04[2:], abs(eps_max_4[2:])/300, "k")
plt.plot(N_05[2:], abs(eps_max_5[2:])/300, "k")
plt.plot(N_06[2:], abs(eps_max_6[2:])/300, "k")
plt.plot(N_07[2:], abs(eps_max_7[2:])/300, "b")
plt.plot(N_08[2:], abs(eps_max_8[2:])/300, "b")
plt.plot(N_09[2:], abs(eps_max_9[2:])/300, "r")
plt.plot(N_10[2:], abs(eps_max_10[2:])/300, "r")
plt.plot(N_11[2:], abs(eps_max_11[2:])/300, "r")

plt.ylim(0.001, 0.0035)
plt.xlim(-0.1 , 3.0 )





#===================================================================================
# comparison

# H-L
# HL_C40_1 = np.array([0.18, 0.22, 0.26, 0.27, 0.45, 1.08, 0.57, 2.67, 2.78])
# HL_C40_2 = np.array([0.15, 0.15, 0.15, 0.15, 0.15, 0.15,  0.305,  0.305, 0.89 ])

HL_C40_1 = np.array([0.18, 0.22, 0.26, 0.27, 0.45, 1.08, 0.14, 0.66, 0.69, 0.92, 0.70])
HL_C40_2 = np.array([0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.16, 0.16 , 0.47, 0.47, 0.47  ])

# L-H
# LH_C40_1 = np.array([0.053, 0.053, 0.053, 0.053, 0.053, 0.053])
# LH_C40_2  = np.array([0.66, 1.02, 0.50, 1.23, 1.67, 2.38])

LH_C40_1 = np.array([0.053, 0.053, 0.053, 0.053, 0.053, 0.053, 0.13, 0.13, 0.18, 0.18, 0.18])
LH_C40_2  = np.array([0.66, 1.02, 0.50, 1.23, 1.67, 2.38, 1.51, 1.39, 0.82, 0.84, 0.35])

# average
PM_1 = np.array([1 , 0])
PM_2 = np.array([0,  1])


plt.subplot(333)
plt.plot(HL_C40_1, HL_C40_2, 'ro', markersize=4, color='r')
plt.plot(LH_C40_1, LH_C40_2, 'ro', markersize=4, color='g')
plt.plot(PM_1, PM_2,  color='r',)

plt.xlim(0, 1.2)
plt.ylim(0, 2.5)
plt.title('sequence effect')
plt.xlabel('L')
plt.ylabel('H')


#===================================================================================================================================================================
"#C80"
#===================================================================================================================================================================


'''
H-L-1
'''

#=============================================
''' H-L (S= 0.85 -0.75) (CT_80-49) '''
#=============================================
N_1 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High-Low\NPY\CT80-49_185821_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_1 = np.arange(1, N_1 +1, 1)
eps_max_1 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High-Low\NPY\CT80-49_185821_Zykl_avg_WA_1_WA_2_WA_3_max.npy') -0.11


#=============================================
''' H-L (S= 0.85 -0.75) (CT_80-50) '''
#=============================================
N_2 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High-Low\NPY\CT80-50_810229_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_2 = np.arange(1, N_2 +1, 1)
eps_max_2 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High-Low\NPY\CT80-50_810229_Zykl_avg_WA_1_WA_2_WA_3_max.npy') -0.071


#=============================================
''' H-L (S= 0.85 -0.75) (CT_80-51) '''
#=============================================
N_3 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High-Low\NPY\CT80-51_1200277_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_3 = np.arange(1, N_3 +1, 1)
eps_max_3 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High-Low\NPY\CT80-51_1200277_Zykl_avg_WA_1_WA_2_WA_3_max.npy') -0.054

#=============================================
''' H-L (S= 0.85 -0.75) (CT_80-41) '''
#=============================================
N_4 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High-Low\NPY\CT80-41_274256_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_4 = np.arange(1, N_4 +1, 1)
eps_max_4 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High-Low\NPY\CT80-41_274256_Zykl_avg_WA_1_WA_2_WA_3_max.npy') -0.08 


'''
H-L-2
'''
#=============================================
''' H-L (S= 0.85 -0.75) (CT_80-53) '''
#=============================================
N_5 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High-Low\NPY\CT80-53_292317_Zykl2_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_5 = np.arange(1, N_5 +1, 1)
eps_max_5 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High-Low\NPY\CT80-53_292317_Zykl2_avg_WA_1_WA_2_WA_3_max.npy') -0.15

#=============================================
''' H-L (S= 0.85 -0.75) (CT_80-62) '''
#=============================================
N_6 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High-Low\NPY\CT80-62_292317_Zykl2_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_6 = np.arange(1, N_6 +1, 1)
eps_max_6 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High-Low\NPY\CT80-62_292317_Zykl2_avg_WA_1_WA_2_WA_3_max.npy') -0.152

#=============================================
''' H-L (S= 0.85 -0.75) (CT_80-63) '''
#=============================================
N_7 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High-Low\NPY\CT80-63_747618_Zykl2_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_7 = np.arange(1, N_7 +1, 1)
eps_max_7 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\High-Low\NPY\CT80-63_747618_Zykl2_avg_WA_1_WA_2_WA_3_max.npy') -0.096


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


plt.subplot(334)
plt.plot(N_01[2:], abs(eps_max_1[2:])/300, "k")
plt.plot(N_02[2:], abs(eps_max_2[2:])/300, "r")
plt.plot(N_03[2:], abs(eps_max_3[2:])/300, "b")
plt.plot(N_04[2:], abs(eps_max_4[2:])/300, "g")
plt.plot(N_05[2:], abs(eps_max_5[2:])/300, "y")
plt.plot(N_06[2:], abs(eps_max_6[2:])/300, "m")
plt.plot(N_07[2:], abs(eps_max_7[2:])/300, "-k")

plt.ylim(0.001, 0.0035)
plt.xlim(-0.1 , 2.5 )
#plt.title('Fatigue creep curve normalized (H-L)')
#plt.xlabel('N/Nf')
#plt.ylabel('Displacement [mm]')


#===================================================================================
'''
L-H-1
'''

#=============================================
''' L-H (S= 0.85 -0.75) (CT_80-54) '''
#=============================================
N_1 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low-High\NPY\CT80-54_63388_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_1 = np.arange(1, N_1 +1, 1)
eps_max_1 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low-High\NPY\CT80-54_63388_Zykl_avg_WA_1_WA_2_WA_3_max.npy') -0.072


#=============================================
''' L-H  (S= 0.85 -0.75) (CT_80-55) '''
#=============================================
N_2 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low-High\NPY\CT80-55_56805_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_2 = np.arange(1, N_2 +1, 1)
eps_max_2 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low-High\NPY\CT80-55_56805_Zykl_avg_WA_1_WA_2_WA_3_max.npy')


#=============================================
''' L-H  (S= 0.85 -0.75) (CT_80-56) '''
#=============================================
N_3 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low-High\NPY\CT80-56_70691_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_3 = np.arange(1, N_3 +1, 1)
eps_max_3 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low-High\NPY\CT80-56_70691_Zykl_avg_WA_1_WA_2_WA_3_max.npy') -0.058




'''
L-H-2
'''
#=============================================
''' L-H  (S= 0.85 -0.75) (CT_80-57) '''
#=============================================
N_4 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low-High\NPY\CT80-57_155538_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_4 = np.arange(1, N_4 +1, 1)
eps_max_4 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low-High\NPY\CT80-57_155538_Zykl_avg_WA_1_WA_2_WA_3_max.npy') - 0.075

#=============================================
''' L-H  (S= 0.85 -0.75) (CT_80-58) '''
#=============================================
N_5 =  len(np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low-High\NPY\CT80-58_162843_Zykl_avg_WA_1_WA_2_WA_3_max.npy'))
N_max_5 = np.arange(1, N_5 +1, 1)
eps_max_5 = np.load( r'D:\Heimarbeit\Rohedata\C80 Charge 2\Low-High\NPY\CT80-58_162843_Zykl_avg_WA_1_WA_2_WA_3_max.npy') - 0.073




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
plt.subplot(335)
plt.plot(N_01[1:], abs(eps_max_1[1:])/300, "k")
plt.plot(N_02[1:], abs(eps_max_2[1:])/300, "r")
plt.plot(N_03[1:], abs(eps_max_3[1:])/300, "b")
plt.plot(N_04[1:], abs(eps_max_4[1:])/300, "g")
plt.plot(N_05[1:], abs(eps_max_5[1:])/300, "y")



plt.ylim(0.0010, 0.0035)
plt.xlim(-0.1 , 3.75 )
#plt.title('Fatigue creep curve normalized (H-L)')
#plt.xlabel('N/Nf')
#plt.ylabel('Displacement [mm]')




#===================================================================================
# comparison
# H-L
HL_C80_1 = np.array([0.32 , 1.41, 2.09, 0.48, 0.51, 0.37, 1.3])
HL_C80_2 = np.array([0.24,  0.24, 0.24, 0.24, 0.47, 0.47, 0.47 ])

# L-H
LH_C80_1 = np.array([0.09, 0.09, 0.09, 0.26, 0.26])
LH_C80_2  = np.array([2.29, 1.2, 3.51, 1.11, 2.33])

# average
PM_1 = np.array([1 , 0])
PM_2 = np.array([0,  1])


plt.subplot(336)
plt.plot(HL_C80_1, HL_C80_2, 'ro', markersize=4, color='r')
plt.plot(LH_C80_1, LH_C80_2, 'ro', markersize=4, color='g')
plt.plot(PM_1, PM_2,  color='r',)

plt.ylim(0,3.7)
plt.xlim(0.0 , 2.15)
plt.title('sequence effect')
plt.xlabel('L')
plt.ylabel('H')







#=================================================================================================================================================================================
"Holmen"
#===================================================================================================================================================================



#creep-fatigue (HL)
N_1 = np.array([ 0.001451371,    0.01741654,    0.05224963,    0.07256893,    0.095791, 0.1248186,    0.1480406,    0.1567489,    0.1552975,    0.1727141, 0.1886792,    0.2162554,    0.2612482,    0.2946299,    0.3323658, 0.3701016,    0.3918723,    0.4150943,    0.4208998,    0.4208998,    0.4223512])
eps_max_1 = np.array([  1.48711,    1.573949,    1.736771,    1.796472,    1.856174,    1.899593, 1.92673,    1.937585,    1.753053,    1.823609,    1.943012,    2.089552, 2.274084,    2.398915,    2.556309,    2.740841,    2.881954,    3.044776, 3.191316,    3.327002,    3.41384])

N_2 = np.array([ 0,    0.01886792,    0.04644411,    0.06676342,    0.095791,    0.1248186,    0.149492,    0.1814223,    0.2089985,    0.2394775,    0.2670537,    0.2859216,    0.3004354,    0.3004354,    0.3222061,    0.3541364,    0.3788099,    0.4092888,    0.436865,    0.467344,    0.4992743,    0.5195936,    0.532656,    0.5428157,    0.5486212,    0.5486212])
eps_max_2 = np.array([1.454545,    1.541384,    1.65536,    1.715061,    1.785617,    1.856174,    1.894166,    1.943012,    1.981004,    2.029851,    2.062415,    2.084125,    2.105834,    1.899593 ,   2.013569,    2.122117,    2.24152,    2.355495,    2.453189,    2.572592,    2.735414,    2.887381,    3.028494,    3.185889,    3.327002,    3.451832])

N_3 = np.array([0.00290275,    0.02467343,    0.04499274,    0.07547169,    0.1030479,    0.137881,    0.1756168,    0.2017417,    0.2220609,    0.2162554,    0.2423802,    0.2714078,   0.298984,    0.3439768,    0.3773585,    0.425254,    0.4629898,    0.4920174,    0.5341074,    0.5645863,    0.5863571,    0.6023222,    0.6081277,    0.6110305 ])
eps_max_3 = np.array([1.530529,    1.693352,    1.796472,    1.883311,    1.953867,    2.024423,    2.07327,    2.105834,    2.127544 ,   1.943012,    2.05156,    2.165536,    2.246947,    2.371778,    2.469471,    2.572592,    2.643148,    2.719132,    2.843962,    2.985075,    3.12076,    3.251018,    3.359566,    3.451832])

N_4 = np.array([ 0,    0.01886792,    0.04208998,    0.07547169,    0.09724237,    0.1306241,    0.1654572,    0.191582,    0.2177068,    0.2380261,    0.2380261,    0.2525399,    0.275762,    0.2902758,    0.3134978,    0.3251089,    0.3381712,    0.352685,    0.3613933,    0.3671988,    0.3671988 ,   0.3701016])
eps_max_4 = np.array([1.519674,    1.698779  ,  1.796472,    1.921303  ,  1.991859 ,   2.05156  ,  2.111262  ,  2.138399 ,   2.192673  ,  2.203528,    1.970149,    2.002714,    2.084125,    2.165536 ,   2.312076 ,   2.415197  ,  2.529172  ,  2.697422,    2.854817,    3.023067,    3.240163,    3.397558])

N_5 = np.array([0.001451371,    0.01596516,    0.04354136,    0.06966618,    0.1059506,    0.137881,    0.1698113,    0.2046444,    0.2409289,    0.2612482,    0.2583454 ,   0.275762 ,   0.3149492 ,   0.3570392,    0.4005806    ,0.448476 ,   0.5036284,    0.5457184,    0.5892598 ,   0.6342525,    0.6734397 ,   0.7155297,    0.7387518 ,   0.754717 ])
eps_max_5 = np.array([1.514247,    1.644505,   1.763908 ,   1.856174  ,  1.959295  ,  2.029851  ,  2.100407 ,   2.160109  ,  2.203528 ,   2.236092 ,   2.018996,    2.105834,    2.252375,    2.36635  ,  2.464044 ,   2.567164 ,   2.68114 ,   2.757123,    2.854817,    2.963365,    3.077341,    3.223881,    3.343284,    3.451832])


plt.subplot(337)
plt.plot(N_1,  eps_max_1 , "r")
plt.plot(N_2,  eps_max_2 , "r")
plt.plot(N_3,  eps_max_3 , "r")
plt.plot(N_4,  eps_max_4 , "r")
plt.plot(N_5,  eps_max_5 , "r")

#plt.ylim(0.001, 0.004)
plt.xlim(-0.1, 2.0)
plt.ylim(1, 3.5)

#=====================================================================
#creep-fatigue (LH)
N_11 = np.array([0,    0.001452426,    0.01161946,    0.03050108,    0.05083514,    0.06826434,    0.08278866,    0.09731299,    0.1263617,    0.1481481,    0.1597676 ,   0.1713871 ,   0.1713871 ,   0.2265795 ,   0.2846768 ,   0.3572985 ,   0.442992,    0.540305,    0.6129267 ,   0.6986202 ,   0.7901235 ,   0.8496732 ,   0.9034132,    0.9469862,    0.9774873,    0.9949165,    1.025418 ,   1.02687])
eps_max_11 = np.array([     1.084011,    1.192412,    1.333333 ,   1.468835,    1.566396  ,  1.636856 ,   1.685637,    1.723577 ,   1.788618,    1.815718 ,   1.842818,    1.853659,    2.059621  ,  2.124661 ,   2.195122,    2.271003 ,   2.368564,    2.487805,    2.569106,    2.666667 ,   2.802168,    2.899729 ,   3.00813 ,   3.138211 ,   3.241192,    3.344173,    3.479675,    3.528455 ])

N_22 = np.array([ 0,    0.008695646,    0.02463767 ,   0.04347825,    0.07101449  ,  0.09999999  ,  0.1434783,    0.1826087 ,   0.2217391 ,   0.2246377,    0.2492754 ,   0.284058,    0.3463768,    0.4101449,    0.457971 ,   0.5144928,    0.5623189 ,   0.6101449,    0.6536232 ,   0.6956522,    0.7188406,    0.7405797,    0.7492754,    0.7478261])
eps_max_22 = np.array([ 1.113514,    1.248649,    1.362162 ,   1.45946,    1.556757,    1.67027 ,   1.724324 ,   1.778378 ,   1.827027 ,   2.048649,    2.059459,    2.108108,    2.162162,    2.221622,    2.286487 ,   2.345946 ,   2.416216,    2.491892,    2.583784 ,   2.713514 ,   2.794595 ,   2.967568 ,   3.140541  ,  3.389189])

N_33 = np.array([ 0,    0.01594202,    0.04057971,    0.07101449,    0.1173913 ,   0.1594203 ,   0.226087 ,   0.2637681 ,   0.2623188,    0.3304348 ,   0.4130435,    0.4884058,    0.5623189,    0.6492754 ,   0.757971 ,   0.8449275 ,   0.9202899,    0.9898551,    1.037681,    1.094203,    1.16087,    1.211594])
eps_max_33= np.array([ 1.108108,    1.259459,    1.394595,    1.535135,    1.648649,    1.713513,    1.794595,    1.854054,    2.081081,    2.135135,    2.232432,    2.297297,    2.351351,    2.427027,    2.52973,    2.610811,    2.708108,    2.810811,    2.908108,    3.048649,    3.281081,    3.481081 ])

N_44 = np.array([ 0,    0.02753623,    0.04492753,    0.07826087,    0.1115942,    0.1652174,    0.2231884,    0.2623188,    0.2913043,    0.2971015,    0.3478261,    0.4217391,    0.5115942,    0.6014493 ,   0.6739131 ,   0.7463768,    0.8217391 ,   0.8565218,    0.8710145,    0.8855073,    0.9028986,    0.9057971 ])
eps_max_44 = np.array([ 1.064865,    1.264865,    1.389189,    1.491892 ,   1.583784,    1.681081 ,   1.735135 ,   1.794595 ,   1.827027 ,   2.064865 ,   2.113513,    2.189189 ,   2.286487 ,   2.372973 ,   2.45946 ,   2.540541 ,   2.697297 ,   2.783784 ,   2.875676,    2.989189,    3.227027 ,   3.454054 ])

N_55 = np.array([ 0,    0.03043478,    0.06521739,    0.09855072,    0.1376812,    0.1753623,    0.2144928,    0.2521739,    0.2855073,    0.2884058,    0.3927536,    0.5405797,    0.7188406,    0.9782609,    1.176812,    1.369565,    1.626087,    1.901449,    2.088406,    2.202899,    2.350725,    2.398551])
eps_max_55 = np.array([ 1.151351,    1.378378,    1.583784 ,   1.67027,    1.713513 ,   1.767568,    1.827027,    1.864865,    1.891892 ,   2.102703 ,   2.124324 ,   2.178378   , 2.243243 ,   2.356757,   2.443243  ,  2.524324 ,   2.65946  ,  2.805405,    2.918919  ,  2.994595  ,  3.216216 ,   3.383784])

plt.subplot(338)
plt.plot(N_11,  eps_max_11 , "r")
plt.plot(N_22,  eps_max_22 , "r")
plt.plot(N_33,  eps_max_33 , "r")
plt.plot(N_44,  eps_max_44 , "r")
plt.plot(N_55,  eps_max_55 , "r")
plt.xlim(-0.1, 3.0)
plt.ylim(1, 3.5)

#=====================================================================================
# comparison
# H-L
HL_C40_1 = np.array([0.266, 0.39, 0.51, 0.11, 0.25])
HL_C40_2 = np.array([0.154, 0.22, 0.24, 0.26, 0.3])



# L-H
LH_C40_1 = np.array([0.17, 0.22, 0.26, 0.29, 0.3])
LH_C40_2  = np.array([0.86, 0.53, 0.64, 0.93, 2.1])



# average
PM_1 = np.array([1 , 0])
PM_2 = np.array([0,  1])


plt.subplot(339)
plt.plot(HL_C40_1, HL_C40_2, 'ro', markersize=4, color='r')
plt.plot(LH_C40_1, LH_C40_2, 'ro', markersize=4, color='g')
plt.plot(PM_1, PM_2,  color='r',)

plt.ylim(0, 2.2)
plt.xlim(0.0, 1.2)
plt.title('sequence effect')
plt.xlabel('L')
plt.ylabel('H')

plt.show()
