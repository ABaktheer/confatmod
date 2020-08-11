'''
Created on 27.03.2019

@author: abaktheer
'''


import matplotlib.pyplot as plt
import numpy as np



#========================================
"LS6 - comparison with assessment rules"
#========================================

#==================================================================================================================================================
"#C40"
#===================================================================================
# comparison
# H-L
# HL_C40_1 = np.array([0.18, 0.22, 0.26, 0.27, 0.45, 1.08, 0.57, 2.67, 2.78])
# HL_C40_2 = np.array([0.15, 0.15, 0.15, 0.15, 0.15, 0.15,  0.305,  0.305, 0.89 ])

HL_C40_2 = np.array([0.18, 0.22, 0.26, 0.27, 0.45, 1.08, 0.14, 0.66,  0.69, 0.93, 0.70])

HL_C40_1 = np.array([0.15, 0.15, 0.15, 0.15, 0.15, 0.15,  0.16, 0.16, 0.47, 0.47, 0.47 ])

# L-H
# LH_C40_1 = np.array([0.053, 0.053, 0.053, 0.053, 0.053, 0.053])
# LH_C40_2  = np.array([0.66, 1.02, 0.50, 1.23, 1.67, 2.38])

LH_C40_1 = np.array([0.053, 0.053, 0.053, 0.053, 0.053, 0.053, 0.13, 0.13, 0.18, 0.18, 0.18, 0.36, 0.36])

LH_C40_2  = np.array([0.66, 1.02, 0.50, 1.23, 1.67, 2.38, 1.51, 1.39, 0.9, 0.84, 0.35, 1.92, 1.10])

#===================================================================================================================================================================
"#C80"
#===================================================================================
# comparison
# H-L
HL_C80_2 = np.array([0.32 , 1.41, 2.09, 0.48, 0.51, 0.37, 1.3])
HL_C80_1 = np.array([0.24,  0.24, 0.24, 0.24, 0.47, 0.47, 0.47 ])

# L-H
LH_C80_1 = np.array([0.09, 0.09, 0.09, 0.26, 0.26])
LH_C80_2  = np.array([2.29, 1.2, 3.51, 1.11, 2.33])

#=================================================================================================================================================================================
"Holmen"
#=====================================================================================
# comparison
# H-L
HL_HOL_2 = np.array([0.266, 0.39, 0.51, 0.11, 0.25])
HL_HOL_1 = np.array([0.154, 0.22, 0.24, 0.26, 0.3])

# L-H
LH_HOL_1 = np.array([0.17, 0.22, 0.26, 0.29, 0.3])
LH_HOL_2  = np.array([0.86, 0.53, 0.64, 0.93, 2.1])


#=====================================================================================
# average
PM_1 = np.array([1 , 0])
PM_2 = np.array([0,  1])



a = HL_C40_1 + HL_C40_2
b = HL_C80_1 + HL_C80_2
c =  HL_HOL_1 + HL_HOL_2

d = LH_C40_1 + LH_C40_2
e = LH_C80_1 + LH_C80_2
f =  LH_HOL_1 + LH_HOL_2


average_HL = np.average( np.hstack((a,b,c) ))

average_LH = np.average( np.hstack((d,e,f) ))

#===================================================================================================================================================================
"Mayer"
#===================================================================================
# H-L
HL_MAY_2 = np.array([1.00, 0.91, 0.82, 0.73, 0.64, 0.55, 0.46, 0.34, 0.23, 0.11, 0.00])
HL_MAY_1 = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])

# L-H
LH_MAY_1 = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
LH_MAY_2 = np.array([1.00, 0.89, 0.78, 0.67, 0.56, 0.45, 0.33, 0.27, 0.18, 0.09, 0.00])

#===================================================================================================================================================================
"Shah"
#===================================================================================
# H-L
HL_SHA_2 = np.array([1.00, 0.70, 0.45, 0.31, 0.23, 0.17, 0.13, 0.10, 0.07, 0.04, 0.00])
HL_SHA_1 = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])

# L-H
LH_SHA_1 = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
LH_SHA_2 = np.array([1.00, 0.70, 0.45, 0.31, 0.23, 0.18, 0.13, 0.10, 0.07, 0.04, 0.00])


#===================================================================================================================================================================
"Baktheer & Chudoba - Enhanced rule"
#===================================================================================
# H-L
HL_BAK_1 = np.array([1.00, 0.92, 0.79, 0.65, 0.50, 0.44, 0.38, 0.33, 0.28, 0.17, 0.11, 0.07, 0.05, 0.03, 0.01, 0.00])
HL_BAK_2 = np.array([0, 0.05, 0.1, 0.15, 0.2, 0.225, 0.25, 0.275, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])

# L-H
LH_BAK_1 = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.925, 0.95, 0.975, 1])
LH_BAK_2 = np.array([1.00, 0.95, 0.90, 0.86, 0.81, 0.76, 0.71, 0.67, 0.62, 0.44, 0.33, 0.22, 0.11, 0])



#=================================================================
plt.subplot(131)
plt.plot(HL_C40_1, HL_C40_1 + HL_C40_2, 'ro', markersize=6, color='r')
plt.plot(LH_C40_1, LH_C40_1 + LH_C40_2, 'go', markersize=6, color='g')

plt.plot(HL_C80_1, HL_C80_1 + HL_C80_2, 'rs', markersize=6, color='r')
plt.plot(LH_C80_1, LH_C80_1 + LH_C80_2, 'gs', markersize=6, color='g')

plt.plot(HL_HOL_1, HL_HOL_1 + HL_HOL_2, 'r^', markersize=6, color='r')
plt.plot(LH_HOL_1, LH_HOL_1 + LH_HOL_2, 'g^', markersize=6, color='g')

plt.plot(HL_MAY_1, HL_MAY_1 + HL_MAY_2, color='r')
plt.plot(LH_MAY_1, LH_MAY_1 + LH_MAY_2, color='g')

plt.plot(PM_1, PM_1 + PM_2,  color='k',)
plt.plot(np.array([0 , 1]), np.array([average_HL , average_HL]) ,color='r')
plt.plot(np.array([0 , 1]), np.array([average_LH , average_LH]) ,color='g')

# plt.ylim(0, 2.5)
# plt.xlim(0, 2.2)
plt.title('sequence effect')
plt.xlabel('L')
plt.ylabel('H')


#====================================================================
plt.subplot(132)
plt.plot(HL_C40_1, HL_C40_1 + HL_C40_2, 'ro', markersize=6, color='r')
plt.plot(LH_C40_1, LH_C40_1 + LH_C40_2, 'go', markersize=6, color='g')

plt.plot(HL_C80_1, HL_C80_1 + HL_C80_2, 'rs', markersize=6, color='r')
plt.plot(LH_C80_1, LH_C80_1 + LH_C80_2, 'gs', markersize=6, color='g')

plt.plot(HL_HOL_1, HL_HOL_1 + HL_HOL_2, 'r^', markersize=6, color='r')
plt.plot(LH_HOL_1, LH_HOL_1 + LH_HOL_2, 'g^', markersize=6, color='g')

plt.plot(HL_SHA_1, HL_SHA_1 + HL_SHA_2, color='r')
plt.plot(LH_SHA_1, LH_SHA_1 + LH_SHA_2, color='g')

plt.plot(PM_1, PM_1 + PM_2,  color='k',)
plt.plot(np.array([0 , 1]), np.array([average_HL , average_HL]) ,color='r')
plt.plot(np.array([0 , 1]), np.array([average_LH , average_LH]) ,color='g')

# plt.ylim(0, 2.5)
# plt.xlim(0, 2.2)
plt.title('sequence effect')
plt.xlabel('L')
plt.ylabel('H')


#====================================================================
plt.subplot(133)
plt.plot(HL_C40_1, HL_C40_1 + HL_C40_2, 'ro', markersize=6, color='r')
plt.plot(LH_C40_1, LH_C40_1 + LH_C40_2, 'go', markersize=6, color='g')

plt.plot(HL_C80_1, HL_C80_1 + HL_C80_2, 'rs', markersize=6, color='r')
plt.plot(LH_C80_1, LH_C80_1 + LH_C80_2, 'gs', markersize=6, color='g')

plt.plot(HL_HOL_1, HL_HOL_1 + HL_HOL_2, 'r^', markersize=6, color='r')
plt.plot(LH_HOL_1, LH_HOL_1 + LH_HOL_2, 'g^', markersize=6, color='g')

plt.plot(HL_BAK_1, HL_BAK_1 + HL_BAK_2, color='r')
plt.plot(LH_BAK_1, LH_BAK_1 + LH_BAK_2, color='g')

plt.plot(PM_1, PM_1 + PM_2,  color='k',)

plt.plot(np.array([0 , 1]), np.array([average_HL , average_HL]) ,color='r')
plt.plot(np.array([0 , 1]), np.array([average_LH , average_LH]) ,color='g')
#plt.plot(0, average_LH, 'bo', markersize=10, color='b')

# plt.ylim(0, 2.5)
# plt.xlim(0, 2.2)
plt.title('sequence effect')
plt.xlabel('L')
plt.ylabel('H')


print(average_HL)
print(average_LH)


plt.show()
