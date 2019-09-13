'''
Created on 28.08.2017

@author: abaktheer
'''


'''
Alliche model- results - two level loading (0.8 - 0.7)
'''
import matplotlib.pyplot as plt
import numpy as np

# #=========================================================================
# # numerical results (creep fatigue),(S 08-07),(H-L)
# #=========================================================================
# ax1 = plt.subplot(221)
#
#
# n_1, eps_1 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_H.txt')
# n_2, eps_2 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_L.txt')
# n_3, eps_3 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_HL_01.txt')
# n_4, eps_4 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_HL_02.txt')
# n_5, eps_5 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_HL_03.txt')
# n_6, eps_6 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_HL_04.txt')
# n_7, eps_7 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_HL_06.txt')
# n_8, eps_8 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_HL_08.txt')
#
#
# ax1.plot(n_1[1:], abs(eps_1[1:] * 1000), "--k", label='H')
# #ax1.plot(n_2[1:], abs(eps_2[1:] * 1000), ":k", label='L')
# ax1.plot(n_3[1:], abs(eps_3[1:] * 1000), "r", label='HL_01')
# ax1.plot(n_4[1:], abs(eps_4[1:] * 1000), "g", label='HL_02')
# ax1.plot(n_5[1:], abs(eps_5[1:] * 1000), "b", label='HL_03')
# ax1.plot(n_6[1:], abs(eps_6[1:] * 1000), "y", label='HL_04')
# ax1.plot(n_7[1:], abs(eps_7[1:] * 1000), "0.6", label='HL_06')
# ax1.plot(n_8[1:], abs(eps_8[1:] * 1000), "0.3", label='HL_08')
#
#
# ax1.set_xlabel('N')
# ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
# #ax1.set_xlim(0, 5)
# #ax1.set_ylim(0.0, 3.0)
# # ax1.legend(loc=2)


# #=========================================================================
# # numerical results (creep fatigue),(S 08-07),(H-L) (damage ratio)
# #=========================================================================
# ax1 = plt.subplot(111)
#
# n_1, eps_1 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_H.txt')
# n_2, eps_2 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_L.txt')
# n_3, eps_3 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_HL_01.txt')
# n_4, eps_4 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_HL_02.txt')
# n_5, eps_5 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_HL_03.txt')
# n_6, eps_6 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_HL_04.txt')
# n_7, eps_7 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_HL_06.txt')
# n_8, eps_8 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_HL_08.txt')
#
#
# # applied number of cycles
# n_01 = 129
# n_02 = 258
# n_03 = 386
# n_04 = 515
# n_06 = 773
# n_08 = 1030
# n_hf = 1288.
# n_lf = 16662.
#
# # normalized accumulative damage
# N_01 = np.hstack(
#     (n_3[0:n_01 + 1] / n_hf, 0.1 + (n_3[n_01 + 1:] - n_3[n_01]) / n_lf))
# N_02 = np.hstack(
#     (n_4[0:n_02 + 1] / n_hf, 0.2 + (n_4[n_02 + 1:] - n_4[n_02]) / n_lf))
# N_03 = np.hstack(
#     (n_5[0:n_03 + 1] / n_hf, 0.3 + (n_5[n_03 + 1:] - n_5[n_03]) / n_lf))
# N_04 = np.hstack(
#     (n_6[0:n_04 + 1] / n_hf, 0.4 + (n_6[n_04 + 1:] - n_6[n_04]) / n_lf))
# N_06 = np.hstack(
#     (n_7[0:n_06 + 1] / n_hf, 0.6 + (n_7[n_06 + 1:] - n_7[n_06]) / n_lf))
# N_08 = np.hstack(
#     (n_8[0:n_08 + 1] / n_hf, 0.8 + (n_8[n_08 + 1:] - n_8[n_08]) / n_lf))
#
# ax1.plot(n_1[1:] / n_1[-1], abs(eps_1[1:] * 1000), "--k", label='H')
# ax1.plot(n_2[1:] / n_2[-1], abs(eps_2[1:] * 1000), ":k", label='L')
# #ax1.plot(N_01[1:], abs(eps_3[1:] * 1000), "r", label='HL_01')
# ax1.plot(N_02[1:], abs(eps_4[1:] * 1000), "r", label='HL_02')
# ax1.plot(N_03[1:], abs(eps_5[1:] * 1000), "g", label='HL_03')
# ax1.plot(N_04[1:], abs(eps_6[1:] * 1000), "b", label='HL_04')
# #ax1.plot(N_06[1:], abs(eps_7[1:] * 1000), "0.6", label='HL_06')
# #ax1.plot(N_08[1:], abs(eps_8[1:] * 1000), "0.3", label='HL_08')
#
#
# ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
# ax1.set_xlabel('cumulative damage ratio $\sum N_i/N^f_i$')
# ax1.set_xlim(0, 2)
# ax1.set_ylim(1.2, 3.4)
# ax1.legend(loc=2)
#


#=========================================================================

# #=========================================================================
# # numerical results (creep fatigue),(S 08-07),(L-H)
# #=========================================================================
# ax1 = plt.subplot(223)
#
#
# n_1, eps_1 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_H.txt')
# n_2, eps_2 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_L.txt')
# n_3, eps_3 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_LH_01.txt')
# n_4, eps_4 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_LH_02.txt')
# n_5, eps_5 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_LH_03.txt')
# n_6, eps_6 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_LH_04.txt')
# n_7, eps_7 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_LH_06.txt')
#
#
# #ax1.plot(n_1[1:], abs(eps_1[1:] * 1000), "--k", label='H')
# ax1.plot(n_2[1:], abs(eps_2[1:] * 1000), "--k", label='L')
# ax1.plot(n_3[1:], abs(eps_3[1:] * 1000), "r", label='HL_01')
# ax1.plot(n_4[1:], abs(eps_4[1:] * 1000), "g", label='HL_02')
# ax1.plot(n_5[1:], abs(eps_5[1:] * 1000), "b", label='HL_03')
# ax1.plot(n_6[1:], abs(eps_6[1:] * 1000), "0.6", label='HL_04')
# ax1.plot(n_7[1:], abs(eps_7[1:] * 1000), "0.3", label='HL_06')
#
# ax1.set_xlabel('N')
# ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
# #ax1.set_xlim(0, 5)
# #ax1.set_ylim(0.0, 3.0)
# ax1.legend(loc=2)

#=========================================================================
# numerical results (creep fatigue),(S 08-07),(L-H) (accumulative damage)
#=========================================================================
ax1 = plt.subplot(111)


n_1, eps_1 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_H.txt')
n_2, eps_2 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_L.txt')
n_3, eps_3 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_LH_01.txt')
n_4, eps_4 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_LH_02.txt')
n_5, eps_5 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_LH_03.txt')
n_6, eps_6 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_LH_04.txt')
n_7, eps_7 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_LH_06.txt')
n_8, eps_8 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_08-07\eps_LH_08.txt')


# applied number of cycles
n_01 = 1662
n_02 = 3332
n_03 = 4999
n_04 = 6665
n_06 = 9997
n_08 = 13330
n_hf = 1288.0
n_lf = 16662.0

# normalized accumulative damage
N_01 = np.hstack(
    (n_3[0:n_01 + 1] / n_lf, 0.1 + (n_3[n_01 + 1:] - n_3[n_01]) / n_hf))
N_02 = np.hstack(
    (n_4[0:n_02 + 1] / n_lf, 0.2 + (n_4[n_02 + 1:] - n_4[n_02]) / n_hf))
N_03 = np.hstack(
    (n_5[0:n_03 + 1] / n_lf, 0.3 + (n_5[n_03 + 1:] - n_5[n_03]) / n_hf))
N_04 = np.hstack(
    (n_6[0:n_04 + 1] / n_lf, 0.4 + (n_6[n_04 + 1:] - n_6[n_04]) / n_hf))
N_06 = np.hstack(
    (n_7[0:n_06 + 1] / n_lf, 0.6 + (n_7[n_06 + 1:] - n_7[n_06]) / n_hf))
N_08 = np.hstack(
    (n_8[0:n_08 + 1] / n_lf, 0.8 + (n_8[n_08 + 1:] - n_8[n_08]) / n_hf))

ax1.plot(n_1[1:] / n_1[-1], abs(eps_1[1:] * 1000), "--k", label='H')
ax1.plot(n_2[1:] / n_2[-1], abs(eps_2[1:] * 1000), ":k", label='L')
#ax1.plot(N_01[1:], abs(eps_3[1:] * 1000), "r", label='HL_01')
ax1.plot(N_02[1:], abs(eps_4[1:] * 1000), "r", label='HL_02')
ax1.plot(N_03[1:], abs(eps_5[1:] * 1000), "g", label='HL_03')
ax1.plot(N_04[1:], abs(eps_6[1:] * 1000), "b", label='HL_04')
#ax1.plot(N_06[1:], abs(eps_7[1:] * 1000), "0.6", label='HL_06')
#ax1.plot(N_08[1:], abs(eps_8[1:] * 1000), "0.3", label='HL_08')

ax1.set_xlabel('cumulative damage ratio $\sum N_i/N^f_i$')
ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')

ax1.set_xlim(0, 2)
ax1.set_ylim(1.2, 3.4)
# ax1.legend(loc=2)


plt.show()
