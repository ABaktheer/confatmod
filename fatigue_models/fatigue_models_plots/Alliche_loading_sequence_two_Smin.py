'''
Created on 28.08.2017

@author: abaktheer
'''


'''
Alliche model- results - two level loading - Smin = const.
'''
import matplotlib.pyplot as plt
import numpy as np


# #=========================================================================
# # numerical results (comparison with P-M rule)
# #=========================================================================
# ax1 = plt.subplot(111)
#
# n_1, m_1 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\Woehler_H-L.txt')
# n_2, m_2 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\Woehler_L-H.txt')
# n_3, m_3 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\Woehler_H-L.txt')
# n_4, m_4 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\Woehler_L-H.txt')
# n_5, m_5 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\Woehler_H-L.txt')
# n_6, m_6 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\Woehler_L-H.txt')
#
# ax1.plot(n_1, m_1, '--b', label='H-L, $\Delta = 0.05$')
# ax1.plot(n_2, m_2, 'b', label='L-H, $\Delta = 0.05$')
# ax1.plot(n_3, m_3, '--r', label='H-L, $\Delta = 0.1$')
# ax1.plot(n_4, m_4, 'r', label='L-H, $\Delta = 0.1$')
# ax1.plot(n_5, m_5, '--g', label='H-L, $\Delta = 0.2$')
# ax1.plot(n_6, m_6, 'g', label='L-H, $\Delta = 0.2$')
#
#
# # Palmgren-Miner rule
# x_2 = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
# y_2 = np.array([1.0, 0.8, 0.6, 0.4, 0.2, 0.0])
#
# ax1.plot(x_2, y_2, "-k", label='Palmgren-Miner rule')
#
# # ax1.fill_between(x_2, y_2, m_1, facecolor='k', alpha=0.1)
# # ax1.fill_between(x_2, y_2, m_2, facecolor='k', alpha=0.1)
# # ax1.fill_between(x_2, y_2, m_3, facecolor='k', alpha=0.1)
# # ax1.fill_between(x_2, y_2, m_4, facecolor='k', alpha=0.1)
# # ax1.fill_between(x_2, y_2, m_5, facecolor='k', alpha=0.1)
# # ax1.fill_between(x_2, y_2, m_6, facecolor='k', alpha=0.1)
#
# ax1.set_xlabel('$N_1/N_f$')
# ax1.set_ylabel('$N_2/N_f$')
# ax1.set_xlim(0, 1.2)
# ax1.set_ylim(0, 1.2)
# ax1.legend(loc=1)


# #=========================================================================
# # numerical results (creep fatigue),(S 0775-0725),(H-L)
# #=========================================================================
# ax1 = plt.subplot(232)
#
# n_1, eps_1 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\eps_H.txt')
# n_2, eps_2 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\eps_L.txt')
# n_3, eps_3 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\eps_HL_02.txt')
# n_4, eps_4 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\eps_HL_04.txt')
# n_5, eps_5 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\eps_HL_06.txt')
# n_6, eps_6 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\eps_HL_08.txt')
#
# ax1.plot(n_1[1:], abs(eps_1[1:] * 1000), "--k", label='H')
# #ax1.plot(n_2[1:], abs(eps_2[1:] * 1000), ":k", label='L')
# ax1.plot(n_3[1:], abs(eps_3[1:] * 1000), "r", label='HL_02')
# ax1.plot(n_4[1:], abs(eps_4[1:] * 1000), "g", label='HL_04')
# ax1.plot(n_5[1:], abs(eps_5[1:] * 1000), "b", label='HL_06')
# ax1.plot(n_6[1:], abs(eps_6[1:] * 1000), "0.6", label='HL_08')
#
# ax1.set_xlabel('N')
# ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
# #ax1.set_xlim(0, 5)
# #ax1.set_ylim(0.0, 3.0)
# # ax1.legend(loc=2)
#
#
#=========================================================================
# numerical results (creep fatigue),(S 0775-0725),(H-L) (damage ratio)
#=========================================================================
ax1 = plt.subplot(111)

n_1, eps_1 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\eps_H.txt')
n_2, eps_2 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\eps_L.txt')
n_3, eps_3 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\eps_HL_02.txt')
n_4, eps_4 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\eps_HL_04.txt')
n_5, eps_5 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\eps_HL_06.txt')
n_6, eps_6 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\eps_HL_08.txt')

# applied number of cycles

n_02 = 479
n_04 = 958
n_06 = 1437
n_08 = 1916
n_hf = 2395.
n_lf = 8550.

# normalized accumulative damage

N_02 = np.hstack(
    (n_3[0:n_02 + 1] / n_hf, 0.2 + (n_3[n_02 + 1:] - n_3[n_02]) / n_lf))
N_04 = np.hstack(
    (n_4[0:n_04 + 1] / n_hf, 0.4 + (n_4[n_04 + 1:] - n_4[n_04]) / n_lf))
N_06 = np.hstack(
    (n_5[0:n_06 + 1] / n_hf, 0.6 + (n_5[n_06 + 1:] - n_5[n_06]) / n_lf))
N_08 = np.hstack(
    (n_6[0:n_08 + 1] / n_hf, 0.8 + (n_6[n_08 + 1:] - n_6[n_08]) / n_lf))

ax1.plot(n_1[1:] / n_1[-1], abs(eps_1[1:] * 1000), "--k", label='H')
ax1.plot(n_2[1:] / n_2[-1], abs(eps_2[1:] * 1000), ":k", label='L')
ax1.plot(N_02[1:], abs(eps_3[1:] * 1000), "r", label='HL_02')
#ax1.plot(N_04[1:], abs(eps_4[1:] * 1000), "b", label='HL_04')
ax1.plot(N_06[1:], abs(eps_5[1:] * 1000), "g", label='HL_06')
#ax1.plot(N_08[1:], abs(eps_6[1:] * 1000), "0.6", label='HL_08')

ax1.set_title('H-L, $\Delta S = 0.05$')
ax1.set_xlabel('N')
ax1.set_xlabel('cumulative damage ratio $\sum N_i/N^f_i$')
ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
ax1.set_xlim(0, 1.2)
ax1.set_ylim(1.2, 3.4)
ax1.legend(loc=2)
#
#
# #=========================================================================
# # numerical results (creep fatigue),(S 0775-0725),(L-H)
# #=========================================================================
# ax1 = plt.subplot(235)
#
# n_1, eps_1 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\eps_H.txt')
# n_2, eps_2 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\eps_L.txt')
# n_3, eps_3 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\eps_LH_02.txt')
# n_4, eps_4 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\eps_LH_04.txt')
# n_5, eps_5 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\eps_LH_06.txt')
# n_6, eps_6 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\eps_LH_08.txt')
#
# #ax1.plot(n_1[1:], abs(eps_1[1:] * 1000), "--k", label='H')
# ax1.plot(n_2[1:], abs(eps_2[1:] * 1000), ":k", label='L')
# ax1.plot(n_3[1:], abs(eps_3[1:] * 1000), "r", label='LH_02')
# ax1.plot(n_4[1:], abs(eps_4[1:] * 1000), "g", label='LH_04')
# ax1.plot(n_5[1:], abs(eps_5[1:] * 1000), "b", label='LH_06')
# ax1.plot(n_6[1:], abs(eps_6[1:] * 1000), "0.6", label='LH_08')
#
# ax1.set_xlabel('N')
# ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
# #ax1.set_xlim(0, 5)
# #ax1.set_ylim(0.0, 3.0)
# # ax1.legend(loc=2)
#
#
# #=========================================================================
# # numerical results (creep fatigue),(S 0775-0725),(H-L) (damage ratio)
# #=========================================================================
# ax1 = plt.subplot(111)
#
# n_1, eps_1 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\eps_H.txt')
# n_2, eps_2 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\eps_L.txt')
# n_3, eps_3 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\eps_LH_02.txt')
# n_4, eps_4 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\eps_LH_04.txt')
# n_5, eps_5 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\eps_LH_06.txt')
# n_6, eps_6 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\eps_LH_08.txt')
#
# # applied number of cycles
#
# n_02 = 1710
# n_04 = 3420
# n_06 = 5130
# n_08 = 6840
# n_hf = 2395.
# n_lf = 8550.
#
# # normalized accumulative damage
#
# N_02 = np.hstack(
#     (n_3[0:n_02 + 1] / n_lf, 0.2 + (n_3[n_02 + 1:] - n_3[n_02]) / n_hf))
# N_04 = np.hstack(
#     (n_4[0:n_04 + 1] / n_lf, 0.4 + (n_4[n_04 + 1:] - n_4[n_04]) / n_hf))
# N_06 = np.hstack(
#     (n_5[0:n_06 + 1] / n_lf, 0.6 + (n_5[n_06 + 1:] - n_5[n_06]) / n_hf))
# N_08 = np.hstack(
#     (n_6[0:n_08 + 1] / n_lf, 0.8 + (n_6[n_08 + 1:] - n_6[n_08]) / n_hf))
#
# ax1.plot(n_1[1:] / n_1[-1], abs(eps_1[1:] * 1000), "--k", label='H')
# ax1.plot(n_2[1:] / n_2[-1], abs(eps_2[1:] * 1000), ":k", label='L')
# ax1.plot(N_02[1:], abs(eps_3[1:] * 1000), "r", label='LH_02')
# #ax1.plot(N_04[1:], abs(eps_4[1:] * 1000), "b", label='LH_04')
# ax1.plot(N_06[1:], abs(eps_5[1:] * 1000), "g", label='LH_06')
# #ax1.plot(N_08[1:], abs(eps_6[1:] * 1000), "0.6", label='LH_08')
#
# ax1.set_title('L-H, $\Delta S = 0.05$')
# ax1.set_xlabel('cumulative damage ratio $\sum N_i/N^f_i$')
# ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
# ax1.set_xlim(0, 1.5)
# ax1.set_ylim(1.2, 3.4)
# ax1.legend(loc=2)


#=========================================================================


# #=========================================================================
# # numerical results (creep fatigue),(S 08-07),(H-L)
# #=========================================================================
# ax1 = plt.subplot(232)
#
# n_1, eps_1 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\eps_H.txt')
# n_2, eps_2 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\eps_L.txt')
# n_3, eps_3 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\eps_HL_02.txt')
# n_4, eps_4 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\eps_HL_04.txt')
# n_5, eps_5 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\eps_HL_06.txt')
# n_6, eps_6 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\eps_HL_08.txt')
#
# ax1.plot(n_1[1:], abs(eps_1[1:] * 1000), "--k", label='H')
# #ax1.plot(n_2[1:], abs(eps_2[1:] * 1000), ":k", label='L')
# ax1.plot(n_3[1:], abs(eps_3[1:] * 1000), "r", label='HL_02')
# ax1.plot(n_4[1:], abs(eps_4[1:] * 1000), "g", label='HL_04')
# ax1.plot(n_5[1:], abs(eps_5[1:] * 1000), "b", label='HL_06')
# ax1.plot(n_6[1:], abs(eps_6[1:] * 1000), "0.6", label='HL_08')
#
# ax1.set_xlabel('N')
# ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
# #ax1.set_xlim(0, 5)
# #ax1.set_ylim(0.0, 3.0)
# # ax1.legend(loc=2)
#
#
#=========================================================================
# numerical results (creep fatigue),(S 08-07),(H-L) (damage ratio)
#=========================================================================
ax1 = plt.subplot(121)

n_1, eps_1 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\eps_H.txt')
n_2, eps_2 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\eps_L.txt')
n_3, eps_3 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\eps_HL_02.txt')
n_4, eps_4 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\eps_HL_04.txt')
n_5, eps_5 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\eps_HL_06.txt')
n_6, eps_6 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\eps_HL_08.txt')

# applied number of cycles

n_02 = 258
n_04 = 515
n_06 = 773
n_08 = 1030
n_hf = 1288.
n_lf = 16662.

# normalized accumulative damage

N_02 = np.hstack(
    (n_3[0:n_02 + 1] / n_hf, 0.2 + (n_3[n_02 + 1:] - n_3[n_02]) / n_lf))
N_04 = np.hstack(
    (n_4[0:n_04 + 1] / n_hf, 0.4 + (n_4[n_04 + 1:] - n_4[n_04]) / n_lf))
N_06 = np.hstack(
    (n_5[0:n_06 + 1] / n_hf, 0.6 + (n_5[n_06 + 1:] - n_5[n_06]) / n_lf))
N_08 = np.hstack(
    (n_6[0:n_08 + 1] / n_hf, 0.8 + (n_6[n_08 + 1:] - n_6[n_08]) / n_lf))

ax1.plot(n_1[1:] / n_1[-1], abs(eps_1[1:] * 1000), "--k", label='H')
ax1.plot(n_2[1:] / n_2[-1], abs(eps_2[1:] * 1000), ":k", label='L')
ax1.plot(N_02[1:], abs(eps_3[1:] * 1000), "r", label='HL_02')
ax1.plot(N_04[1:], abs(eps_4[1:] * 1000), "g", label='HL_04')
ax1.plot(N_06[1:], abs(eps_5[1:] * 1000), "b", label='HL_06')
ax1.plot(N_08[1:], abs(eps_6[1:] * 1000), "0.6", label='HL_08')


ax1.set_xlabel('N')
ax1.set_xlabel('cumulative damage ratio $\sum N_i/N^f_i$')
ax1.set_xlim(0, 2)
ax1.set_ylim(1.2, 3.4)
# ax1.legend(loc=2)
#
#
# #=========================================================================
# # numerical results (creep fatigue),(S 08-07),(L-H)
# #=========================================================================
# ax1 = plt.subplot(235)
#
# n_1, eps_1 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\eps_H.txt')
# n_2, eps_2 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\eps_L.txt')
# n_3, eps_3 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\eps_LH_02.txt')
# n_4, eps_4 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\eps_LH_04.txt')
# n_5, eps_5 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\eps_LH_06.txt')
# n_6, eps_6 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\eps_LH_08.txt')
#
# #ax1.plot(n_1[1:], abs(eps_1[1:] * 1000), "--k", label='H')
# ax1.plot(n_2[1:], abs(eps_2[1:] * 1000), ":k", label='L')
# ax1.plot(n_3[1:], abs(eps_3[1:] * 1000), "r", label='LH_02')
# ax1.plot(n_4[1:], abs(eps_4[1:] * 1000), "g", label='LH_04')
# ax1.plot(n_5[1:], abs(eps_5[1:] * 1000), "b", label='LH_06')
# ax1.plot(n_6[1:], abs(eps_6[1:] * 1000), "0.6", label='LH_08')
#
# ax1.set_xlabel('N')
# ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
# #ax1.set_xlim(0, 5)
# #ax1.set_ylim(0.0, 3.0)
# # ax1.legend(loc=2)
#
#
#=========================================================================
# numerical results (creep fatigue),(S 08-07),(H-L) (damage ratio)
#=========================================================================
ax1 = plt.subplot(122)

n_1, eps_1 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\eps_H.txt')
n_2, eps_2 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\eps_L.txt')
n_3, eps_3 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\eps_LH_02.txt')
n_4, eps_4 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\eps_LH_04.txt')
n_5, eps_5 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\eps_LH_06.txt')
n_6, eps_6 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\eps_LH_08.txt')

# applied number of cycles

n_02 = 3332
n_04 = 6665
n_06 = 9997
n_08 = 13330
n_hf = 1288.
n_lf = 16662.

# normalized accumulative damage

N_02 = np.hstack(
    (n_3[0:n_02 + 1] / n_lf, 0.2 + (n_3[n_02 + 1:] - n_3[n_02]) / n_hf))
N_04 = np.hstack(
    (n_4[0:n_04 + 1] / n_lf, 0.4 + (n_4[n_04 + 1:] - n_4[n_04]) / n_hf))
N_06 = np.hstack(
    (n_5[0:n_06 + 1] / n_lf, 0.6 + (n_5[n_06 + 1:] - n_5[n_06]) / n_hf))
N_08 = np.hstack(
    (n_6[0:n_08 + 1] / n_lf, 0.8 + (n_6[n_08 + 1:] - n_6[n_08]) / n_hf))

ax1.plot(n_1[1:] / n_1[-1], abs(eps_1[1:] * 1000), "--k", label='H')
ax1.plot(n_2[1:] / n_2[-1], abs(eps_2[1:] * 1000), ":k", label='L')
ax1.plot(N_02[1:], abs(eps_3[1:] * 1000), "r", label='LH_02')
ax1.plot(N_04[1:], abs(eps_4[1:] * 1000), "g", label='LH_04')
ax1.plot(N_06[1:], abs(eps_5[1:] * 1000), "b", label='LH_06')
ax1.plot(N_08[1:], abs(eps_6[1:] * 1000), "0.6", label='LH_08')


ax1.set_xlabel('N')
ax1.set_xlabel('cumulative damage ratio $\sum N_i/N^f_i$')
ax1.set_xlim(0, 2)
ax1.set_ylim(1.2, 3.4)
# ax1.legend(loc=2)

#=========================================================================

# #=========================================================================
# # numerical results (creep fatigue),(S 085-065),(H-L)
# #=========================================================================
# ax1 = plt.subplot(232)
#
# n_1, eps_1 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\eps_H.txt')
# n_2, eps_2 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\eps_L.txt')
# n_3, eps_3 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\eps_HL_02.txt')
# n_4, eps_4 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\eps_HL_04.txt')
# n_5, eps_5 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\eps_HL_06.txt')
# n_6, eps_6 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\eps_HL_08.txt')
#
# ax1.plot(n_1[1:], abs(eps_1[1:] * 1000), "--k", label='H')
# #ax1.plot(n_2[1:], abs(eps_2[1:] * 1000), ":k", label='L')
# ax1.plot(n_3[1:], abs(eps_3[1:] * 1000), "r", label='HL_02')
# #ax1.plot(n_4[1:], abs(eps_4[1:] * 1000), "g", label='HL_04')
# ax1.plot(n_5[1:], abs(eps_5[1:] * 1000), "b", label='HL_06')
# #ax1.plot(n_6[1:], abs(eps_6[1:] * 1000), "0.6", label='HL_08')
#
# ax1.set_xlabel('N')
# ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
# #ax1.set_xlim(0, 5)
# ax1.set_ylim(1.2, 3.3)
# # ax1.legend(loc=2)
#
#
# #=========================================================================
# # numerical results (creep fatigue),(S 085-065),(H-L) (damage ratio)
# #=========================================================================
# ax1 = plt.subplot(111)
#
# n_1, eps_1 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\eps_H.txt')
# n_2, eps_2 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\eps_L.txt')
# n_3, eps_3 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\eps_HL_02.txt')
# n_4, eps_4 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\eps_HL_04.txt')
# n_5, eps_5 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\eps_HL_06.txt')
# n_6, eps_6 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\eps_HL_08.txt')
#
# # applied number of cycles
#
# n_02 = 72
# n_04 = 144
# n_06 = 216
# n_08 = 288
# n_hf = 360.
# n_lf = 68848.
#
# # normalized accumulative damage
#
# N_02 = np.hstack(
#     (n_3[0:n_02 + 1] / n_hf, 0.2 + (n_3[n_02 + 1:] - n_3[n_02]) / n_lf))
# N_04 = np.hstack(
#     (n_4[0:n_04 + 1] / n_hf, 0.4 + (n_4[n_04 + 1:] - n_4[n_04]) / n_lf))
# N_06 = np.hstack(
#     (n_5[0:n_06 + 1] / n_hf, 0.6 + (n_5[n_06 + 1:] - n_5[n_06]) / n_lf))
# N_08 = np.hstack(
#     (n_6[0:n_08 + 1] / n_hf, 0.8 + (n_6[n_08 + 1:] - n_6[n_08]) / n_lf))
#
# ax1.plot(n_1[1:] / n_1[-1], abs(eps_1[1:] * 1000), "--k", label='H')
# ax1.plot(n_2[1:] / n_2[-1], abs(eps_2[1:] * 1000), ":k", label='L')
# ax1.plot(N_02[1:], abs(eps_3[1:] * 1000), "r", label='HL_02')
# #ax1.plot(N_04[1:], abs(eps_4[1:] * 1000), "b", label='HL_04')
# ax1.plot(N_06[1:], abs(eps_5[1:] * 1000), "g", label='HL_06')
# #ax1.plot(N_08[1:], abs(eps_6[1:] * 1000), "0.6", label='HL_08')
#
#
# ax1.set_title('H-L, $\Delta S = 0.2$')
# ax1.set_xlabel('cumulative damage ratio $\sum N_i/N^f_i$')
# ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
# ax1.set_xlim(0, 1.1)
# ax1.set_ylim(1.2, 3.4)
# ax1.legend(loc=2)

#
# #=========================================================================
# # numerical results (creep fatigue),(S 085-065),(L-H)
# #=========================================================================
# ax1 = plt.subplot(235)
#
# n_1, eps_1 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\eps_H.txt')
# n_2, eps_2 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\eps_L.txt')
# n_3, eps_3 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\eps_LH_02.txt')
# n_4, eps_4 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\eps_LH_04.txt')
# n_5, eps_5 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\eps_LH_06.txt')
# n_6, eps_6 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\eps_LH_08.txt')
#
# #ax1.plot(n_1[1:], abs(eps_1[1:] * 1000), "--k", label='H')
# ax1.plot(n_2[1:], abs(eps_2[1:] * 1000), ":k", label='L')
# ax1.plot(n_3[1:], abs(eps_3[1:] * 1000), "r", label='LH_02')
# #ax1.plot(n_4[1:], abs(eps_4[1:] * 1000), "g", label='LH_04')
# ax1.plot(n_5[1:], abs(eps_5[1:] * 1000), "b", label='LH_06')
# #ax1.plot(n_6[1:], abs(eps_6[1:] * 1000), "0.6", label='LH_08')
#
# ax1.set_xlabel('N')
# ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
# #ax1.set_xlim(0, 5)
# ax1.set_ylim(1.2, 3.3)
# ax1.legend(loc=2)
#
#
# #=========================================================================
# # numerical results (creep fatigue),(S 085-065),(H-L) (damage ratio)
# #=========================================================================
# ax1 = plt.subplot(111)
#
# n_1, eps_1 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\eps_H.txt')
# n_2, eps_2 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\eps_L.txt')
# n_3, eps_3 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\eps_LH_02.txt')
# n_4, eps_4 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\eps_LH_04.txt')
# n_5, eps_5 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\eps_LH_06.txt')
# n_6, eps_6 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\eps_LH_08.txt')
#
# # applied number of cycles
#
# n_02 = 13770
# n_04 = 27539
# n_06 = 41309
# n_08 = 55078
# n_hf = 360.
# n_lf = 68848.
#
# # normalized accumulative damage
#
# N_02 = np.hstack(
#     (n_3[0:n_02 + 1] / n_lf, 0.2 + (n_3[n_02 + 1:] - n_3[n_02]) / n_hf))
# N_04 = np.hstack(
#     (n_4[0:n_04 + 1] / n_lf, 0.4 + (n_4[n_04 + 1:] - n_4[n_04]) / n_hf))
# N_06 = np.hstack(
#     (n_5[0:n_06 + 1] / n_lf, 0.6 + (n_5[n_06 + 1:] - n_5[n_06]) / n_hf))
# N_08 = np.hstack(
#     (n_6[0:n_08 + 1] / n_lf, 0.8 + (n_6[n_08 + 1:] - n_6[n_08]) / n_hf))
#
# ax1.plot(n_1[1:] / n_1[-1], abs(eps_1[1:] * 1000), "--k", label='H')
# ax1.plot(n_2[1:] / n_2[-1], abs(eps_2[1:] * 1000), ":k", label='L')
# ax1.plot(N_02[1:], abs(eps_3[1:] * 1000), "r", label='LH_02')
# #ax1.plot(N_04[1:], abs(eps_4[1:] * 1000), "b", label='LH_04')
# ax1.plot(N_06[1:], abs(eps_5[1:] * 1000), "g", label='LH_06')
# #ax1.plot(N_08[1:], abs(eps_6[1:] * 1000), "0.6", label='LH_08')
#
#
# ax1.set_title('H-L, $\Delta S = 0.2$')
# ax1.set_xlabel('cumulative damage ratio $\sum N_i/N^f_i$')
# ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
# ax1.set_xlim(0, 1.5)
# ax1.set_ylim(1.2, 3.4)
# ax1.legend(loc=2)


plt.show()
