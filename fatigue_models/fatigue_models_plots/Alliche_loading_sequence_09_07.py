'''
Created on 28.08.2017

@author: abaktheer
'''


'''
Alliche model- results - two level loading (0.9 -0.7)
'''
import matplotlib.pyplot as plt
import numpy as np

#=========================================================================
# numerical results (creep fatigue),(S 09-08),(H-L)
#=========================================================================
# ax1 = plt.subplot(231)
#
#
# n_1, eps_1 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_H.txt')
# n_2, eps_2 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_L.txt')
# n_3, eps_3 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_HL_02.txt')
# n_4, eps_4 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_HL_04.txt')
# n_5, eps_5 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_HL_06.txt')
# n_6, eps_6 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_HL_08.txt')
#
#
# ax1.plot(n_1[1:], abs(eps_1[1:] * 1000), "--k", label='H')
# ax1.plot(n_2[1:], abs(eps_2[1:] * 1000), ":k", label='L')
# ax1.plot(n_3[1:], abs(eps_3[1:] * 1000), "r", label='HL_02')
# ax1.plot(n_4[1:], abs(eps_4[1:] * 1000), "g", label='HL_04')
# ax1.plot(n_5[1:], abs(eps_5[1:] * 1000), "b", label='HL_06')
# ax1.plot(n_6[1:], abs(eps_6[1:] * 1000), "y", label='HL_08')
#
#
# ax1.set_xlabel('N')
# ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
# #ax1.set_xlim(0, 5)
# #ax1.set_ylim(0.0, 3.0)
# # ax1.legend(loc=2)


# #=========================================================================
# # numerical results (creep fatigue),(S 09-08),(H-L) (damage ratio)
# #=========================================================================
# ax1 = plt.subplot(231)
#
#
# n_1, eps_1 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_H.txt')
# n_2, eps_2 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_L.txt')
# n_3, eps_3 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_HL_02.txt')
# n_4, eps_4 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_HL_04.txt')
# n_5, eps_5 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_HL_06.txt')
# n_6, eps_6 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_HL_08.txt')
#
#
# # applied number of cycles
# n_02 = 19
# n_04 = 38
# n_06 = 57
# n_08 = 76
# n_hf = 95.0
# n_lf = 1288.0
#
# # normalized accumulative damage
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
# ax1.plot(N_04[1:], abs(eps_4[1:] * 1000), "g", label='HL_04')
# ax1.plot(N_06[1:], abs(eps_5[1:] * 1000), "b", label='HL_06')
# ax1.plot(N_08[1:], abs(eps_6[1:] * 1000), "y", label='HL_08')
#
#
# ax1.set_xlabel('N')
# ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
# #ax1.set_xlim(0, 5)
# #ax1.set_ylim(0.0, 3.0)
# # ax1.legend(loc=2)
#
# #=========================================================================
# # numerical results (creep fatigue),(S 09-08),(H-L) (log)
# #=========================================================================
# ax1 = plt.subplot(232)
#
#
# n_1, eps_1 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_H.txt')
# n_2, eps_2 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_L.txt')
# n_3, eps_3 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_HL_02.txt')
# n_4, eps_4 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_HL_04.txt')
# n_5, eps_5 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_HL_06.txt')
# n_6, eps_6 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_HL_08.txt')
#
#
# ax1.plot(np.log10(n_1[1:]), abs(eps_1[1:] * 1000), "--k", label='H')
# ax1.plot(np.log10(n_2[1:]), abs(eps_2[1:] * 1000), ":k", label='L')
# ax1.plot(np.log10(n_3[1:]), abs(eps_3[1:] * 1000), "r", label='HL_02')
# ax1.plot(np.log10(n_4[1:]), abs(eps_4[1:] * 1000), "g", label='HL_04')
# ax1.plot(np.log10(n_5[1:]), abs(eps_5[1:] * 1000), "b", label='HL_06')
# ax1.plot(np.log10(n_6[1:]), abs(eps_6[1:] * 1000), "y", label='HL_08')
#
#
# ax1.set_xlabel('log(N)')
# ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
# #ax1.set_xlim(0, 5)
# #ax1.set_ylim(0.0, 3.0)
# ax1.legend(loc=2)
#
# #=========================================================================
# # numerical results (creep fatigue),(S 09-08),(H-L), (Normalized N)
# #=========================================================================
# ax3 = plt.subplot(233)
#
# n_1, eps_1 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_H.txt')
# n_2, eps_2 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_L.txt')
# n_3, eps_3 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_HL_02.txt')
# n_4, eps_4 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_HL_04.txt')
# n_5, eps_5 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_HL_06.txt')
# n_6, eps_6 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_HL_08.txt')
#
#
# ax3.plot(n_1[1:] / n_1[-1], abs(eps_1[1:] * 1000), "--k", label='H')
# ax3.plot(n_2[1:] / n_2[-1], abs(eps_2[1:] * 1000), ":k", label='L')
# ax3.plot(n_3[1:] / n_3[-1], abs(eps_3[1:] * 1000), "r", label='HL_02')
# ax3.plot(n_4[1:] / n_4[-1], abs(eps_4[1:] * 1000), "g", label='HL_04')
# ax3.plot(n_5[1:] / n_5[-1], abs(eps_5[1:] * 1000), "b", label='HL_06')
# ax3.plot(n_6[1:] / n_6[-1], abs(eps_6[1:] * 1000), "y", label='HL_08')
#
#
# ax3.set_xlabel('N/Nf')
# ax3.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
# #ax1.set_xlim(0, 5)
# #ax1.set_ylim(0.0, 3.0)
# # ax3.legend(loc=2)
#
#
# # #=========================================================================
# # # numerical results (creep fatigue),(S 09-08),(L-H)
# # #=========================================================================
# ax1 = plt.subplot(234)
#
#
# n_1, eps_1 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_H.txt')
# n_2, eps_2 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_L.txt')
# n_3, eps_3 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_LH_02.txt')
# n_4, eps_4 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_LH_04.txt')
# n_5, eps_5 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_LH_06.txt')
# n_6, eps_6 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_LH_08.txt')
#
# ax1.plot(n_1[1:], abs(eps_1[1:] * 1000), "--k", label='H')
# ax1.plot(n_2[1:], abs(eps_2[1:] * 1000), ":k", label='L')
# ax1.plot(n_3[1:], abs(eps_3[1:] * 1000), "r", label='LH_02')
# ax1.plot(n_4[1:], abs(eps_4[1:] * 1000), "g", label='LH_04')
# ax1.plot(n_5[1:], abs(eps_5[1:] * 1000), "b", label='LH_06')
# ax1.plot(n_6[1:], abs(eps_6[1:] * 1000), "y", label='LH_08')
#
#
# ax1.set_xlabel('N')
# ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
# #ax1.set_xlim(0, 5)
# #ax1.set_ylim(0.0, 3.0)
# # ax1.legend(loc=2)
#
#
# #=========================================================================
# # # numerical results (creep fatigue),(S 09-08),(L-H) (log)
# # #=========================================================================
# ax1 = plt.subplot(235)
#
#
# n_1, eps_1 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_H.txt')
# n_2, eps_2 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_L.txt')
# n_3, eps_3 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_LH_02.txt')
# n_4, eps_4 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_LH_04.txt')
# n_5, eps_5 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_LH_06.txt')
# n_6, eps_6 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_LH_08.txt')
#
# ax1.plot(np.log10(n_1[1:]), abs(eps_1[1:] * 1000), "--k", label='H')
# ax1.plot(np.log10(n_2[1:]), abs(eps_2[1:] * 1000), ":k", label='L')
# ax1.plot(np.log10(n_3[1:]), abs(eps_3[1:] * 1000), "r", label='LH_02')
# ax1.plot(np.log10(n_4[1:]), abs(eps_4[1:] * 1000), "g", label='LH_04')
# ax1.plot(np.log10(n_5[1:]), abs(eps_5[1:] * 1000), "b", label='LH_06')
# ax1.plot(np.log10(n_6[1:]), abs(eps_6[1:] * 1000), "y", label='LH_08')
#
#
# ax1.set_xlabel('N')
# ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
# #ax1.set_xlim(0, 5)
# #ax1.set_ylim(0.0, 3.0)
# ax1.legend(loc=2)
#
# # #=========================================================================
# # # numerical results (creep fatigue),(S 09-08),(L-H), (normalized N)
# # #=========================================================================
# ax4 = plt.subplot(236)
#
#
# n_1, eps_1 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_H.txt')
# n_2, eps_2 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_L.txt')
# n_3, eps_3 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_LH_02.txt')
# n_4, eps_4 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_LH_04.txt')
# n_5, eps_5 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_LH_06.txt')
# n_6, eps_6 = np.loadtxt(
#     'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-08\eps_LH_08.txt')
#
# ax4.plot(n_1[1:] / n_1[-1], abs(eps_1[1:] * 1000), "--k", label='H')
# ax4.plot(n_2[1:] / n_2[-1], abs(eps_2[1:] * 1000), ":k", label='L')
# ax4.plot(n_3[1:] / n_3[-1], abs(eps_3[1:] * 1000), "r", label='LH_02')
# ax4.plot(n_4[1:] / n_4[-1], abs(eps_4[1:] * 1000), "g", label='LH_04')
# ax4.plot(n_5[1:] / n_5[-1], abs(eps_5[1:] * 1000), "b", label='LH_06')
# ax4.plot(n_6[1:] / n_6[-1], abs(eps_6[1:] * 1000), "y", label='LH_08')
#
#
# ax4.set_xlabel('N/Nf')
# ax4.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
# #ax1.set_xlim(0, 5)
# #ax1.set_ylim(0.0, 3.0)
# # ax4.legend(loc=2)

#=========================================================================

#=========================================================================
# numerical results (creep fatigue),(S 09-07),(H-L)
#=========================================================================
ax1 = plt.subplot(241)


n_1, eps_1 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_H.txt')
n_2, eps_2 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_L.txt')
n_3, eps_3 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_HL_02.txt')
n_4, eps_4 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_HL_04.txt')
n_5, eps_5 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_HL_06.txt')
n_6, eps_6 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_HL_08.txt')


ax1.plot(n_1[1:], abs(eps_1[1:] * 1000), "--k", label='H')
#ax1.plot(n_2[1:], abs(eps_2[1:] * 1000), "k", label='L')
ax1.plot(n_3[1:], abs(eps_3[1:] * 1000), "r", label='HL_02')
ax1.plot(n_4[1:], abs(eps_4[1:] * 1000), "g", label='HL_04')
ax1.plot(n_5[1:], abs(eps_5[1:] * 1000), "b", label='HL_06')
ax1.plot(n_6[1:], abs(eps_6[1:] * 1000), "y", label='HL_08')

ax1.set_xlabel('N')
ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
#ax1.set_xlim(0, 5)
#ax1.set_ylim(0.0, 3.0)
ax1.legend(loc=2)

#=========================================================================
# numerical results (creep fatigue),(S 09-07),(H-L) (accumulative damage)
#=========================================================================
ax1 = plt.subplot(242)


n_1, eps_1 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_H.txt')
n_2, eps_2 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_L.txt')
n_3, eps_3 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_HL_02.txt')
n_4, eps_4 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_HL_04.txt')
n_5, eps_5 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_HL_06.txt')
n_6, eps_6 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_HL_08.txt')


# applied number of cycles
n_02 = 19
n_04 = 38
n_06 = 57
n_08 = 76
n_hf = 95.0
n_lf = 16662.0

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
#ax1.plot(n_2[1:] / n_2[-1], abs(eps_2[1:] * 1000), ":k", label='L')
ax1.plot(N_02[1:], abs(eps_3[1:] * 1000), "r", label='HL_02')
ax1.plot(N_04[1:], abs(eps_4[1:] * 1000), "g", label='HL_04')
ax1.plot(N_06[1:], abs(eps_5[1:] * 1000), "b", label='HL_06')
ax1.plot(N_08[1:], abs(eps_6[1:] * 1000), "y", label='HL_08')

ax1.set_xlabel('cumulative damage ratio $\sum N_i/N^f_i$')
ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')

# #ax1.set_xlim(0, 5)
ax1.set_ylim(1.8, 3.3)
# ax1.legend(loc=2)

#=========================================================================
# numerical results (creep fatigue),(S 09-07),(H-L)
#=========================================================================
ax1 = plt.subplot(243)


n_1, eps_1 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_H.txt')
n_2, eps_2 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_L.txt')
n_3, eps_3 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_HL_02.txt')
n_4, eps_4 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_HL_04.txt')
n_5, eps_5 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_HL_06.txt')
n_6, eps_6 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_HL_08.txt')

ax1.plot(np.log10(n_1[1:]), abs(eps_1[1:] * 1000), "--k", label='H')
#ax1.plot(np.log10(n_2[1:]), abs(eps_2[1:] * 1000), "k", label='L')
ax1.plot(np.log10(n_3[1:]), abs(eps_3[1:] * 1000), "r", label='HL_02')
ax1.plot(np.log10(n_4[1:]), abs(eps_4[1:] * 1000), "g", label='HL_04')
ax1.plot(np.log10(n_5[1:]), abs(eps_5[1:] * 1000), "b", label='HL_06')
ax1.plot(np.log10(n_6[1:]), abs(eps_6[1:] * 1000), "y", label='HL_08')

ax1.set_xlabel('log(N)')
ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
#ax1.set_xlim(0, 5)
#ax1.set_ylim(0.0, 3.0)
ax1.legend(loc=2)

#=========================================================================
# numerical results (creep fatigue),(S 09-07),(H-L), (Normalized N)
#=========================================================================
ax3 = plt.subplot(244)


n_1, eps_1 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_H.txt')
n_2, eps_2 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_L.txt')
n_3, eps_3 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_HL_02.txt')
n_4, eps_4 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_HL_04.txt')
n_5, eps_5 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_HL_06.txt')
n_6, eps_6 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_HL_08.txt')


ax3.plot(n_1[1:] / n_1[-1], abs(eps_1[1:] * 1000), "--k", label='H')
#ax3.plot(n_2[1:] / n_2[-1], abs(eps_2[1:] * 1000), "k", label='L')
ax3.plot(n_3[1:] / n_3[-1], abs(eps_3[1:] * 1000), "r", label='HL_02')
ax3.plot(n_4[1:] / n_4[-1], abs(eps_4[1:] * 1000), "g", label='HL_04')
ax3.plot(n_5[1:] / n_5[-1], abs(eps_5[1:] * 1000), "b", label='HL_06')
ax3.plot(n_6[1:] / n_6[-1], abs(eps_6[1:] * 1000), "y", label='HL_08')


ax3.set_xlabel('N/Nf')
ax3.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
#ax1.set_xlim(0, 5)
#ax1.set_ylim(0.0, 3.0)
ax3.legend(loc=2)

#=========================================================================
# numerical results (creep fatigue),(S 09-07),(L-H)
#=========================================================================
ax1 = plt.subplot(245)
n_1, eps_1 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_H.txt')
n_2, eps_2 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_L.txt')
n_3, eps_3 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_LH_01.txt')
n_4, eps_4 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_LH_02.txt')
n_5, eps_5 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_LH_04.txt')
n_6, eps_6 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_LH_06.txt')

#ax1.plot(n_1[1:], abs(eps_1[1:] * 1000), "k", label='H')
ax1.plot(n_2[1:], abs(eps_2[1:] * 1000), "--k", label='L')
ax1.plot(n_3[1:], abs(eps_3[1:] * 1000), "r", label='LH_01')
ax1.plot(n_4[1:], abs(eps_4[1:] * 1000), "g", label='LH_02')
ax1.plot(n_5[1:], abs(eps_5[1:] * 1000), "b", label='LH_04')
ax1.plot(n_6[1:], abs(eps_6[1:] * 1000), "y", label='LH_06')


ax1.set_xlabel('N')
ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
#ax1.set_xlim(0, 5)
#ax1.set_ylim(0.0, 3.0)
ax1.legend(loc=2)

#=========================================================================
# numerical results (creep fatigue),(S 09-07),(L-H) ( accumulated damage)
#=========================================================================
ax1 = plt.subplot(246)
n_1, eps_1 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_H.txt')
n_2, eps_2 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_L.txt')
n_3, eps_3 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_LH_02.txt')
n_4, eps_4 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_LH_04.txt')
n_5, eps_5 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_LH_06.txt')
n_6, eps_6 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_LH_08.txt')

# applied number of cycles

n_02 = 3332
n_04 = 6665
n_06 = 9997
n_08 = 13330
n_hf = 95.0
n_lf = 16662.0

# normalized accumulative damage
N_02 = np.hstack(
    (n_3[0:n_02 + 1] / n_lf, 0.2 + (n_3[n_02 + 1:] - n_3[n_02]) / n_hf))
N_04 = np.hstack(
    (n_4[0:n_04 + 1] / n_lf, 0.4 + (n_4[n_04 + 1:] - n_4[n_04]) / n_hf))
N_06 = np.hstack(
    (n_5[0:n_06 + 1] / n_lf, 0.6 + (n_5[n_06 + 1:] - n_5[n_06]) / n_hf))
N_08 = np.hstack(
    (n_6[0:n_08 + 1] / n_lf, 0.8 + (n_6[n_08 + 1:] - n_6[n_08]) / n_hf))

#ax1.plot(n_1[1:] / n_1[-1], abs(eps_1[1:] * 1000), "--k", label='H')
ax1.plot(n_2[1:] / n_2[-1], abs(eps_2[1:] * 1000), "--k", label='L')
ax1.plot(N_02[1:], abs(eps_3[1:] * 1000), "r", label='HL_02')
ax1.plot(N_04[1:], abs(eps_4[1:] * 1000), "g", label='HL_04')
ax1.plot(N_06[1:], abs(eps_5[1:] * 1000), "b", label='HL_06')
ax1.plot(N_08[1:], abs(eps_6[1:] * 1000), "y", label='HL_08')

ax1.set_xlabel('N')
ax1.set_xlabel('cumulative damage ratio $\sum N_i/N^f_i$')

# #ax1.set_xlim(0, 5)
ax1.set_ylim(1.5, 3.2)
# ax1.legend(loc=2)
#

# #=========================================================================
# # numerical results (creep fatigue),(S 09-07),(L-H)
# #=========================================================================
ax1 = plt.subplot(247)
n_1, eps_1 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_H.txt')
n_2, eps_2 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_L.txt')
n_3, eps_3 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_LH_01.txt')
n_4, eps_4 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_LH_02.txt')
n_5, eps_5 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_LH_04.txt')
n_6, eps_6 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_LH_06.txt')

#ax1.plot(np.log10(n_1[1:]), abs(eps_1[1:] * 1000), "k", label='H')
ax1.plot(np.log10(n_2[1:]), abs(eps_2[1:] * 1000), "--k", label='L')
ax1.plot(np.log10(n_3[1:]), abs(eps_3[1:] * 1000), "r", label='LH_01')
ax1.plot(np.log10(n_4[1:]), abs(eps_4[1:] * 1000), "g", label='LH_02')
ax1.plot(np.log10(n_5[1:]), abs(eps_5[1:] * 1000), "b", label='LH_04')
ax1.plot(np.log10(n_6[1:]), abs(eps_6[1:] * 1000), "y", label='LH_06')


ax1.set_xlabel('log(N)')
ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
#ax1.set_xlim(0, 5)
#ax1.set_ylim(0.0, 3.0)
ax1.legend(loc=2)


# #=========================================================================
# # numerical results (creep fatigue),(S 09-07),(L-H), (normalized N)
# #=========================================================================
ax4 = plt.subplot(248)

n_1, eps_1 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_H.txt')
n_2, eps_2 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_L.txt')
n_3, eps_3 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_LH_01.txt')
n_4, eps_4 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_LH_02.txt')
n_5, eps_5 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_LH_04.txt')
n_6, eps_6 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\S_09-07\eps_LH_06.txt')

#ax4.plot(n_1[1:] / n_1[-1], abs(eps_1[1:] * 1000), "k", label='H')
ax4.plot(n_2[1:] / n_2[-1], abs(eps_2[1:] * 1000), "--k", label='L')
ax4.plot(n_3[1:] / n_3[-1], abs(eps_3[1:] * 1000), "r", label='LH_01')
ax4.plot(n_4[1:] / n_4[-1], abs(eps_4[1:] * 1000), "g", label='LH_02')
ax4.plot(n_5[1:] / n_5[-1], abs(eps_5[1:] * 1000), "b", label='LH_04')
ax4.plot(n_6[1:] / n_6[-1], abs(eps_6[1:] * 1000), "y", label='LH_06')


ax4.set_xlabel('N/Nf')
ax4.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
#ax1.set_xlim(0, 5)
#ax1.set_ylim(0.0, 3.0)
ax4.legend(loc=2)


plt.show()
