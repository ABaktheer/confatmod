'''
Created on 07.09.2017

@author: abaktheer
'''
import matplotlib.pyplot as plt
import numpy as np


#============================================================
# plot damage d=0.01 , In _ De _ InDe , Alliche vs. Miner
#============================================================
ax1 = plt.subplot(221)


n_1, d_1 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\001\In_d001_damage_Miner.txt')
n_2, d_2 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\001\De_d001_damage_Miner.txt')
n_3, d_3 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\001\InDe_d001_damage_Miner.txt')
n_4, d_4 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\001\In_d001_damage.txt')
n_5, d_5 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\001\De_d001_damage.txt')
n_6, d_6 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\001\InDe_d001_damage.txt')


ax1.plot(n_1, d_1, "--k", lw=2, label='In_Miner')
ax1.plot(n_2, d_2, ":k", lw=2, label='De_Miner')
ax1.plot(n_3, d_3, "k", lw=0.5, label='InDe_Miner')
ax1.plot(n_4[0:], d_4[0:] / d_4[-1], "--r", label='In_Alliche')
ax1.plot(n_5[0:], d_5[0:] / d_5[-1], ":r", label='De_Alliche')
ax1.plot(n_6[0:], d_6[0:] / d_6[-1], "r", label='InDe_Alliche')


ax1.set_xlabel('N')
ax1.set_ylabel('damage')
ax1.legend()


#============================================================
# plot creep_fatigue d=0.01 , In _ De _ InDe , Alliche vs. Miner
#============================================================
ax1 = plt.subplot(223)

n_4, eps_4 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\001\In_d001_creep_fatigue.txt')
n_5, eps_5 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\001\De_d001_creep_fatigue.txt')
n_6, eps_6 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\001\InDe_d001_creep_fatigue.txt')

ax1.plot(n_4[0:], abs(eps_4[0:] * 1000), "k", label='In')
ax1.plot(n_5[0:], abs(eps_5[0:] * 1000), "r", label='De')
ax1.plot(n_6[0:], abs(eps_6[0:] * 1000), "g", label='InDe')

ax1.set_xlabel('N')
ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
ax1.legend()


# #============================================================
# # plot damage d=0.05 , In _ De _ InDe , Alliche vs. Miner
# #============================================================
# ax1 = plt.subplot(232)
#
# n_1, d_1 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\005\In_d005_damage_Miner.txt')
# n_2, d_2 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\005\De_d005_damage_Miner.txt')
# n_3, d_3 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\005\InDe_d005_damage_Miner.txt')
# n_4, d_4 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\005\In_d005_damage.txt')
# n_5, d_5 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\005\De_d005_damage.txt')
# n_6, d_6 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\005\InDe_d005_damage.txt')
#
#
# ax1.plot(n_1, d_1, "--k", lw=2, label='In_Miner')
# ax1.plot(n_2, d_2, ":k", lw=2, label='De_Miner')
# ax1.plot(n_3, d_3, "k", lw=0.5, label='InDe_Miner')
# ax1.plot(n_4[0:], d_4[0:] / d_4[-1], "--r", label='In_Alliche')
# ax1.plot(n_5[0:], d_5[0:] / d_5[-1], ":r", label='De_Alliche')
# ax1.plot(n_6[0:], d_6[0:] / d_6[-1], "r", label='InDe_Alliche')
#
#
# ax1.set_xlabel('N')
# ax1.set_ylabel('damage')
# ax1.legend()


# #===================================================================
# # plot creep_fatigue d=0.1 , In _ De _ InDe , Alliche vs. Miner
# #===================================================================
# ax1 = plt.subplot(235)
#
# n_4, eps_4 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\005\In_d005_creep_fatigue.txt')
# n_5, eps_5 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\005\De_d005_creep_fatigue.txt')
# n_6, eps_6 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\005\InDe_d005_creep_fatigue.txt')
#
# ax1.plot(n_4[0:], abs(eps_4[0:] * 1000), "k", label='In')
# ax1.plot(n_5[0:], abs(eps_5[0:] * 1000), "r", label='De')
# ax1.plot(n_6[0:], abs(eps_6[0:] * 1000), "g", label='InDe')
#
# ax1.set_xlabel('N')
# ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
# ax1.legend()


# #============================================================
# # plot damage d=0.1 , In _ De _ InDe , Alliche vs. Miner
# #============================================================
# ax1 = plt.subplot(222)
#
# n_1, d_1 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\01\In_d01_damage_Miner.txt')
# n_2, d_2 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\01\De_d01_damage_Miner.txt')
# n_3, d_3 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\01\InDe_d01_damage_Miner.txt')
# n_4, d_4 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\01\In_d01_damage.txt')
# n_5, d_5 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\01\De_d01_damage.txt')
# n_6, d_6 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\01\InDe_d01_damage.txt')
#
#
# ax1.plot(n_1, d_1, "--k", lw=2, label='In_Miner')
# ax1.plot(n_2, d_2, ":k", lw=2, label='De_Miner')
# ax1.plot(n_3, d_3, "k", lw=0.5, label='InDe_Miner')
# ax1.plot(n_4[0:], d_4[0:] / d_4[-1], "--r", label='In_Alliche')
# ax1.plot(n_5[0:], d_5[0:] / d_5[-1], ":r", label='De_Alliche')
# ax1.plot(n_6[0:], d_6[0:] / d_6[-1], "r", label='InDe_Alliche')
#
#
# ax1.set_xlabel('N')
# ax1.set_ylabel('damage')
# ax1.legend()
#
#
# #===================================================================
# # plot creep_fatigue d=0.1 , In _ De _ InDe , Alliche vs. Miner
# #===================================================================
# ax1 = plt.subplot(224)
#
# n_4, eps_4 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\01\In_d01_creep_fatigue.txt')
# n_5, eps_5 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\01\De_d01_creep_fatigue.txt')
# n_6, eps_6 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_damage\01\InDe_d01_creep_fatigue.txt')
#
# ax1.plot(n_4[0:], abs(eps_4[0:] * 1000), "k", label='In')
# ax1.plot(n_5[0:], abs(eps_5[0:] * 1000), "r", label='De')
# ax1.plot(n_6[0:], abs(eps_6[0:] * 1000), "g", label='InDe')
#
# ax1.set_xlabel('N')
# ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
# ax1.legend()

#=========================================================================

#============================================================
# plot damage N=10 , In _ De _ InDe , Alliche vs. Miner
#============================================================
ax1 = plt.subplot(222)


# n_1, d_1 = np.loadtxt(
# r'E:\Models_Implementation\Concrete Fatigue
# models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\10\In_N10_damage_Miner22222.txt')
n_2, d_2 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\10\De_N10_damage_Miner.txt')
n_3, d_3 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\10\InDe_N10_damage_Miner.txt')
n_4, d_4 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\10\In_N10_damage.txt')
n_5, d_5 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\10\De_N10_damage.txt')
n_6, d_6 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\10\InDe_N10_damage.txt')


#ax1.plot(n_1, d_1, "--k", lw=1, label='In_Miner')
ax1.plot(n_2, d_2, ":k", lw=1, label='De')
ax1.plot(n_3, d_3, "k", lw=1, label='InDe')
ax1.plot(n_4[0:], d_4[0:] / d_4[-1], "--r", label='In')
ax1.plot(n_5[0:], d_5[0:] / d_5[-1], ":r", label='De')
ax1.plot(n_6[0:], d_6[0:] / d_6[-1], "r", label='InDe')

ax1.set_ylim(0, 1)
ax1.set_xlabel('N')
ax1.set_ylabel('damage')
ax1.legend()


#============================================================
# plot creep_fatigue N=10 , In _ De _ InDe , Alliche vs. Miner
#============================================================
ax1 = plt.subplot(224)

n_4, eps_4 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\10\In_N10_creep_fatigue.txt')
n_5, eps_5 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\10\De_N10_creep_fatigue.txt')
n_6, eps_6 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\10\InDe_N10_creep_fatigue.txt')

ax1.plot(n_4[0:], abs(eps_4[0:] * 1000), "k", label='In')
ax1.plot(n_5[0:], abs(eps_5[0:] * 1000), "r", label='De')
ax1.plot(n_6[0:], abs(eps_6[0:] * 1000), "g", label='InDe')

ax1.set_xlabel('N')
ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
ax1.legend()

#
#
# #============================================================
# # plot damage N=50 , In _ De _ InDe , Alliche vs. Miner
# #============================================================
# ax1 = plt.subplot(232)
#
# n_1, d_1 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\50\In_N50_damage_Miner.txt')
# n_2, d_2 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\50\De_N50_damage_Miner.txt')
# n_3, d_3 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\50\InDe_N50_damage_Miner.txt')
# n_4, d_4 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\50\In_N50_damage.txt')
# n_5, d_5 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\50\De_N50_damage.txt')
# n_6, d_6 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\50\InDe_N50_damage.txt')
#
#
# ax1.plot(n_1, d_1, "--k", lw=1, label='In')
# ax1.plot(n_2, d_2, ":k", lw=1, label='De')
# ax1.plot(n_3, d_3, "k", lw=1, label='InDe')
# ax1.plot(n_4[0:], d_4[0:] / d_4[-1], "--r", label='In')
# ax1.plot(n_5[0:], d_5[0:] / d_5[-1], ":r", label='De')
# ax1.plot(n_6[0:], d_6[0:] / d_6[-1], "r", label='InDe')
#
# ax1.set_ylim(0, 1)
# ax1.set_xlabel('N')
# ax1.set_ylabel('damage')
# ax1.legend()
#
#
# #===================================================================
# # plot creep_fatigue N=50 , In _ De _ InDe , Alliche vs. Miner
# #===================================================================
# ax1 = plt.subplot(235)
#
# n_4, eps_4 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\50\In_N50_creep_fatigue.txt')
# n_5, eps_5 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\50\De_N50_creep_fatigue.txt')
# n_6, eps_6 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\50\InDe_N50_creep_fatigue.txt')
#
# ax1.plot(n_4[0:], abs(eps_4[0:] * 1000), "k", label='In')
# ax1.plot(n_5[0:], abs(eps_5[0:] * 1000), "r", label='De')
# ax1.plot(n_6[0:], abs(eps_6[0:] * 1000), "g", label='InDe')
#
# ax1.set_xlabel('N')
# ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
# ax1.legend()
#
#
# #============================================================
# # plot damage N=100 , In _ De _ InDe , Alliche vs. Miner
# #============================================================
# ax1 = plt.subplot(233)
#
# n_1, d_1 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\100\In_N100_damage_Miner.txt')
# n_2, d_2 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\100\De_N100_damage_Miner.txt')
# n_3, d_3 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\100\InDe_N100_damage_Miner.txt')
# n_4, d_4 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\100\In_N100_damage.txt')
# n_5, d_5 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\100\De_N100_damage.txt')
# n_6, d_6 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\100\InDe_N100_damage.txt')
#
#
# ax1.plot(n_1, d_1, "--k", lw=1, label='In')
# ax1.plot(n_2, d_2, ":k", lw=1, label='De')
# ax1.plot(n_3, d_3, "k", lw=1, label='InDe')
# ax1.plot(n_4[0:], d_4[0:] / d_4[-1], "--r", label='In')
# ax1.plot(n_5[0:], d_5[0:] / d_5[-1], ":r", label='De')
# ax1.plot(n_6[0:], d_6[0:] / d_6[-1], "r", label='InDe')
#
# ax1.set_ylim(0, 1)
# ax1.set_xlabel('N')
# ax1.set_ylabel('damage')
# ax1.legend()
#
#
# #===================================================================
# # plot creep_fatigue N=100 , In _ De _ InDe , Alliche vs. Miner
# #===================================================================
# ax1 = plt.subplot(236)
#
# n_4, eps_4 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\100\In_N100_creep_fatigue.txt')
# n_5, eps_5 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\100\De_N100_creep_fatigue.txt')
# n_6, eps_6 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\5_levels_deltaS_N\100\InDe_N100_creep_fatigue.txt')
#
# ax1.plot(n_4[0:], abs(eps_4[0:] * 1000), "k", label='In')
# ax1.plot(n_5[0:], abs(eps_5[0:] * 1000), "r", label='De')
# ax1.plot(n_6[0:], abs(eps_6[0:] * 1000), "g", label='InDe')
#
#
# ax1.set_xlabel('N')
# ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
# ax1.legend()


plt.show()
