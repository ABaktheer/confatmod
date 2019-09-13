'''
Created on 07.09.2017

@author: abaktheer
'''
import matplotlib.pyplot as plt
import numpy as np


#============================================================
# plot damage 5_Smin_d , In _ De _ InDe , Alliche vs. Miner
#============================================================
ax1 = plt.subplot(221)

n_0, d_0 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\damage_comparison\d_Miner.txt')
n_1, d_1 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\damage_comparison\5_Smin_d_In.txt')
n_2, d_2 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\damage_comparison\5_Smin_d_De.txt')
n_3, d_3 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\damage_comparison\5_Smin_d_InDe.txt')


ax1.plot(n_0, d_0, "--k", lw=2, label='P-M rule')
ax1.plot(n_1, d_1, marker='o', markersize=5, color="green", label='In')
ax1.plot(n_2, d_2, marker='o', markersize=5, color="red", label='De')
ax1.plot(n_3, d_3,  marker='o', markersize=5, color="blue", label='In-De')

ax1.set_ylim(0, 1.2)
ax1.set_xlabel('$N_i/N_i^f$')
ax1.set_ylabel('$\sum N_i/N_i^f$')
ax1.legend()


#============================================================
# plot damage 5_Smin_N , In _ De _ InDe , Alliche vs. Miner
#============================================================
ax1 = plt.subplot(222)

n_0, d_0 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\damage_comparison\N_Miner.txt')
n_1, d_1 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\damage_comparison\5_Smin_N_In.txt')
n_2, d_2 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\damage_comparison\5_Smin_N_De.txt')
n_3, d_3 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\damage_comparison\5_Smin_N_InDe.txt')

ax1.plot(n_0, d_0, "--k", lw=2, label='P-M rule')
ax1.plot(n_1, d_1, marker='o', markersize=5, color="green", label='In')
ax1.plot(n_2, d_2, marker='o', markersize=5, color="red", label='De')
ax1.plot(n_3, d_3,  marker='o', markersize=5, color="blue", label='In-De')

ax1.set_ylim(0, 1.2)
ax1.set_xlabel('$N_i$')
ax1.set_ylabel('$\sum N_i/N_i^f$')
ax1.legend()


#============================================================
# plot damage 5_DS_d , In _ De _ InDe , Alliche vs. Miner
#============================================================
ax1 = plt.subplot(223)

n_0, d_0 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\damage_comparison\d_Miner.txt')
n_1, d_1 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\damage_comparison\5_DS_d_In.txt')
n_2, d_2 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\damage_comparison\5_DS_d_De.txt')
n_3, d_3 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\damage_comparison\5_DS_d_InDe.txt')

ax1.plot(n_0, d_0, "--k", lw=2, label='P-M rule')
ax1.plot(n_1, d_1, marker='o', markersize=5, color="green", label='In')
ax1.plot(n_2, d_2, marker='o', markersize=5, color="red", label='De')
ax1.plot(n_3, d_3,  marker='o', markersize=5, color="blue", label='In-De')

ax1.set_ylim(0, 1.2)
ax1.set_xlabel('$N_i/N_i^f$')
ax1.set_ylabel('$\sum N_i/N_i^f$')
ax1.legend()


#============================================================
# plot damage 5_DS_N , In _ De _ InDe , Alliche vs. Miner
#============================================================
ax1 = plt.subplot(224)

n_0, d_0 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\damage_comparison\N_Miner.txt')
n_1, d_1 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\damage_comparison\5_DS_N_In.txt')
n_2, d_2 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\damage_comparison\5_DS_N_De.txt')
n_3, d_3 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\damage_comparison\5_DS_N_InDe.txt')

ax1.plot(n_0, d_0, "--k", lw=2, label='P-M rule')
ax1.plot(n_1, d_1, marker='o', markersize=5, color="green", label='In')
ax1.plot(n_2, d_2, marker='o', markersize=5, color="red", label='De')
ax1.plot(n_3, d_3,  marker='o', markersize=5, color="blue", label='In-De')

ax1.set_ylim(0, 1.2)
ax1.set_xlabel('$N_i$')
ax1.set_ylabel('$\sum N_i/N_i^f$')
ax1.legend()


plt.show()
