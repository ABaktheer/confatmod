'''
Created on 07.09.2017

@author: abaktheer
'''
import matplotlib.pyplot as plt
import numpy as np


# #=========================================================================================
# # plot damage comparison depend on the shape of the block (Linear-quadratic)(d=0.01, N=10)
# #=========================================================================================
# ax1 = plt.subplot(111)
#
# n_1, d_1 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\amplitude_shape_(Linear, Qudratic)\Smin_d001_In.txt')
# n_2, d_2 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\amplitude_shape_(Linear, Qudratic)\Smin_d001_De.txt')
# n_3, d_3 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\amplitude_shape_(Linear, Qudratic)\Smin_d001_InDe.txt')
#
# n_4, d_4 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\amplitude_shape_(Linear, Qudratic)\DS_d001_In.txt')
# n_5, d_5 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\amplitude_shape_(Linear, Qudratic)\DS_d001_De.txt')
# n_6, d_6 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\amplitude_shape_(Linear, Qudratic)\DS_d001_InDe.txt')
#
#
# ax1.plot(n_1, d_1, marker='o', markersize=5, color="green", label='In')
# ax1.plot(n_2, d_2, marker='o', markersize=5, color="red", label='De')
# ax1.plot(n_3, d_3,  marker='o', markersize=5, color="blue", label='In-De')
# ax1.plot(n_4, d_4, marker='o', markersize=5, color="green")
# ax1.plot(n_5, d_5, marker='o', markersize=5, color="red")
# ax1.plot(n_6, d_6,  marker='o', markersize=5, color="blue")
#
# ax1.set_ylim(0, 1.2)
# ax1.set_xlabel('$N_i/N_i^f$')
# ax1.set_ylabel('$\sum N_i/N_i^f$')
# ax1.legend()


#=========================================================================
# plot damage comparison depend on the shape of the block (Linear-quadratic)(d=0.05, N=50)
#=========================================================================
ax1 = plt.subplot(111)

n_1, d_1 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\amplitude_shape_(Linear, Qudratic)\Smin_d05_In.txt')
n_2, d_2 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\amplitude_shape_(Linear, Qudratic)\Smin_d05_De.txt')
n_3, d_3 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\amplitude_shape_(Linear, Qudratic)\Smin_d05_InDe.txt')

n_4, d_4 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\amplitude_shape_(Linear, Qudratic)\DS_d05_In.txt')
n_5, d_5 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\amplitude_shape_(Linear, Qudratic)\DS_d05_De.txt')
n_6, d_6 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_loading_blocks\amplitude_shape_(Linear, Qudratic)\DS_d05_InDe.txt')


ax1.plot(n_1, d_1, marker='o', markersize=5, color="green", label='In')
ax1.plot(n_2, d_2, marker='o', markersize=5, color="red", label='De')
ax1.plot(n_3, d_3,  marker='o', markersize=5, color="blue", label='In-De')
ax1.plot(n_4, d_4, marker='o', markersize=5, color="green")
ax1.plot(n_5, d_5, marker='o', markersize=5, color="red")
ax1.plot(n_6, d_6,  marker='o', markersize=5, color="blue")

ax1.set_ylim(0, 1.2)
ax1.set_xlabel('$N_i/N_i^f$')
ax1.set_ylabel('$\sum N_i/N_i^f$')
ax1.legend()


plt.show()
