'''
Created on 28.08.2017

@author: abaktheer
'''


'''
Alliche model- results - two level loading - Smin = const.
'''
import matplotlib.pyplot as plt
import numpy as np


#=========================================================================
# numerical results (comparison with P-M rule) (S_min)
#=========================================================================
ax1 = plt.subplot(111)

n_1, m_1 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\Woehler_H-L.txt')
n_2, m_2 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_0775-0725\Woehler_L-H.txt')
n_3, m_3 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\Woehler_H-L.txt')
n_4, m_4 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_08-07\Woehler_L-H.txt')
n_5, m_5 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\Woehler_H-L.txt')
n_6, m_6 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\Smin_085-065\Woehler_L-H.txt')

ax1.plot(n_1, m_1, 'b', label='H-L, $\Delta = 0.05$')
ax1.plot(n_2, m_2, 'b', label='L-H, $\Delta = 0.05$')
ax1.plot(n_3, m_3, 'r', label='H-L, $\Delta = 0.1$')
ax1.plot(n_4, m_4, 'r', label='L-H, $\Delta = 0.1$')
ax1.plot(n_5, m_5, 'g', label='H-L, $\Delta = 0.2$')
ax1.plot(n_6, m_6, 'g', label='L-H, $\Delta = 0.2$')


# Palmgren-Miner rule
x_2 = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
y_2 = np.array([1.0, 0.8, 0.6, 0.4, 0.2, 0.0])

ax1.plot(x_2, y_2, "-k", label='Palmgren-Miner rule')

# ax1.fill_between(x_2, y_2, m_1, facecolor='k', alpha=0.1)
# ax1.fill_between(x_2, y_2, m_2, facecolor='k', alpha=0.1)
# ax1.fill_between(x_2, y_2, m_3, facecolor='k', alpha=0.1)
# ax1.fill_between(x_2, y_2, m_4, facecolor='k', alpha=0.1)
# ax1.fill_between(x_2, y_2, m_5, facecolor='k', alpha=0.1)
# ax1.fill_between(x_2, y_2, m_6, facecolor='k', alpha=0.1)

ax1.set_xlabel('$N_1/N_f$')
ax1.set_ylabel('$N_2/N_f$')
ax1.set_xlim(0, 1.2)
ax1.set_ylim(0, 1.2)
ax1.legend(loc=1)

#=========================================================================
# numerical results (comparison with P-M rule) (DS)
#=========================================================================
ax1 = plt.subplot(111)

n_11, m_11 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\DS_0775-0725\Woehler_H-L.txt')
n_22, m_22 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\DS_0775-0725\Woehler_L-H.txt')
n_33, m_33 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\DS_08-07\Woehler_H-L.txt')
n_44, m_44 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\DS_08-07\Woehler_L-H.txt')
n_55, m_55 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\DS_085-065\Woehler_H-L.txt')
n_66, m_66 = np.loadtxt(
    r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\LS4_two_levels\DS_085-065\Woehler_L-H.txt')

ax1.plot(n_11, m_11, '--b', label='H-L, $\Delta = 0.05$')
ax1.plot(n_22, m_22, '--b', label='L-H, $\Delta = 0.05$')
ax1.plot(n_33, m_33, '--r', label='H-L, $\Delta = 0.1$')
ax1.plot(n_44, m_44, '--r', label='L-H, $\Delta = 0.1$')
ax1.plot(n_55, m_55, '--g', label='H-L, $\Delta = 0.2$')
ax1.plot(n_66, m_66, '--g', label='L-H, $\Delta = 0.2$')


# Palmgren-Miner rule
x_2 = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
y_2 = np.array([1.0, 0.8, 0.6, 0.4, 0.2, 0.0])

ax1.plot(x_2, y_2, "-k", label='Palmgren-Miner rule')

# ax1.fill_between(x_2, y_2, m_1, facecolor='k', alpha=0.1)
# ax1.fill_between(x_2, y_2, m_2, facecolor='k', alpha=0.1)
# ax1.fill_between(x_2, y_2, m_3, facecolor='k', alpha=0.1)
# ax1.fill_between(x_2, y_2, m_4, facecolor='k', alpha=0.1)
# ax1.fill_between(x_2, y_2, m_5, facecolor='k', alpha=0.1)
# ax1.fill_between(x_2, y_2, m_6, facecolor='k', alpha=0.1)

ax1.set_xlabel('$N_1/N_f$')
ax1.set_ylabel('$N_2/N_f$')
ax1.set_xlim(0, 1.2)
ax1.set_ylim(0, 1.2)
ax1.legend(loc=1)


plt.show()
