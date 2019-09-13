'''
Created on 28.08.2017

@author: abaktheer
'''


'''
Created on 28.08.2017

@author: abaktheer
'''
import matplotlib.pyplot as plt
import numpy as np

# #=========================================================================
# # numerical results (creep fatigue),(eps_1),(s06-S07-S08-S09)
# #=========================================================================
ax1 = plt.subplot(221)


n_1, eps_1 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\constant_loading\creep_fatigue\S06\eps_1_S06.txt')
n_2, eps_2 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\constant_loading\creep_fatigue\S07\eps_1_S07.txt')
n_3, eps_3 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\constant_loading\creep_fatigue\S08\eps_1_S08.txt')
n_4, eps_4 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\constant_loading\creep_fatigue\S09\eps_1_S09.txt')


ax1.plot(n_1[1:] / n_1[-1], eps_1[1:] * 1000, "--k", label='S=0.6')
ax1.plot(n_2[1:] / n_2[-1], eps_2[1:] * 1000, "--r", label='S=0.7')
ax1.plot(n_3[1:] / n_3[-1], eps_3[1:] * 1000, "--g", label='S=0.8')
ax1.plot(n_4[1:] / n_4[-1], eps_4[1:] * 1000, "--b", label='S=0.9')


ax1.set_xlabel('N/Nf')
ax1.set_ylabel('max strain $\epsilon_1 \;.10^{-3}$')
#ax1.set_xlim(0, 5)
ax1.set_ylim(-4.5, -1.0)
ax1.legend(loc=3)
#
#
# #=========================================================================
# # numerical results (creep fatigue)(eps_2),(s06-S07-S08-S09)
# #=========================================================================
ax2 = plt.subplot(222)

n_1, eps_1 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\constant_loading\creep_fatigue\S06\eps_2_S06.txt')
n_2, eps_2 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\constant_loading\creep_fatigue\S07\eps_2_S07.txt')
n_3, eps_3 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\constant_loading\creep_fatigue\S08\eps_2_S08.txt')
n_4, eps_4 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\constant_loading\creep_fatigue\S09\eps_2_S09.txt')


ax2.plot(n_1[1:] / n_1[-1], eps_1[1:] * 1000, "--k", label='S=0.6')
ax2.plot(n_2[1:] / n_2[-1], eps_2[1:] * 1000, "--r", label='S=0.7')
ax2.plot(n_3[1:] / n_3[-1], eps_3[1:] * 1000, "--g", label='S=0.8')
ax2.plot(n_4[1:] / n_4[-1], eps_4[1:] * 1000, "--b", label='S=0.9')


ax2.set_xlabel('N/Nf')
ax2.set_ylabel('max strain $\epsilon_2 \;.10^{-3}$')

ax2.legend(loc=2)

# #=========================================================================
# # numerical results (creep fatigue)(damage),(s06-S07-S08-S09)
# #=========================================================================
ax3 = plt.subplot(223)

n_1, w_1 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\constant_loading\creep_fatigue\S06\damage_S06.txt')
n_2, w_2 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\constant_loading\creep_fatigue\S07\damage_S07.txt')
n_3, w_3 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\constant_loading\creep_fatigue\S08\damage_S08.txt')
n_4, w_4 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\constant_loading\creep_fatigue\S09\damage_S09.txt')


ax3.plot(n_1[1:] / n_1[-1], w_1[1:] / w_1[-1], "--k", label='S=0.6')
ax3.plot(n_2[1:] / n_2[-1], w_2[1:] / w_2[-1], "--r", label='S=0.7')
ax3.plot(n_3[1:] / n_3[-1], w_3[1:] / w_3[-1], "--g", label='S=0.8')
ax3.plot(n_4[1:] / n_4[-1], w_4[1:] / w_4[-1], "--b", label='S=0.9')

ax3.set_xlabel('N/Nf')
ax3.set_ylabel('normalized damage')

ax3.legend(loc=2)


# #=========================================================================
# # numerical results (creep fatigue)(damage),(s06-S07-S08-S09)
# #=========================================================================
ax4 = plt.subplot(224)

n_1, K_1 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\constant_loading\creep_fatigue\S06\K_S06.txt')
n_2, K_2 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\constant_loading\creep_fatigue\S07\K_S07.txt')
n_3, K_3 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\constant_loading\creep_fatigue\S08\K_S08.txt')
n_4, K_4 = np.loadtxt(
    'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\constant_loading\creep_fatigue\S09\K_S09.txt')


ax4.plot(n_1[1:] / n_1[-1], K_1[1:] / K_1[1], "--k", label='S=0.6')
ax4.plot(n_2[1:] / n_2[-1], K_2[1:] / K_2[1], "--r", label='S=0.7')
ax4.plot(n_3[1:] / n_3[-1], K_3[1:] / K_3[1], "--g", label='S=0.8')
ax4.plot(n_4[1:] / n_4[-1], K_4[1:] / K_4[1], "--b", label='S=0.9')

ax4.set_xlabel('N/Nf')
ax4.set_ylabel('stiffness[MPa]')

ax4.legend(loc=3)


# # #=========================================================================
# # # numerical results (wohler)
# # #=========================================================================
# ax4 = plt.subplot(111)
#
# n_1, s_1 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\constant_loading\woehler\exp.txt')
# n_2, s_2 = np.loadtxt(
#     r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\constant_loading\woehler\num.txt')
#
#
# ax4.plot(n_1[:], s_1[:], "--k", label='experiment')
# ax4.plot(n_2[:], s_2[:], "k", label='model')
#
#
# ax4.set_xlabel('log(N)')
# ax4.set_ylabel('$S = F_{max}/F_{u}$')
# ax4.set_ylim(0.5, 1.0)
# ax4.legend(loc=1)


plt.show()
