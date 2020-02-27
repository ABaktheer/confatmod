'''
Created on 14.06.2017

@author: ABaktheer
'''

'''
Implementation of the fatigue model for plain concrete [A.Alliche, 2004]under uniaxialpha compressive loading
(stress driven alphagorithm) - implicit code
'''

from scipy.optimize import newton
import matplotlib.pyplot as plt
import numpy as np


def get_stress_strain(sigma_1_arr, lamda, mu, alpha, beta, g, C0, C1, K, n):

    #-----------------------------------------------------------------------
    # arrays to store the values
    #-----------------------------------------------------------------------

    # normal strain
    eps_1_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)
    # lateral strain
    eps_2_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)
    # damage factor
    w_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)

    #-----------------------------------------------------------------------
    # material parameters
    #-----------------------------------------------------------------------

    # lame constants [MPa]
    lamda = lamda
    mu = mu
    # fatigue model material parameter
    alpha = alpha
    beta = beta
    g = g
    C0 = C0
    C1 = C1
    K = K
    n = n

    #-----------------------------------------------------------------------
    # state variables
    #-----------------------------------------------------------------------
    eps_1_i = 0.0
    eps_2_i = 0.0
    w_n = 0.0
    w_i = 0.0

    for i in range(1, len(sigma_1_arr)):
        print ('increment', i)

        sigma_1_i = sigma_1_arr[i]
        d_sigma_1 = sigma_1_arr[i] - sigma_1_arr[i - 1]

        if d_sigma_1 < 0:
            f_w_n = lambda w_n:  w_i - w_n - abs(g) / (2.0 * C1) * (((abs(g) * (-1.0 * ((lamda + alpha * w_i) * sigma_1_i + g * w_i * (lamda + 2.0 * mu)) /
                                                                                ((lamda + 2.0 * mu) * (2.0 * (lamda + mu) + 4.0 *
                                                                                                       w_i * (alpha + beta)) - 2.0 * (lamda + alpha * w_i) ** 2)) -
                                                                      (C0 + 2 * C1 * w_i)) / K)**n) * d_sigma_1 * ((lamda + alpha * w_i) /
                                                                                                                   ((lamda + 2.0 * mu) * (2.0 * (lamda + mu) + 4.0 * w_i * (alpha + beta) -
                                                                                                                                          alpha * (g / (2.0 * C1)) * (2.0 * (-1.0 * ((lamda + alpha * w_i) * sigma_1_i + g * w_i * (lamda + 2.0 * mu)) /
                                                                                                                                                                             ((lamda + 2.0 * mu) * (2.0 * (lamda + mu) + 4.0 * w_i * (alpha + beta)) - 2.0 * (lamda + alpha * w_i) ** 2)) + (sigma_1_i /
                                                                                                                                                                                                                                                                                             (lamda + 2.0 * mu) - 2.0 * (-1.0 * ((lamda + alpha * w_i) * sigma_1_i + g * w_i * (lamda + 2.0 * mu)) /
                                                                                                                                                                                                                                                                                                                         ((lamda + 2.0 * mu) * (2.0 * (lamda + mu) + 4.0 * w_i * (alpha + beta)) - 2.0 * (lamda + alpha * w_i) ** 2)) *
                                                                                                                                                                                                                                                                                             (lamda + alpha * w_i) / (lamda + 2.0 * mu))) - (g**2.0 / (2.0 * C1))) - 2.0 * (lamda + alpha * w_i)**2))
            w_i = newton(f_w_n, 0.0, tol=1e-8, maxiter=10)

            w_n = w_i
        else:
            w_n = w_n

        eps_2_i = -1.0 * ((lamda + alpha * w_n) * sigma_1_i + g * w_n * (lamda + 2.0 * mu)) / \
            ((lamda + 2.0 * mu) * (2.0 * (lamda + mu) + 4.0 *
                                   w_n * (alpha + beta)) - 2.0 * (lamda + alpha * w_n) ** 2)
        eps_1_i = sigma_1_i / \
            (lamda + 2.0 * mu) - 2.0 * eps_2_i  *  \
            (lamda + alpha * w_n) / (lamda + 2.0 * mu)

        eps_1_arr[i] = eps_1_i
        eps_2_arr[i] = eps_2_i
        w_arr[i] = w_n

        if w_i > 5.0:
            print( ' ----------> No Convergence any more')
            break

        if abs(eps_1_i) > 0.0035:
            print (' ----------> No Convergence any more')
            break

    return sigma_1_arr, eps_1_arr, eps_2_arr, w_arr, i


if __name__ == '__main__':

    maximum_stress = -80.6  # [MPa]
    number_of_cycles = 5000  # 2000  # 4410  # 4410  # 4650
    #unloading_ratio = 0.333
    minimum_stress = -21
    m = 50  # number of increments in each cycle
    n = number_of_cycles

    s_levels = np.linspace(0, maximum_stress, n * 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] = minimum_stress
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 1] = maximum_stress
    s_history = s_levels.flatten()

    sig_1_arr = np.hstack([np.linspace(s_history[i], s_history[i + 1], m, dtype=np.float_)
                           for i in range(len(s_levels) - 1)])

    sigma_1_arr, eps_1_arr, eps_2_arr, w_arr, inc = get_stress_strain(
        sig_1_arr, lamda=10555.5, mu=15833.33, alpha=2237.5, beta=-2216.5, g=-9.788, C0=0.0, C1=0.002033, K=0.003345, n=10)

    t_arr = np.linspace(0, 1, len(sigma_1_arr))
    #-----------------------------------------------------------------------
    # plot 1
    #-----------------------------------------------------------------------
    ax1 = plt.subplot(231)
    ax1.plot(-eps_1_arr[0:inc], - sigma_1_arr[0:inc],
             'k', linewidth=1, alpha=1.0)
    ax1.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax1.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.xlabel('$\epsilon_{11}$')
    plt.ylabel('$\sigma_{11}$[MPa]')

    #-----------------------------------------------------------------------
    # plot 2
    #-----------------------------------------------------------------------
    ax2 = plt.subplot(232)
    ax2.plot(t_arr[0:inc], w_arr[0:inc] / w_arr[inc],
             'k', linewidth=1, alpha=1.0)
    ax2.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax2.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('Damage evolution')
    plt.xlabel('Time')
    plt.ylabel('Damage')

    #-----------------------------------------------------------------------
    # plot 3
    #-----------------------------------------------------------------------
    plt.subplot(233)
    plt.plot(-sigma_1_arr[0:inc], w_arr[0:inc] /
             w_arr[inc], 'k', linewidth=1, alpha=1)
    plt.xlabel('$\sigma_{11}$')
    plt.ylabel('Damage')

    #-----------------------------------------------------------------
    # plot 4 (creep_fatigue)
    #-----------------------------------------------------------------
    plt.subplot(234)

    eps_1 = np.zeros(n)
    cycle = np.zeros(n)
    for i in range(0, n, 1):
        idx = m + 2 * i * m - 1
        if idx <= len(eps_1_arr[0:inc]):
            idx = idx
        else:
            idx = m + 2 * (i - 1.0) * m - 1
            break

        eps_1[i] = eps_1_arr[idx]
        cycle[i] = i + 1

    plt.plot(cycle[0:i], eps_1[0:i] * -1000.0, 'k', linewidth=1,
             alpha=1, label='Alliche model')

    N_exp, e1_exp = np.loadtxt(
        r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\Exp_creepfatigue_S075.txt')

    plt.plot(N_exp, e1_exp * -1000.0, "--k", label='experiment')

    plt.xlabel('number of cycles')
    plt.ylabel('max $\epsilon_{11} \;.10^{-3}$')
    plt.legend()

    #-----------------------------------------------------------------
    # plot 5
    #-----------------------------------------------------------------
    plt.subplot(235)

    eps_2 = np.zeros(n)
    cycle_2 = np.zeros(n)
    for i in range(0, n, 1):

        idx = m + 2 * i * m - 1
        if idx <= len(eps_1_arr[0:inc]):
            idx = idx
        else:
            idx = m + 2 * (i - 1.0) * m - 1
            break

        eps_2[i] = eps_2_arr[idx]
        cycle_2[i] = i + 1


#         idx = m + 2 * i * m - 1
#         if idx < len(eps_2_arr[0:inc]):
#             eps_2[i] = eps_2_arr[idx]
#             cycle[i] = i + 1
#         else:
#             break
    plt.plot(cycle_2[0:i], eps_2[0:i], 'k', linewidth=1, alpha=1)
    plt.plot(cycle[0:i], eps_1[0:i], 'k', linewidth=1,
             alpha=1, label='Alliche model')
    plt.xlabel('number of cycles')
    plt.ylabel('max $\epsilon_{22}$')

    #-----------------------------------------------------------------
    # plot 6
    #-----------------------------------------------------------------
    plt.subplot(236)
    w = np.zeros(n)
    cycle_3 = np.zeros(n)
    for i in range(0, n, 1):

        idx = m + 2 * i * m - 1
        if idx <= len(eps_1_arr[0:inc]):
            idx = idx
        else:
            idx = m + 2 * (i - 1.0) * m - 1
            break

        w[i] = w_arr[idx]
        cycle_3[i] = i + 1
#
#     plt.plot(cycle_3[0:i], w[0:i], 'k', linewidth=1, alpha=1)
#     plt.xlabel('number of cycles')
#     plt.ylabel('Damage')
# #
#
#     #-----------------------------------------------------------------
#     # plot 4 (Stiffness)
#     #-----------------------------------------------------------------
#     plt.subplot(236)
#
#     eps_1 = np.zeros(n)
#     cycle = np.zeros(n)
#     for i in range(0, n, 1):
#         idx = m + 2 * i * m - 1
#         if idx <= len(eps_1_arr[0:inc]):
#             idx = idx
#         else:
#             idx = m + 2 * (i - 1.0) * m - 1
#             break
#
#         eps_1[i] = eps_1_arr[idx]
#         cycle[i] = i + 1
#
#     plt.plot(cycle[0:i], maximum_stress / eps_1[0:i],
#              'k', linewidth=1, alpha=1)
#
# #     N_exp, e1_exp = np.loadtxt(
# #         r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\Exp_creepfatigue_S075.txt')
# #
# #     plt.plot(N_exp, e1_exp, "--k", label='Experiment')
#
#     plt.xlabel('number of cycles')
#     plt.ylabel('Stiffness[MPa]')
#
#     #=========================================================================
#     # saving results
#     #=========================================================================
#     eps_record_1 = np.zeros(1)
#     eps_record_2 = np.zeros(1)
#     N_record = np.zeros(1)
#     w_record = np.zeros(1)
#     stiffness_record = np.zeros(1)
#
#     for i in range(0, i, 1):
#         eps_record_1 = np.vstack((eps_record_1, eps_1[i]))
#         eps_record_2 = np.vstack((eps_record_2, eps_2[i]))
#         N_record = np.vstack((N_record, cycle[i]))
#         w_record = np.vstack((w_record, w[i]))
#         stiffness_record = np.vstack(
#             (stiffness_record, maximum_stress / eps_1[i]))
#
#     np.savetxt(r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\constant_loading\creep_fatigue\eps_1.txt',
#                np.transpose(eps_record_1), delimiter=" ", fmt="%s")
#     np.savetxt(r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\constant_loading\creep_fatigue\eps_2.txt',
#                np.transpose(eps_record_2), delimiter=" ", fmt="%s")
#     np.savetxt(r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\constant_loading\creep_fatigue\N.txt',
#                np.transpose(N_record), delimiter=" ", fmt="%s")
#     np.savetxt(r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\constant_loading\creep_fatigue\w.txt',
#                np.transpose(w_record), delimiter=" ", fmt="%s")
#     np.savetxt(r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\constant_loading\creep_fatigue\stiffness.txt',
#                np.transpose(stiffness_record), delimiter=" ", fmt="%s")

    plt.show()
