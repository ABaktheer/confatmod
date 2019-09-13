'''
Created on 14.06.2017

@author: abaktheer
'''

'''
Implementation of the fatigue model for plain concrete [A.Alliche, 2004] under uniaxial compressive loading
(stress driven algorithm)
'''

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
    f_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)

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
    #sigma_1_arr[0] = 0
    eps_1_i = 0.0
    eps_2_i = 0.0
    w_i = 0.0

    for i in range(1, len(sigma_1_arr)):
        # print 'increment', i

        sigma_1_i = sigma_1_arr[i]

        eps_2_i = -1.0 * ((lamda + alpha * w_i) * sigma_1_i + g * w_i * (lamda + 2.0 * mu)) / \
                         ((lamda + 2.0 * mu) * (2.0 * (lamda + mu) + 4.0 *
                                                w_i * (alpha + beta)) - 2.0 * (lamda + alpha * w_i) ** 2)

        eps_1_i = sigma_1_i / \
            (lamda + 2.0 * mu) - 2.0 * eps_2_i *  \
            (lamda + alpha * w_i) / (lamda + 2.0 * mu)

        f_i = abs(g) * eps_2_i - (C0 + 2 * C1 * w_i)

        kappa_i = (lamda + 2.0 * mu) * (2.0 * (lamda + mu) + 4.0 * w_i * (alpha + beta) -
                                        alpha * (g / (2.0 * C1)) * (2.0 * eps_2_i + eps_1_i) -
                                        (g**2.0 / (2.0 * C1))) - 2.0 * (lamda + alpha * w_i)**2

        # print 'kappa_i', kappa_i

        d_sigma_1 = sigma_1_arr[i] - sigma_1_arr[i - 1]
        m = -1.0 * ((lamda + alpha * w_i) / kappa_i) * d_sigma_1

        # loading stage (damage is active) based on (Marigo.85) model)\
        if m > 0:
            d_w = m * abs(g) / (2.0 * C1) * (f_i / K)**n
        else:  # unloading stage (no fatigue damage)
            d_w = 0

        w_i = w_i + d_w

        if w_i > 5.0:
            break

        if abs(eps_1_i) > 0.0035:
            break

        eps_1_arr[i] = eps_1_i
        eps_2_arr[i] = eps_2_i
        w_arr[i] = w_i
        f_arr[i] = f_i

    return sigma_1_arr, eps_1_arr, eps_2_arr, w_arr, f_arr


if __name__ == '__main__':

    maximum_stress = -120  # [MPa]
    number_of_cycles = 10000  # 4554
    minimum_stress = -21
    #unloading_ratio = 0.33
    m = 100  # number of increments in each cycle
    n = number_of_cycles

    s_levels = np.linspace(0, maximum_stress, n * 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] = minimum_stress
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 1] = maximum_stress
    s_history = s_levels.flatten()

    sigma_1_arr = np.hstack([np.linspace(s_history[i], s_history[i + 1], m, dtype=np.float_)
                             for i in range(len(s_levels) - 1)])

    t_arr = np.linspace(0, 1, len(sigma_1_arr))

    sigma_1_arr, eps_1_arr, eps_2_arr, w_arr, f_arr = get_stress_strain(
        sigma_1_arr, lamda=13972.2, mu=20958.3, alpha=2237.5, beta=-2216.5, g=-10, C0=0.00, C1=0.00188, K=0.00334, n=10)

    #-----------------------------------------------------------------------
    # plot 1
    #-----------------------------------------------------------------------
    plt.subplot(231)
    plt.plot(t_arr, sigma_1_arr, 'k', linewidth=0.5, alpha=1.0)
    plt.title('Loading history')
    plt.xlabel('Time')
    plt.ylabel('$\sigma_{11}$')

    #-----------------------------------------------------------------------
    # plot 2
    #-----------------------------------------------------------------------
    plt.subplot(232)
    plt.plot(eps_1_arr, sigma_1_arr, 'k', linewidth=0.5, alpha=1.0)
    plt.title('$ \epsilon_{11} - \sigma_{11}$')
    plt.xlabel('$\epsilon_{11}$')
    plt.ylabel('$\sigma_{11}$[MPa]')

    #-----------------------------------------------------------------------
    # plot 3
    #-----------------------------------------------------------------------
    plt.subplot(233)
    plt.plot(sigma_1_arr, w_arr, 'k', linewidth=0.5, alpha=1)
    plt.xlabel('$\sigma_{11}$')
    plt.ylabel('Damage')
    #plt.ylim(0, 1)
    # plt.legend()

    #-----------------------------------------------------------------------
    # plot 4
    #-----------------------------------------------------------------------
    plt.subplot(234)
    eps_1 = np.zeros(n)
    cycle = np.zeros(n)
    for i in range(0, n, 1):
        idx = m + 2 * i * m - 1
        eps_1[i] = eps_1_arr[idx]
        cycle[i] = i
    plt.plot(cycle, eps_1, 'k', linewidth=1, alpha=1)
    plt.xlabel('number of cycles')
    plt.ylabel('max $\epsilon_{11}$')
    #plt.ylim(0, 1)

    #-----------------------------------------------------------------------
    # plot 5
    #-----------------------------------------------------------------------
#     plt.subplot(235)
#
#     eps_2 = np.zeros(n)
#     cycle = np.zeros(n)
#     for i in range(0, n, 1):
#         idx = m + 2 * i * m - 1
#         eps_2[i] = eps_2_arr[idx]
#         cycle[i] = i
#     plt.plot(cycle, eps_2, 'b', linewidth=1, alpha=1)
#     plt.xlabel('number of cycles')
#     plt.ylabel('max $\epsilon_{22}$')
#     # plt.legend()

#     #-----------------------------------------------------------------------
#     # plot 6
#     #-----------------------------------------------------------------------
    plt.subplot(235)
    w = np.zeros(n)
    cycle = np.zeros(n)
    for i in range(0, n, 1):
        idx = m + 2 * i * m - 1

        w[i] = w_arr[idx]
        cycle[i] = i
    plt.plot(cycle, w, 'k', linewidth=1, alpha=1)
    plt.xlabel('number of cycles')
    plt.ylabel('Damage')
    #plt.ylim(0, 1)
    # plt.legend()

    #-----------------------------------------------------------------------
    # plot 7
    #-----------------------------------------------------------------------
    plt.subplot(236)
    plt.plot(eps_2_arr, sigma_1_arr, 'k', linewidth=0.5, alpha=1.0)
    plt.title('$ \epsilon_{22} - \sigma_{11}$')
    plt.xlabel('$\epsilon_{22}$')
    plt.ylabel('$\sigma_{11}$[MPa]')
    # plt.legend(loc=4)

    plt.show()
