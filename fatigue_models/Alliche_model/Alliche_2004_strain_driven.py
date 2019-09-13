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


def get_stress_strain(eps_1_arr, lamda, mu, alpha, beta, g, C0, C1, K, n):

    #-----------------------------------------------------------------------
    # arrays to store the values
    #-----------------------------------------------------------------------
    # normal strain
    sig_1_arr = np.zeros_like(eps_1_arr, dtype=np.float_)
    # lateral strain
    eps_2_arr = np.zeros_like(eps_1_arr, dtype=np.float_)
    # damage factor
    w_arr = np.zeros_like(eps_1_arr, dtype=np.float_)
    f_arr = np.zeros_like(eps_1_arr, dtype=np.float_)
    D_arr = np.zeros_like(eps_1_arr, dtype=np.float_)

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
    sig_1_i = 0.0
    eps_2_i = 0.0
    w_i = 0.0
    D_i = 0.0

    for i in range(1, len(eps_1_arr)):
        # print 'increment', i

        eps_1_i = eps_1_arr[i]

        eps_2_i = -1.0 * ((lamda + alpha * w_i) * eps_1_i + g * w_i) / \
            (2.0 * (lamda + mu) + 4.0 * w_i * (alpha + beta))

        f_i = abs(g) * eps_2_i - (C0 + 2.0 * C1 * w_i)

        m = -1.0 * (eps_2_arr[i] - eps_2_arr[i - 1])
        # loading stage (evolve of the fatigue damage based on (Marigo.85)
        # model)
        if m > 0:
            d_w = m * abs(g) / (2.0 * C1) * (f_i / K)**n
        else:  # unloading stage (no fatigue damage)
            d_w = 0

        w_i = w_i + d_w

        # Energy release rate
        Y_norm = np.sqrt((-g * eps_1_i - alpha * (eps_1_i + 2. * eps_2_i) * eps_1_i
                          - 2. * beta * (eps_1_i**2.0))**2.0 + 2.0 * (-g * eps_2_i - alpha * (eps_1_i + 2. * eps_2_i) * eps_1_i
                                                                      - 2 * beta * (eps_1_i**2.0))**2.0)
        # dissipation
        d_D = Y_norm * d_w
        D_i += d_D
        print('D=', D_i)

        if w_i >= 5.0:
            print(' ----------> No Convergence any more')
            print(i)
            break

        sig_1_i = (lamda + 2.0 * mu) * eps_1_i + 2.0 * \
            (lamda + alpha * w_i) * eps_2_i

        sig_1_arr[i] = sig_1_i
        eps_2_arr[i] = eps_2_i
        w_arr[i] = w_i
        f_arr[i] = f_i
        D_arr[i] = D_i

    return sig_1_arr, eps_2_arr, w_arr, f_arr, D_arr, i


if __name__ == '__main__':

    eps_max = - 0.001
    number_of_cycles = 51
    eps_min = - 0.0005

    m = 2000  # number of increments in each cycle
    n = number_of_cycles

    s_levels = np.linspace(0, eps_max, n * 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] = eps_max
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 1] = eps_min
    s_history = s_levels.flatten()

#     eps_1_arr = np.hstack([np.linspace(s_history[i], s_history[i + 1], m, dtype=np.float_)
#                            for i in range(len(s_history) - 1)])

    eps_1_arr = np.zeros(1)
    for i in range(len(s_levels) - 1):
        eps_part = np.linspace(s_history[i], s_history[i + 1], m)
        eps_1_arr = np.hstack((eps_1_arr, eps_part[:-1]))

    t_arr = np.linspace(0, 1, len(eps_1_arr))

#     sig_1_arr, eps_2_arr, w_arr, f_arr, D_arr, inc = get_stress_strain(
# eps_1_arr, lamda=10555.5, mu=15833.33, alpha=2237.5, beta=-2216.5,
# g=-9.788, C0=0.00, C1=0.002033, K=0.00389, n=10)

    sig_1_arr, eps_2_arr, w_arr, f_arr, D_arr, inc = get_stress_strain(
        eps_1_arr, lamda=10277.778, mu=15416.667, alpha=4237.5, beta=-2216.5, g=-10.0, C0=0.00, C1=0.0013, K=0.00295, n=9)

#
    idx_1 = np.where(eps_1_arr == eps_max)

    sig_arr_max_1 = sig_1_arr[idx_1]
    D_tot_arr_max_1 = D_arr[idx_1]
    idx_2 = np.where(sig_arr_max_1 < 0.0)

    sig_arr_max = sig_arr_max_1[idx_2]
    D_tot_arr_max = D_tot_arr_max_1[idx_2]

    n_1 = sig_arr_max.size
    n_arr = np.arange(1, n_1)

    #-----------------------------------------------------------------------
    # plot 1
    #-----------------------------------------------------------------------
    plt.subplot(231)
    plt.plot(t_arr, eps_1_arr, 'k', linewidth=0.5, alpha=1.0)
    plt.title('Loading history')
    plt.xlabel('Time')
    plt.ylabel('$\epsilon_{11}$')

    #-----------------------------------------------------------------------
    # plot 2
    #-----------------------------------------------------------------------
    plt.subplot(232)
    plt.plot(eps_1_arr, sig_1_arr, 'k', linewidth=0.5, alpha=1.0)
    plt.title('$ \epsilon_{11} - \sigma_{11}$')
    plt.xlabel('$\epsilon_{11}$')
    plt.ylabel('$\sigma_{11}$[MPa]')

    #-----------------------------------------------------------------------
    # plot 3
    #-----------------------------------------------------------------------
    plt.subplot(233)
    plt.plot(eps_1_arr, w_arr, 'k', linewidth=0.5, alpha=1)
    plt.xlabel('$\epsilon_{11}$')
    plt.ylabel('Damage')
    #plt.ylim(0, 1)
    # plt.legend()

    #-----------------------------------------------------------------------
    # plot 4
    #-----------------------------------------------------------------------
    plt.subplot(234)
#     sig_1 = np.zeros(n)
#     cycle = np.zeros(n)
#     for i in range(0, n, 1):
#         idx = m + 2 * i * m - 1
#         sig_1[i] = sig_1_arr[idx]
#         cycle[i] = i
#     plt.plot(cycle, -1 * sig_1, 'k', linewidth=1, alpha=1)
    plt.plot(n_arr, -sig_arr_max[:-1], 'r', linewidth=1, alpha=1)
    plt.xlabel('number of cycles')
    plt.ylabel('max $\epsilon_{11}$')
    #plt.ylim(0, 1)

    #-----------------------------------------------------------------------
    # plot 5
    #-----------------------------------------------------------------------
    plt.subplot(235)
#
#     D = np.zeros(n)
#     cycle = np.zeros(n)
#     for i in range(0, n, 1):
#         idx = m + 2 * i * m   - 1
#         D[i] = D_arr[idx]
#         cycle[i] = i
#     plt.plot(cycle[1:], D[1:], 'b', linewidth=1, alpha=1)
    plt.plot(n_arr, D_tot_arr_max[:-1], 'k', linewidth=1.0, alpha=1.0)
    plt.xlabel('number of cycles')
    plt.ylabel('max $\epsilon_{22}$')
    # plt.legend()

    #-----------------------------------------------------------------------
    # plot 7
    #-----------------------------------------------------------------------
    plt.subplot(236)
    plt.plot(eps_2_arr, sig_1_arr, 'k', linewidth=1.0, alpha=1.0)
    plt.title('$ \epsilon_{22} - \sigma_{11}$')
    plt.xlabel('$\epsilon_{22}$')
    plt.ylabel('$\sigma_{11}$[MPa]')
    # plt.legend(loc=4)

    plt.show()
