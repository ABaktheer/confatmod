'''
Created on 23.08.2018

@author: abaktheer

Uniaxial fatigue model for concrete - stress driven
'''
import matplotlib.pyplot as plt
import numpy as np


def get_stress_strain(sig_arr, sig_0, E, K, gamma, S, c, r):

    # arrays to store the values
    eps_arr = np.zeros_like(sig_arr)
    eps_p_arr = np.zeros_like(sig_arr)
    w_arr = np.zeros_like(sig_arr)
    eps_p_cum = np.zeros_like(sig_arr)

    # material parameters
    E = E
    K = K
    gamma = gamma
    sig_0 = sig_0
    S = S
    c = c
    r = r

    # state variables
    eps_i = 0.
    alpha_i = 0.
    eps_p_i = 0.
    z_i = 0.
    w_i = 0.

    delta_lamda = 0.
    eps_p_cum_i = 0.
    Z = K * z_i
    X_i = gamma * alpha_i

    for i in range(1, len(sig_arr)):
        print 'increment', i

        sig_i = sig_arr[i]
        sig_i_1 = sig_arr[i] / (1.0 - w_i)

        # Threshold
        f_pi_i = np.fabs(sig_i_1 - X_i) - sig_0 - Z

        if f_pi_i > 1e-8:

            delta_lamda = f_pi_i / \
                (E / (1.0 - w_i) + gamma + K)

            # update all the state variables
            eps_p_i = eps_p_i + delta_lamda * \
                np.sign(sig_i_1 - X_i) / (1. - w_i)

            eps_i = sig_i / (E * (1.0 - w_i)) + eps_p_i

            Y_i = 0.5 * E * (eps_i - eps_p_i) ** 2.0

#             w_i = w_i + ((1.0 - w_i) ** c) * (delta_lamda *
#                                               (Y_i / S) ** r)

            w_i = w_i + ((1.0 - w_i * 0.0) ** c) * (delta_lamda * (Y_i / S) ** r) * \
                (np.exp(-5.0 * w_i) + np.exp(-50. * (1. - w_i)))

            eps_i = sig_i / (E * (1.0 - w_i)) + eps_p_i

            print 'eps_i', eps_i
            if eps_i >= 0.0035:
                print '***************************'
                break

            alpha_i = alpha_i + delta_lamda * np.sign(sig_i - X_i)
            z_i = z_i + delta_lamda
            X_i = gamma * alpha_i
            eps_p_cum_i = eps_p_cum_i + delta_lamda

        else:
            eps_p_i = eps_p_i
            Y_i = 0.5 * E * (eps_i - eps_p_i) ** 2
            w_i = w_i
            eps_i = sig_i / (E * (1.0 - w_i)) + eps_p_i
            alpha_i = alpha_i
            z_i = z_i
            eps_p_cum_i = eps_p_cum_i

            print 'eps_i', eps_i
            if eps_i >= 0.0035:
                print '******************************'
                break

        print 'damage: ', w_i

        sig_arr[i] = sig_i
        eps_arr[i] = eps_i
        w_arr[i] = w_i
        eps_p_arr[i] = eps_p_i
        eps_p_cum[i] = eps_p_cum_i

    return eps_arr, sig_arr, w_arr, eps_p_arr, eps_p_cum,  i


if __name__ == '__main__':

    m = 100  # number of increments in each cycle

    n1 = 100
    n2 = 649
    n3 = 1
    n4 = 100
    n5 = 100
    n6 = 10000

    b = 1  # number_of_repeted_blocks

    sigma_u = 100

    stress_level_1_max = 0.6 * sigma_u
    stress_level_2_max = 0.75 * sigma_u
    stress_level_3_max = 0.55 * sigma_u
    stress_level_4_max = 0.55 * sigma_u
    stress_level_5_max = 0.55 * sigma_u
    stress_level_6_max = 0.55 * sigma_u

    stress_level_1_min = 0.2 * sigma_u
    stress_level_2_min = 0.1 * sigma_u
    stress_level_3_min = 0.1 * sigma_u
    stress_level_4_min = 0.1 * sigma_u
    stress_level_5_min = 0.1 * sigma_u
    stress_level_6_min = 0.1 * sigma_u
    #unloading_ratio = 0.0

    d_0 = np.zeros(1)

    d_1 = np.linspace(0, stress_level_1_max, n1 * 2)
    d_1.reshape(-1, 2)[:, 0] = stress_level_1_max
    d_1.reshape(-1, 2)[:, 1] = stress_level_1_min
    d_history_1 = d_1.flatten()
    sig_1_arr = np.hstack([np.linspace(d_history_1[i], d_history_1[i + 1], m, dtype=np.float_)
                           for i in range(len(d_1) - 1)])

    d_2 = np.linspace(0, stress_level_2_max, n2 * 2)
    d_2.reshape(-1, 2)[:, 0] = stress_level_2_max
    d_2.reshape(-1, 2)[:, 1] = stress_level_2_min
    d_history_2 = d_2.flatten()
    sig_2_arr = np.hstack([np.linspace(d_history_2[i], d_history_2[i + 1], m, dtype=np.float_)
                           for i in range(len(d_2) - 1)])

    d_3 = np.linspace(0, stress_level_3_max, n3 * 2)
    d_3.reshape(-1, 2)[:, 0] = stress_level_3_max
    d_3.reshape(-1, 2)[:, 1] = stress_level_3_min
    d_history_3 = d_3.flatten()
    sig_3_arr = np.hstack([np.linspace(d_history_3[i], d_history_3[i + 1], m, dtype=np.float_)
                           for i in range(len(d_3) - 1)])

    d_4 = np.linspace(0, stress_level_4_max, n4 * 2)
    d_4.reshape(-1, 2)[:, 0] = stress_level_4_max
    d_4.reshape(-1, 2)[:, 1] = stress_level_4_min
    d_history_4 = d_4.flatten()
    sig_4_arr = np.hstack([np.linspace(d_history_4[i], d_history_4[i + 1], m, dtype=np.float_)
                           for i in range(len(d_4) - 1)])

    d_5 = np.linspace(0, stress_level_5_max, n5 * 2)
    d_5.reshape(-1, 2)[:, 0] = stress_level_5_max
    d_5.reshape(-1, 2)[:, 1] = stress_level_5_min
    d_history_5 = d_5.flatten()
    sig_5_arr = np.hstack([np.linspace(d_history_5[i], d_history_5[i + 1], m, dtype=np.float_)
                           for i in range(len(d_5) - 1)])

    d_6 = np.linspace(0, stress_level_6_max, n6 * 2)
    d_6.reshape(-1, 2)[:, 0] = stress_level_6_max
    d_6.reshape(-1, 2)[:, 1] = stress_level_6_min
    d_history_6 = d_6.flatten()
    sig_6_arr = np.hstack([np.linspace(d_history_6[i], d_history_6[i + 1], m, dtype=np.float_)
                           for i in range(len(d_6) - 1)])

    sig_0_1_arr = np.linspace(
        d_0[-1], d_history_1[0], m, dtype=np.float_)
    sig_1_2_arr = np.linspace(
        d_history_1[-1], d_history_2[0], m, dtype=np.float_)
    sig_2_3_arr = np.linspace(
        d_history_2[-1], d_history_3[0], m, dtype=np.float_)
    sig_3_4_arr = np.linspace(
        d_history_3[-1], d_history_4[0], m, dtype=np.float_)
    sig_4_5_arr = np.linspace(
        d_history_4[-1], d_history_5[0], m, dtype=np.float_)
    sig_5_6_arr = np.linspace(
        d_history_5[-1], d_history_6[0], m, dtype=np.float_)

    sig_1_1_arr = np.linspace(
        sig_1_arr[-1], sig_1_arr[0], m, dtype=np.float_)
    sig_6_1_arr = np.linspace(
        sig_6_arr[-1], sig_1_arr[0], m, dtype=np.float_)
    sig_6_6_arr = np.linspace(
        sig_6_arr[-1], sig_6_arr[0], m, dtype=np.float_)

    sigma_1_arr = np.hstack(
        (d_0, sig_0_1_arr, sig_1_arr, sig_1_2_arr, sig_2_arr, sig_2_3_arr, sig_3_arr, sig_3_4_arr, sig_4_arr, sig_4_5_arr, sig_5_arr, sig_5_6_arr, sig_6_arr))

    sigma_arr = sigma_1_arr

    for i in range(1, b):

        # for constant repeating
        sigma_arr = np.hstack((sigma_arr, sig_6_1_arr, sig_1_arr, sig_1_2_arr, sig_2_arr, sig_2_3_arr,
                               sig_3_arr, sig_3_4_arr, sig_4_arr, sig_4_5_arr, sig_5_arr,
                               sig_5_6_arr, sig_6_arr))

#         # for repeating (H-L-H) or (L-H-L)
#         if int(i) % 2 == 0:
#             # print 'even'
#             sigma_arr = np.hstack(
#                 (sigma_arr, sig_6_1_arr, sig_1_arr, sig_1_2_arr, sig_2_arr, sig_2_3_arr, sig_3_arr, sig_3_4_arr, sig_4_arr, sig_4_5_arr, sig_5_arr, sig_5_6_arr, sig_6_arr))
#         else:
#             # print 'odd'
#             sigma_arr = np.hstack(
#                 (sigma_arr, sig_6_6_arr, sig_6_arr, sig_5_6_arr, sig_5_arr, sig_4_5_arr, sig_4_arr, sig_3_4_arr, sig_3_arr, sig_2_3_arr, sig_2_arr, sig_1_2_arr, sig_1_arr))

    # print sigma_arr.shape

    t_arr = np.linspace(0, 1, len(sigma_arr))

    eps_arr, sig_arr, w_arr, eps_p_arr, eps_p_cum,  inc = get_stress_strain(
        sigma_arr, sig_0=25, E=45000, K=600000, gamma=120000, S=10e-7, c=4.5, r=0.5)

    print 'inc', inc

    #-----------------------------------------------------------------------
    # plot 1
    #-----------------------------------------------------------------------
    plt.subplot(221)
    plt.plot(t_arr[0:inc], abs(sigma_arr[0:inc]), 'k', linewidth=1, alpha=1.0)
    plt.title('loading history')

    plt.xlabel('Time')
    plt.ylabel('$\sigma_{1}$')
    # plt.legend(loc=4)

    #-----------------------------------------------------------------------
    # plot 2
    #-----------------------------------------------------------------------
    plt.subplot(222)
    plt.plot(abs(eps_arr[0:inc]), abs(
        sigma_arr[0:inc]), 'k', linewidth=1, alpha=1.0)
    plt.title('$ \epsilon_{11} - \sigma_{11}$')
    plt.xlabel('$\epsilon_{11}$')
    plt.ylabel('$\sigma_{11}$[MPa]')
    # plt.legend(loc=4)

    #-----------------------------------------------------------------------
    # plot 3
    #-----------------------------------------------------------------------
#     plt.subplot(222)
#     plt.plot(abs(sigma_arr), w_arr, 'k', linewidth=1, alpha=1)
#     plt.xlabel('$\sigma_{11}$')
#     plt.ylabel('Damage')
#     #plt.ylim(0, 1)
#     # plt.legend()

    #-----------------------------------------------------------------------
    # plot 4
    #-----------------------------------------------------------------------
#     plt.subplot(223)
#     plt.plot(t_arr, w_arr, 'b', linewidth=1, alpha=1)
#     plt.xlabel('Time')
#     plt.ylabel('Damage')
#     #plt.ylim(0, 1)
#     # plt.legend()

    #-----------------------------------------------------------------------
    # plot 5
    #-----------------------------------------------------------------------
#     plt.subplot(225)
#
    n = b * (n1 + n2 + n3 + n4 + n5 + n6)
#     eps_1 = np.zeros(n)
#     cycle = np.zeros(n)
#     sig_1 = np.zeros(n)
#
#     for i in range(0, n, 1):
#         idx = m + 2 * i * m  # - 1
#         eps_1[i] = eps_1_arr[idx]
#         cycle[i] = i
#         #sig_1[i] = sigma_1_arr[idx]
#         # print sig_1
#     plt.plot(cycle, eps_1, 'b', linewidth=1, alpha=1)
#     plt.xlabel('number of cycles')
#     plt.ylabel('max $\epsilon_{11}$')
    # plt.legend()

#     #-----------------------------------------------------------------------
#     # plot 5
#     #-----------------------------------------------------------------------
    plt.subplot(224)

    eps_1 = np.zeros(n)
    cycle = np.zeros(n)
    for i in range(0, n, 1):
        idx = m + 2 * i * m - 1
        if idx <= len(eps_arr[0:inc]):
            idx = idx
        else:
            idx = m + 2 * (i - 1.0) * m - 1
            break

        eps_1[i] = eps_arr[idx]
        cycle[i] = i + 1

    plt.plot(cycle[0:i], abs(eps_1[0:i]), 'k', linewidth=1, alpha=1)

    #plt.plot(cycle, eps_2 / 1.5, color='gray', linewidth=1.0, alpha=0.5)
    #plt.fill_between(cycle, eps_2, eps_2 / 1.5, facecolor='gray', alpha=0.5)
    plt.xlabel('number of cycles')
    plt.ylabel('max $\epsilon_{11}$')
    #plt.ylim(0, 1)
    # plt.legend()

#     #-----------------------------------------------------------------------
#     # plot 6
#     #-----------------------------------------------------------------------
    plt.subplot(223)
    w = np.zeros(n)
    cycle = np.zeros(n)
    for i in range(0, n, 1):
        idx = m + 2 * i * m - 1
        if idx <= len(w_arr[0:inc]):
            idx = idx
        else:
            idx = m + 2 * (i - 1.0) * m - 1
            break

        w[i] = w_arr[idx]
        cycle[i] = i + 1

    plt.plot(cycle[0:i], w[0:i], 'k', linewidth=1, alpha=1)
    plt.xlabel('number of cycles')
    plt.ylabel('Damage')
    #plt.ylim(0, 1)
    # plt.legend()

    plt.show()
