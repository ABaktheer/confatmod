'''
Created on 23.08.2018

@author: abaktheer

Uniaxial fatigue model for concrete - stress driven
'''
import matplotlib.pyplot as plt
import numpy as np


def get_stress_strain(sig_arr, sig_0, E1, E2, K, gamma, S, c, r):

    # arrays to store the values
    eps_arr = np.zeros_like(sig_arr)
    eps_pi_arr = np.zeros_like(sig_arr)
    sig_pi_arr = np.zeros_like(sig_arr)
    w_arr = np.zeros_like(sig_arr)
    eps_p_cum = np.zeros_like(sig_arr)
    phi_arr = np.zeros_like(sig_arr)

    # state variables
    eps_i = 0.
    alpha_i = 0.
    eps_pi_i = 0.
    z_i = 0.
    w_i = 0.

    delta_pi = 0.
    eps_p_cum_i = 0.
    Z = K * z_i
    X_i = gamma * alpha_i

    for i in range(1, len(sig_arr)):
        print('increment', i)

        sig_i = sig_arr[i]
        #sig_i_1 = sig_arr[i] / (1.0 - w_i)

        eps_i = sig_i / ((E1 + E2) * (1.0 - w_i)) + \
            eps_pi_i * (E2 / (E1 + E2))

        sig_pi_i = E2 * (eps_i - eps_pi_i)

        # Threshold
        f_pi_i = np.fabs(sig_pi_i - X_i) - sig_0 - Z

        if f_pi_i > 1e-8:

            delta_pi = f_pi_i / \
                (E2 + (gamma + K) * (1.0 - w_i))

            # update all the state variables
            eps_pi_i = eps_pi_i + delta_pi * \
                np.sign(sig_pi_i - X_i)

            eps_i = sig_i / ((E1 + E2) * (1.0 - w_i)) + \
                eps_pi_i * (E2 / (E1 + E2))

            Y_i = 0.5 * E2 * eps_i ** 2.0 + 0.5 * \
                E2 * (eps_i - eps_pi_i) ** 2.0

            delta_lamda = delta_pi * (1. - w_i)

            w_i = w_i + ((1.0 - w_i) ** c) * (delta_lamda *
                                              (Y_i / S) ** r)

#             w_i = w_i + delta_lamda * (Y_i / S) * (1. - w_i)**(-1)

            eps_i = sig_i / ((E1 + E2) * (1.0 - w_i)) + \
                eps_pi_i * (E2 / (E1 + E2))

            if w_i >= 0.99:
                break

            alpha_i = alpha_i + delta_pi * np.sign(sig_pi_i - X_i)
            z_i = z_i + delta_pi
            X_i = gamma * alpha_i
            eps_p_cum_i = eps_p_cum_i + delta_pi

        else:
            eps_p_i = eps_pi_i

#             Y_i = 0.5 * E2 * eps_i ** 2.0 + 0.5 * \
#                 E2 * (eps_i - eps_pi_i) ** 2.0
            w_i = w_i
            eps_i = sig_i / ((E1 + E2) * (1.0 - w_i)) + \
                eps_p_i * (E2 / (E1 + E2))
            alpha_i = alpha_i
            z_i = z_i
            eps_p_cum_i = eps_p_cum_i

            if w_i >= 1.0:
                break

        # Helholtz free energy
        phi_i = 0.5 * (1.0 - w_i) * E1 * eps_i**2.0 + 0.5 * (1.0 - w_i) * E2 * \
            (eps_i - eps_pi_i)**2.0 + 0.5 * K * \
            z_i ** 2.0 + 0.5 * gamma * alpha_i**2.0

#         print 'damage: ', w_i

        sig_arr[i] = sig_i
        sig_pi_arr[i] = sig_pi_i
        eps_arr[i] = eps_i
        w_arr[i] = w_i
        eps_pi_arr[i] = eps_p_i
        eps_p_cum[i] = eps_p_cum_i
        phi_arr[i] = phi_i

    return eps_arr, sig_arr, w_arr, eps_pi_arr, eps_p_cum,  i, phi_arr


if __name__ == '__main__':

    # monotonic loading
    s_levels = np.linspace(0, 60, 5000)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= 0
    #s_levels.reshape(-1, 2)[:, 1] = 40
    s_history = s_levels.flatten()

#     sig_arr_1 = np.hstack([np.linspace(s_history[i], s_history[i + 1], 100)
#                            for i in range(len(s_levels) - 1)])

    sig_arr_1 = np.zeros(1)
    for i in range(len(s_levels) - 1):
        sig_part_1 = np.linspace(s_history[i], s_history[i + 1], 100)
        sig_arr_1 = np.hstack((sig_arr_1, sig_part_1[:-1]))

    # fatigue loading
    fc = 30
    sig_max = 0.9 * fc
    sig_min = 0.0 * fc
    cycles = 100
    inc = 100.
    s_levels_2 = np.linspace(0, 100, cycles)
    s_levels_2.reshape(-1, 2)[:, 0] = sig_min
    s_levels_2.reshape(-1, 2)[:, 1] = sig_max
    s_levels_2[0] = 0.0
    s_history_2 = s_levels_2.flatten()

    sig_arr_2 = np.zeros(1)
    for i in range(len(s_levels_2) - 1):
        sig_part = np.linspace(s_history_2[i], s_history_2[i + 1], inc)
        sig_arr_2 = np.hstack((sig_arr_2, sig_part[:-1]))

    print(sig_arr_1)

    sig_0 = 1.0
    E1 = 20000.
    E2 = 15000.
    K = 0.0
    gamma = 0.0
    S = 0.000324

    # material parameters input
    eps_arr_1, sig_arr_1, w_arr_1, eps_pi_arr_1, eps_pi_cum_1, i1, phi_arr = get_stress_strain(
        sig_arr_1, sig_0=sig_0, E1=E1, E2=E2, K=K, gamma=gamma, S=S, c=-1.0, r=1.0)
    eps_arr_2, sig_arr_2, w_arr_2, eps_pi_arr_2, eps_pi_cum_2, i2, phi_arr = get_stress_strain(
        sig_arr_1, sig_0=sig_0, E1=E1, E2=E2, K=K, gamma=gamma, S=S, c=1.0, r=1.0)
    eps_arr_3, sig_arr_3, w_arr_3, eps_pi_arr_3, eps_pi_cum_3, i3, phi_arr = get_stress_strain(
        sig_arr_1, sig_0=sig_0, E1=E1, E2=E2, K=K, gamma=gamma, S=S, c=1.25, r=1.0)
    eps_arr_4, sig_arr_4, w_arr_4, eps_pi_arr_4, eps_pi_cum_4, i4, phi_arr = get_stress_strain(
        sig_arr_1, sig_0=sig_0, E1=E1, E2=E2, K=K, gamma=gamma, S=S, c=1.5, r=1.0)
    eps_arr_5, sig_arr_5, w_arr_5, eps_pi_arr_5, eps_pi_cum_5, i5, phi_arr = get_stress_strain(
        sig_arr_1, sig_0=sig_0, E1=E1, E2=E2, K=K, gamma=gamma, S=S, c=1.75, r=1.0)
    eps_arr_6, sig_arr_6, w_arr_6, eps_pi_arr_6, eps_pi_cum_6, i6, phi_arr = get_stress_strain(
        sig_arr_1, sig_0=sig_0, E1=E1, E2=E2, K=K, gamma=gamma, S=S, c=2.0, r=1.0)
    eps_arr_7, sig_arr_7, w_arr_7, eps_pi_arr_7, eps_pi_cum_7, i7, phi_arr = get_stress_strain(
        sig_arr_1, sig_0=sig_0, E1=E1, E2=E2, K=K, gamma=gamma, S=S, c=2.25, r=1.0)
    eps_arr_8, sig_arr_8, w_arr_8, eps_pi_arr_8, eps_pi_cum_8, i8, phi_arr = get_stress_strain(
        sig_arr_1, sig_0=sig_0, E1=E1, E2=E2, K=K, gamma=gamma, S=S, c=2.5, r=1.0)
    eps_arr_9, sig_arr_9, w_arr_9, eps_pi_arr_9, eps_pi_cum_9, i9, phi_arr = get_stress_strain(
        sig_arr_1, sig_0=sig_0, E1=E1, E2=E2, K=K, gamma=gamma, S=S, c=2.75, r=1.0)
    eps_arr_10, sig_arr_10, w_arr_10, eps_pi_arr_10, eps_pi_cum_10, i10, phi_arr = get_stress_strain(
        sig_arr_1, sig_0=sig_0, E1=E1, E2=E2, K=K, gamma=gamma, S=S, c=3.0, r=1.0)


#     idx = np.where(sig_arr_2 == sig_max)
#     eps_arr_max_1 = eps_arr_2[idx]
#     idx_2 = np.where(eps_arr_max_1 > 0.0)
#     eps_arr_max = eps_arr_max_1[idx_2]
#     phi_arr_n = phi_arr[idx_2]
#     n = eps_arr_max.size
#     n_arr = np.arange(1, n)
#
#     idx_min = np.where(sig_arr_2 == sig_min)
#     eps_arr_min_1 = eps_arr_2[idx_min]
#
#     idx_min_2 = np.where(eps_arr_min_1 >= 0.0)
#     eps_arr_min = eps_arr_min_1[idx_min_2]


#     np.savetxt(r'E:\Publishing\Educative_fatigue_Article\results\Desmorat\saved\N.txt',
#                n_arr, delimiter=" ", fmt="%s")
#     np.savetxt(r'E:\Publishing\Educative_fatigue_Article\results\Desmorat\saved\eps_max.txt',
#                eps_arr_max, delimiter=" ", fmt="%s")
#     np.savetxt(r'E:\Publishing\Educative_fatigue_Article\results\Desmorat\saved\w.txt',
#                w_arr_2, delimiter=" ", fmt="%s")
#     np.savetxt(r'E:\Publishing\Educative_fatigue_Article\results\Desmorat\saved\eps_pi_cum.txt',
#                eps_pi_cum_2, delimiter=" ", fmt="%s")

    idx_1 = np.where(sig_arr_1 == sig_max)
    eps_arr_max = eps_arr_2[idx_1]
    w_arr_max = w_arr_2[idx_1]
    eps_pi_cum_max = eps_pi_cum_1[idx_1]
    idx_11 = np.where(eps_arr_max > 0.0)
    eps_arr_max_n_1 = eps_arr_max[idx_11]
    phi_arr_n_1 = phi_arr[idx_11]
    w_arr_n_1 = w_arr_max[idx_11]
    eps_pi_cum_n_1 = eps_pi_cum_max[idx_11]
    n_1 = eps_arr_max_n_1.size
    n_arr_1 = np.arange(1, n_1)

    idx_10 = np.where(sig_arr_10 == sig_max)
    eps_arr_max = eps_arr_2[idx_10]
    w_arr_max = w_arr_2[idx_10]
    eps_pi_cum_max = eps_pi_cum_10[idx_10]
    idx_110 = np.where(eps_arr_max > 0.0)
    eps_arr_max_n_10 = eps_arr_max[idx_110]
    phi_arr_n_10 = phi_arr[idx_110]
    w_arr_n_10 = w_arr_max[idx_110]
    eps_pi_cum_n_10 = eps_pi_cum_max[idx_110]
    n_10 = eps_arr_max_n_10.size
    n_arr_10 = np.arange(1, n_10)


#===========================================
# plot slip-tau
#===========================================
    ax1 = plt.subplot(231)

    ax1.plot(eps_arr_1[0:i1], sig_arr_1[0:i1], 'k', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax1.plot(eps_arr_2[0:i2], sig_arr_2[0:i2], 'r', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax1.plot(eps_arr_3[0:i3], sig_arr_3[0:i3], 'g', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax1.plot(eps_arr_4[0:i4], sig_arr_4[0:i4], 'b', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax1.plot(eps_arr_5[0:i5], sig_arr_5[0:i5], 'y', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax1.plot(eps_arr_6[0:i6], sig_arr_6[0:i6], 'm', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax1.plot(eps_arr_7[0:i7], sig_arr_7[0:i7], 'm', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)

    ax1.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax1.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.xlabel('Strain')
    plt.ylabel('Stress(MPa)')
    # plt.legend(loc=4)

#============================================
# plot slip-damage
#============================================
    ax2 = plt.subplot(232)
    ax2.plot(eps_arr_1, w_arr_1, 'k', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax2.plot(eps_arr_2[0:i2], w_arr_2[0:i2], 'r', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)

    ax2.axhline(y=0, color='k', linewidth=1, alpha=1)
    ax2.axvline(x=0, color='k', linewidth=1, alpha=1)
    plt.title('Damage evolution')
    #plt.ylim(0, 1)
    plt.xlabel('Slip(mm)')
    plt.ylabel('Damage')
    # plt.legend(loc=4)


#================================================
# plot cum_sliding-damage
#================================================
    ax4 = plt.subplot(233)

    ax4.plot(eps_pi_cum_1[0:i1], w_arr_1[0:i1], 'k', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax4.plot(eps_pi_cum_2[0:i2], w_arr_2[0:i2], 'r', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax4.plot(eps_pi_cum_3[0:i3], w_arr_3[0:i3], 'g', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax4.plot(eps_pi_cum_4[0:i4], w_arr_4[0:i4], 'b', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax4.plot(eps_pi_cum_5[0:i5], w_arr_5[0:i5], 'y', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax4.plot(eps_pi_cum_6[0:i6], w_arr_6[0:i6], 'm', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax4.plot(eps_pi_cum_7[0:i7], w_arr_7[0:i7], '--r', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax4.plot(eps_pi_cum_8[0:i8], w_arr_8[0:i8], 'r', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax4.plot(eps_pi_cum_9[0:i9], w_arr_9[0:i9], '--k', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax4.plot(eps_pi_cum_10[0:i10], w_arr_10[0:i10], '--k', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)

    plt.xlabel('Cumulative sliding(mm)')
    plt.ylabel('Damage')
    plt.ylim(0, 1)
    plt.xlim(0, 1)
    # plt.legend()


# #================================================
# # plot HFE
# #================================================
#     ax4 = plt.subplot(235)
#
#     ax4.plot(n_arr, phi_arr_n[:-1], '--k', linewidth=1,
#              label='$ \sigma_N = 0$ MPa', alpha=1)
#     ax4.plot(n_arr, phi_arr_n[:-1], 'r', linewidth=1,
#              label='$ \sigma_N = 0$ MPa', alpha=1)
#
#     plt.xlabel('number of cycles')
#     plt.ylabel('Damage')
#     #plt.ylim(0, 1)
#     plt.legend()


#================================================
# plot damage with cycles
#================================================
    ax4 = plt.subplot(235)

    ax4.plot(n_arr_1, w_arr_n_1[:-1], '--k', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax4.plot(n_arr_10, w_arr_n_10[:-1], 'r', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)

    plt.xlabel('number of cycles')
    plt.ylabel('Damage')
    #plt.ylim(0, 1)
    plt.legend()


#================================================
# plot cumulative sliding with cycles
#================================================
    ax4 = plt.subplot(236)

    ax4.plot(n_arr_1, eps_pi_cum_n_1[:-1], '--k', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax4.plot(n_arr_10, eps_pi_cum_n_10[:-1], 'r', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)

    plt.xlabel('number of cycles')
    plt.ylabel('cumulative sliding')
    #plt.ylim(0, 1)
    plt.legend()


# #===========================================
# # plot fatigue creep
# #===========================================
#     ax1 = plt.subplot(234)
#
#     ax1.plot(n_arr, eps_arr_max[:-1], 'r', linewidth=1, alpha=1)
#
#     #ax1.axhline(y=0, color='k', linewidth=1, alpha=0.5)
#     #ax1.axvline(x=0, color='k', linewidth=1, alpha=0.5)
#     plt.xlabel('N')
#     plt.ylabel('strain')
#     plt.legend(loc=4)
#
#
# #================================================
# # plot HFE
# #================================================
#     ax4 = plt.subplot(235)
#
#     ax4.plot(n_arr, phi_arr_n[:-1], 'k', linewidth=1,
#              label='$ \sigma_N = 0$ MPa', alpha=1)
#     ax4.plot(n_arr, phi_arr_n[:-1], 'r', linewidth=1,
#              label='$ \sigma_N = 0$ MPa', alpha=1)
#
#     plt.xlabel('Cumulative sliding(mm)')
#     plt.ylabel('Damage')
#     #plt.ylim(0, 1)
#     plt.legend()

# #===========================================
# # plot stiffness
# #===========================================
#     ax1 = plt.subplot(235)
#
#     stifness_0 = sig_max / eps_arr_max[1]
#
#     ax1.plot(
#         n_arr / (1.0 * n_arr[-1]),  (sig_max / (stifness_0 * eps_arr_max[:-1])), 'r', linewidth=1, alpha=1)
#
#     #ax1.axhline(y=0, color='k', linewidth=1, alpha=0.5)
#     #ax1.axvline(x=0, color='k', linewidth=1, alpha=0.5)
#     plt.xlabel('N')
#     plt.ylabel('stiffness')
#     plt.legend(loc=4)

# #===========================================
# # plot fatigue creep ( normalized)
# #===========================================
#     ax1 = plt.subplot(236)
#
#     ax1.plot(n_arr / (1.0 * n_arr[-1]),
#              eps_arr_max[:-1], 'r', linewidth=1, alpha=1)
#
#     n_1, s_1 = np.loadtxt(
#         r'E:\Publishing\Uniaxial_fatigue_model_Article\results\Do_1993\data\concrete_A\creep_fatigue\A_creep_fatigue_075_1.txt')
#     n_2, s_2 = np.loadtxt(
#         r'E:\Publishing\Uniaxial_fatigue_model_Article\results\Do_1993\data\concrete_A\creep_fatigue\A_creep_fatigue_075_2.txt')
#
#     ax1.plot(n_1, s_1, '--k', label='S=0.75')
#     ax1.plot(n_2, s_2, '--k', label='S=0.75')
#
#     plt.xlabel('N')
#     plt.ylabel('strain')
#     plt.legend(loc=4)

    plt.show()
