'''
Created on 23.08.2018

@author: abaktheer

Uniaxial fatigue model for concrete - stress driven
'''
import matplotlib.pyplot as plt
import numpy as np


def get_stress_strain(sig_arr, sig_0, E, K, gamma, S, c, r, m):

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
    m = m

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
        print('increment', i)

        sig_i = sig_arr[i]
        sig_i_1 = sig_arr[i] / (1.0 - w_i)

        # Threshold
        f_pi_i = np.fabs(sig_i_1 - X_i) - sig_0 - Z

        if f_pi_i > 1e-8:

            delta_lamda = f_pi_i / \
                (E / (1.0 - w_i) + gamma * (1.0 + m * X_i * np.fabs(sig_i_1 - X_i)) + K)

            # update all the state variables
            eps_p_i = eps_p_i + delta_lamda * \
                np.sign(sig_i_1 - X_i) / (1. - w_i)

            Y_i = 0.5 * E * (eps_i - eps_p_i) ** 2.0

#             w_i = w_i + ((1.0 - w_i) ** c) * (delta_lamda *
#                                               (Y_i / S) ** r)

            w_i = w_i + ((1.0 - w_i) ** c) * (delta_lamda * (Y_i / S) ** r) * \
                (np.exp(-35.0 * w_i) + np.exp(-10.0 * (1. - w_i)))

            eps_i = sig_i / (E * (1.0 - w_i)) + eps_p_i

            if eps_i >= 0.005:
                break

            alpha_i = alpha_i + delta_lamda * (np.sign(sig_i - X_i) + m * X_i)
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

            if eps_i >= 0.005:
                break

        print('damage: ', w_i)

        sig_arr[i] = sig_i
        eps_arr[i] = eps_i
        w_arr[i] = w_i
        eps_p_arr[i] = eps_p_i
        eps_p_cum[i] = eps_p_cum_i

    return eps_arr, sig_arr, w_arr, eps_p_arr, eps_p_cum,  i


if __name__ == '__main__':

    # monotonic loading
    s_levels = np.linspace(0, 130, 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= 0
    # s_levels.reshape(-1, 2)[:, 1] = 2
    s_history = s_levels.flatten()

    sig_arr_1 = np.hstack([np.linspace(s_history[i], s_history[i + 1], 500)
                           for i in range(len(s_levels) - 1)])

    # fatigue loading
    sig_max = 85.5
    sig_min = 10
    cycles = 10000
    s_levels_2 = np.linspace(0, 70, cycles)
    s_levels_2.reshape(-1, 2)[:, 0] = sig_min
    s_levels_2.reshape(-1, 2)[:, 1] = sig_max
    s_levels_2[0] = 0
    s_history_2 = s_levels_2.flatten()
    sig_arr_2 = np.hstack([np.linspace(s_history_2[i], s_history_2[i + 1], 100)
                           for i in range(len(s_levels_2) - 1)])

    # material parameters input
    eps_arr_1, sig_arr_1, w_arr_1, eps_pi_arr_1, eps_p_cum_1, i1 = get_stress_strain(
        sig_arr_1, sig_0=40, K=5000, gamma=200000, E=50000, S=0.00001, c=1.5, r=0.8, m=0.0)
    eps_arr_2, sig_arr_2, w_arr_2, eps_pi_arr_2, eps_pi_cum_2, i2 = get_stress_strain(
        sig_arr_2, sig_0=40, K=5000, gamma=200000, E=50000, S=0.00001,
        c=1.5, r=0.8, m=0.0)

    idx = np.where(sig_arr_2 == sig_max)
    eps_arr_max_1 = eps_arr_2[idx]
    idx_2 = np.where(eps_arr_max_1 > 0.0)
    eps_arr_max = eps_arr_max_1[idx_2]
    n = eps_arr_max.size
    n_arr = np.arange(1, n)

    print(eps_arr_max)
    print(n)
    print(n_arr)


#===========================================
# plot slip-tau
#===========================================
    ax1 = plt.subplot(221)

    ax1.plot(eps_arr_1, sig_arr_1, '--k', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax1.plot(eps_arr_2[0:i2], sig_arr_2[0:i2], 'r', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)

    ax1.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax1.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.xlabel('Strain')
    plt.ylabel('Stress(MPa)')
    plt.legend(loc=4)

#============================================
# plot slip-damage
#============================================
    ax2 = plt.subplot(222)
    ax2.plot(eps_arr_1, w_arr_1, '--k', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax2.plot(eps_arr_2[0:i2], w_arr_2[0:i2], 'r', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)

    ax2.axhline(y=0, color='k', linewidth=1, alpha=1)
    ax2.axvline(x=0, color='k', linewidth=1, alpha=1)
    plt.title('Damage evolution')
    plt.ylim(0, 1)
    plt.xlabel('Slip(mm)')
    plt.ylabel('Damage')
    plt.legend(loc=4)


#===========================================
# plot fatigue creep
#===========================================
    ax1 = plt.subplot(223)

    ax1.plot(n_arr, eps_arr_max[:-1], 'r', linewidth=1, alpha=1)

    #ax1.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    #ax1.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.xlabel('N')
    plt.ylabel('strain')
    plt.legend(loc=4)


#================================================
# plot cum_sliding-damage
#================================================
    ax4 = plt.subplot(224)

    ax4.plot(eps_p_cum_1[0:i1], w_arr_1[0:i1], '--k', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax4.plot(eps_pi_cum_2[0:i2], w_arr_2[0:i2], 'r', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)

    plt.xlabel('Cumulative sliding(mm)')
    plt.ylabel('Damage')
    plt.ylim(0, 1)
    plt.legend()

    plt.show()
