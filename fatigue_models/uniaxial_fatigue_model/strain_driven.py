'''
Created on 14.11.2016

@author: abaktheer

Uniaxial fatigue model for concrete - strain driven
'''

from scipy.optimize import newton

import matplotlib.pyplot as plt
import numpy as np


def get_stress_strain(eps_arr, sig_pi_bar, K, gamma, E, S, c, r):

    # arrays to store the values
    # nominal stress
    sig_arr = np.zeros_like(eps_arr)
    # sliding stress
    sig_pi_arr = np.zeros_like(eps_arr)
    # damage factor
    w_arr = np.zeros_like(eps_arr)
    # sliding slip
    eps_pi_arr = np.zeros_like(eps_arr)
    # max sliding
    eps_max = np.zeros_like(eps_arr)
    # max stress
    sig_max = np.zeros_like(eps_arr)
    # cumulative sliding
    eps_pi_cum = np.zeros_like(eps_arr)
    diss = np.zeros_like(eps_arr)
    diss_pi = np.zeros_like(eps_arr)
    diss_w = np.zeros_like(eps_arr)

    # material parameters
    E = E
    K = K
    gamma = gamma
    sig_pi_bar = sig_pi_bar
    S = S
    c = c
    r = r

    # state variables
    sig_i = 0
    alpha_i = 0.
    eps_pi_i = 0
    z_i = 0.
    w_i = 0.  # damage
    X_i = gamma * alpha_i
    delta_lamda = 0
    Z = K * z_i
    eps_pi_cum_i = 0
    diss_i = 0
    diss_i_pi = 0
    diss_i_w = 0

    for i in range(1, len(eps_arr)):
        print 'increment', i
        eps_i = eps_arr[i]
        eps_max_i = np.fabs(eps_i)

        sig_i = (1. - w_i) * E * (eps_i - eps_pi_i)

        sig_i_1 = E * (eps_i - eps_pi_i)

        Y_i = 0.5 * E * (eps_i - eps_pi_i) ** 2

        # Threshold
        f_pi_i = np.fabs(sig_i_1 - X_i) - sig_pi_bar - Z

        if f_pi_i > 1e-8:
            # Return mapping
            delta_lamda = f_pi_i / (E / (1. - w_i) + gamma + K)
            # update all the state variables

            ''' 
            f_w_n = lambda w_n :  w_n - w_i - ((1 - w_n) ** c) * (delta_lamda * (Y_i / S) ** r)
            f_w_n2 = lambda w_n : 1 + c * ((1 - w_n) ** (c - 1)) * (delta_lamda * (Y_i / S) ** r)
            w_n = newton(f_w_n, 0., fprime=f_w_n2 , tol=1e-6, maxiter=10) 
            w_i = w_n
            '''
            #w_i = w_i + ((1 - w_i) ** c) * (delta_lamda * (Y_i / S) ** r)

            eps_pi_i = eps_pi_i + delta_lamda * \
                np.sign(sig_i_1 - X_i) / (1. - w_i)

#             delta_eps_pi_i = delta_lamda * \
#                 np.sign(sig_i_1 - X_i) / (1. - w_i)

            Y_i = 0.5 * E * (eps_i - eps_pi_i) ** 2

            delta_w_i = ((1 - w_i) ** c) * (delta_lamda *
                                            (Y_i / S) ** r)

#             w_i = w_i + ((1 - w_i) ** c) * (delta_lamda *
#                                             (Y_i / S) ** r)

            w_i = w_i + ((1.0 - w_i) ** c) * \
                (delta_lamda * (Y_i / S) ** r) * \
                (np.exp(-30.0 * w_i) + np.exp(-20.0 * (1. - w_i)))

            sig_i = E * (1. - w_i) * (eps_i - eps_pi_i)
            #X_i = X_i + gamma * delta_lamda * np.sign(tau_i_1 - X_i)
            alpha_i = alpha_i + delta_lamda * np.sign(sig_i_1 - X_i)
            #delta_alpha_i = delta_lamda * np.sign(sig_i_1 - X_i)

            z_i = z_i + delta_lamda
            #Z_i = K * z_i

            X_i = gamma * alpha_i
            eps_pi_cum_i = eps_pi_cum_i + delta_lamda / (1. - w_i)

            diss_i += (sig_pi_bar) * \
                delta_lamda + Y_i * delta_w_i

            diss_i_pi += (sig_pi_bar) * delta_lamda
            diss_i_w += Y_i * delta_w_i

        sig_max_i = np.fabs(sig_i)
        sig_arr[i] = sig_i
        w_arr[i] = w_i
        eps_pi_arr[i] = eps_pi_i
        eps_max[i] = eps_max_i
        sig_max[i] = sig_max_i
        eps_pi_cum[i] = eps_pi_cum_i
        diss[i] = diss_i
        diss_pi[i] = diss_i_pi
        diss_w[i] = diss_i_w

    return eps_arr, sig_arr, w_arr, eps_pi_arr, eps_max, sig_max, eps_pi_cum, diss, diss_pi, diss_w


if __name__ == '__main__':

    s_levels = np.linspace(0, 0.01, 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= 0
    # s_levels.reshape(-1, 2)[:, 1] = 2
    s_history = s_levels.flatten()

    # slip array as input
    eps_arr_1 = np.hstack([np.linspace(s_history[i], s_history[i + 1], 200)
                           for i in range(len(s_levels) - 1)])

    s_levels_2 = np.linspace(0, 0.01, 2)
    #s_levels_2.reshape(-1, 2)[:, 0] = 0.0002
    s_levels_2.reshape(-1, 2)[:, 1] = 0.002
    s_levels_2[0] = 0
    s_history_2 = s_levels_2.flatten()
#
#     s_history_2 = [0, 0.00169, 0.0003,
#                    0.0021, 0.00035, 0.0032, 0.0004,
#                    0.0038, 0.00042, 0.0055, ]

    # slip array as input
    eps_arr_2 = np.hstack([np.linspace(s_history_2[i], s_history_2[i + 1], 200)
                           for i in range(len(s_history_2) - 1)])

    eps_arr_1, sig_arr_1, w_arr_1, eps_pi_arr_1, eps_max_1, sig_max_1, eps_pi_cum_1, diss_1,  diss_pi_1, diss_w_1 = get_stress_strain(
        eps_arr_1, sig_pi_bar=20, K=100000, gamma=20000, E=40000, S=0.00025, c=2.5, r=1.4)
    eps_arr_2, sig_arr_2, w_arr_2, eps_pi_arr_2, eps_max_2, sig_max_2, eps_pi_cum_2, diss_2,  diss_pi_2, diss_w_2 = get_stress_strain(
        eps_arr_2, sig_pi_bar=20, K=100000, gamma=20000, E=40000, S=0.00025, c=2.5, r=1.4)


#===========================================
# plot slip-tau
#===========================================
    ax1 = plt.subplot(221)

    ax1.plot(eps_arr_1, sig_arr_1, '--k', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax1.plot(eps_arr_2, sig_arr_2, 'k', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)

    ax1.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax1.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('stress-strain')
    plt.xlabel('strain')
    plt.ylabel('Stress(MPa)')
    plt.legend(loc=4)

#============================================
# plot slip-damage
#============================================
    ax2 = plt.subplot(222)
    ax2.plot(eps_arr_1, w_arr_1, '--k', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax2.plot(eps_arr_2, w_arr_2, 'k', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)

    ax2.axhline(y=0, color='k', linewidth=1, alpha=1)
    ax2.axvline(x=0, color='k', linewidth=1, alpha=1)
    plt.title('Damage evolution')
    plt.ylim(0, 1)
    plt.xlabel('Slip(mm)')
    plt.ylabel('Damage')
    plt.legend(loc=4)

#==============================================
# plot slip-dissip
#==============================================
    ax3 = plt.subplot(223)
    ax3.plot(eps_arr_1, diss_1, '--k', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax3.plot(eps_arr_2, diss_2, 'k', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)

    plt.xlabel('Slip(mm)')
    plt.ylabel('Dissipation')
    plt.legend()


# #================================================
# # plot cum_sliding-damage
# #================================================
#     ax4 = plt.subplot(224)
#
#     ax4.plot(eps_pi_cum_1, w_arr_1, '--k', linewidth=1,
#              label='$ \sigma_N = 0$ MPa', alpha=1)
#     ax4.plot(eps_pi_cum_2, w_arr_2, 'k', linewidth=1,
#              label='$ \sigma_N = 0$ MPa', alpha=1)
#
#     plt.xlabel('Cumulative sliding(mm)')
#     plt.ylabel('Damage')
#     plt.ylim(0, 1)
#     plt.legend()


#================================================
# plot cum_sliding-damage
#================================================
    ax4 = plt.subplot(224)

    ax4.plot(eps_arr_1, eps_pi_arr_1, '--k', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax4.plot(eps_arr_2, eps_pi_arr_2, 'k', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)

    plt.xlabel('Cumulative sliding(mm)')
    plt.ylabel('Damage')
    #plt.ylim(0, 1)
    plt.legend()

    plt.show()
