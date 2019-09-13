'''
Created on 14.11.2016

@author: abaktheer

Uniaxial fatigue model for concrete - strain driven
'''

from scipy.optimize import newton

import matplotlib.pyplot as plt
import numpy as np


def get_stress_strain(eps_arr, sig_0, E, K, gamma):

    # arrays to store the values
    # nominal stress
    sig_arr = np.zeros_like(eps_arr)

    # sliding slip
    eps_p_arr = np.zeros_like(eps_arr)

    # material parameters
    E = E
    K = K
    gamma = gamma
    sig_0 = sig_0

    # state variables
    sig_i = 0
    alpha_i = 0.
    eps_p_i = 0
    z_i = 0.

    X_i = gamma * alpha_i
    delta_lamda = 0
    Z = K * z_i

    for i in range(1, len(eps_arr)):
        print 'increment', i
        eps_i = eps_arr[i]

        sig_i = E * (eps_i - eps_p_i)

        # Threshold
        f_pi_i = np.fabs(sig_i - X_i) - sig_0 - Z

        if f_pi_i > 1e-8:
            # Return mapping
            delta_lamda = f_pi_i / (E + gamma + K)
            # update all the state variables

            eps_p_i = eps_p_i + delta_lamda * \
                np.sign(sig_i - X_i)

            sig_i = E * (eps_i - eps_p_i)
            alpha_i = alpha_i + delta_lamda * np.sign(sig_i - X_i)

            z_i = z_i + delta_lamda

            X_i = gamma * alpha_i

        sig_arr[i] = sig_i

        eps_p_arr[i] = eps_p_i

    return eps_arr, sig_arr, eps_p_arr


if __name__ == '__main__':

    s_levels = np.linspace(0, 0.01, 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= 0
    # s_levels.reshape(-1, 2)[:, 1] = 2
    s_history = s_levels.flatten()

    # slip array as input
    eps_arr_1 = np.hstack([np.linspace(s_history[i], s_history[i + 1], 2000)
                           for i in range(len(s_levels) - 1)])

    s_levels_2 = np.linspace(0, 0.01, 10)
    s_levels_2.reshape(-1, 2)[:, 0] = 0.0002
    s_levels_2.reshape(-1, 2)[:, 1] = 0.002
    s_levels_2[0] = 0
    s_history_2 = s_levels_2.flatten()

    # slip array as input
    eps_arr_2 = np.hstack([np.linspace(s_history_2[i], s_history_2[i + 1], 2000)
                           for i in range(len(s_history_2) - 1)])

    eps_arr_1, sig_arr_1, eps_p_arr_1 = get_stress_strain(
        eps_arr_1, sig_0=10, K=10000, gamma=2000, E=30000)
    eps_arr_2, sig_arr_2, eps_p_arr_2 = get_stress_strain(
        eps_arr_2, sig_0=10, K=10000, gamma=2000, E=30000)


#===========================================
# plot slip-tau
#===========================================
    ax1 = plt.subplot(211)

    ax1.plot(eps_arr_1, sig_arr_1, '--k', linewidth=2,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax1.plot(eps_arr_2, sig_arr_2, 'k', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)

    ax1.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax1.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('stress-strain')
    plt.xlabel('strain')
    plt.ylabel('Stress(MPa)')
    plt.legend(loc=4)


#================================================
# plot cum_sliding-damage
#================================================
    ax4 = plt.subplot(212)

    ax4.plot(eps_arr_1, eps_p_arr_1, '--k', linewidth=2,
             label='$ \sigma_N = 0$ MPa', alpha=1)
    ax4.plot(eps_arr_2, eps_p_arr_2, 'k', linewidth=1,
             label='$ \sigma_N = 0$ MPa', alpha=1)

    plt.xlabel('Cumulative sliding(mm)')
    plt.ylabel('Damage')
    #plt.ylim(0, 1)
    plt.legend()

    plt.show()
