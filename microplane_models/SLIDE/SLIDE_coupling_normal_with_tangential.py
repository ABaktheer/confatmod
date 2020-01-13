'''
Created on 14.11.2016

@author: abaktheer
'''

#from scipy.optimize import newton

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np


def get_bond_slip(s_arr, tau_pi_bar, K, gamma, E_b, S, c, r, m, sigma_n):
    '''for plotting the bond slip fatigue - Initial version modified modified threshold with cumulation-2 implicit
    '''
    # arrays to store the values
    # nominal stress
    tau_arr = np.zeros_like(s_arr)
    # sliding stress
    tau_pi_arr = np.zeros_like(s_arr)
    # damage factor
    w_arr = np.zeros_like(s_arr)
    # sliding slip
    s_pi_arr = np.zeros_like(s_arr)
    # max sliding
    s_max = np.zeros_like(s_arr)
    # max stress
    tau_max = np.zeros_like(s_arr)
    # cumulative sliding
    s_pi_cum = np.zeros_like(s_arr)
    diss = np.zeros_like(s_arr)
    # material parameters
    # shear modulus [MPa]
    E_b = E_b
    # damage - brittleness [MPa^-1]
    K = K
    # Kinematic hardening modulus [MPa]
    gamma = gamma
    # constant in the sliding threshold function
    tau_pi_bar = tau_pi_bar
    # material parameters
    S = S
    c = c
    r = r
    m = m
    sigma_n = sigma_n
    # state variables
    tau_i = 0
    alpha_i = 0.
    s_pi_i = 0
    z_i = 0.
    w_i = 0.  # damage
    X_i = gamma * alpha_i
    delta_lamda = 0
    Z = K * z_i
    s_pi_cum_i = 0
    diss_i = 0

    for i in range(1, len(s_arr)):
        print('increment', i)
        s_i = s_arr[i]
        s_max_i = np.fabs(s_i)

        tau_i = (1 - w_i) * E_b * (s_i - s_pi_i)

        tau_i_1 = E_b * (s_i - s_pi_i)

        Y_i = 0.5 * E_b * (s_i - s_pi_i) ** 2

        # Threshold
        f_pi_i = np.fabs(tau_i_1 - gamma * alpha_i) - \
            tau_pi_bar - Z + m * sigma_n / 3

        if f_pi_i > 1e-6:
            # Return mapping
            delta_lamda = f_pi_i / (E_b / (1. - w_i) + gamma + K)
            # update all the state variables

            s_pi_i = s_pi_i + delta_lamda * \
                np.sign(tau_i_1 - gamma * alpha_i) / (1 - w_i)

            Y_i = 0.5 * E_b * (s_i - s_pi_i) ** 2

            diss_i += (tau_pi_bar - m * sigma_n / 3) * delta_lamda + \
                Y_i * ((1 - w_i) ** c) * (delta_lamda * (Y_i / S) ** r)

            w_i = w_i + ((1 - w_i) ** c) * (delta_lamda *
                                            (Y_i / S) ** r)  * (1 - sigma_n / (0.5 * tau_pi_bar))

            tau_i = E_b * (1 - w_i) * (s_i - s_pi_i)

            alpha_i = alpha_i + delta_lamda * \
                np.sign(tau_i_1 - gamma * alpha_i)
            z_i = z_i + delta_lamda
            s_pi_cum_i = s_pi_cum_i + delta_lamda / (1 - w_i)

        tau_max_i = np.fabs(tau_i)
        tau_arr[i] = tau_i
        w_arr[i] = w_i
        s_pi_arr[i] = s_pi_i
        s_max[i] = s_max_i
        tau_max[i] = tau_max_i
        s_pi_cum[i] = s_pi_cum_i
        diss[i] = diss_i

    return s_arr, tau_arr, w_arr, s_pi_arr, s_max, tau_max, s_pi_cum, diss


if __name__ == '__main__':
    s_levels = np.linspace(0, 8, 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= -1
    # s_levels.reshape(-1, 2)[:, 1] = 2
    s_history = s_levels.flatten()

    # slip array as input
    s_arr_1 = np.hstack([np.linspace(s_history[i], s_history[i + 1], 100)
                         for i in range(len(s_levels) - 1)])



    tau_pi_bar=1
    K=0.1
    gamma=0.2
    E_b=1
    S=0.0001
    c=1
    r=0.001
    m=0
    sigma_n=0


    s_arr_1, tau_arr_1, w_arr_1, s_pi_arr_1, s_max_1, tau_max_1, s_pi_cum_1, diss_1 = get_bond_slip(
        s_arr_1, tau_pi_bar=tau_pi_bar, K=K, gamma=gamma, E_b=E_b, S=S, c=c, r=r, m=m, sigma_n=sigma_n)

    

    ax1 = plt.subplot(221)
    ax1.plot(s_arr_1, tau_arr_1, 'b', linewidth=2,
             label='increments = 10')

    
    ax1.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax1.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('Bond_slip')
    plt.xlabel('Slip(mm)')
    plt.ylabel('Stress(MPa)')
    plt.legend(loc=4)

    ax2 = plt.subplot(222)

    ax2.plot(s_arr_1, w_arr_1, 'b', linewidth=2,
             label='increments = 10')

    ax2.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax2.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('Damage evolution')
    plt.ylim(0, 1)
    plt.xlabel('Slip(mm)')
    plt.ylabel('Damage')
    plt.legend(loc=4)

    #'''
    plt.subplot(223)

    plt.plot(s_pi_cum_1, w_arr_1, 'b', linewidth=2,
             label='increments = 10')


    plt.xlabel('Cumulative sliding(mm)')
    plt.ylabel('Damage')
    plt.ylim(0, 1)
    plt.legend()

    plt.show()
