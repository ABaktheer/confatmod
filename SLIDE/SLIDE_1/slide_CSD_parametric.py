'''
Created on 07.04.2017

@author: abaktheer
'''


from numpy import \
    array, zeros, dot, trace, \
    tensordot, einsum, zeros_like,\
    identity, sign, linspace, hstack, \
    sqrt, copy, fabs
from numpy.linalg import norm
from scipy.linalg import \
    eigh, inv
from scipy.optimize import newton
from traits.api import \
    Constant, \
    Bool, Enum, Float, HasTraits, \
    Int, Property, cached_property
from traitsui.api import \
    Item, View, Group, Spring, Include

import matplotlib.pyplot as plt


#===================================================================
# Implementation of the (SLIDE) - Cumulative Sliding Damage (CSD)
#===================================================================
def get_CSD_Law(eps, tau_pi_bar, K, gamma, G, S, c, r, m, sigma_n, p):

    # nominal stress
    tau_arr = zeros_like(eps)

    # damage factor
    w_arr = zeros_like(eps)
    # sliding slip
    xs_pi_arr = zeros_like(eps)
    xs_pi_cum_arr = zeros_like(eps)

    # material parameters
    # shear modulus [MPa]
#     G = G
#     # damage - brittleness [MPa^-1]
#     K = K
#     # Kinematic hardening modulus [MPa]
#     gamma = gamma
#     # constant in the sliding threshold function
#     tau_pi_bar = tau_pi_bar
#     # material parameters
#     S = S
#     c = c
#     r = r
#     a = a
#     sigma_n = sigma_n

    # state variables
    tau_i = 0
    alpha_i = 0.
    xs_pi_i = 0
    z_i = 0.
    w_i = 0.  # damage
    X_i = gamma * alpha_i
    delta_lamda = 0
    Z = K * z_i
    xs_pi_cum_i = 0

    for i in range(1, len(eps)):

        s_i = eps[i]

        tau_i = (1 - w_i) * G * (s_i - xs_pi_i)

        tau_i_1 = G * (s_i - xs_pi_i)

        Y_i = 0.5 * G * (s_i - xs_pi_i) ** 2

        # Threshold
        f_pi_i = fabs(tau_i_1 - X_i) - tau_pi_bar - Z + m * sigma_n / 3

        if f_pi_i > 1e-6:
            # Return mapping
            delta_lamda = f_pi_i / (G / (1 - w_i) + gamma + K)
            # update all the state variables

            #w_i = w_i + ((1 - w_i) ** c) * (delta_lamda * (Y_i / S) ** r)

            xs_pi_i = xs_pi_i + delta_lamda * \
                sign(tau_i_1 - X_i) / (1 - w_i)

            Y_i = 0.5 * G * (s_i - xs_pi_i) ** 2

            w_i = w_i + ((1 - w_i) ** c) * (delta_lamda * (Y_i / S) ** r)  *(tau_pi_bar / (tau_pi_bar + m * sigma_n))**p

            if w_i >= 0.99:
                break

            tau_i = G * (1 - w_i) * (s_i - xs_pi_i)
            X_i = X_i + gamma * delta_lamda * sign(tau_i_1 - X_i)
            alpha_i = alpha_i + delta_lamda * sign(tau_i_1 - X_i)
            z_i = z_i + delta_lamda
            xs_pi_cum_i = xs_pi_cum_i + delta_lamda / (1 - w_i)

        tau_arr[i] = tau_i
        w_arr[i] = w_i
        xs_pi_arr[i] = xs_pi_i
        xs_pi_cum_arr[i] = xs_pi_cum_i

    return eps, tau_arr, w_arr, xs_pi_arr, xs_pi_cum_arr



if __name__ == '__main__':

    # Check the model behavior
    n = 1000

    s_levels = linspace(0,  0.003, 500)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= 0
    s_history_1 = s_levels.flatten()

    s_levels = linspace(0,  0.003, 500)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= 0
    s_history_2 = s_levels.flatten()

    s_levels = linspace(0,  -0.005, 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= -1
    s_history_3 = s_levels.flatten()

    s_arr_1 = hstack([linspace(s_history_1[i], s_history_1[i + 1], 100)
                      for i in range(len(s_history_1) - 1)])

    s_arr_2 = hstack([linspace(s_history_2[i], s_history_2[i + 1], 100)
                      for i in range(len(s_history_2) - 1)])

    s_arr_3 = hstack([linspace(s_history_3[i], s_history_3[i + 1], 1000)
                      for i in range(len(s_history_3) - 1)])

    s_arr_1, tau_arr_1, w_arr_1, xs_pi_arr_1, xs_pi_cum_arr_1 = get_CSD_Law(
        s_arr_1, tau_pi_bar=4, K=50, gamma=10000, G=14166, S=.00004, c=-1.0, r=1.0, m=0.0, sigma_n=0, p=0)

    s_arr_2, tau_arr_2, w_arr_2, xs_pi_arr_2, xs_pi_cum_arr_2 = get_CSD_Law(
        s_arr_2, tau_pi_bar=4, K=50, gamma=10000, G=14166, S=.00004, c=1.0, r=1.0, m=0.0, sigma_n=0, p=0)
#     s_arr_3, tau_arr_3, w_arr_3, xs_pi_arr_3, xs_pi_cum_arr_3 = get_CSD_Law(
#         s_arr_2, tau_pi_bar=4, K=50, gamma=10000, G=14166, S=.00004, c=2.0, r=1.0, a=0.0, sigma_n=0)
#     s_arr_4, tau_arr_4, w_arr_4, xs_pi_arr_4, xs_pi_cum_arr_4 = get_CSD_Law(
#         s_arr_2, tau_pi_bar=4, K=50, gamma=10000, G=14166, S=.00004, c=3.0, r=1.0, a=0.0, sigma_n=0)
#     s_arr_5, tau_arr_5, w_arr_5, xs_pi_arr_5, xs_pi_cum_arr_5 = get_CSD_Law(
#         s_arr_2, tau_pi_bar=4, K=50, gamma=10000, G=14166, S=.00004, c=4.0, r=1.0, a=0.0, sigma_n=0)
#     s_arr_6, tau_arr_6, w_arr_6, xs_pi_arr_6, xs_pi_cum_arr_6 = get_CSD_Law(
#         s_arr_2, tau_pi_bar=4, K=50, gamma=10000, G=14166, S=.00004, c=5.0, r=1.0, a=0.0, sigma_n=0)
    #

#     s_arr_8, tau_arr_8, w_arr_8, xs_pi_arr_8 = get_CSD_Law(
# s_arr_8, tau_pi_bar=4, K=50, gamma=10000, G=14166, S=.00004, c=1.0,
# r=1.2, a=0.0, sigma_n=0)

    ax1 = plt.subplot(221)
    ax1.plot(s_arr_1, tau_arr_1, 'k')
    ax1.plot(s_arr_2, tau_arr_2, 'r')

    plt.title('Tangential behavior (CSD)')
    plt.xlabel('strain')
    plt.ylabel('stress(MPa)')
    plt.legend(loc=2)

    ax2 = plt.subplot(222)

    plt.plot(s_arr_1, w_arr_1, 'k')
    plt.plot(s_arr_2, w_arr_2, 'r')

    plt.title('Damage evolution')
    plt.ylim(0, 1)
    plt.xlabel('Slip(mm)')
    plt.ylabel('Damage')
    plt.legend(loc=4)

    ax2 = plt.subplot(223)

    plt.plot(xs_pi_cum_arr_1[:-2], w_arr_1[:-2], 'k')
    plt.plot(xs_pi_cum_arr_2[:-2], w_arr_2[:-2], 'r')
#     plt.plot(xs_pi_cum_arr_3[:-2], w_arr_3[:-2], 'r')
#     plt.plot(xs_pi_cum_arr_4[:-2], w_arr_4[:-2], 'r')
#     plt.plot(xs_pi_cum_arr_5[:-2], w_arr_5[:-2], 'r')
#     plt.plot(xs_pi_cum_arr_6[:-2], w_arr_6[:-2], 'r')

    ax2.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax2.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('Damage evolution')
    plt.ylim(0, 1)
    plt.xlabel('Slip(mm)')
    plt.ylabel('Damage')
    plt.legend(loc=4)

    plt.show()
