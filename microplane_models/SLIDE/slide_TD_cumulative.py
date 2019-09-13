'''
Created on 07.04.2017

@author: abaktheer
'''


from ibvpy.mats.mats3D.mats3D_eval import MATS3DEval
from ibvpy.mats.mats_eval import \
    IMATSEval
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
    Constant, implements,\
    Bool, Enum, Float, HasTraits, \
    Int, Property, cached_property
from traitsui.api import \
    Item, View, Group, Spring, Include

import matplotlib.pyplot as plt


#===================================================================
# Implementation of the (SLIDE) - Cumulative Sliding Damage (CSD)
#===================================================================
def get_TD_Law(eps, tau_pi_bar, K, gamma, G, S, c, r, a, sigma_n):

    # nominal stress
    tau_arr = zeros_like(eps)

    # damage factor
    w_arr = zeros_like(eps)
    # sliding slip
    xs_pi_arr = zeros_like(eps)

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
        f_pi_i = fabs(tau_i_1 - X_i) - tau_pi_bar - Z + a * sigma_n / 3

        if f_pi_i > 1e-6:
            # Return mapping
            delta_lamda = f_pi_i / (G / (1 - w_i) + gamma + K)
            # update all the state variables

            #w_i = w_i + ((1 - w_i) ** c) * (delta_lamda * (Y_i / S) ** r)

            xs_pi_i = xs_pi_i + delta_lamda * \
                sign(tau_i_1 - X_i) / (1 - w_i)

            Y_i = 0.5 * G * (s_i - xs_pi_i) ** 2

            w_i = w_i + ((1 - w_i) ** c) * (delta_lamda * (Y_i / S) ** r)

            tau_i = G * (1 - w_i) * (s_i - xs_pi_i)
            X_i = X_i + gamma * delta_lamda * sign(tau_i_1 - X_i)
            alpha_i = alpha_i + delta_lamda * sign(tau_i_1 - X_i)
            z_i = z_i + delta_lamda
            xs_pi_cum_i = xs_pi_cum_i + delta_lamda / (1 - w_i)

        tau_arr[i] = tau_i
        w_arr[i] = w_i
        xs_pi_arr[i] = xs_pi_i

    return eps, tau_arr, w_arr, xs_pi_arr

#===================================================================
# Heviside function
#===================================================================


if __name__ == '__main__':

    # Check the model behavior
    n = 100

    s_levels = linspace(0,  0.001, 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= 0
    s_history_1 = s_levels.flatten()

    s_levels = linspace(0,  0.0005, 50)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= 0
    s_history_2 = s_levels.flatten()

    s_arr_1 = hstack([linspace(s_history_1[i], s_history_1[i + 1], 100)
                      for i in range(len(s_history_1) - 1)])

    s_arr_2 = hstack([linspace(s_history_2[i], s_history_2[i + 1], 100)
                      for i in range(len(s_history_2) - 1)])

    s_arr_1, tau_arr_1, w_arr_1, xs_pi_arr_1 = get_TD_Law(
        s_arr_1, tau_pi_bar=3, K=0, gamma=10000, G=30000, S=.000005, c=1.0, r=2.8, a=0.0, sigma_n=0)

    s_arr_2, tau_arr_2, w_arr_2, xs_pi_arr_2 = get_TD_Law(
        s_arr_2, tau_pi_bar=3, K=0, gamma=10000, G=30000, S=.000005, c=1.0, r=2.8, a=0.0, sigma_n=0)

    ax1 = plt.subplot(221)
    ax1.plot(s_arr_1, tau_arr_1, 'r', linewidth=2,
             alpha=1)
    ax1.plot(s_arr_2, tau_arr_2, 'k', linewidth=2,
             alpha=1)

    ax1.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax1.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('Tangential behavior (CSD)')
    plt.xlabel('strain')
    plt.ylabel('stress(MPa)')
    plt.legend(loc=2)

    ax2 = plt.subplot(222)

    plt.plot(s_arr_1, w_arr_1, 'r', linewidth=2,
             label='pressure = 0 MPa', alpha=1.0)
    plt.plot(s_arr_2, w_arr_2, 'k', linewidth=2,
             label='pressure = 10 MPa', alpha=1.0)

    ax2.axhline(y=0, color='r', linewidth=1, alpha=0.5)
    ax2.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('Damage evolution')
    plt.ylim(0, 1)
    plt.xlabel('Slip(mm)')
    plt.ylabel('Damage')
    plt.legend(loc=4)

    plt.show()
