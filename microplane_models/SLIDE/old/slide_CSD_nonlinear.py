'''
Created on 07.04.2017

@author: abaktheer
'''


from numpy.linalg import norm

from ibvpy.mats.mats3D.mats3D_eval import MATS3DEval
from ibvpy.mats.mats_eval import \
    IMATSEval
from numpy import \
    array, zeros, dot, trace, \
    tensordot, einsum, zeros_like,\
    identity, sign, linspace, hstack, \
    sqrt, copy, fabs
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
def get_CSD_Law(eps, tau_pi_bar, K, gamma, G, S, c, r, a, sigma_n, m):

    # nominal stress
    tau_arr = zeros_like(eps)

    # damage factor
    w_arr = zeros_like(eps)
    # sliding slip
    xs_pi_arr = zeros_like(eps)

    # state variables
    tau_i = 0
    alpha_i = 0.
    xs_pi_i = 0
    z_i = 0.
    w_i = 0.
    X_i = gamma * alpha_i
    delta_lamda = 0.
    Z = K * z_i
    xs_pi_cum_i = 0.

    for i in range(1, len(eps)):

        s_i = eps[i]

        tau_i = (1.0 - w_i) * G * (s_i - xs_pi_i)

        tau_i_1 = G * (s_i - xs_pi_i)

        Y_i = 0.5 * G * (s_i - xs_pi_i) ** 2

        # Threshold
        f_pi_i = fabs(tau_i_1 - X_i) - tau_pi_bar - Z + a * sigma_n

        if f_pi_i > 1e-6:
            # Return mapping
            delta_lamda = f_pi_i / \
                (G / (1.0 - w_i) + gamma +
                 K)
            # update all the state variables
            #w_i = w_i + ((1.0 - w_i) ** c) * (delta_lamda * (Y_i / S) ** r)
            xs_pi_i = xs_pi_i + delta_lamda * sign(tau_i_1 - X_i) / (1.0 - w_i)

            z_i = z_i + delta_lamda

            alpha_i = alpha_i + delta_lamda * sign(tau_i_1 - X_i)

            Y_i = 0.5 * G * (s_i - xs_pi_i) ** 2.0
            w_i = w_i + ((1.0 - w_i) ** c) * (delta_lamda * (Y_i / S) ** r)
            tau_i = G * (1.0 - w_i) * (s_i - xs_pi_i)

            xs_pi_cum_i = xs_pi_cum_i + delta_lamda / (1.0 - w_i)

        tau_arr[i] = tau_i
        w_arr[i] = w_i
        xs_pi_arr[i] = xs_pi_i

    return eps, tau_arr, w_arr, xs_pi_arr


if __name__ == '__main__':

    # Check the model behavior
    n = 1000
    s_levels = linspace(0,  0.01, 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= -1
    s_history_1 = s_levels.flatten()

    s_levels = linspace(0,  0.01, 20)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= -1
    s_history_2 = s_levels.flatten()

    s_levels = linspace(0,  -0.01, 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= -1
    s_history_3 = s_levels.flatten()
    #s_history_3 = [0, 0.0002, -0.001, -0.00001]

    s_arr_1 = hstack([linspace(s_history_1[i], s_history_1[i + 1], n)
                      for i in range(len(s_history_1) - 1)])
    s_arr_2 = hstack([linspace(s_history_2[i], s_history_2[i + 1], n)
                      for i in range(len(s_history_2) - 1)])
    s_arr_3 = hstack([linspace(s_history_3[i], s_history_3[i + 1], n)
                      for i in range(len(s_history_3) - 1)])

    s_arr_1, tau_arr_1, w_arr_1, xs_pi_arr_1 = get_CSD_Law(
        s_arr_1, tau_pi_bar=1., K=0, gamma=1000000., G=3000., S=0.00001, c=1.2, r=1.2, a=0.0, sigma_n=0.0, m=0.0)

    s_arr_2, tau_arr_2, w_arr_2, xs_pi_arr_2 = get_CSD_Law(
        s_arr_2, tau_pi_bar=1., K=0, gamma=1000000., G=3000., S=0.00001, c=1.2, r=1.2, a=0.0, sigma_n=0.0, m=0.0)

    s_arr_3, tau_arr_3, w_arr_3, xs_pi_arr_3 = get_CSD_Law(
        s_arr_3, tau_pi_bar=1., K=0, gamma=1000000., G=3000., S=0.00001, c=1.2, r=1.2, a=0.0, sigma_n=0.0, m=0.0)

    ax2 = plt.subplot(221)
    ax2.plot(s_arr_1, tau_arr_1, color='k', label='a')
    ax2.plot(s_arr_2, tau_arr_2, color='r', label='b')
    ax2.plot(s_arr_3, tau_arr_3, color='g')
    ax2.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax2.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('Normal behavior (CP,TD) - monotonic')
    plt.xlabel('strain')
    plt.ylabel('stress(MPa)')
    plt.legend(loc=2)

    ax2 = plt.subplot(222)
    ax2.plot(s_arr_1, w_arr_1, color='k', label='a')
    ax2.plot(s_arr_2, w_arr_2, color='r', label='b')
    ax2.plot(s_arr_3, w_arr_3, color='g')
    ax2.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax2.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('Normal behavior (CP,TD) - monotonic')
    plt.xlabel('strain')
    plt.ylabel('stress(MPa)')
    plt.legend(loc=2)

    plt.show()
