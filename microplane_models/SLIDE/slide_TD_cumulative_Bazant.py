'''
Created on 07.04.2017

@author: abaktheer
'''

from numpy import \
    zeros_like,\
    sign, linspace, hstack, \
    fabs

import matplotlib.pyplot as plt
import numpy as np

#===================================================
# Cumulative damage function [Kirane & Bazant 2015]
#===================================================


def get_TD_cum_bazant(eps, A, E, r, p):

    q = p + 1.0
    s = r * q

    sigma_arr = zeros_like(eps)

    eps_N_p_arr = zeros_like(eps)

    w_N_arr = zeros_like(eps)
    kappa_arr = zeros_like(eps)

    sigma_i = 0.0
    w_N_i = 0.0
    kappa_i = 0.0
    kappa_cum = 0.0

    for i in range(1, len(eps)):

        eps_i = eps[i]

        #kappa_i = np.abs(eps_i)

        kappa_i = max(eps_i, kappa_i)
        kappa_cum += kappa_i
        w_N_i = 1. - (np.exp(-A * kappa_i)) * (1.0 / (1.0 + r *
                                                      (kappa_cum ** q) + (r**2) * (kappa_cum ** (2 * q))))
        sigma_i = (1. - w_N_i) * E * eps_i

        sigma_arr[i] = sigma_i
        kappa_arr[i] = kappa_i
        w_N_arr[i] = w_N_i

    return eps, sigma_arr, w_N_arr, kappa_arr


if __name__ == '__main__':

    # Check the model behavior
    n = 100
    s_levels = linspace(0,  0.002, 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= 0
    s_history_1 = s_levels.flatten()

    s_levels = linspace(0,  0.002, 100)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] = 0.000005
    s_levels.reshape(-1, 2)[:, 1] = 0.00015
    s_levels[0] = 0
    s_history_2 = s_levels.flatten()

    s_arr_1 = hstack([linspace(s_history_1[i], s_history_1[i + 1], n)
                      for i in range(len(s_history_1) - 1)])
    s_arr_2 = hstack([linspace(s_history_2[i], s_history_2[i + 1], n)
                      for i in range(len(s_history_2) - 1)])

    s_arr_1, sigma_arr_1,  w_N_arr_1, kappa_arr_1 = get_TD_cum_bazant(
        s_arr_1, A=4500., E=30000., r=1., p=1.0)
    s_arr_2, sigma_arr_2,  w_N_arr_2, kappa_arr_1 = get_TD_cum_bazant(
        s_arr_2, A=4500., E=30000., r=1., p=1.0)

    ax2 = plt.subplot(221)
    ax2.plot(s_arr_1, sigma_arr_1, color='k', label='mon.')
    ax2.plot(s_arr_2, sigma_arr_2, color='r', label='cyc.')

    ax2.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax2.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('Normal behavior (CP,TD) - monotonic')
    plt.xlabel('strain')
    plt.ylabel('stress(MPa)')
    plt.legend(loc=2)

    ax2 = plt.subplot(222)
    ax2.plot(s_arr_1, w_N_arr_1, color='k', label='mon.')
    ax2.plot(s_arr_2, w_N_arr_2, color='r', label='cyc.')

    ax2.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax2.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('Normal behavior (CP,TD) - monotonic')
    plt.xlabel('strain')
    plt.ylabel('damage')
    plt.legend(loc=4)

    plt.show()
