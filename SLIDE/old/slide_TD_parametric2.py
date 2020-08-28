'''
Created on 07.04.2017

@author: abaktheer
'''


#from scipy.optimize import newton

from numpy import \
    array, zeros, dot, trace, \
    tensordot, einsum, zeros_like,\
    identity, sign, linspace, hstack, \
    sqrt, copy, fabs
from numpy.linalg import norm
from scipy.linalg import \
    eigh, inv
from traits.api import \
    Constant,\
    Bool, Enum, Float, HasTraits, \
    Int, Property, cached_property
from traitsui.api import \
    Item, View, Group, Spring, Include

import matplotlib.pyplot as plt

import numpy as np


#=========================================================================
# Implementation of the (SLIDE) - Compression plasticity(CP), and Tensile Damage (TD)
#=========================================================================


def get_CP_TD_Law(eps, E, eps_0, Ad):

    sigma_arr = zeros_like(eps)

    eps_N_p_arr = zeros_like(eps)

    w_N_arr = zeros_like(eps)

    sigma_i = 0.0
    alpha_i = 0.0
    r_i = 0.0
    eps_N_p_i = 0.0

    w_i = 0.0
    z_i = 0.0

    for i in range(1, len(eps)):

        eps_i = eps[i]

        def Z_N(z_i): return 1. / Ad * (-z_i) / (1 + z_i)
        Y_N = 0.5 * E * eps_i ** 2.
        Y_0 = 0.5 * E * eps_0 ** 2
        f = Y_N - (Y_0 + Z_N(z_i))

        if f > 1e-6:
            def f_w(Y): return 1. - 1. / (1. + Ad * (Y - Y_0))

            w_i = f_w(Y_N)
            z_i = - w_i

        sigma_i = (1. - w_i) * E * (eps_i - eps_N_p_i)

        sigma_arr[i] = sigma_i
        w_N_arr[i] = w_i

    return eps, sigma_arr, w_N_arr


if __name__ == '__main__':

    # Check the model behavior
    n = 500
    s_levels = linspace(0,  0.003, 2)
    s_levels[0] = 0
    #s_levels.reshape(-1, 2)[:, 0] *= -1
    s_history_1 = s_levels.flatten()


    s_arr_1 = hstack([linspace(s_history_1[i], s_history_1[i + 1], n)
                      for i in range(len(s_history_1) - 1)])



    s_arr_1, sigma_arr_1, w_N_arr_1 = get_CP_TD_Law(s_arr_1, E=40000, eps_0=0.00001, Ad=4800)
    s_arr_2, sigma_arr_2, w_N_arr_2 = get_CP_TD_Law(s_arr_1, E=40000, eps_0=0.000005, Ad=4800)
    s_arr_3, sigma_arr_3, w_N_arr_3 = get_CP_TD_Law(s_arr_1, E=40000, eps_0=0.0000001, Ad=4800)
    s_arr_4, sigma_arr_4, w_N_arr_4 = get_CP_TD_Law(s_arr_1, E=40000, eps_0=0.00005, Ad=4800)
    
    eps_1 = np.array([0,    4.43E-05,    0.000229705,    0.000412362,    0.000636531,
                          0.000824723,    0.001087638,    0.001295203,    0.001663284,
                              0.001942805,    0.002197417,    0.002515683,    0.002997233])
    sig_1 = np.array([0,    2.008949,    1.668904,    1.378076,    1.096197,    0.8993289,
                          0.6935124,    0.5548099,    0.3803133,    0.2863536,    0.2192395,    0.1610739,    0.1029084])

    ax2 = plt.subplot(111)
    ax2.plot(eps_1, sig_1, color='g', linewidth=2, label='Mars')
    
    ax2.plot(s_arr_1, sigma_arr_1, color='b', linewidth=2, label='eps_0=0.00001')
    ax2.plot(s_arr_2, sigma_arr_2, color='r', linewidth=2, label='eps_0=0.000005')
    ax2.plot(s_arr_3, sigma_arr_3, color='k', linewidth=2, label='eps_0=0.0000001')
    ax2.plot(s_arr_4, sigma_arr_4, color='y', linewidth=2, label='eps_0=0.00005')
    

    ax2.fill_between(s_arr_1, sigma_arr_1, facecolor='gray', alpha=0.2)
#     ax2.fill_between(s_arr_2, sigma_arr_2, facecolor='gray', alpha=0.2)
#     ax2.fill_between(s_arr_3, sigma_arr_3, facecolor='gray', alpha=0.2)
#     ax2.fill_between(s_arr_4, sigma_arr_4, facecolor='gray', alpha=0.2)
#     ax2.fill_between(s_arr_5, sigma_arr_5, facecolor='gray', alpha=0.2)

    ax2.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax2.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('Normal behavior (CP,TD) - monotonic')
    plt.xlabel('strain')
    plt.ylabel('stress(MPa)')
    plt.legend(loc=1)

    plt.show()
