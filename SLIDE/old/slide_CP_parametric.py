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
def get_CSD_Law(eps, tau_pi_bar, K, gamma, G, S, c, r, a, sigma_n):

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
            #X_i = X_i + gamma * delta_lamda * sign(tau_i_1 - X_i)
            alpha_i = alpha_i + delta_lamda * sign(tau_i_1 - X_i)
            z_i = z_i + delta_lamda
            xs_pi_cum_i = xs_pi_cum_i + delta_lamda / (1.0 - w_i)

        tau_arr[i] = tau_i
        w_arr[i] = w_i
        xs_pi_arr[i] = xs_pi_i

    return eps, tau_arr, w_arr, xs_pi_arr

#===================================================================
# Heviside function
#===================================================================


def get_heviside(eps):
    if eps > 0:
        return 1.0
    else:
        return 0.0

#=========================================================================
# Implementation of the (SLIDE) - Compression plasticity(CP), and Tensile Damage (TD)
#=========================================================================


def get_CP_TD_Law(eps, sigma_0, K, gamma, E, eps_0, Ad, m, a):

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
        H = get_heviside(sigma_i)
        sigma_i = (1 - H * w_i) * E * (eps_i - eps_N_p_i)

        h = max(0., (sigma_0 + K * r_i ** (m + 1.0)))
        f_trial = abs(sigma_i - gamma * alpha_i) - h

        # plasticity yield function
        if f_trial > 1e-6:

            def f_dw_n(delta_lamda_p): return delta_lamda_p - f_trial / \
                (E + abs(K) + gamma + m * gamma * gamma *
                 alpha_i * sign(sigma_i - gamma * alpha_i))
            #             f_dw_n2 = lambda dw_n: 1 + (f_trial * abs((m + 1.0) * K ** m)) /\
            #                 (E + abs((m + 1.0) * delta_lamda_p * K ** m) + gamma)**2.0
            delta_lamda_p = newton(
                f_dw_n, 0., tol=1e-6, maxiter=10)

#             delta_lamda_p = f_trial / \
#                 (E + abs(K) + gamma + m * gamma * gamma *
#                  alpha_i * sign(sigma_i - gamma * alpha_i))
            eps_N_p_i = eps_N_p_i + delta_lamda_p * \
                sign(sigma_i - gamma * alpha_i)
            r_i = r_i + delta_lamda_p
            alpha_i = alpha_i + \
                (delta_lamda_p * sign(sigma_i - gamma * alpha_i))

            alpha_i = alpha_i + \
                (delta_lamda_p / (1.0 + a * gamma * delta_lamda_p)) * \
                (sign(sigma_i - gamma * alpha_i) + a * gamma * alpha_i)

#         Y_0 = 0.5 * E * eps_0**2.0
#         Y_N = 0.5 * H * E * (eps_i - eps_N_p_i)**2.0
#         Z_N = (1.0 / Ad) * (- z_i / (1 + z_i))
#         f_w_trial = Y_N - (Y_0 + Z_N)
#
#         # damage threshold
#         if f_w_trial > 1e-6:
#
#             #             delta_lamda_w = E * eps_i * Ad * \
#             #                 (1 + z_i)**2.0 * (eps[i] - eps[i - 1])
#
#             f_dw_n = lambda dw_n:  dw_n - E * \
#                 (eps_i) * (eps_i - eps[i - 1]) * Ad * (1 + z_i - dw_n) ** 2
#             f_dw_n2 = lambda dw_n: 1 + 2 * E * \
#                 (eps_i) * (eps_i - eps[i - 1]) * Ad * (1 + z_i - dw_n)
#             dw_n = newton(f_dw_n, 0., fprime=f_dw_n2, tol=1e-6, maxiter=50)
#
#             w_i = w_i + dw_n
#             z_i = z_i - dw_n

        def Z_N(z_N): return 1. / Ad * (-z_i) / (1 + z_i)
        Y_N = 0.5 * H * E * eps_i ** 2.
        Y_0 = 0.5 * E * eps_0 ** 2
        f = Y_N - (Y_0 + Z_N(z_i))

        if f > 1e-6:
            def f_w(Y): return 1. - 1. / (1. + Ad * (Y - Y_0))

            w_i = f_w(Y_N)
            z_i = - w_i

        sigma_i = (1. - H * w_i) * E * (eps_i - eps_N_p_i)

        sigma_arr[i] = sigma_i
        eps_N_p_arr[i] = eps_N_p_i
        w_N_arr[i] = w_i

    return eps, sigma_arr, eps_N_p_arr, w_N_arr


if __name__ == '__main__':

    # Check the model behavior
    n = 500
    s_levels = linspace(0,  -0.005, 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= 0
    s_history_1 = s_levels.flatten()

    s_levels = linspace(0,  -0.005, 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= -1
    s_history_2 = s_levels.flatten()

    s_levels = linspace(0,  -0.005, 2)
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

    s_arr_1, sigma_arr_1, eps_N_p_arr_1, w_N_arr_1 = get_CP_TD_Law(
        s_arr_1, sigma_0=20, K=0.0, gamma=0.0, E=34000., eps_0=1.0e-4, Ad=20000, m=0.0, a=0.0)
    s_arr_2, sigma_arr_2, eps_N_p_arr_2, w_N_arr_2 = get_CP_TD_Law(
        s_arr_2, sigma_0=20, K=1000.0, gamma=10000.0, E=34000., eps_0=1.0e-4, Ad=10000, m=0.0, a=-0.05)
    s_arr_3, sigma_arr_3, eps_N_p_arr_3, w_N_arr_3 = get_CP_TD_Law(
        s_arr_3, sigma_0=20, K=1000.0, gamma=10000.0, E=34000., eps_0=1.0e-4, Ad=5000, m=0.0, a=-0.1)
    s_arr_4, sigma_arr_4, eps_N_p_arr_4, w_N_arr_4 = get_CP_TD_Law(
        s_arr_2, sigma_0=20, K=1000.0, gamma=10000.0, E=34000., eps_0=1.0e-4, Ad=5000, m=0.0, a=-0.2)
    s_arr_5, sigma_arr_5, eps_N_p_arr_5, w_N_arr_5 = get_CP_TD_Law(
        s_arr_2, sigma_0=20, K=1000.0, gamma=10000.0, E=34000., eps_0=1.0e-4, Ad=5000, m=0.0, a=0.0)

    ax2 = plt.subplot(111)
    ax2.plot(-s_arr_1, -sigma_arr_1, color='k', label='no hardening')
    ax2.plot(-s_arr_2, -sigma_arr_2, color='r', label='a=-0.05')
    ax2.plot(-s_arr_3, -sigma_arr_3, color='g', label='a=-0.1')
    ax2.plot(-s_arr_4, -sigma_arr_4, color='b', label='a=-0.2')
    ax2.plot(-s_arr_5, -sigma_arr_5, color='y', label='a=0.0')

    ax2.fill_between(-s_arr_1, -sigma_arr_1, facecolor='gray', alpha=0.2)
    ax2.fill_between(-s_arr_2, -sigma_arr_2, facecolor='gray', alpha=0.2)
    ax2.fill_between(-s_arr_3, -sigma_arr_3, facecolor='gray', alpha=0.2)
    ax2.fill_between(-s_arr_4, -sigma_arr_4, facecolor='gray', alpha=0.2)
    ax2.fill_between(-s_arr_5, -sigma_arr_5, facecolor='gray', alpha=0.2)

    ax2.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax2.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('Normal behavior (CP,TD) - monotonic')
    plt.xlabel('strain')
    plt.ylabel('stress(MPa)')
    plt.legend(loc=2)

    plt.show()
