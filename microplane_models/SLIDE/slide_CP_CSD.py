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


def get_heviside(eps):
    if eps > 0:
        return 1.0
    else:
        return 0.0

#=========================================================================
# Implementation of the (SLIDE) - Compression plasticity(CP), and Tensile Damage (TD)
#=========================================================================


def get_CP_TD_Law(eps, sigma_0, K, gamma, E, eps_0, Ad):

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

        h = max(0., (sigma_0 + K * r_i))
        f_trial = abs(sigma_i - gamma * alpha_i) - h

        # plasticity yield function
        if f_trial > 1e-6:

            delta_lamda_p = f_trial / (E + abs(K) + gamma)
            eps_N_p_i = eps_N_p_i + delta_lamda_p * \
                sign(sigma_i - gamma * alpha_i)
            r_i = r_i + delta_lamda_p
            alpha_i = alpha_i + delta_lamda_p * \
                gamma * sign(sigma_i - gamma * alpha_i)

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
        Y_N = 0.5 * H * E * eps_i ** 2
        Y_0 = 0.5 * E * eps_0 ** 2
        f = Y_N - (Y_0 + Z_N(z_i))

        if f > 1e-6:
            def f_w(Y): return 1 - 1. / (1 + Ad * (Y - Y_0))

            w_i = f_w(Y_N)
            z_i = - w_i

        sigma_i = (1 - H * w_i) * E * (eps_i - eps_N_p_i)

        sigma_arr[i] = sigma_i
        eps_N_p_arr[i] = eps_N_p_i
        w_N_arr[i] = w_i

    return eps, sigma_arr, eps_N_p_arr, w_N_arr


def get_CP_Law(eps, sigma_0, K, gamma, E):

    sigma_arr = zeros_like(eps)

    eps_N_p_arr = zeros_like(eps)

    sigma_i = 0
    alpha_i = 0.
    r_i = 0.
    eps_N_p_i = 0

    for i in range(0, len(eps)):

        eps_i = eps[i]

        sigma_i = E * (eps_i - eps_N_p_i)

        h = max(0., (sigma_0 + K * r_i))
        f_trial = abs(sigma_i - gamma * alpha_i) - h

        if f_trial > 1e-6:

            delta_lamda = f_trial / (E + abs(K) + gamma)
            eps_N_p_i = eps_N_p_i + delta_lamda * \
                sign(sigma_i - gamma * alpha_i)
            r_i = r_i + delta_lamda
            alpha_i = alpha_i + delta_lamda * \
                gamma * sign(sigma_i - gamma * alpha_i)

        sigma_i = E * (eps_i - eps_N_p_i)

        sigma_arr[i] = sigma_i
        eps_N_p_arr[i] = eps_N_p_i

    return eps, sigma_arr, eps_N_p_arr


if __name__ == '__main__':

    # Check the model behavior
    n = 1000
    s_levels = linspace(0, - 0.005, 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= -1
    s_history_1 = s_levels.flatten()

    s_levels = linspace(0,  0.005, 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= -1
    s_history_2 = s_levels.flatten()

    s_levels = linspace(0,  0.0010, 20)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= -1
    s_history_3 = s_levels.flatten()
    #s_history_3 = [0, 0.0002, -0.001, -0.00001, -0.001]
    # cyclic loading
    s_history_4 = [-0, -0.01, -0.007, -0.015, -0.011, -0.023, -0.017, -0.031, -.023, -
                   0.04, -0.032, -0.05, -0.041, -0.06, -0.050, -.07, -.058, -.08, -0.066, -.09, -0.075, - 0.1]
#
    s_history_5 = [0, 0.01, 0.007, 0.015, 0.011, 0.023, 0.017, 0.031, .023,
                   0.04, 0.032, 0.05, 0.041, 0.06, 0.050, .07, .058, .08, 0.066, .09,
                   0.075,  0.1]

    s_levels = linspace(0,  0.005, 50)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= -1
    s_history_6 = s_levels.flatten()

    s_levels = linspace(0,  0.005, 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= -1
    s_history_7 = s_levels.flatten()

    s_levels = linspace(0,  -0.005, 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= -1
    s_history_8 = s_levels.flatten()

    s_arr_1 = hstack([linspace(s_history_1[i], s_history_1[i + 1], n)
                      for i in range(len(s_history_1) - 1)])
    s_arr_2 = hstack([linspace(s_history_2[i], s_history_2[i + 1], n)
                      for i in range(len(s_history_2) - 1)])
    s_arr_3 = hstack([linspace(s_history_3[i], s_history_3[i + 1], n)
                      for i in range(len(s_history_3) - 1)])

    s_arr_4 = hstack([linspace(s_history_4[i], s_history_4[i + 1], 50)
                      for i in range(len(s_history_4) - 1)])
    s_arr_5 = hstack([linspace(s_history_5[i], s_history_5[i + 1], 50)
                      for i in range(len(s_history_5) - 1)])

    s_arr_6 = hstack([linspace(s_history_6[i], s_history_6[i + 1], 50)
                      for i in range(len(s_history_6) - 1)])

    s_arr_7 = hstack([linspace(s_history_7[i], s_history_7[i + 1], 50)
                      for i in range(len(s_history_7) - 1)])

    s_arr_8 = hstack([linspace(s_history_8[i], s_history_8[i + 1], 50)
                      for i in range(len(s_history_8) - 1)])


# #
#     s_arr_4, tau_arr_4, w_arr_4, xs_pi_arr_4 = get_CSD_Law(
#         s_arr_4, tau_pi_bar=5, K=1000, gamma=10000, G=14166, S=.00003, c=1.0, r=1.2, a=0.0, sigma_n=0)
# #
#     s_arr_5, tau_arr_5, w_arr_5, xs_pi_arr_5 = get_CSD_Law(
# s_arr_5, tau_pi_bar=5, K=1000, gamma=10000, G=14166, S=.00003, c=1.0,
# r=1.2, a=0.0, sigma_n=0)

    s_arr_6, tau_arr_6, w_arr_6, xs_pi_arr_6 = get_CSD_Law(
        s_arr_6, tau_pi_bar=4, K=100, gamma=10000, G=14166, S=.00003, c=1.0, r=1.2, a=0.0, sigma_n=0)

    s_arr_7, tau_arr_7, w_arr_7, xs_pi_arr_7 = get_CSD_Law(
        s_arr_7, tau_pi_bar=4, K=100, gamma=10000, G=14166, S=.00003, c=1.0, r=1.2, a=0.0, sigma_n=0)

    s_arr_8, tau_arr_8, w_arr_8, xs_pi_arr_8 = get_CSD_Law(
        s_arr_8, tau_pi_bar=4, K=100, gamma=10000, G=14166, S=.00003, c=1.0, r=1.2, a=0.0, sigma_n=0)

    s_arr_1, sigma_arr_1, eps_N_p_arr_1, w_N_arr_1 = get_CP_TD_Law(
        s_arr_1, sigma_0=20, K=100., gamma=100., E=34000, eps_0=1.0e-4, Ad=5000)
    s_arr_2, sigma_arr_2, eps_N_p_arr_2, w_N_arr_2 = get_CP_TD_Law(
        s_arr_2, sigma_0=20, K=100., gamma=100., E=34000, eps_0=1.0e-4, Ad=5000)
    s_arr_3, sigma_arr_3, eps_N_p_arr_3, w_N_arr_3 = get_CP_TD_Law(
        s_arr_3, sigma_0=20, K=100., gamma=100., E=34000, eps_0=1.0e-4, Ad=5000)

    # print'w_N_arr', w_N_arr
    # print 'Max_slip', np.amax(s_max)
    # print 'Max_stress', np.amax(tau_max)

    ax1 = plt.subplot(221)

#     ax1.plot(s_arr_4, tau_arr_4, 'k', linewidth=1,
#              alpha=1)
#     ax1.plot(s_arr_5, tau_arr_5, 'k', linewidth=1,
#              alpha=1)
    ax1.plot(s_arr_6, tau_arr_6, 'k', linewidth=1,
             alpha=1)
    ax1.plot(s_arr_7, tau_arr_7, 'r', linewidth=1,
             alpha=1)
    ax1.plot(s_arr_8, tau_arr_8, 'r', linewidth=1,
             alpha=1)

    ax1.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax1.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('Tangential behavior (CSD)')
    plt.xlabel('strain')
    plt.ylabel('stress(MPa)')
    plt.legend(loc=2)

#
#     ax2 = plt.subplot(222)
#     ax2.plot(s_arr_1, sigma_arr_1, 'r', linewidth=1,
#              label='Compression plasticity (CP)', alpha=1)
#     ax2.plot(s_arr_2, sigma_arr_2, 'w', linewidth=1,
#              alpha=1)
#
#     ax2.axhline(y=0, color='k', linewidth=1, alpha=0.5)
#     ax2.axvline(x=0, color='k', linewidth=1, alpha=0.5)
#     plt.title('Normal behavior (CP,TD..)')
#     plt.xlabel('strain')
#     plt.ylabel('stress(MPa)')
#     plt.legend(loc=2)

#     ax2 = plt.subplot(221)
#     ax2.plot(s_arr_1, sigma_arr_1, color='r')
#     ax2.plot(s_arr_2, sigma_arr_2, color='k')
#     ax2.axhline(y=0, color='k', linewidth=1, alpha=0.5)
#     ax2.axvline(x=0, color='k', linewidth=1, alpha=0.5)
#     plt.title('Normal behavior (CP,TD) - monotonic')
#     plt.xlabel('strain')
#     plt.ylabel('stress(MPa)')
#     plt.legend(loc=2)
#
#     ax2 = plt.subplot(222)
#     ax2.plot(s_arr_3, sigma_arr_3, color='r')
#     ax2.axhline(y=0, color='k', linewidth=1, alpha=0.5)
#     ax2.axvline(x=0, color='k', linewidth=1, alpha=0.5)
#     plt.title('Normal behavior (CP,TD) - cyclic')
#     plt.xlabel('strain')
#     plt.ylabel('stress(MPa)')
#     plt.legend(loc=2)

    # print s_arr_2

#     ax2 = plt.subplot(222)
#     # plt.plot(s_arr_2, tau_arr_2, '--')#, label='c = 1')
#     ax2.plot(s_arr_1, w_arr_1, '--r', linewidth=1,
#              label='monotonic loading', alpha=1)
#
#     ax2.plot(s_arr_2, w_arr_2, '--r', linewidth=1,
#              alpha=1)
#
#     ax2.plot(s_arr_3, w_arr_3, 'k', linewidth=1,
#              label='cyclic loading', alpha=1)
#
#     ax2.plot(s_arr_4, w_arr_4, 'k', linewidth=1,
#              alpha=1)
# #     ax1.plot(s_arr_1, sigma_arr_2, 'r', linewidth=1,
# #              label='$ \sigma_N = 10$ MPa', alpha=0.5)
#
#     ax2.axhline(y=0, color='k', linewidth=1, alpha=0.5)
#     ax2.axvline(x=0, color='k', linewidth=1, alpha=0.5)
#     plt.title('Cumulative sliding damage (CSD)')
#     plt.xlabel('strain')
#     plt.ylabel('damage')
#     plt.legend(loc=3)

    ax2 = plt.subplot(224)
    # plt.plot(s_arr_2 , w_arr_2, '--')
    plt.plot(s_arr_6, w_arr_6, 'b', linewidth=2,
             label='pressure = 0 MPa', alpha=0.5)
    plt.plot(s_arr_7, w_arr_7, 'r', linewidth=2,
             label='pressure = 10 MPa', alpha=0.5)
    plt.plot(s_arr_8, w_arr_8, 'g', linewidth=2,
             label='pressure = 20 MPa', alpha=0.5)

    #ax2.plot(s_arr_3, w_N_arr_3, 'b')
    ax2.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax2.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('Damage evolution')
    plt.ylim(0, 1)
    plt.xlabel('Slip(mm)')
    plt.ylabel('Damage')
    plt.legend(loc=4)
#
#     #'''
#     plt.subplot(223)
#     #gs = gridspec.GridSpec(2, 2)
#     #plt.subplot(gs[-1, :])
#     # plt.plot(xs_pi_cum_2, w_arr_2 , '--')
#     plt.plot(xs_pi_cum_1, w_arr_1, 'b', linewidth=1,
#              label='$ \sigma_N = 0$ MPa', alpha=1)
#     plt.plot(xs_pi_cum_2, w_arr_2, 'r', linewidth=1,
#              label='$ \sigma_N = 10$ MPa', alpha=1)
#     plt.plot(xs_pi_cum_3, w_arr_3, 'g', linewidth=1,
#              label='$ \sigma_N = 20$ MPa', alpha=1)
#
#     plt.xlabel('Cumulative sliding(mm)')
#     plt.ylabel('Damage')
#     plt.ylim(0, 1)
#     plt.legend()

    plt.show()
