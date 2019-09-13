'''
Created on 07.04.2017

@author: abaktheer
'''

from numpy import \
    zeros_like, sign, linspace, hstack, fabs

from scipy.optimize import newton

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


def get_CP_TDP_Law(eps, sigma_0_1, sigma_0_2, K_1, K_2, gamma_1, gamma_2, E, m, c, S, r):

    sigma_arr = zeros_like(eps)

    eps_N_p_arr = zeros_like(eps)

    w_N_arr = zeros_like(eps)

    sigma_i = 0.0
    alpha_i = 0.0
    r_i = 0.0
    eps_N_p_i = 0.0
    w_i = 0.0
    eps_N_p_i_cum = 0.0

    for i in range(1, len(eps)):

        eps_i = eps[i]
        H = get_heviside(sigma_i)
        sigma_i = (1 - H * w_i) * E * (eps_i - eps_N_p_i)

        if H == 1.0:  # tension
            sigma_0 = sigma_0_1
            K = K_1
            gamma = gamma_1

            sigma_i_1 = E * (eps_i - eps_N_p_i)
            Y_i = 0.5 * E * (eps_i - eps_N_p_i) ** 2.0

            f_trial = abs(sigma_i_1 - gamma * alpha_i) - sigma_0

            # plasticity-damage yield function
            if f_trial > 1e-6:

                delta_lamda_p = f_trial / ((E / (1 - w_i)) + K + gamma)

                eps_N_p_i = eps_N_p_i + delta_lamda_p * \
                    sign(sigma_i_1 - gamma * alpha_i)
                r_i = r_i + delta_lamda_p
                alpha_i = alpha_i + \
                    (delta_lamda_p * sign(sigma_i_1 - gamma * alpha_i))
                w_i = w_i + ((1 - w_i)**c) * ((Y_i / S)**r) * delta_lamda_p

                sigma_i = (1 - w_i) * E * (eps_i - eps_N_p_i)
                alpha_i = alpha_i + delta_lamda_p * \
                    sign(sigma_i_1 - gamma * alpha_i)
                eps_N_p_i_cum = eps_N_p_i_cum + delta_lamda_p / (1.0 - w_i)

        else:  # compression
            sigma_0 = sigma_0_2
            K = K_2
            gamma = gamma_2

            f_trial = abs(sigma_i - gamma * alpha_i) - \
                (sigma_0 + K * r_i ** (m + 1.0))

            # plasticity yield function
            if f_trial > 1e-6:

                def f_dw_n(delta_lamda_p): return delta_lamda_p - f_trial / \
                    (E + abs(K) + gamma + m * gamma * gamma *
                     alpha_i * sign(sigma_i - gamma * alpha_i))

                delta_lamda_p = newton(
                    f_dw_n, 0., tol=1e-6, maxiter=10)

                eps_N_p_i = eps_N_p_i + delta_lamda_p * \
                    sign(sigma_i - gamma * alpha_i)
                r_i = r_i + delta_lamda_p
                alpha_i = alpha_i + \
                    (delta_lamda_p *
                     (sign(sigma_i - gamma * alpha_i) + m * gamma * alpha_i))

                sigma_i = E * (eps_i - eps_N_p_i)

#                 alpha_i = alpha_i + \
#                     (delta_lamda_p / (1.0 +  gamma * delta_lamda_p)) * \
#                     (sign(sigma_i - gamma * alpha_i) +  m* gamma * alpha_i)

        sigma_arr[i] = sigma_i
        eps_N_p_arr[i] = eps_N_p_i
        w_N_arr[i] = w_i

    return eps, sigma_arr, eps_N_p_arr, w_N_arr


if __name__ == '__main__':

    # Check the model behavior
    n = 1000
    s_levels = linspace(0,  0.002, 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= 0
    s_history_1 = s_levels.flatten()

    s_levels = linspace(0,  0.002, 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= -1
    s_history_2 = s_levels.flatten()

    s_levels = linspace(0,  0.002, 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= -1
    s_history_3 = s_levels.flatten()
    #s_history_3 = [0, 0.0002, 0.0, 0.0005, 0.0, 0.005]

    s_levels = linspace(0,  0.0008, 4)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= -1
    s_history_4 = s_levels.flatten()
    s_history_4 = [0, 0.00025, -0.0015]

    s_levels = linspace(0,  -0.01, 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= -1
    s_history_5 = s_levels.flatten()

    s_arr_1 = hstack([linspace(s_history_1[i], s_history_1[i + 1], n)
                      for i in range(len(s_history_1) - 1)])
    s_arr_2 = hstack([linspace(s_history_2[i], s_history_2[i + 1], n)
                      for i in range(len(s_history_2) - 1)])
    s_arr_3 = hstack([linspace(s_history_3[i], s_history_3[i + 1], n)
                      for i in range(len(s_history_3) - 1)])
    s_arr_4 = hstack([linspace(s_history_4[i], s_history_4[i + 1], n)
                      for i in range(len(s_history_4) - 1)])
    s_arr_5 = hstack([linspace(s_history_5[i], s_history_5[i + 1], n)
                      for i in range(len(s_history_5) - 1)])

    #eps, sigma_0_1, sigma_0_2, K_1, K_2, gamma_1, gamma_2, E, m, c, S, r

    s_arr_1, sigma_arr_1, eps_N_p_arr_1, w_N_arr_1 = get_CP_TDP_Law(
        s_arr_1, sigma_0_1=2.0, sigma_0_2=30, K_1=10.0, K_2=10.0, gamma_1=100.0, gamma_2=10.0, E=34000,  m=-0.0, c=1.0, S=2.5e-8, r=1.0)
    s_arr_2, sigma_arr_2, eps_N_p_arr_2, w_N_arr_2 = get_CP_TDP_Law(
        s_arr_2, sigma_0_1=3.0, sigma_0_2=30, K_1=10.0, K_2=10.0, gamma_1=100.0, gamma_2=10.0, E=34000,  m=-0.0, c=1.0, S=2.5e-8, r=1.0)
    s_arr_3, sigma_arr_3, eps_N_p_arr_3, w_N_arr_3 = get_CP_TDP_Law(
        s_arr_3, sigma_0_1=4.0, sigma_0_2=30, K_1=10.0, K_2=0.0, gamma_1=100.0, gamma_2=10.0, E=34000,  m=-0.0, c=1.0, S=2.5e-8, r=1.0)
    s_arr_4, sigma_arr_4, eps_N_p_arr_4, w_N_arr_4 = get_CP_TDP_Law(
        s_arr_4, sigma_0_1=5.0, sigma_0_2=30, K_1=10.0, K_2=10000.0, gamma_1=10000.0, gamma_2=10000.0, E=34000,  m=-0.0, c=1.0,  S=2.5e-8, r=1.0)
    s_arr_5, sigma_arr_5, eps_N_p_arr_5, w_N_arr_5 = get_CP_TDP_Law(
        s_arr_5, sigma_0_1=2.0, sigma_0_2=30, K_1=10.0, K_2=10.0, gamma_1=10000.0, gamma_2=10000.0, E=34000,  m=-0.1, c=1.0, S=0.0000005, r=1.1)

    ax2 = plt.subplot(221)
    #ax2.plot(s_arr_1, sigma_arr_1, color='k', label='sigma_0_1=2.0')
    #ax2.plot(s_arr_2, sigma_arr_2, color='r', label='sigma_0_1=3.0')
    #ax2.plot(s_arr_3, sigma_arr_3, color='g', label='sigma_0_1=4.0')
    ax2.plot(s_arr_4, sigma_arr_4, color='b', label='sigma_0_1=5.0')

    ax2.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax2.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('Normal behavior (CP,TD) - monotonic')
    plt.xlabel('strain')
    plt.ylabel('stress(MPa)')
    # plt.legend(loc=1)

    ax2 = plt.subplot(222)
    ax2.plot(s_arr_1, w_N_arr_1, color='k', label='no hardening')
    ax2.plot(s_arr_2, w_N_arr_2, color='r', label='a=-0.05')
    ax2.plot(s_arr_3, w_N_arr_3, color='g', label='a=-0.1')
    ax2.plot(s_arr_4, w_N_arr_4, color='b', label='a=-0.2')

    ax2.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax2.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('Normal behavior (CP,TD) - monotonic')
    plt.xlabel('strain')
    plt.ylabel('stress(MPa)')
    plt.legend(loc=4)

    ax2 = plt.subplot(223)

    ax2.plot(s_arr_5, sigma_arr_5, color='y', label='a=0.0')
    ax2.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax2.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('Normal behavior (CP,TD) - monotonic')
    plt.xlabel('strain')
    plt.ylabel('stress(MPa)')
    plt.legend(loc=4)

    plt.show()
