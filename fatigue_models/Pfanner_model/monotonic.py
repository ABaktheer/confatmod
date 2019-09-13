'''
Created on 24.08.2017

@author: abaktheer
'''
'''
Implementation of the fatigue model for plain concrete [D.Pffanner, 2003] under uniaxial compressive loading
'''
import math

from scipy.optimize import newton, bisect

import matplotlib.pyplot as plt
import numpy as np


def get_monotonic_response(d, leq, Gcl, Ec, fc, eps_c, b):

    # monotonic strain
    eps_cd = np.zeros_like(d, dtype=np.float_)
    # damage factor
    w_arr = np.zeros_like(d, dtype=np.float_)

    #======================================================================
    # calculation of material parameters (Eci, Y_0, eps_c0, gamma_c, g_c12, Y_g)
    #======================================================================
    Eci = (1.0 / (2.0 * Ec)) * (fc / eps_c)**2.0 - \
        (fc / eps_c) + (3.0 / 2.0) * Ec
    print('Eci', Eci)

    gamma_c = ((math.pi ** 2.0) * fc * eps_c) /\
        (2.0 * (((Gcl / leq) - 0.5 * fc *
                 (eps_c * (b - 1.0)) + b * (fc / Ec)))**2.0)
    print('gamma_c ', gamma_c)

    g_c12 = math.pi * np.sqrt(fc * eps_c / (2.0 * gamma_c))
    print('g_c12 ', g_c12)

    #Y_g = 2.0 * Eci / (eps_c * fc)
    Y_g = Eci / fc - 2.0 / eps_c
    print('Y_g', Y_g)

    # calculation of Dc
    D_c = 1.0 - (fc / (eps_c * Ec * (1.0 - b) + b * fc))
    print('D_c', D_c)

    #======================================================================
    # calculation of the strain eps_cd(D)
    #======================================================================
    def get_eps_cd(D):
        if D <= D_c:
            eps_cd = (Ec * (1.0 - D / (1.0 - b * (1.0 - D))) - Eci) /\
                (Ec * (1.0 - D / (1.0 - b * (1.0 - D)))
                 * ((Eci / fc) - (2.0 / eps_c)) + (fc / eps_c ** 2.0))
        else:
            eps_cd_1 = -(eps_c +
                         np.sqrt(2.0 * eps_c *
                                 (1.0 + 1.0 / 0.01) / (gamma_c * fc)))
            eps_cd_2 = -eps_c

        #====================================================================
        # bisection method
        #====================================================================
            def f_1(eps_cd):
                return (gamma_c / (2.0 * eps_c)) * (eps_cd ** 3.0) +\
                    gamma_c * eps_cd**2.0 + ((2.0 + gamma_c * fc * eps_c) /
                                             (2.0 * fc)) * eps_cd + (1.0 - b * (1.0 - D)) / \
                    (Ec * (1.0 - b) * (1.0 - D))
            eps_cd = bisect(f_1, eps_cd_1, eps_cd_2, xtol=1e-10, maxiter=100)

        #====================================================================
        # newton method
        #====================================================================
#             f_dw_n = lambda eps_cd:  (gamma_c / (2.0 * eps_c)) * (eps_cd ** 3.0)   +\
#                 gamma_c * eps_cd**2.0 + ((2.0 + gamma_c * fc * eps_c) /
#                                          (2.0 * fc)) * eps_cd + (1.0 - b * (1.0 - D)) / \
#                 (Ec * (1.0 - b) * (1.0 - D))
#             #             f_dw_n2 = lambda dw_n: 1 + (f_trial * abs((m + 1.0) * K ** m)) /\
#             #                 (E + abs((m + 1.0) * delta_lamda_p * K ** m) + gamma)**2.0
#             eps_cd = newton(
#                 f_dw_n, 0.0, tol=1e-4, maxiter=100)

        return eps_cd

    #======================================================================
    # calculation of dissipated energy due to damage
    #======================================================================
    def get_g_tot_da(D):

        eps_cd = get_eps_cd(D)

        if D <= D_c and -eps_cd <= (fc / (3.0 * Ec)):
            g_casc_da = 0.5 * Ec * (eps_cd**2.0)
            g_cdesc_da = 0.0

        elif D <= D_c and (fc / (3.0 * Ec)) < -eps_cd <= eps_c:
            g_casc_da = (1.0 / Y_g) * ((fc / (2.0 * eps_c**2.0) *
                                        ((fc / (3.0 * Ec))**2.0 - eps_cd**2.0)) +
                                       (Eci + (fc / (Y_g * eps_c**2.0))) *
                                       ((1.0 / Y_g) * np.log((1.0 + Y_g * fc / (3.0 * Ec)) /
                                                             (1.0 - Y_g * eps_cd)) - eps_cd - (fc / (3.0 * Ec)))) +\
                (fc**2.0 / (18.0 * Ec))
            g_cdesc_da = 0.0

#         elif D > D_c:
        else:
            g_casc_da = get_g_casc_da_eps(-eps_c)

            g_cdesc_da = - np.sqrt((2.0 * fc * eps_c) / (gamma_c)) *\
                math.atan(
                np.sqrt(0.5 * gamma_c * fc * eps_c) * (1.0 + (eps_cd / eps_c)))

        get_g_tot_da = g_casc_da + g_cdesc_da
        return get_g_tot_da

    #======================================================================
    # calculation of the energy g_casc(eps_c0)
    #======================================================================
    def get_g_casc_da_eps(eps_cd):
        if -eps_cd <= (fc / (3.0 * Ec)):
            g_casc_da = 0.5 * (eps_cd**2.0) * Ec
        elif (fc / (3.0 * Ec)) < -eps_cd <= eps_c:
            g_casc_da = (1.0 / Y_g) * ((fc / (2.0 * eps_c**2.0) *
                                        ((fc / (3.0 * Ec))**2.0 - eps_cd**2.0)) +
                                       (Eci + (fc / (Y_g * eps_c**2.0))) *
                                       ((1.0 / Y_g) * np.log((1.0 + Y_g * fc / (3.0 * Ec)) /
                                                             (1.0 - Y_g * eps_cd)) - eps_cd - (fc / (3.0 * Ec)))) +\
                (fc**2.0 / (18.0 * Ec))
        return g_casc_da

    #======================================================================
    # calculation of the stress sig_cd(D)
    #======================================================================
    def get_sig_cd(D):
        eps_cd = get_eps_cd(D)

        if -eps_cd <= (fc / (3.0 * Ec)):
            sig_cd = Ec * eps_cd
        elif (fc / (3.0 * Ec)) < -eps_cd <= eps_c:
            sig_cd = (Eci * (eps_cd / fc) + (eps_cd / eps_c)**2.0) / \
                (1.0 - (Eci * (eps_c / fc) - 2.0)
                 * (eps_cd / eps_c)) * fc
        else:
            sig_cd = -1.0 / ((2.0 + gamma_c * fc * eps_c) / (2.0 * fc) +
                             gamma_c * eps_cd + (gamma_c / (2.0 * eps_c)) * eps_cd**2.0)

        return sig_cd

#     #======================================================================
#     # calculation of the stress sig_cd(eps_cd)
#     #======================================================================
#     def get_sig_cd_2(eps_cd):
#
#         if -eps_cd <= (fc / (3.0 * Ec)):
#             sig_cd = Ec * eps_cd
#
#         elif (fc / (3.0 * Ec)) > -eps_cd >= eps_c:
#             sig_cd = (Eci * (eps_cd / fc) + (eps_cd / eps_c)**2.0 ) / \
#                 (1.0 - (Eci * (eps_c / fc) - 2.0)
#                  * (eps_cd / eps_c)) * fc
#         else:
#             sig_cd = -1.0 / ((2.0 + gamma_c * fc * eps_c) / (2.0 * fc) +
#                              gamma_c * eps_cd + (gamma_c / (2.0 * eps_c)) * eps_cd**2.0)
#         print '***************sig_cd', sig_cd
#         return sig_cd

    #======================================================================
    # calculation of the energy g_c12-da(D)
    #======================================================================
    def get_g_c12_da(D):

        eps_cd = get_eps_cd(D)
        g_c12_da = g_c12 + np.sqrt((2.0 * fc * eps_c) / (gamma_c)) *\
            math.atan(
            np.sqrt(0.5 * gamma_c * fc * eps_c) * (1.0 + (eps_cd / eps_c)))

        return g_c12_da

    def get_g_p_da(D):
        eps_cd = get_eps_cd(D)
        sig_da = get_sig_cd(D)
        return 0.5 * sig_da ** 2.0 / (2.0 * Ec * (1.0 - D))

    def get_g_da(D):

        g_tot_da = get_g_tot_da(D)
        g_p_da = get_g_p_da(D)
        return g_tot_da - g_p_da

    eps = get_eps_cd(d)
    sig = get_sig_cd(d)
    g_c12_da = get_g_c12_da(d)
    g_p_da = get_g_p_da(d)
    g_tot_da = get_g_tot_da(d)
    g_da = get_g_da(d)

    return eps, sig, g_c12_da, g_p_da, g_tot_da,  g_da


if __name__ == '__main__':

    n = 1000
    s_levels = np.linspace(0, 0.99, 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= -1
    s_history_1 = s_levels.flatten()
    d_arr = np.hstack([np.linspace(s_history_1[i], s_history_1[i + 1], n)
                       for i in range(len(s_levels) - 1)])

    arr_1 = np.zeros(len(d_arr) + 1)
    arr_2 = np.zeros(len(d_arr) + 1)
    arr_3 = np.zeros(len(d_arr) + 1)
    arr_4 = np.zeros(len(d_arr) + 1)
    arr_5 = np.zeros(len(d_arr) + 1)
    arr_6 = np.zeros(len(d_arr) + 1)

    for i in range(0, len(d_arr)):
        model = get_monotonic_response(
            d_arr[i], 0.25, 0.03, 50300.0, 120.0, 0.003, 0.2)  # (N, leq, Gcl, Ec, fc, eps_c, b):

        arr_1[i + 1] = model[0]
        arr_2[i + 1] = model[1]
        arr_3[i + 1] = model[2]
        arr_4[i + 1] = model[3]
        arr_5[i + 1] = model[4]
        arr_6[i + 1] = model[5]

    plt.subplot(231)
    plt.plot(arr_1[0:-1], arr_2[0:-1])
    plt.xlabel('strain')
    plt.ylabel('stress')

    plt.subplot(232)
    plt.plot(arr_1[0:-1], d_arr)
    plt.xlabel('strain')
    plt.ylabel('damage')

    plt.subplot(233)
    plt.plot(arr_1[0:-1], arr_3[0:-1])
    plt.xlabel('strain')
    plt.ylabel('g_c12')

    plt.subplot(234)
    plt.plot(arr_1[0:-1], arr_4[0:-1])
    plt.xlabel('strain')
    plt.ylabel('g_p')

    plt.subplot(235)
    plt.plot(arr_1[0:-1], arr_5[0:-1])
    plt.xlabel('strain')
    plt.ylabel('g_tot')

    plt.subplot(236)
    plt.plot(arr_1[0:-1], arr_6[0:-1])
    plt.xlabel('strain')
    plt.ylabel('g_da')

    # print e_N

    plt.show()
