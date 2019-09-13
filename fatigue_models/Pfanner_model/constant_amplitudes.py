'''
Created on 24.08.2017

@author: abaktheer
'''
'''
Implementation of a fatigue model for plain concrete [D.Pfanner, 2003] under uniaxial compressive loading
'''


import math
from scipy.optimize import newton, bisect
import matplotlib.pyplot as plt
import numpy as np


def get_fatigue_response(N, sig_max, sig_min, T, leq, Gcl, Ec, fc, eps_c, b, xci_c, inc):

    #======================================================================
    # calculation of material parameters (Eci, Y_0, eps_c0, gamma_c, g_c12, Y_g)
    #======================================================================
    Eci = (1.0 / (2.0 * Ec)) * (fc / eps_c)**2.0 - \
        (fc / eps_c) + (3.0 / 2.0) * Ec

    Y_0 = (Eci / fc) * \
        (((1.0 / fc) - (2.0 / (Eci * eps_c))) * sig_max + 1)

    if np.abs(sig_max) <= (fc / 3.0):
        eps_c0 = sig_max / Ec
    else:
        eps_c0 = ((eps_c ** 2.0) / 2.0) * \
            (np.sqrt((Y_0**2.0) + ((4.0 * sig_max) / (fc * eps_c**2.0))) - Y_0)

    gamma_c = ((math.pi ** 2.0) * fc * eps_c) /\
        (2.0 * (((Gcl / leq) - 0.5 * fc *
                 (eps_c * (b - 1.0)) + b * (fc / Ec)))**2.0)

    g_c12 = math.pi * np.sqrt(fc * eps_c / (2.0 * gamma_c))

    Y_g = Eci / fc - 2.0 / eps_c

    R = sig_min / sig_max
    kappa_r = 1.0 + R * (sig_max / fc)**10.0

    # calculation of Dc
    D_c = 1.0 - (fc / (eps_c * Ec * (1.0 - b) + b * fc))

    # calculation of Dc_fat
    D_c_da_fat = 1.0 + 1.0 / ((1.0 / sig_max) * Ec * (1.0 - b) * (
        eps_c + np.sqrt(-2.0 * eps_c * (1.0 + fc / sig_max) / (gamma_c * fc))) - b)

    #print('D_c_da_fat ', D_c_da_fat)

    #======================================================================
    # calculation of the fatigue life (Nf) (based on the given Woehler curves)
    #======================================================================

    def get_Nf(sig_max, sig_min, T):

        #         Nf = 10.0**((np.abs(sig_max) / fc + 0.053 * (1.0 - 0.445 * R) *
        # np.log10(T) + 0.2 * R - 1.2) / (-0.133 * (1.0 - 0.779 * R)))

        # fib model code 2010
        S_max = np.abs(sig_max) / fc
        S_min = np.abs(sig_min) / fc
        Y = (0.45 + 1.8 * S_min) / (1.0 + 1.8 * S_min - 0.3 * S_min**2.0)
        Nf = 10.0**((8.0 / (Y - 1.0)) * (S_max - 1.0))
        return Nf

    '''====================================================================
    part 1 - damage dissipation
    ===================================================================='''

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
            eps_cd = bisect(f_1, eps_cd_1, eps_cd_2, xtol=1e-6, maxiter=100)

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
    # calculation of the strain eps_cd_3(D)
    #======================================================================
    def get_eps_cd_3(D):

        #         eps_cd_1 = -(eps_c +
        #                      np.sqrt(2.0 * eps_c *
        #                              (1.0 + 1.0 / 0.01) / (gamma_c * fc)))
        #         eps_cd_2 = -eps_c
        #
        #         #==================================================================
        #         # bisection method
        #         #==================================================================
        #         def f_1(eps_cd):
        #             return (gamma_c / (2.0 * eps_c)) * (eps_cd ** 3.0)   +\
        #                 gamma_c * eps_cd**2.0 + ((2.0 + gamma_c * fc * eps_c) /
        #                                          (2.0 * fc)) * eps_cd + (1.0 - b * (1.0 - D)) / \
        #                 (Ec * (1.0 - b) * (1.0 - D))
        #         eps_cd_3 = bisect(
        #             f_1, eps_cd_1, eps_cd_2, xtol=1e-5, maxiter=100)

        #==================================================================
        # newton method
        #==================================================================
        def f_dw_n(eps_cd_3): return (gamma_c / (2.0 * eps_c)) * (eps_cd_3 ** 3.0) +\
            gamma_c * eps_cd_3**2.0 + ((2.0 + gamma_c * fc * eps_c) /
                                       (2.0 * fc)) * eps_cd_3 + (1.0 - b * (1.0 - D)) / \
            (Ec * (1.0 - b) * (1.0 - D))
        #             f_dw_n2 = lambda dw_n: 1 + (f_trial * abs((m + 1.0) * K ** m)) /\
        #                 (E + abs((m + 1.0) * delta_lamda_p * K ** m) + gamma)**2.0
        eps_cd_3 = newton(
            f_dw_n, 0.0, tol=1e-6, maxiter=500)
        return eps_cd_3

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
        else:
            g_casc_da = get_g_casc_da_eps(-eps_c)
        return g_casc_da

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

        elif D > D_c:
            g_casc_da = get_g_casc_da_eps(-eps_c)
            g_cdesc_da = - np.sqrt((2.0 * fc * eps_c) / (gamma_c)) *\
                math.atan(
                np.sqrt(0.5 * gamma_c * fc * eps_c) * (1.0 + (eps_cd / eps_c)))

        get_g_tot_da = g_casc_da + g_cdesc_da
        return get_g_tot_da

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

    #======================================================================
    # calculation of the stress sig_cd(D)
    #======================================================================
    def get_sig_cd_3(D):
        eps_cd = get_eps_cd(D)

        sig_cd_3 = -1.0 / ((2.0 + gamma_c * fc * eps_c) / (2.0 * fc) +
                           gamma_c * eps_cd + (gamma_c / (2.0 * eps_c)) * eps_cd**2.0)

        return sig_cd_3

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
        sig_da = get_sig_cd(D)
        return 0.5 * sig_da ** 2.0 / (2.0 * Ec * (1.0 - D))

    def get_g_da(D):
        g_tot_da = get_g_tot_da(D)
        g_p_da = get_g_p_da(D)
        return g_tot_da - g_p_da

    '''======================================================================
    part 2 - fatigue damage dissipation
    ======================================================================'''
    #======================================================================
    # calculate the damage using bisection method (implicit function)
    #======================================================================
    def get_dcf(sig_max):

        d_1 = 0.01
        d_2 = 0.99

        def f(D): return 0.5 * ((sig_max ** 2.0 - get_sig_cd_3(D) ** 2.0) /
                                ((1.0 - D) * Ec)) - (1.0 - (1.0 - D) ** xci_c) * get_g_c12_da(D)

        #print('D_c', D_c)

        d_cf = bisect(f, D_c, d_2, xtol=1e-6, maxiter=100)


#         f_dw_n = lambda d_cf:  0.5 * ((sig_max ** 2.0 - get_sig_cd_3(d_cf) ** 2.0) /
#                                       ((1.0 - d_cf) * Ec)) - (1.0 - (1.0 - d_cf) ** xci_c) * get_g_c12_da(d_cf)
#
#         d_cf = newton(
#             f_dw_n, D_c, tol=1e-4, maxiter=100)
        # print 'd_cf ', d_cf
        # print 'get_sig_cd_3(d_cf) ', get_sig_cd_3(d_cf)
        # print 'get_g_c12_da(d_cf) ', get_g_c12_da(d_cf)

        return d_cf

    #======================================================================
    # calculation of fatigue strain eps_f_fat
    #======================================================================
    def get_eps_f_fat(sig_max):

        D_cf = get_dcf(sig_max)
        #print('D_cf', D_cf)

        g_tot_da = get_g_tot_da(D_c_da_fat)

        g_casc_eps_c0 = get_g_casc_da_eps(eps_c0)

        eps_cd = get_eps_cd_3(D_cf)
        #print('eps_cd', eps_cd)
        eps_c3_fat = get_eps_cd(D_c_da_fat)
        #print('eps_c3_fat', eps_c3_fat)
        sig_c3 = get_sig_cd_3(D_cf)
        #print('sig_c3', sig_c3)
        eps_f_fat = (1.0 / sig_max) * (g_tot_da - g_casc_eps_c0) +\
            eps_c0 + eps_cd - eps_c3_fat + \
            (sig_max - sig_c3) / ((1.0 - D_cf) * Ec)

        return eps_f_fat

    #======================================================================
    # calculation of eps_2_fat, eps_3_fat
    #======================================================================
    def get_delta_eps_2_fat(sig_max):
        delta_eps_2_fat = kappa_r * (0.3 * (sig_max / fc)**3.0 + 0.46)
        return delta_eps_2_fat

    def get_delta_eps_3_fat(sig_max):
        delta_eps_3_fat = kappa_r * (-0.3 * (sig_max / fc)**2.0 + 0.8)
        return delta_eps_3_fat

    #======================================================================
    # calculation of J(N)
    #======================================================================
    def get_delta_J(N, sig_max, sig_min, T):

        delta_eps_2_fat = get_delta_eps_2_fat(sig_max)
        delta_eps_3_fat = get_delta_eps_3_fat(sig_max)
        eps_f_fat = get_eps_f_fat(sig_max)

        Nf = get_Nf(sig_max, sig_min, T)

        eps_0 = eps_c0
        delta_N = 0.2

        delta_J_1 = (eps_f_fat / eps_0 - 1.0) *\
            ((delta_eps_2_fat * ((delta_N**2.0 - 1.0) *
                                 delta_eps_2_fat + (1.0 - delta_N) * delta_eps_3_fat * delta_N)) /
             (delta_eps_2_fat - delta_N * delta_eps_3_fat)) *\
            ((((delta_N * delta_eps_3_fat - delta_eps_2_fat) /
                ((delta_N**2.0) * (delta_eps_2_fat - delta_eps_3_fat))) * (N / Nf) + 1.0)**(-1.0) - 1.0)

        delta_J_2 = (eps_f_fat / eps_0 - 1.0) * \
            (delta_eps_3_fat - delta_eps_2_fat) * (N / Nf)

        delta_J_3 = (eps_f_fat / eps_0 - 1.0) * (1.0 - delta_eps_3_fat) * \
            (N / Nf) ** (math.log((delta_N * (delta_eps_2_fat - delta_eps_3_fat)) /
                                  (delta_eps_3_fat - 1.0), (1.0 - delta_N)))

        delta_J = 1.0 + delta_J_1 + delta_J_2 + delta_J_3

        return delta_J

    #======================================================================
    # calculation of the fatigue strain eps_f
    #======================================================================
    def get_eps_f(N, sig_max, sig_min, T):
        eps_0 = eps_c0
        delta_J = get_delta_J(N, sig_max, sig_min, T)
        eps_f = delta_J * eps_0
        return eps_f

    #======================================================================
    # calculation of the energy g_tot_fat(D)
    #======================================================================
    def get_g_tot_fat(N, sig_max, sig_min, T):

        eps_f = get_eps_f(N, sig_max, sig_min, T)
        g_casc_eps_c0 = get_g_casc_da_eps(eps_c0)
        g_tot_fat = g_casc_eps_c0 + sig_max * (eps_f - eps_c0)
        return g_tot_fat

    #======================================================================
    # calculation of the energy g_p(D)
    #======================================================================
    def get_g_p_fat(sig_max, D):
        g_p_fat = 0.5 * ((sig_max) ** 2.0) / (Ec * (1.0 - D))
        return g_p_fat

    #======================================================================
    # calculation of the energy g_fat(D)
    #======================================================================
    def get_g_fat(N, sig_max, sig_min, T, D):
        g_tot_fat = get_g_tot_fat(N, sig_max, sig_min, T)
        g_p_fat = get_g_p_fat(sig_max, D)
        g_fat = g_tot_fat - g_p_fat
        return g_fat

    #======================================================================
    # check the energy equality (g_da = g_fat)
    #======================================================================
    def get_fatigue_state(N, sig_max, sig_min, T):

        #         D_1 = 0.01
        #         D_2 = 0.99
        #
        #         def f_2(D):
        #             return get_g_da(D) - get_g_fat(N, sig_max, sig_min, T, D)
        #
        #         D = bisect(f_2, D_1, D_2, xtol=1e-8, maxiter=200)

        def f_dw_n(D): return get_g_da(
            D) - get_g_fat(N, sig_max, sig_min, T, D)

        D = newton(
            f_dw_n, 0.0, tol=1e-3, maxiter=200)

        return D

    # fatigue strain
    eps_f = np.zeros(1, dtype=np.float_)
    # damage factor
    w_arr = np.zeros(1, dtype=np.float_)
    # number of cycles
    N_arr = np.zeros(1, dtype=np.float_)

    for i in range(1, N, inc):

        w_arr = np.vstack((w_arr, get_fatigue_state(i, sig_max, sig_min, T)))
        eps_f = np.vstack((eps_f, get_eps_f(i, sig_max, sig_min, T)))
        N_arr = np.vstack((N_arr, i))
        eps_f_fat = get_eps_f_fat(sig_max)
        #print('eps_f_fat', eps_f_fat)
        print('cycle: ', i)
        # print('---------------------------------------------')
        if get_eps_f(i, sig_max, sig_min, T) < eps_f_fat * 1.0:
            print('maximum fatigue strain reached')

            break

        if get_eps_f(i, sig_max, sig_min, T) < eps_f_fat * 1.0:
            print('maximum fatigue strain reached')

            break

    return N_arr, w_arr, eps_f


if __name__ == '__main__':

    N = 1000000

    S_max_1 = 0.65
    S_max_2 = 0.7
    S_max_3 = 0.7
    S_max_4 = 0.7
    S_max_5 = 0.7

    S_min = 0.2

    fc = 120.0

    sig_max_1 = - S_max_1 * fc
    sig_max_2 = - S_max_2 * fc
    sig_max_3 = - S_max_3 * fc
    sig_max_4 = - S_max_4 * fc
    sig_max_5 = - S_max_5 * fc

    sig_min = - S_min * fc

    D_1 = get_fatigue_response(
        N, sig_max_1, sig_min, 0.1, 0.25, 0.03, 50000, fc, 0.003, 0.2, 0.2, 10)
#     D_2 = get_fatigue_response(
#         N, sig_max_2, sig_min, 0.1, 0.25, 0.03, 50000, fc, 0.003, 0.2, 0.2, 100)
#     D_3 = get_fatigue_response(
#         N, sig_max_3, sig_min, 0.1, 0.25, 0.03, 50000.0, fc, 0.003, 0.2, 0.2, 100)
#     D_4 = get_fatigue_response(
#         N, sig_max_4, sig_min, 0.1, 0.25, 0.03, 50000.0, fc, 0.003, 0.2, 0.2, 100)
#     D_5 = get_fatigue_response(
# N, sig_max_5, sig_min, 0.1, 0.25, 0.03, 50000.0, fc, 0.003, 0.2, 0.2,
# 100)

    plt.subplot(231)
    plt.plot(D_1[0][1:], D_1[2][1:], 'k', label='S=0.65')
#     plt.plot(D_2[0][1:], D_2[2][1:], 'r', label='S=0.70')
#     plt.plot(D_3[0][1:], D_3[2][1:], 'g', label='S=0.75')
#     plt.plot(D_4[0][1:], D_4[2][1:], 'b', label='S=0.80')
#     plt.plot(D_5[0][1:], D_5[2][1:], 'y', label='S=0.85')
    plt.xlabel('N')
    plt.ylabel('strain')
    plt.legend(loc=3)

    plt.subplot(232)
    plt.plot(D_1[0][1:], ((1.0 - D_1[1][1:]) / (1.0 - D_1[1][1])) * 100.0, 'k')
#     plt.plot(D_2[0][1:], ((1.0 - D_2[1][1:]) / (1.0 - D_2[1][1])) * 100.0, 'r')
#     plt.plot(D_3[0][1:], ((1.0 - D_3[1][1:]) / (1.0 - D_3[1][1])) * 100.0, 'g')
#     plt.plot(D_4[0][1:], ((1.0 - D_4[1][1:]) / (1.0 - D_4[1][1])) * 100.0, 'b')
#     plt.plot(D_5[0][1:], ((1.0 - D_5[1][1:]) / (1.0 - D_5[1][1])) * 100.0, 'y')
    plt.xlabel('N')
    plt.ylabel('stiffness %')

    plt.subplot(233)
    plt.plot(D_1[0][1:], D_1[1][1:], 'k')
#     plt.plot(D_2[0][1:], D_2[1][1:], 'r')
#     plt.plot(D_3[0][1:], D_3[1][1:], 'g')
#     plt.plot(D_4[0][1:], D_4[1][1:], 'b')
#     plt.plot(D_5[0][1:], D_5[1][1:], 'y')
    plt.xlabel('N')
    plt.ylabel('damage')

    # for calibration
    n_1, s_1 = np.loadtxt(
        r'E:\Publishing\Uniaxial_fatigue_model_Article\results\Kim_1996\data\C84\creep_fatigue\S075.txt')

    plt.subplot(234)
    plt.plot(D_1[0][1:] / D_1[0][-1], -D_1[2][1:], 'k')
    plt.plot(n_1, s_1, 'r')
#     plt.plot(D_2[0][1:] / D_2[0][-1], D_2[2][1:], 'r')
#     plt.plot(D_3[0][1:] / D_3[0][-1], D_3[2][1:], 'g')
#     plt.plot(D_4[0][1:] / D_4[0][-1], D_4[2][1:], 'b')
#     plt.plot(D_5[0][1:] / D_5[0][-1], D_5[2][1:], 'y')
    plt.xlabel('N/Nf')
    plt.ylabel('strain')

    plt.subplot(235)
    plt.plot(D_1[0][1:] / D_1[0][-1],
             ((1.0 - D_1[1][1:]) / (1.0 - D_1[1][1])) * 100.0, 'k')
#     plt.plot(D_2[0][1:] / D_2[0][-1],
#              ((1.0 - D_2[1][1:]) / (1.0 - D_2[1][1])) * 100.0, 'r')
#     plt.plot(D_3[0][1:] / D_3[0][-1],
#              ((1.0 - D_3[1][1:]) / (1.0 - D_3[1][1])) * 100.0, 'g')
#     plt.plot(D_4[0][1:] / D_4[0][-1],
#              ((1.0 - D_4[1][1:]) / (1.0 - D_4[1][1])) * 100.0, 'b')
#     plt.plot(D_5[0][1:] / D_5[0][-1],
#              ((1.0 - D_5[1][1:]) / (1.0 - D_5[1][1])) * 100.0, 'y')
    plt.xlabel('N/Nf')
    plt.ylabel('stiffness %')

    plt.subplot(236)
    plt.plot(D_1[0][1:] / D_1[0][-1], D_1[1][1:], 'k')
#     plt.plot(D_2[0][1:] / D_2[0][-1], D_2[1][1:], 'r')
#     plt.plot(D_3[0][1:] / D_3[0][-1], D_3[1][1:], 'g')
#     plt.plot(D_4[0][1:] / D_4[0][-1], D_4[1][1:], 'b')
#     plt.plot(D_5[0][1:] / D_5[0][-1], D_5[1][1:], 'y')
    plt.xlabel('N/Nf')
    plt.ylabel('damage')

#     np.savetxt(r'E:\Publishing\Educative_fatigue_Article\results\Pfanner\saving\N.txt',
#                D_1[0], delimiter=" ", fmt="%s")
#     np.savetxt(r'E:\Publishing\Educative_fatigue_Article\results\Pfanner\saving\w.txt',
#                D_1[1], delimiter=" ", fmt="%s")
#     np.savetxt(r'E:\Publishing\Educative_fatigue_Article\results\Pfanner\saving\eps.txt',
#                D_1[2], delimiter=" ", fmt="%s")
#     np.savetxt(r'E:\Publishing\Educative_fatigue_Article\results\Pfanner\saving\stiffness.txt',
#                ((1.0 - D_1[1][1:]) / (1.0 - D_1[1][1])) * 100.0, delimiter=" ", fmt="%s")

    plt.show()
