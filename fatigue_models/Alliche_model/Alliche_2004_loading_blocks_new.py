'''
Created on 14.06.2017

@author: abaktheer
'''

'''
Implementation of the fatigue model for plain concrete [A.Alliche, 2004] under uniaxial compressive loading
(stress driven algorithm)
'''


'''
To do#

1. model class
2. loading scenario class (reduce repetition)
3. improve printing
'''

import matplotlib.pyplot as plt
import numpy as np


def get_stress_strain(sigma_1_arr, lamda, mu, alpha, beta, g, C0, C1, K, n):

    #-----------------------------------------------------------------------
    # arrays to store the values
    #-----------------------------------------------------------------------
    # normal strain
    eps_1_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)
    # lateral strain
    eps_2_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)
    # damage factor
    w_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)
    f_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)
    D_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)
    phi_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)

    #-----------------------------------------------------------------------
    # material parameters
    #-----------------------------------------------------------------------
    # lame constants [MPa]
    lamda = lamda
    mu = mu
    # fatigue model material parameter
    alpha = alpha
    beta = beta
    g = g
    C0 = C0
    C1 = C1
    K = K
    n = n

    #-----------------------------------------------------------------------
    # state variables
    #-----------------------------------------------------------------------
    #sigma_1_arr[0] = 0
    eps_1_i = 0.0
    eps_2_i = 0.0
    w_i = 0.0
    D_i = 0.0

    for i in range(1, len(sigma_1_arr)):

        sigma_1_i = sigma_1_arr[i]

        eps_2_i = -1.0 * ((lamda + alpha * w_i) * sigma_1_i + g * w_i * (lamda + 2.0 * mu)) / \
            ((lamda + 2.0 * mu) * (2.0 * (lamda + mu) + 4.0 *
                                   w_i * (alpha + beta)) - 2.0 * (lamda + alpha * w_i) ** 2)

        eps_1_i = sigma_1_i / \
            (lamda + 2.0 * mu) - 2.0 * eps_2_i *  \
            (lamda + alpha * w_i) / (lamda + 2.0 * mu)

        f_i = abs(g) * eps_2_i - (C0 + 2 * C1 * w_i)

        kappa_i = (lamda + 2.0 * mu) * (2.0 * (lamda + mu) +
                                        4.0 * w_i * (alpha + beta) -
                                        alpha * (g / (2.0 * C1)) *
                                        (2.0 * eps_2_i + eps_1_i) -
                                        (g**2.0 / (2.0 * C1))) - 2.0 * (lamda + alpha * w_i)**2

        d_sigma_1 = sigma_1_arr[i] - sigma_1_arr[i - 1]
        m = -1.0 * ((lamda + alpha * w_i) / kappa_i) * d_sigma_1

        # loading stage (evolve of the fatigue damage based on (Marigo.85)
        # model)
        if m > 0:
            d_w = m * abs(g) / (2.0 * C1) * (f_i / K)**n
        else:  # unloading stage (no fatigue damage)
            d_w = 0

        w_i = w_i + d_w

        # Energy release rate
        Y_norm = np.sqrt((-g * eps_1_i - alpha * (eps_1_i + 2. * eps_2_i) * eps_1_i
                          - 2. * beta * (eps_1_i**2.0))**2.0 + 2.0 * (-g * eps_2_i - alpha * (eps_1_i + 2. * eps_2_i) * eps_1_i
                                                                      - 2 * beta * (eps_1_i**2.0))**2.0)
        d_D = Y_norm  # * d_w
        D_i += d_D
        # print 'Y=', Y_norm
        print('D=', D_i)

        # Helmholtz free energy
        phi_i = 0.5 * lamda * (eps_1_i + 2.0 * eps_2_i)**2.0 + mu * ((eps_1_i)**2.0 + 2.0 * eps_2_i**2.0) + 2.0 * g * w_i * eps_2_i + alpha * \
            (2.0 * w_i * eps_1_i * eps_2_i + 4.0 * w_i *
             eps_2_i**2.0) + 4.0 * beta * w_i * eps_2_i**2.0

        if w_i > 5.0:
            print(' ----------> No Convergence any more')
            print(i)
            break

        if abs(eps_1_i) > 0.005:
            print(' ----------> No Convergence any more')
            print(i)
            break

        eps_1_arr[i] = eps_1_i
        eps_2_arr[i] = eps_2_i
        w_arr[i] = w_i
        f_arr[i] = f_i
        D_arr[i] = D_i
        phi_arr[i] = phi_i

    return sigma_1_arr, eps_1_arr, eps_2_arr, w_arr, f_arr, D_arr, i, phi_arr


if __name__ == '__main__':

    m = 200  # number of increments in each cycle

#     n1 = 3
#     n2 = 13
#     n3 = 44
#     n4 = 164
#     n5 = 78
#     n6 = 600

    n1 = 10
    n2 = 10
    n3 = 10
    n4 = 50
    n5 = 50
    n6 = 50

    b = 3  # number_of_repeted_blocks

    sigma_u = - 120

    stress_level_1_max = 0.85 * sigma_u
    stress_level_2_max = 0.85 * sigma_u
    stress_level_3_max = 0.85 * sigma_u
    stress_level_4_max = 0.75 * sigma_u
    stress_level_5_max = 0.75 * sigma_u
    stress_level_6_max = 0.75 * sigma_u

    stress_level_1_min = 0.2 * sigma_u
    stress_level_2_min = 0.2 * sigma_u
    stress_level_3_min = 0.2 * sigma_u
    stress_level_4_min = 0.2 * sigma_u
    stress_level_5_min = 0.2 * sigma_u
    stress_level_6_min = 0.2 * sigma_u
    #unloading_ratio = 0.0

    d_0 = np.zeros(1)

    d_1 = np.linspace(0, stress_level_1_max, n1 * 2)
    d_1.reshape(-1, 2)[:, 0] = stress_level_1_max
    #d_1[0] = 0.0
    d_1.reshape(-1, 2)[:, 1] = stress_level_1_min
    d_history_1 = d_1.flatten()
    sig_1_arr = np.hstack([np.linspace(d_history_1[i], d_history_1[i + 1], m, dtype=np.float_)
                           for i in range(len(d_1) - 1)])

    d_2 = np.linspace(0, stress_level_2_max, n2 * 2)
    d_2.reshape(-1, 2)[:, 0] = stress_level_2_max
    d_2.reshape(-1, 2)[:, 1] = stress_level_2_min
    d_history_2 = d_2.flatten()
    sig_2_arr = np.hstack([np.linspace(d_history_2[i], d_history_2[i + 1], m, dtype=np.float_)
                           for i in range(len(d_2) - 1)])

    d_3 = np.linspace(0, stress_level_3_max, n3 * 2)
    d_3.reshape(-1, 2)[:, 0] = stress_level_3_max
    d_3.reshape(-1, 2)[:, 1] = stress_level_3_min
    d_history_3 = d_3.flatten()
    sig_3_arr = np.hstack([np.linspace(d_history_3[i], d_history_3[i + 1], m, dtype=np.float_)
                           for i in range(len(d_3) - 1)])

    d_4 = np.linspace(0, stress_level_4_max, n4 * 2)
    d_4.reshape(-1, 2)[:, 0] = stress_level_4_max
    d_4.reshape(-1, 2)[:, 1] = stress_level_4_min
    d_history_4 = d_4.flatten()
    sig_4_arr = np.hstack([np.linspace(d_history_4[i], d_history_4[i + 1], m, dtype=np.float_)
                           for i in range(len(d_4) - 1)])

    d_5 = np.linspace(0, stress_level_5_max, n5 * 2)
    d_5.reshape(-1, 2)[:, 0] = stress_level_5_max
    d_5.reshape(-1, 2)[:, 1] = stress_level_5_min
    d_history_5 = d_5.flatten()
    sig_5_arr = np.hstack([np.linspace(d_history_5[i], d_history_5[i + 1], m, dtype=np.float_)
                           for i in range(len(d_5) - 1)])

    d_6 = np.linspace(0, stress_level_6_max, n6 * 2)
    d_6.reshape(-1, 2)[:, 0] = stress_level_6_max
    d_6.reshape(-1, 2)[:, 1] = stress_level_6_min
    d_history_6 = d_6.flatten()
    sig_6_arr = np.hstack([np.linspace(d_history_6[i], d_history_6[i + 1], m, dtype=np.float_)
                           for i in range(len(d_6) - 1)])

    sig_0_1_arr = np.linspace(
        d_0[-1], d_history_1[0], m, dtype=np.float_)
    sig_1_2_arr = np.linspace(
        d_history_1[-1], d_history_2[0], m, dtype=np.float_)
    sig_2_3_arr = np.linspace(
        d_history_2[-1], d_history_3[0], m, dtype=np.float_)
    sig_3_4_arr = np.linspace(
        d_history_3[-1], d_history_4[0], m, dtype=np.float_)
    sig_4_5_arr = np.linspace(
        d_history_4[-1], d_history_5[0], m, dtype=np.float_)
    sig_5_6_arr = np.linspace(
        d_history_5[-1], d_history_6[0], m, dtype=np.float_)

    sig_1_1_arr = np.linspace(
        sig_1_arr[-1], sig_1_arr[0], m, dtype=np.float_)
    sig_6_1_arr = np.linspace(
        sig_6_arr[-1], sig_1_arr[0], m, dtype=np.float_)
    sig_6_6_arr = np.linspace(
        sig_6_arr[-1], sig_6_arr[0], m, dtype=np.float_)

    sigma_1_arr = np.hstack(
        (d_0, sig_0_1_arr, sig_1_arr, sig_1_2_arr, sig_2_arr, sig_2_3_arr, sig_3_arr, sig_3_4_arr, sig_4_arr, sig_4_5_arr, sig_5_arr, sig_5_6_arr, sig_6_arr))

    sigma_arr = sigma_1_arr

    for i in range(1, b):

        # for constant repeating
        sigma_arr = np.hstack((sigma_arr, sig_6_1_arr, sig_1_arr, sig_1_2_arr, sig_2_arr, sig_2_3_arr,
                               sig_3_arr, sig_3_4_arr, sig_4_arr, sig_4_5_arr, sig_5_arr,
                               sig_5_6_arr, sig_6_arr))

#         # for repeating (H-L-H) or (L-H-L)
#         if int(i) % 2 == 0:
#             # print 'even'
#             sigma_arr = np.hstack(
#                 (sigma_arr, sig_6_1_arr, sig_1_arr, sig_1_2_arr, sig_2_arr, sig_2_3_arr, sig_3_arr, sig_3_4_arr, sig_4_arr, sig_4_5_arr, sig_5_arr, sig_5_6_arr, sig_6_arr))
#         else:
#             # print 'odd'
#             sigma_arr = np.hstack(
#                 (sigma_arr, sig_6_6_arr, sig_6_arr, sig_5_6_arr, sig_5_arr, sig_4_5_arr, sig_4_arr, sig_3_4_arr, sig_3_arr, sig_2_3_arr, sig_2_arr, sig_1_2_arr, sig_1_arr))

    print(sigma_arr)

    t_arr = np.linspace(0, 1, len(sigma_arr))

    # C120
    sigma_arr, eps_1_arr, eps_2_arr, w_arr, f_arr, D_arr, inc, phi_arr = get_stress_strain(
        sigma_arr, lamda=12500, mu=18750, alpha=2237.5, beta=-2216.5, g=-10.0,
        C0=0.00, C1=0.0019, K=0.00485, n=10)

#     sigma_arr, eps_1_arr, eps_2_arr, w_arr, f_arr, inc = get_stress_strain(
# sigma_arr, lamda=13972.2, mu=20958.3, alpha=2237.5, beta=-2216.5,
# g=-10.0, C0=0.00, C1=0.00188, K=0.003345, n=10)


#     # C120
#     sigma_arr, eps_1_arr, eps_2_arr, w_arr, f_arr, inc = get_stress_strain(
#         sigma_arr, lamda=12500, mu=18750, alpha=2237.5, beta=-2216.5, g=-10.0,
#         C0=0.00, C1=0.0019, K=0.00485, n=10)

#     # C80 - alliche paper
#     sigma_arr, eps_1_arr, eps_2_arr, w_arr, f_arr, inc = get_stress_strain(
# sigma_arr, lamda=10555.55, mu=15833.33, alpha=2237.5, beta=-2216.5,
# g=-9.788, C0=0.00, C1=0.002033, K=0.003345, n=10)

    #-----------------------------------------------------------------------
    # plot 1
    #-----------------------------------------------------------------------
    plt.subplot(231)
    plt.plot(t_arr[0:inc], abs(sigma_arr[0:inc]), 'k', linewidth=1, alpha=1.0)
    plt.title('loading history')

    plt.xlabel('Time')
    plt.ylabel('$\sigma_{1}$')
    # plt.legend(loc=4)

    #-----------------------------------------------------------------------
    # plot 2
    #-----------------------------------------------------------------------
    plt.subplot(232)
    plt.plot(abs(eps_1_arr[0:inc]), abs(
        sigma_arr[0:inc]), 'k', linewidth=1, alpha=1.0)
    plt.title('$ \epsilon_{11} - \sigma_{11}$')
    plt.xlabel('$\epsilon_{11}$')
    plt.ylabel('$\sigma_{11}$[MPa]')
    # plt.legend(loc=4)

    #-----------------------------------------------------------------------
    # plot 2
    #-----------------------------------------------------------------------
    plt.subplot(236)
    plt.plot(abs(sigma_arr[0:inc]), w_arr[0:inc], 'k', linewidth=1, alpha=1.0)
    plt.title('$ \epsilon_{11} - \sigma_{11}$')
    plt.xlabel('$\epsilon_{11}$')
    plt.ylabel('$\sigma_{11}$[MPa]')
    # plt.legend(loc=4)
    #-----------------------------------------------------------------------
    # plot 3
    #-----------------------------------------------------------------------
#     plt.subplot(222)
#     plt.plot(abs(sigma_arr), w_arr, 'k', linewidth=1, alpha=1)
#     plt.xlabel('$\sigma_{11}$')
#     plt.ylabel('Damage')
#     #plt.ylim(0, 1)
#     # plt.legend()

    #-----------------------------------------------------------------------
    # plot 4
    #-----------------------------------------------------------------------
#     plt.subplot(223)
#     plt.plot(t_arr, w_arr, 'b', linewidth=1, alpha=1)
#     plt.xlabel('Time')
#     plt.ylabel('Damage')
#     #plt.ylim(0, 1)
#     # plt.legend()

    #-----------------------------------------------------------------------
    # plot 5
    #-----------------------------------------------------------------------
#     plt.subplot(225)
#
    n = b * (n1 + n2 + n3 + n4 + n5 + n6)
#     eps_1 = np.zeros(n)
#     cycle = np.zeros(n)
#     sig_1 = np.zeros(n)
#
#     for i in range(0, n, 1):
#         idx = m + 2 * i * m  # - 1
#         eps_1[i] = eps_1_arr[idx]
#         cycle[i] = i
#         #sig_1[i] = sigma_1_arr[idx]
#         # print sig_1
#     plt.plot(cycle, eps_1, 'b', linewidth=1, alpha=1)
#     plt.xlabel('number of cycles')
#     plt.ylabel('max $\epsilon_{11}$')
    # plt.legend()

#     #-----------------------------------------------------------------------
#     # plot 5
#     #-----------------------------------------------------------------------
    plt.subplot(234)

    eps_1_max = np.zeros(n)
    eps_1_min = np.zeros(n)
    cycle = np.zeros(n)
    for i in range(0, n, 1):
        idx_1 = m + 2 * i * m - 1
        idx_2 = 2 * i * m
        if idx_1 <= len(eps_1_arr[0:inc]):
            idx_1 = idx_1
        else:
            idx_1 = m + 2 * (i - 1.0) * m - 1
            break

        if idx_2 <= len(eps_1_arr[0:inc]):
            idx_2 = idx_2
        else:
            idx_2 = 1 * (i - 1.0) * m
            break

        eps_1_max[i] = eps_1_arr[int(idx_1)]
        eps_1_min[i] = eps_1_arr[int(idx_2)]
        cycle[i] = i + 1

    plt.plot(cycle[0:i], abs(eps_1_max[0:i]), 'k', linewidth=1, alpha=1)

    #plt.plot(cycle, eps_2 / 1.5, color='gray', linewidth=1.0, alpha=0.5)
    #plt.fill_between(cycle, eps_2, eps_2 / 1.5, facecolor='gray', alpha=0.5)
    plt.xlabel('number of cycles')
    plt.ylabel('max $\epsilon_{11}$')
    #plt.ylim(0, 1)
    # plt.legend()

#     #-----------------------------------------------------------------------
#     # plot 6
#     #-----------------------------------------------------------------------
    plt.subplot(235)
    w = np.zeros(n)
    cycle = np.zeros(n)
    for i in range(0, n, 1):
        idx = m + 2 * i * m - 1
        if idx <= len(w_arr[0:inc]):
            idx = idx
        else:
            idx = m + 2 * (i - 1.0) * m - 1
            break

        w[i] = w_arr[idx]
        cycle[i] = i + 1

    plt.plot(cycle[0:i], w[0:i], 'k', linewidth=1, alpha=1)
    plt.xlabel('number of cycles')
    plt.ylabel('Damage')
    #plt.ylim(0, 1)
    # plt.legend()


#     #-----------------------------------------------------------------------
#     # plot 6
#     #-----------------------------------------------------------------------
    plt.subplot(233)
    D = np.zeros(n)
    cycle = np.zeros(n)
    for i in range(0, n, 1):
        idx = m + 2 * i * m - 1
        if idx <= len(D_arr[0:inc]):
            idx = idx
        else:
            idx = m + 2 * (i - 1.0) * m - 1
            break

        D[i] = D_arr[idx]
        cycle[i] = i + 1

    plt.plot(cycle[0:i], D[0:i], 'k', linewidth=1, alpha=1)
    plt.xlabel('number of cycles')
    plt.ylabel('$Dissipation$')


# #     #-----------------------------------------------------------------
# #     # plot 6
# #     #-----------------------------------------------------------------
#     plt.subplot(236)
#     phi = np.zeros(n)
#     cycle = np.zeros(n)
#     for i in range(0, n, 1):
#         idx = m + 2 * i * m - 1
#         if idx <= len(phi_arr[0:inc]):
#             idx = idx
#         else:
#             idx = m + 2 * (i - 1.0) * m - 1
#             break
#
#         phi[i] = phi_arr[idx]
#         cycle[i] = i + 1
#
#     plt.plot(cycle[0:i], phi[0:i], 'k', linewidth=1, alpha=1)
#     plt.xlabel('number of cycles')
#     plt.ylabel('free energy')

    #-----------------------------------------------------------------------
    # plot 7
    #-----------------------------------------------------------------------
#     plt.subplot(236)
#     plt.plot(eps_2_arr, sigma_1_arr, 'b', linewidth=1, alpha=1.0)
#     plt.title('$ \epsilon_{22} - \sigma_{11}$')
#     plt.xlabel('$\epsilon_{22}$')
#     plt.ylabel('$\sigma_{11}$[MPa]')
#     # plt.legend(loc=4)

    #=========================================================================
    # saving results
    #=========================================================================
    eps_max_record = np.zeros(1)
    eps_min_record = np.zeros(1)
    N_record = np.zeros(1)
    w_record = np.zeros(1)
    stiffness_record = np.zeros(1)

    for i in range(0, i, 1):
        eps_max_record = np.vstack((eps_max_record, eps_1_max[i]))
        eps_min_record = np.vstack((eps_min_record, eps_1_min[i]))
        N_record = np.vstack((N_record, cycle[i]))
        w_record = np.vstack((w_record, w[i]))
        stiffness_record = np.vstack(
            (stiffness_record, stress_level_1_max / eps_1_max[i]))

#     np.savetxt(r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\eps_max.txt',
#                np.transpose(eps_max_record), delimiter=" ", fmt="%s")
#     np.savetxt(r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\eps_2.txt',
#                np.transpose(eps_record_2), delimiter=" ", fmt="%s")
#     np.savetxt(r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\N.txt',
#                np.transpose(N_record), delimiter=" ", fmt="%s")
#     np.savetxt(r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\w.txt',
#                np.transpose(w_record), delimiter=" ", fmt="%s")
#     np.savetxt(r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\stiffness.txt',
#                np.transpose(stiffness_record), delimiter=" ", fmt="%s")

# # saving for 3rd journal paper
#     np.savetxt(r'E:\Publishing\Educative_fatigue_Article\results\Alliche\saving\eps_arr.txt',
#                eps_1_arr[0:inc], delimiter=" ", fmt="%s")
#     np.savetxt(r'E:\Publishing\Educative_fatigue_Article\results\Alliche\saving\sig_arr.txt',
#                sigma_arr[0:inc], delimiter=" ", fmt="%s")
#     np.savetxt(r'E:\Publishing\Educative_fatigue_Article\results\Alliche\saving\w_arr.txt',
#                w_arr[0:inc], delimiter=" ", fmt="%s")
#
#     np.savetxt(r'E:\Publishing\Educative_fatigue_Article\results\Alliche\saving\N.txt',
#                N_record, delimiter=" ", fmt="%s")
#     np.savetxt(r'E:\Publishing\Educative_fatigue_Article\results\Alliche\saving\w.txt',
#                w_record, delimiter=" ", fmt="%s")
#     np.savetxt(r'E:\Publishing\Educative_fatigue_Article\results\Alliche\saving\eps_max.txt',
#                eps_max_record, delimiter=" ", fmt="%s")
#
#     np.savetxt(r'E:\Publishing\Educative_fatigue_papre\results\Alliche\saving\eps_min.txt',
#                eps_min_record, delimiter=" ", fmt="%s")
#     np.savetxt(r'E:\Publishing\Educative_fatigue_papre\results\Alliche\saving\stiffness.txt',
#                stiffness_record, delimiter=" ", fmt="%s")

    plt.show()
