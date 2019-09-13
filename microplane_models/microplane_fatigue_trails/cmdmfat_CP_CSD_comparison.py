'''
Created on 06.04.2017

@author: abaktheer
'''

from numpy import \
    array, zeros, trace, zeros_like, linspace, hstack

from cmdmfat_CP_CSD_Jir import MATS3DMicroplaneDamageJir
from cmdmfat_CP_CSD_Wu import MATS3DMicroplaneDamageWu
import matplotlib.pyplot as plt

import numpy as np


def map3d_tns4_to_tns2(tns4):
    '''
    Map a fourth order tensor to a matrix assuming minor and major symmetry,
    e.g. D_tns(3,3,3,3) to D_mtx (6x6) in engineering notation.
    (Note: Explicit assignment of components used for speed-up.)
    '''
    n_eng = 6
    tns2 = zeros([n_eng, n_eng])

    tns2[0, 0] = tns4[0, 0, 0, 0]
    tns2[0, 1] = tns2[1, 0] = tns4[0, 0, 1, 1]
    tns2[0, 2] = tns2[2, 0] = tns4[0, 0, 2, 2]
    tns2[0, 3] = tns2[3, 0] = tns4[0, 0, 1, 2]
    tns2[0, 4] = tns2[4, 0] = tns4[0, 0, 0, 2]
    tns2[0, 5] = tns2[5, 0] = tns4[0, 0, 0, 1]

    tns2[1, 1] = tns4[1, 1, 1, 1]
    tns2[1, 2] = tns2[2, 1] = tns4[1, 1, 2, 2]
    tns2[1, 3] = tns2[3, 1] = tns4[1, 1, 1, 2]
    tns2[1, 4] = tns2[4, 1] = tns4[1, 1, 0, 2]
    tns2[1, 5] = tns2[5, 1] = tns4[1, 1, 0, 1]

    tns2[2, 2] = tns4[2, 2, 2, 2]
    tns2[2, 3] = tns2[3, 2] = tns4[2, 2, 1, 2]
    tns2[2, 4] = tns2[4, 2] = tns4[2, 2, 0, 2]
    tns2[2, 5] = tns2[5, 2] = tns4[2, 2, 0, 1]

    tns2[3, 3] = tns4[1, 2, 1, 2]
    tns2[3, 4] = tns2[4, 3] = tns4[1, 2, 0, 2]
    tns2[3, 5] = tns2[5, 3] = tns4[1, 2, 0, 1]

    tns2[4, 4] = tns4[0, 2, 0, 2]
    tns2[4, 5] = tns2[5, 4] = tns4[0, 2, 0, 1]

    tns2[5, 5] = tns4[0, 1, 0, 1]

    return tns2


if __name__ == '__main__':

    # Check the model behavior
    n = 10
    s_levels = linspace(0, -0.015, 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= -1
    s_history = s_levels.flatten()

    # cyclic loading
#     s_history_2 = [-0, -0.01, -0.007, -0.015, -0.011, -0.023, -0.017, -0.031, -.023, -
# 0.04, -0.032, -0.05, -0.041, -0.06, -0.050, -.07, -.058, -.08, -0.066,
# -.09, -0.075, - 0.1]

    s_history_2 = [0, -0.01]

    s_arr_1 = hstack([linspace(s_history[i], s_history[i + 1], n)
                      for i in range(len(s_levels) - 1)])
    s_arr_2 = hstack([linspace(s_history_2[i], s_history_2[i + 1], 50)
                      for i in range(len(s_history_2) - 1)])

    eps_1 = array([array([[s_arr_1[i], 0, 0],
                          [0, 0, 0],
                          [0, 0, 0]]) for i in range(0, len(s_arr_1))])

    eps_2 = array([array([[s_arr_2[i], 0, 0],
                          [0, 0, 0],
                          [0, 0, 0]]) for i in range(0, len(s_arr_2))])

    m1 = MATS3DMicroplaneDamageJir()

    m2 = MATS3DMicroplaneDamageWu()
    m2.zeta_G = 0

    m3 = MATS3DMicroplaneDamageWu()
    m3.zeta_G = 0.5

    m4 = MATS3DMicroplaneDamageWu()
    m4.zeta_G = 1.0

    sigma_1 = zeros_like(eps_1)
    sigma_kk_1 = zeros(len(s_arr_1) + 1)
    w_1 = zeros((len(eps_1[:, 0, 0]), 28))
    eps_P_N_1 = zeros((len(eps_1[:, 0, 0]), 28))
    eps_Pi_T_1 = zeros((len(eps_1[:, 0, 0]), 28, 3))
    e_1 = zeros((len(eps_1[:, 0, 0]), 28, 3))
    e_T_1 = zeros((len(eps_1[:, 0, 0]), 28, 3))
    e_N_1 = zeros((len(eps_1[:, 0, 0]), 28))
    sctx_1 = zeros((len(eps_1[:, 0, 0]) + 1, 28, 13))

    sigma_2 = zeros_like(eps_1)
    sigma_kk_2 = zeros(len(s_arr_1) + 1)
    w_2 = zeros((len(eps_1[:, 0, 0]), 28))
    eps_P_N_2 = zeros((len(eps_1[:, 0, 0]), 28))
    eps_Pi_T_2 = zeros((len(eps_1[:, 0, 0]), 28, 3))
    e_2 = zeros((len(eps_1[:, 0, 0]), 28, 3))
    e_T_2 = zeros((len(eps_1[:, 0, 0]), 28, 3))
    e_N_2 = zeros((len(eps_1[:, 0, 0]), 28))
    sctx_2 = zeros((len(eps_1[:, 0, 0]) + 1, 28, 13))

    sigma_3 = zeros_like(eps_1)
    sigma_kk_3 = zeros(len(s_arr_1) + 1)
    sctx_3 = zeros((len(eps_1[:, 0, 0]) + 1, 28, 13))

    sigma_4 = zeros_like(eps_1)
    sigma_kk_4 = zeros(len(s_arr_1) + 1)
    sctx_4 = zeros((len(eps_1[:, 0, 0]) + 1, 28, 13))

    for i in range(0, len(eps_1[:, 0, 0])):

        sigma_1[i, :] = m1.get_corr_pred(
            sctx_1[i, :], eps_1[i, :], sigma_kk_1[i])[0]
        S_1 = m4.get_corr_pred(
            sctx_1[i, :], eps_1[i, :], sigma_kk_1[i])[1]
        print(S_1)
        print(S_1.shape)

        S_new = map3d_tns4_to_tns2(S_1)
        print(S_new)
        print(S_new.shape)
        plt.subplot(221)
        plt.imshow(S_new, cmap='gray')
        plt.colorbar()
        # plt.show()

        sigma_kk_1[i + 1] = trace(sigma_1[i, :])
        sctx_1[
            i + 1] = m1._get_state_variables(sctx_1[i, :], eps_1[i, :], sigma_kk_1[i])

        sigma_2[i, :] = m2.get_corr_pred(
            sctx_2[i, :], eps_1[i, :], sigma_kk_2[i])[0]

        S_2 = m2.get_corr_pred(
            sctx_2[i, :], eps_1[i, :], sigma_kk_2[i])[1]
        print(S_2)
        print(S_2.shape)

        S_new_2 = map3d_tns4_to_tns2(S_2)
        print(S_new_2)
        print(S_new_2.shape)
        plt.subplot(222)
        plt.imshow(S_new_2, cmap='gray')
        plt.colorbar()
        plt.show()

        sigma_kk_2[i + 1] = trace(sigma_1[i, :])
        sctx_2[
            i + 1] = m2._get_state_variables(sctx_2[i, :], eps_1[i, :], sigma_kk_2[i])

        sigma_3[i, :] = m3.get_corr_pred(
            sctx_3[i, :], eps_1[i, :], sigma_kk_3[i])[0]
        sigma_kk_3[i + 1] = trace(sigma_3[i, :])
        sctx_3[
            i + 1] = m3._get_state_variables(sctx_3[i, :], eps_1[i, :], sigma_kk_3[i])

        sigma_4[i, :] = m4.get_corr_pred(
            sctx_4[i, :], eps_1[i, :], sigma_kk_4[i])[0]
        sigma_kk_4[i + 1] = trace(sigma_4[i, :])
        sctx_4[
            i + 1] = m4._get_state_variables(sctx_4[i, :], eps_1[i, :], sigma_kk_4[i])

        w_1[i, :] = sctx_1[i, :, 5]
        eps_P_N_1[i, :] = sctx_1[i, :, 4]
        eps_Pi_T_1[i, :, :] = sctx_1[i, :, 10:13]
        e_1[i, :] = m2._get_e_vct_arr(eps_1[i, :])
        e_T_1[i, :] = m2._get_e_T_vct_arr_2(eps_1[i, :])
        e_N_1[i, :] = m2._get_e_N_arr(e_1[i, :])

#     for i in range(0, len(eps_2[:, 0, 0])):
#
#         sigma_1[i, :] = m1.get_corr_pred(
#             sctx_1[i, :], eps_2[i, :], sigma_kk_1[i])[0]
#         sigma_kk_1[i + 1] = trace(sigma_1[i, :])
#         sctx_1[
#             i + 1] = m1._get_state_variables(sctx_1[i, :], eps_2[i, :], sigma_kk_1[i])
#
#
#         w_2[i, :] = sctx_2[i, :, 5]
#         eps_P_N_2[i, :] = sctx_2[i, :, 4]
#         eps_Pi_T_2[i, :, :] = sctx_2[i, :, 10:13]
#
#         e_2[i, :] = m2._get_e_vct_arr(eps_2[i, :])
#         e_T_2[i, :] = m2._get_e_T_vct_arr_2(eps_2[i, :])
#         e_N_2[i, :] = m2._get_e_N_arr(e_2[i, :])

    plt.subplot(221)
    plt.plot(eps_1[:, 0, 0], sigma_1[:, 0, 0], 'r',
             linewidth=1, label='(Jirasek)')
    plt.plot(eps_1[:, 0, 0], sigma_2[:, 0, 0], 'g--',
             linewidth=1, label='(Wu, zata=0 )')
    plt.plot(eps_1[:, 0, 0], sigma_3[:, 0, 0], 'b--',
             linewidth=1, label='(Wu, zeta=0.5)')
    plt.plot(eps_1[:, 0, 0], sigma_4[:, 0, 0], 'k--',
             linewidth=1, label='(Wu, zeta=1 )')
    #plt.plot(eps_1[:, 0, 0], sigma_1[:, 1, 1], linewidth=1, label='sigma_22')
    #plt.plot(eps_2[:, 0, 0], sigma_2[:, 1, 1], linewidth=1, label='sigma_22')

    plt.xlabel('strain')
    plt.ylabel('stress(MPa)')
    plt.legend()

    plt.subplot(222)
    for i in range(0, 28):
        plt.plot(
            eps_1[:, 0, 0], w_1[:, i], 'r', linewidth=1.0,  label='cyclic', alpha=1)
        plt.plot(
            eps_1[:, 0, 0], w_2[:, i], linewidth=1.0, label='monotonic', alpha=1)

        plt.xlabel('strain')
        plt.ylabel('damage')
        # plt.legend()

    plt.subplot(223)
    for i in range(0, 28):
        plt.plot(
            eps_P_N_1[:, i], w_1[:, i], 'r', linewidth=1,  label='plastic_strain')
        plt.plot(
            eps_P_N_2[:, i], w_2[:, i], linewidth=1, label='plastic_strain')

        plt.xlabel('Strain')
        plt.ylabel('normal_strain')
        # plt.legend()

    plt.subplot(224)
    for i in range(0, 28):
        plt.plot(
            eps_Pi_T_1[:, i, 1], w_1[:, i],  'r', linewidth=1, label='sliding strain')
        plt.plot(
            eps_Pi_T_2[:, i, 1], w_2[:, i],  linewidth=1, label='sliding strain')

        plt.xlabel('sliding strain')
        plt.ylabel('damage')

    plt.show()
