'''
Created on 05.02.2017

@author: abaktheer

M7 Microplane model [F.Caner & Z.P.Bazant 2013] 
'''

from numpy.linalg import norm

from ibvpy.mats.mats3D.mats3D_eval import MATS3DEval
from ibvpy.mats.mats_eval import \
    IMATSEval
from numpy import \
    array, zeros, dot, trace, einsum, zeros_like,\
    identity, sign, linspace, hstack
from traits.api import Constant, implements,\
    Float, HasTraits, Property, cached_property
from traitsui.api import View, Include

import matplotlib.pyplot as plt
import numpy as np


class MATSEvalMicroplaneM7(HasTraits):
    #--------------------------
    # Elasticity material parameters
    #--------------------------
    E = Float(30000.0,
              label="E",
              desc="Young modulus",
              enter_set=True,
              auto_set=False)

    nu = Float(0.18,
               label="nu",
               desc="poission ratio",
               enter_set=True,
               auto_set=False)

    #---------------------------------------
    # material model parameters
    #---------------------------------------
    c1 = Float(0.089,
               label="c1",
               desc="",
               enter_set=True,
               auto_set=False)

    c2 = Float(0.176,
               label="c2",
               desc="",
               enter_set=True,
               auto_set=False)

    c3 = Float(4.0,
               label="c3",
               desc="",
               enter_set=True,
               auto_set=False)

    c4 = Float(50.0,
               label="c4",
               desc="",
               enter_set=True,
               auto_set=False)

    c5 = Float(3500.0,
               label="c5",
               desc="",
               enter_set=True,
               auto_set=False)

    c6 = Float(20.0,
               label="c6",
               desc="",
               enter_set=True,
               auto_set=False)

    c7 = Float(1.0,
               label="c7",
               desc="",
               enter_set=True,
               auto_set=False)

    c8 = Float(8.0,
               label="c8",
               desc="",
               enter_set=True,
               auto_set=False)

    c9 = Float(0.012,
               label="c9",
               desc="",
               enter_set=True,
               auto_set=False)

    c10 = Float(0.33,
                label="c10",
                desc="",
                enter_set=True,
                auto_set=False)

    c11 = Float(0.5,
                label="c11",
                desc="",
                enter_set=True,
                auto_set=False)

    c12 = Float(2.36,
                label="c12",
                desc="",
                enter_set=True,
                auto_set=False)

    c13 = Float(4500.0,
                label="c13",
                desc="",
                enter_set=True,
                auto_set=False)

    c14 = Float(300.0,
                label="c14",
                desc="",
                enter_set=True,
                auto_set=False)

    c15 = Float(4000.0,
                label="c15",
                desc="",
                enter_set=True,
                auto_set=False)

    c16 = Float(60.0,
                label="c16",
                desc="",
                enter_set=True,
                auto_set=False)

    c17 = Float(1.4,
                label="c17",
                desc="",
                enter_set=True,
                auto_set=False)

    c18 = Float(1.6e-3,
                label="c18",
                desc="",
                enter_set=True,
                auto_set=False)

    c19 = Float(1000.0,
                label="c19",
                desc="",
                enter_set=True,
                auto_set=False)

    c20 = Float(2.0,
                label="c20",
                desc="",
                enter_set=True,
                auto_set=False)

    c21 = Float(250.0,
                label="c21",
                desc="",
                enter_set=True,
                auto_set=False)

    K1 = Float(0.00010,
               label="K1",
               desc="",
               enter_set=True,
               auto_set=False)

    K2 = Float(110.0,
               label="K2",
               desc="",
               enter_set=True,
               auto_set=False)

    K3 = Float(20.0,
               label="K3",
               desc="",
               enter_set=True,
               auto_set=False)

    K4 = Float(40.0,
               label="K4",
               desc="",
               enter_set=True,
               auto_set=False)

    K5 = Float(0.0001,
               label="K5",
               desc="",
               enter_set=True,
               auto_set=False)

    E0 = Float(20000.0,
               label="E0",
               desc="",
               enter_set=True,
               auto_set=False)

    fc0 = Float(15.08,
                label="fc0",
                desc="",
                enter_set=True,
                auto_set=False)

    fc = Float(40.0,
               label="fc",
               desc="",
               enter_set=True,
               auto_set=False)

    a = Float(0.1,
              label="a",
              desc="",
              enter_set=True,
              auto_set=False)

    q = Float(2.0,
              label="a",
              desc="",
              enter_set=True,
              auto_set=False)


class MATSXDMicroplaneM7(MATSEvalMicroplaneM7):
    '''
    Microplane model M7.
    '''

    #-------------------------------------------------------------------------
    # MICROPLANE-Kinematic constraints
    #-------------------------------------------------------------------------
    def _get_e_Emna(self, eps_Emab):
        # Projection of apparent strain onto the individual microplanes
        e_ni = einsum('nb,...ba->...na', self._MPN, eps_Emab)
        return e_ni

    def _get_e_N_Emn(self, e_Emna):
        # get the normal strain array for each microplane
        e_N_Emn = einsum('...na, na->...n', e_Emna, self._MPN)
        return e_N_Emn

    def _get_e_T_Emna(self, e_Emna):
        # get the tangential strain vector array for each microplane
        e_N_Emn = self._get_e_N_Emn(e_Emna)
        e_N_Emna = einsum('...n,na->...na', e_N_Emn, self._MPN)
        return e_Emna - e_N_Emna

    #-------------------------------------------------
    # Alternative methods for the kinematic constraint
    #-------------------------------------------------
    # get the dyadic product of the microplane normals
    _MPNN = Property(depends_on='n_mp')

    @cached_property
    def _get__MPNN(self):
        MPNN_nij = einsum('ni,nj->nij', self._MPN, self._MPN)
        return MPNN_nij

    # get the third order tangential tensor (operator) for each microplane
    _MPTT = Property(depends_on='n_mp')

    @cached_property
    def _get__MPTT(self):
        delta = identity(2)
        MPTT_nijr = 0.5 * (einsum('ni,jr -> nijr', self._MPN, delta) +
                           einsum('nj,ir -> njir', self._MPN, delta) - 2 *
                           einsum('ni,nj,nr -> nijr', self._MPN, self._MPN, self._MPN))
        return MPTT_nijr

    def _get_e_N_Emn_2(self, eps_Emab):
        # Projection of apparent strain onto the individual microplanes
        return einsum('nij,...ij->...n', self._MPNN, eps_Emab)

    def _get_e_T_Emna_2(self, eps_Emab):
        # get the normal strain array for each microplane
        MPTT_ijr = self._get__MPTT()
        return einsum('nija,...ij->...na', MPTT_ijr, eps_Emab)

    def _get_e_Emna_2(self, eps_Emab):
        # get the tangential strain vector array for each microplane
        return self._get_e_N_Emn_2(eps_Emab) * self._MPN +\
            self._get_e_T_Emna_2(eps_Emab)

    #============================================================
    # normal stress - micro
    #============================================================
    def _get_normal_stress(self, eps_Emab, deps_Emab, sigma_v_Em, ev_0_Em, max_en_pos_Emn, max_en_neg_Emn, sigma_N0_Emn, ev_cum_Em):

        E_N0 = self.E / (1.0 - 2.0 * self.nu)

        ev_0_Em = einsum('...aa->...', 1. / 3. * (eps_Emab - deps_Emab))
        d_ev_Em = einsum('...aa->...', 1. / 3. * deps_Emab)
        ev_Em = ev_0_Em + d_ev_Em

        ev_cum_Em = ev_cum_Em + d_ev_Em

        en_0_Emn = self._get_e_N_Emn_2(eps_Emab - deps_Emab)
        d_en_Emn = self._get_e_N_Emn_2(deps_Emab)
        en_Emn = self._get_e_N_Emn_2(eps_Emab)

        ed_0_Emn = en_0_Emn - ev_0_Em
        d_ed_Emn = d_en_Emn - d_ev_Em
        ed_Emn = ed_0_Emn + d_ed_Emn

        eps_0_Emab = eps_Emab - deps_Emab

        eps_I_Em = einsum('...aa->...', eps_0_Emab)
        eps_III_Em = np.linalg.det(eps_0_Emab)

        eps_e_Em = np.minimum(self.c21, np.maximum(-sigma_v_Em, 0.0))

        alpha = (self.K5 / (1.0 + (eps_e_Em / E_N0))) * \
            ((eps_I_Em - eps_III_Em) / self.K1)**self.c20 + self.K4

        sigma_v_b_Em = -1.0 * self.E * self.K1 * self.K3 * \
            np.exp(-1.0 * ev_Em / (self.K1 * alpha))

#         print 'alpha', alpha
#         print 'sigma_v_b', sigma_v_b
#         print '************************'

        gamma_0 = self.fc0 / self.E - self.fc / self.E
        gamma_1 = np.exp(gamma_0) * np.tanh(self.c9 *
                                            np.maximum(0.0, -1.0 * ev_Em) / self.K1)
        beta_2 = self.c5 * gamma_1 + self.c7
        beta_3 = self.c6 * gamma_1 + self.c8

        sigma_d_b_Emn = -1.0 * (self.E * self.K1 * beta_3) / \
            (1.0 + (np.maximum(0.0, -1.0 * ed_Emn) / (self.K1 * beta_2))**2.0)
#         print 'sigma_d_b', sigma_d_b
#         print '************************'

        # save the maximum normal strain during loading history
        idx_1 = np.where(np.any(en_Emn) >= 0.0 and np.any(
            en_Emn) >= np.any(max_en_pos_Emn))
        max_en_pos[idx_1] = en_Emn[idx_1]
#         print 'max_en_pos', max_en_pos
#         print '************************'

        idx_2 = np.where(np.any(en_Emn) < 0.0 and np.any(
            en_Emn) < np.any(max_en_neg_Emn))
        print idx_2
        print idx_1
        max_en_neg[idx_2] = en_Emn[idx_2]
        # print 'max_en_neg', max_en_neg
        # print '************************'
        f_zi_Em = (1.0 + self.a * ev_cum_Em**2.0)**(-1.0)

        #  Fatigue extension
        f_zi_Em = 1.0 / (1.0 + self.a * ev_cum_Em**self.q +
                         (self.a**2.0) * ev_cum_Em**(2.0 * self.q))
        # print 'f_zi', f_zi
        # print '************************'
        eps_e_2_Em = max(-1.0 * sigma_v_Em / E_N0, 0)
#         print 'eps_e_2', eps_e_2
#         print '************************'
        idx_3 = np.where(sigma_N0_Emn >= 0.0)
        EN_Emn = np.zeros_like(sigma_N0_Emn)
        EN_Emn[idx_3] = E_N0 * \
            np.exp(-1.0 * self.c13 * max_en_pos_Emn * f_zi_Em)

        idx_4 = np.where(sigma_N0_Emn > E_N0 *
                         en_Emn and sigma_N0_Emn * d_en_Emn < 0.0)
        EN_Emn[idx_4] = E_N0

        idx_5 = np.where(sigma_N0_Emn < 0.0)

        EN_Emn[idx_5] = E_N0 * (np.exp(-1.0 * self.c14 * np.abs(max_en_neg_Emn) /
                                       (1.0 + self.c15 * eps_e_2_Em)) + self.c16 * eps_e_2_Em)
#             print '333333333333333'
#         print 'EN', EN
#         print '************************'
        sigma_N_e_Emn = sigma_N0_Emn + EN_Emn * d_en_Emn
#         print 'sigma_N_e', sigma_N_e

        beta_1_Em = -1.0 * self.c1 + self.c17 * \
            np.exp(-1.0 * self.c19 * max(0.0, -1.0 * sigma_v_Em - self.c18) / E_N0)

        sigma_N_b_Emn = self.E * self.K1 * beta_1_Em * np.exp(-1.0 * max(0.0, en_Emn - beta_1_Em * self.c2 * self.K1) /
                                                              (self.c4 * eps_e_2_Em + self.K1 * self.c3))

#         print 'sigma_N_b', sigma_N_b
#         print '************************'

        sigma_N_Emn = np.maximum(
            np.minimum(sigma_N_e_Emn, sigma_N_b_Emn), sigma_v_b_Em + sigma_d_b_Emn)

#         print 'sigma_N', sigma_N
#         print '************************'

        return sigma_N_Emn, max_en_pos_Emn, max_en_neg_Emn, ev_cum_Em, EN_Emn, ev_Em

    #============================================================
    # tangential stress - micro
    #============================================================
    def _get_tangential_stress(self, eps_Emab, deps_Emab, sigma_N_Emn, sigma_T0_Emn, EN_Emn, ev_Em):

        eT_0 = np.linalg.norm(self._get_e_T_vct_arr_2(eps_Emab - deps_Emab))
        d_eT = np.linalg.norm(self._get_e_T_vct_arr_2(deps_Emab))
        eT = eT_0 + d_eT

        # print 'd_eT', d_eT.shape

        ET_Emn = EN_Emn * (1.0 - 4.0 * self.nu) / (1.0 + self.nu)

        sigma_N_0_hat_Emn = ET_Emn * np.maximum(
            0.0,  self.K1 * self.c11 - self.c21 * np.maximum(0.0, ev_Em))

        # print 'sigma_N_0_hat', sigma_N_0_hat
        idx_1 = np.where(sigma_N_Emn <= 0.0)
        sigma_T_b_Emn = np.zeros_like(sigma_N_Emn)
        sigma_T_b_Emn[idx_1] = 1.0 / \
            (1.0 / (self.c10 * np.maximum(0.0, sigma_N_0_hat_Emn - sigma_N_Emn)) +
             1.0 / (ET_Emn * self.K1 * self.K2))

        idx_2 = np.where(sigma_N_Emn > 0.0)
        sigma_T_b_Emn[idx_2] = 1.0 / \
            (1.0 / (self.c10 * sigma_N_0_hat_Emn) +
             1.0 / (ET_Emn * self.K1 * self.K2))

        sigma_T_e_Emn = sigma_T0_Emn + ET_Emn * d_eT

        # print 'sigma_T_b ', sigma_T_b
        # print 'sigma_T_e ', sigma_T_e

        sigma_T_Emn = np.maximum(sigma_T_b_Emn, np.abs(sigma_T_e_Emn))

        return sigma_T_Emn

    #============================================================
    # stress tensor - macro
    #============================================================
    def _get_stress_tns(self, eps_Emab, sigma_N_Emn, sigma_T_Emna):

        eT_Emna = self._get_e_T_vct_arr_2(eps_Emab)
        # print 'eT ', eT.shape
        # print 'sigma_T ', sigma_T.shape
        sigma_T_Emna = einsum('...n,...na -> ...na', sigma_T, eT_Emna)

        # print 'sigma_T ', sigma_T.shape

        MPNN_nij = self._get__MPNN()
        MPTT_nijk = self._get__MPTT()
        MPW_n = self._MPW

#         sigma_N = zeros(28)
#         sigma_T = zeros(( 28, 3))
#
#         for i in range(1,28):
#             sigma_N [i]= self._get

        sigma_Emij = einsum('n,nij,...n -> ...ij', MPW_n, MPNN_nij, sigma_N_Emn) + \
            einsum('n,nkij,...nk -> ...ij', MPW_n, MPTT_nijk, sigma_T_Emna)

        # print 'sigma ', sigma_ij

        return sigma_Emij

#     #-------------------------------------------------------------------------
#     # Evaluation - get the corrector and predictor
#     #-------------------------------------------------------------------------
#     def get_corr_pred(self, sctx, eps_app_eng, sigma_kk):
#
#         #----------------------------------------------------------------------
#         # Return stresses (corrector) and damaged secant stiffness matrix (predictor)
#         #----------------------------------------------------------------------
#
#         # secant stiffness tensor
#         S_ijkl = self._get_S_4_tns(sctx, eps_app_eng, sigma_kk)
#
#         # plastic strain tensor
#         eps_p_ij = self._get_eps_p_mtx(sctx, eps_app_eng, sigma_kk)
#
#         # elastic strain tensor
#         eps_e_mtx = eps_app_eng - eps_p_ij
#
#         # calculation of the stress tensor
#         sig_eng = einsum('ijmn,mn -> ij', S_ijkl, eps_e_mtx)
#
#         return sig_eng, S_ijkl


class MATS3DMicroplaneM7(MATSXDMicroplaneM7, MATS3DEval):

    implements(IMATSEval)

    #-----------------------------------------------
    # number of microplanes - currently fixed for 3D
    #-----------------------------------------------
    n_mp = Constant(28)

    #-----------------------------------------------
    # get the normal vectors of the microplanes
    #-----------------------------------------------
    _MPN = Property(depends_on='n_mp')

    @cached_property
    def _get__MPN(self):
        # microplane normals:
        return array([[.577350259, .577350259, .577350259],
                      [.577350259, .577350259, -.577350259],
                      [.577350259, -.577350259, .577350259],
                      [.577350259, -.577350259, -.577350259],
                      [.935113132, .250562787, .250562787],
                      [.935113132, .250562787, -.250562787],
                      [.935113132, -.250562787, .250562787],
                      [.935113132, -.250562787, -.250562787],
                      [.250562787, .935113132, .250562787],
                      [.250562787, .935113132, -.250562787],
                      [.250562787, -.935113132, .250562787],
                      [.250562787, -.935113132, -.250562787],
                      [.250562787, .250562787, .935113132],
                      [.250562787, .250562787, -.935113132],
                      [.250562787, -.250562787, .935113132],
                      [.250562787, -.250562787, -.935113132],
                      [.186156720, .694746614, .694746614],
                      [.186156720, .694746614, -.694746614],
                      [.186156720, -.694746614, .694746614],
                      [.186156720, -.694746614, -.694746614],
                      [.694746614, .186156720, .694746614],
                      [.694746614, .186156720, -.694746614],
                      [.694746614, -.186156720, .694746614],
                      [.694746614, -.186156720, -.694746614],
                      [.694746614, .694746614, .186156720],
                      [.694746614, .694746614, -.186156720],
                      [.694746614, -.694746614, .186156720],
                      [.694746614, -.694746614, -.186156720]])

    #-------------------------------------
    # get the weights of the microplanes
    #-------------------------------------
    _MPW = Property(depends_on='n_mp')

    @cached_property
    def _get__MPW(self):
        # Note that the values in the array must be multiplied by 6 (cf. [Baz05])!
        # The sum of of the array equals 0.5. (cf. [BazLuz04]))
        # The values are given for an Gaussian integration over the unit
        # hemisphere.
        return array([.0160714276, .0160714276, .0160714276, .0160714276, .0204744730,
                      .0204744730, .0204744730, .0204744730, .0204744730, .0204744730,
                      .0204744730, .0204744730, .0204744730, .0204744730, .0204744730,
                      .0204744730, .0158350505, .0158350505, .0158350505, .0158350505,
                      .0158350505, .0158350505, .0158350505, .0158350505, .0158350505,
                      .0158350505, .0158350505, .0158350505]) * 6.0

    #-------------------------------------------------------------------------
    # Cached elasticity tensors
    #-------------------------------------------------------------------------

    @cached_property
    def _get_elasticity_tensors(self):
        '''
        Intialize the fourth order elasticity tensor for 3D or 2D plane strain or 2D plane stress
        '''
        # ----------------------------------------------------------------------------
        # Lame constants calculated from E and nu
        # ----------------------------------------------------------------------------

        # first Lame paramter
        la = self.E * self.nu / ((1 + self.nu) * (1 - 2 * self.nu))
        # second Lame parameter (shear modulus)
        mu = self.E / (2 + 2 * self.nu)

        # -----------------------------------------------------------------------------------------------------
        # Get the fourth order elasticity and compliance tensors for the 3D-case
        # -----------------------------------------------------------------------------------------------------

        # construct the elasticity tensor (using Numpy - einsum function)
        delta = identity(3)
        D_ijkl = (einsum(',ij,kl->ijkl', la, delta, delta) +
                  einsum(',ik,jl->ijkl', mu, delta, delta) +
                  einsum(',il,jk->ijkl', mu, delta, delta))

        return D_ijkl

    #-------------------------------------------------------------------------
    # Dock-based view with its own id
    #-------------------------------------------------------------------------

    traits_view = View(Include('polar_fn_group'),
                       dock='tab',
                       id='ibvpy.mats.mats3D.mats_3D_cmdm.MATS3D_cmdm',
                       kind='modal',
                       resizable=True,
                       scrollable=True,
                       width=0.6, height=0.8,
                       buttons=['OK', 'Cancel']
                       )


if __name__ == '__main__':
    #==========================
    # Check the model behavior
    #==========================

    # monotonic loading
    n = 10
    s_levels = linspace(0, -0.01, 2)
    #s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] = 0.0001
    s_levels[0] = 0
    s_history_1 = s_levels.flatten()
    s_arr_1 = hstack([linspace(s_history_1[i], s_history_1[i + 1], n)
                      for i in range(len(s_history_1) - 1)])

    eps_1 = array([array([[s_arr_1[i], 0, 0],
                          [0, 0, 0],
                          [0, 0, 0]]) for i in range(0, len(s_arr_1))])

    m = MATS3DMicroplaneM7()

    # print 'eps_1', eps_1.shape

    # construct the arrays
    sigma = zeros_like(eps_1)
    sigma_kk = zeros(len(s_arr_1) + 1)

    ev = zeros((len(eps_1[:, 0, 0]) + 1, 28))
    sigma_N = zeros((len(eps_1[:, 0, 0]) + 1, 28))
    max_en_pos = zeros((len(eps_1[:, 0, 0]) + 1, 28))
    max_en_neg = zeros((len(eps_1[:, 0, 0]) + 1, 28))
    ev_cum = zeros((len(eps_1[:, 0, 0]) + 1, 28))
    EN = zeros((len(eps_1[:, 0, 0]) + 1, 28))

    sigma_T = zeros((len(eps_1[:, 0, 0]) + 1, 28))
    #sigma_T = zeros((len(eps_1[:, 0, 0]), 28, 3))
    #w_1_N = zeros((len(eps_1[:, 0, 0]), 28))
    #w_1_T = zeros((len(eps_1[:, 0, 0]), 28))
    #eps_P_N_1 = zeros((len(eps_1[:, 0, 0]), 28))
    #eps_Pi_T_1 = zeros((len(eps_1[:, 0, 0]), 28, 3))
    #sctx_1 = zeros((len(eps_1[:, 0, 0]) + 1, 28, 13))

    # monotonic loading
    for i in range(1, len(eps_1[:, 0, 0])):
        print '-----------------------------------------------------------------> i', i
        sigma_N_0, max_en_pos_0, max_en_neg_0, ev_cum_0,\
            EN_0, ev_0 = m._get_normal_stress(eps_1[i, :],
                                              eps_1[i, :] -
                                              eps_1[i - 1, :],
                                              sigma_kk[i] / 3.0,
                                              ev[i],
                                              max_en_pos[i, :],
                                              max_en_neg[i, :],
                                              sigma_N[i, :],
                                              ev_cum[i])
        # print 'ev', ev
        # print 'sigma_N', sigma_N

        ev[i + 1] = ev_0
        max_en_pos[i + 1] = max_en_pos_0
        max_en_neg[i + 1] = max_en_neg_0
        ev_cum[i + 1] = ev_cum_0
        EN[i + 1] = EN_0
        sigma_N[i + 1] = sigma_N_0

        sigma_T_0 = m._get_tangential_stress(eps_1[i, :],
                                             eps_1[i, :] - eps_1[i - 1, :],
                                             sigma_N[
                                                 i + 1], sigma_T[i],
                                             EN[i + 1], ev[i + 1])
        # print 'sigma_T_0', sigma_T_0.shape
        sigma_T[i + 1] = sigma_T_0

        sigma[i, :] = m._get_stress_tns(
            eps_1[i, :], sigma_N[i, :], sigma_T[i, :])

        sigma_kk[i + 1] = trace(sigma[i, :])

    #======================
    # plotting
    #======================
    # stress -strain
    plt.subplot(221)
    plt.plot(eps_1[1:, 0, 0], sigma[1:, 0, 0], color='k',
             linewidth=1, label='sigma_11_(monotonic)')

    plt.title('$\sigma - \epsilon$')
    plt.xlabel('Strain')
    plt.ylabel('Stress(MPa)')
    plt.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    plt.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.legend()

    plt.subplot(222)
    plt.plot(eps_1[:, 0, 0], EN[1:, 0], color='k',
             linewidth=1, label='sigma_11_(monotonic)')

    plt.title('$\sigma - \epsilon$')
    plt.xlabel('Strain')
    plt.ylabel('Stress(MPa)')
    plt.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    plt.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.legend()
    plt.show()
