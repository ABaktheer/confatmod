'''
Created on 30.03.2017

@author: abaktheer
'''

from ibvpy.core.rtrace_eval import \
    RTraceEval
from ibvpy.mats.mats3D.mats3D_eval import MATS3DEval
from ibvpy.mats.mats_eval import \
    IMATSEval
from numpy import \
    array, zeros, outer, inner, transpose, dot, trace, \
    fabs, identity, tensordot, einsum, zeros_like,\
    float_, identity, sign, fabs, linspace, hstack, \
    sqrt, copy
from numpy.linalg import norm
from scipy.linalg import \
    eigh, inv
from traits.api import \
    Constant,  Property, cached_property, implements,\
    Bool, Callable, Enum, Float, HasTraits, \
    Int, Trait, on_trait_change, \
    Dict, Property, cached_property
from traitsui.api import \
    Item, View, Group, Spring, Include

import matplotlib.pyplot as plt
import numpy as np


class MATSEvalMicroplaneFatigue(HasTraits):

     #--------------------------
    # material model parameters
    #--------------------------

    E = Float(34000.,
              label="E",
              desc="Young modulus",
              enter_set=True,
              auto_set=False)

    nu = Float(0.2,
               label="nu",
               desc="poission ratio",
               enter_set=True,
               auto_set=False)

    #---------------------------------------
    # Tangential constitutive law parameters
    #---------------------------------------
    gamma_T = Float(5000.,
                    label="Gamma",
                    desc=" Tangential Kinematic hardening modulus",
                    enter_set=True,
                    auto_set=False)

    K_T = Float(5000.0,
                label="K",
                desc="Tangential Isotropic harening",
                enter_set=True,
                auto_set=False)

    S = Float(0.00001,
              label="S",
              desc="Damage strength",
              enter_set=True,
              auto_set=False)

    r = Float(1.2,
              label="r",
              desc="Damage cumulation parameter",
              enter_set=True,
              auto_set=False)

    c = Float(1.0,
              label="c",
              desc="Damage cumulation parameter",
              enter_set=True,
              auto_set=False)

    tau_pi_bar = Float(10.0,
                       label="Tau_bar",
                       desc="Reversibility limit",
                       enter_set=True,
                       auto_set=False)

    a = Float(0.0,
              label="a",
              desc="Lateral pressure coefficient",
              enter_set=True,
              auto_set=False)

    #-------------------------------------------
    # Normal_Tension constitutive law parameters
    #-------------------------------------------
    Ad = Float(10000.0,
               label="a",
               desc="brittleness coefficient",
               enter_set=True,
               auto_set=False)

    eps_0 = Float(0.5e-4,
                  label="a",
                  desc="threshold strain",
                  enter_set=True,
                  auto_set=False)

    #-----------------------------------------------
    # Normal_Compression constitutive law parameters
    #-----------------------------------------------
    K_N = Float(200.,
                label="K_N",
                desc=" Normal isotropic harening",
                enter_set=True,
                auto_set=False)

    gamma_N = Float(200.,
                    label="gamma_N",
                    desc="Normal kinematic hardening",
                    enter_set=True,
                    auto_set=False)

    sigma_0 = Float(20.,
                    label="sigma_0",
                    desc="Yielding stress",
                    enter_set=True,
                    auto_set=False)

    #--------------------------------------------------------------
    # Heviside function for the unilateral effect (CP - TD)
    #--------------------------------------------------------------
    def get_heviside(self, eps):
        if eps > 0:
            return 1.0
        else:
            return 0.0

    #--------------------------------------------------------------
    # microplane constitutive law (normal behavior CP + TD)
    #--------------------------------------------------------------
    def get_normal_Law(self, eps, sctx):
        E_N = self.E / (1. - 2. * self.nu)
        w_N = sctx[0]
        z_N = sctx[1]
        alpha_N = sctx[2]
        r_N = sctx[3]
        eps_N_p = sctx[4]
        sigma_N = sctx[5]

        H = self.get_heviside(sigma_N)

        sigma_n_trial = (1 - H * w_N) * E_N * (eps - eps_N_p)

        Z = self.K_N * r_N
        X = self.gamma_N * alpha_N
        h = max(0., (self.sigma_0 + Z))
        f_trial = abs(sigma_n_trial - X) - h

        if f_trial > 1e-6:
            delta_lamda = f_trial / (E_N + abs(self.K_N) + self.gamma_N)
            eps_N_p = eps_N_p + delta_lamda * sign(sigma_n_trial - X)
            r_N = r_N + delta_lamda
            alpha_N = alpha_N + delta_lamda * sign(sigma_n_trial - X)

        def Z_N(z_N): return 1. / self.Ad * (-z_N) / (1 + z_N)
        Y_N = 0.5 * H * E_N * eps ** 2
        Y_0 = 0.5 * E_N * self.eps_0 ** 2
        f = Y_N - (Y_0 + Z_N(z_N))

        if f > 1e-6:
            def f_w(Y): return 1 - 1. / (1 + self.Ad * (Y - Y_0))

            w_N = f_w(Y_N)
            z_N = - w_N

        sigma_N = (1 - H * w_N) * E_N * (eps - eps_N_p)

        new_sctx = zeros(6)

        new_sctx[0] = w_N
        new_sctx[1] = z_N
        new_sctx[2] = alpha_N
        new_sctx[3] = r_N
        new_sctx[4] = eps_N_p
        new_sctx[5] = sigma_N
        return new_sctx

    #-------------------------------------------------------------------------
    # microplane constitutive law (Tangential CSD)-(Pressure sensitive cumulative damage)
    #-------------------------------------------------------------------------
    def get_tangential_Law(self, e_T, sctx, sigma_kk):

        G = self.E / (1 + 2.0 * self.nu)

        w_T = sctx[6]
        z_T = sctx[7]
        alpha_T = sctx[8:11]
        eps_T_pi = sctx[11:14]
        sig_T = sctx[14:17]

        sig_pi_trial = G * (e_T - eps_T_pi)
        Z = self.K_T * z_T
        X = self.gamma_T * alpha_T
        f = norm(sig_pi_trial - X) - self.tau_pi_bar - \
            Z + self.a * sigma_kk / 3

        if f > 1e-6:
            delta_lamda = f / \
                (G / (1 - w_T) + self.gamma_T + self.K_T)
            eps_T_pi = eps_T_pi + delta_lamda * \
                ((sig_pi_trial - X) / (1 - w_T)) / norm(sig_pi_trial - X)
            Y = 0.5 * G * dot((e_T - eps_T_pi), (e_T - eps_T_pi))
            w_T += ((1 - w_T) ** self.c) * \
                (delta_lamda * (Y / self.S) ** self.r)

            alpha_T = alpha_T + delta_lamda * \
                (sig_pi_trial - X) / norm(sig_pi_trial - X)
            z_T = z_T + delta_lamda

        sig_T = (1 - w_T) * G * (e_T - eps_T_pi)

        new_sctx = zeros(11)
        new_sctx[0:2] = w_T, z_T
        new_sctx[2:5] = alpha_T
        new_sctx[5:8] = eps_T_pi
        new_sctx[8:11] = sig_T
        return new_sctx

#     def get_normal_compression_Law(self, eps, sctx):
#         # microplane constitutive law (normal-compression)-(Plasticity)
#         E_N = self.E / (1. - 2. * self.nu)
#
#         alpha_N = sctx[3]
#         r_N = sctx[4]
#         eps_N_p = sctx[5]
#         #sigma_N = sctx[6]
#
#         sigma_n_trial = E_N * (eps - eps_N_p)
#         Z = self.K * r_N
#         X = self.gamma_N * alpha_N
#
#         f_trial = abs(sigma_n_trial - X) - self.sigma_0 - Z
#
#         if f_trial > 1e-6:
#             delta_lamda = f_trial / (E_N + self.K + self.gamma_N)
#             eps_N_p = eps_N_p + delta_lamda * sign(sigma_n_trial - X)
#             r_N = r_N + delta_lamda
#             alpha_N = alpha_N + delta_lamda * sign(sigma_n_trial - X)
#
#         sigma_N = E_N * (eps - eps_N_p)
#
# #         print 'eps',  eps
# #         print 'eps_N_p',  eps_N_p
#         # print 'sigma_N',  sigma_N
#
#         new_sctx = zeros(4)
#         new_sctx[0] = alpha_N
#         new_sctx[1] = r_N
#         new_sctx[2] = eps_N_p
#         new_sctx[3] = sigma_N
#         return new_sctx
#
#     def get_tangential_Law(self, e_T, sctx, sigma_kk):
#         # microplane constitutive law (Tangential)-(Pressure sensitive
#         # cumulative damage)
#         w_T = sctx[7]
#         z_T = sctx[8]
#         alpha_T = sctx[9:12]
#         eps_T_pi = sctx[12:15]
#         #sig_T_pi = sctx[15:18]
#
#         K0 = self.E / 3. * (1. - 2. * self.nu)
#         G0 = self.E / 2. * (1. + self.nu)
#
# #         E_T =  (10. / 3.) * G0 - 2 * K0
#         E_T = self.E * (1. - 4. * self.nu) / \
#             ((1. + self.nu) * (1. - 2. * self.nu))
#
#         sig_pi_trial = E_T * (e_T - eps_T_pi)
#         Z = self.K * z_T
#         X = self.gamma * alpha_T
#         f = norm(sig_pi_trial - X) - self.tau_pi_bar - \
#             Z + self.a * sigma_kk / 3
#
#         if f > 1e-6:
#             delta_lamda = f / \
#                 (E_T / (1 - w_T) + self.gamma + self.K)
#             eps_T_pi = eps_T_pi + delta_lamda * \
#                 ((sig_pi_trial - X) / (1 - w_T)) / norm(sig_pi_trial - X)
#             Y = 0.5 * E_T * dot((e_T - eps_T_pi), (e_T - eps_T_pi))
#             w_T += ((1 - w_T) ** self.c) * \
#                 (delta_lamda * (Y / self.S) ** self.r)
#             # print 'w', w
#             X = X + self.gamma * delta_lamda * \
#                 (sig_pi_trial - X) / norm(sig_pi_trial - X)
#             alpha_T = alpha_T + delta_lamda * \
#                 (sig_pi_trial - X) / norm(sig_pi_trial - X)
#             z_T = z_T + delta_lamda
#
#         sig_T_pi = (1 - w_T) * E_T * (e_T - eps_T_pi)
#
#         # print 'sig_T_pi',  sig_T_pi
#
#         new_sctx = zeros(11)
#         new_sctx[0:2] = w_T, z_T
#         new_sctx[2:5] = alpha_T
#         new_sctx[5:8] = eps_T_pi
#         new_sctx[8:11] = sig_T_pi
#
#         return new_sctx


class MATSXDMicroplaneDamageFatigue(MATSEvalMicroplaneFatigue):

    '''
    Microplane Damage Fatigue Model.
    '''
    #-------------------------------------------------------------------------
    # Classification traits (supplied by the dimensional subclasses)
    #-------------------------------------------------------------------------

    # specification of the model dimension (2D, 3D)
    n_dim = Int

    # specification of number of engineering strain and stress components
    n_eng = Int

    #-------------------------------------------------------------------------
    # Configuration parameters
    #-------------------------------------------------------------------------

    model_version = Enum("compliance", "stiffness")

    symmetrization = Enum("product-type", "sum-type")

    regularization = Bool(False,
                          desc='Flag to use the element length projection'
                          ' in the direction of principle strains',
                          enter_set=True,
                          auto_set=False)

    elastic_debug = Bool(False,
                         desc='Switch to elastic behavior - used for debugging',
                         auto_set=False)

    double_constraint = Bool(False,
                             desc='Use double constraint to evaluate microplane elastic and fracture energy (Option effects only the response tracers)',
                             auto_set=False)

    #-------------------------------------------------------------------------
    # View specification
    #-------------------------------------------------------------------------

    config_param_vgroup = Group(Item('model_version', style='custom'),
                                #     Item('stress_state', style='custom'),
                                Item('symmetrization', style='custom'),
                                Item('elastic_debug@'),
                                Item('double_constraint@'),
                                Spring(resizable=True),
                                label='Configuration parameters',
                                show_border=True,
                                dock='tab',
                                id='ibvpy.mats.matsXD.MATSXD_cmdm.config',
                                )

    traits_view = View(Include('polar_fn_group'),
                       dock='tab',
                       id='ibvpy.mats.matsXD.MATSXD_cmdm',
                       kind='modal',
                       resizable=True,
                       scrollable=True,
                       width=0.6, height=0.8,
                       buttons=['OK', 'Cancel']
                       )

    #-------------------------------------------------------------------------
    # Setup for computation within a supplied spatial context
    #-------------------------------------------------------------------------

    D4_e = Property

    def _get_D4_e(self):
        # Return the elasticity tensor
        return self.elasticity_tensors[0]

    #-------------------------------------------------------------------------
    # MICROPLANE-PROJECTION (Kinematic constraint)
    #-------------------------------------------------------------------------
    # get the dyadic product of the microplane normals
    _MPNN = Property(depends_on='n_mp')

    @cached_property
    def _get__MPNN(self):
        # dyadic product of the microplane normals
        # return array([outer(mpn, mpn) for mpn in self._MPN]) # old
        # implementation
        # n identify the microplane
        MPNN_nij = einsum('ni,nj->nij', self._MPN, self._MPN)
        return MPNN_nij

    # get the third order tangential tensor (operator) for each microplane

    _MPTT = Property(depends_on='n_mp')

    @cached_property
    def _get__MPTT(self):
        # Third order tangential tensor for each microplane
        delta = identity(3)
        MPTT_nijr = 0.5 * (einsum('ni,jr -> nijr', self._MPN, delta) +
                           einsum('nj,ir -> njir', self._MPN, delta) - 2 *
                           einsum('ni,nj,nr -> nijr', self._MPN, self._MPN, self._MPN))
        return MPTT_nijr

    def _get_e_vct_arr(self, eps_eng):
        # Projection of apparent strain onto the individual microplanes
        e_ni = einsum('nj,ji->ni', self._MPN, eps_eng)
        return e_ni

    def _get_e_N_arr(self, e_vct_arr):
        # get the normal strain array for each microplane
        eN_n = einsum('ni,ni->n', e_vct_arr, self._MPN)
        return eN_n

    def _get_e_T_vct_arr(self, e_vct_arr):
        # get the tangential strain vector array for each microplane
        eN_n = self._get_e_N_arr(e_vct_arr)

        eN_vct_ni = einsum('n,ni->ni', eN_n, self._MPN)

        return e_vct_arr - eN_vct_ni

    #-------------------------------------------------
    # Alternative methods for the kinematic constraint
    #-------------------------------------------------

    def _get_e_N_arr_2(self, eps_eng):

        #eps_mtx = self.map_eps_eng_to_mtx(eps_eng)
        return einsum('nij,ij->n', self._MPNN, eps_eng)

    def _get_e_T_vct_arr_2(self, eps_eng):

        #eps_mtx = self.map_eps_eng_to_mtx(eps_eng)
        MPTT_ijr = self._get__MPTT()
        return einsum('nijr,ij->nr', MPTT_ijr, eps_eng)

    def _get_e_vct_arr_2(self, eps_eng):

        return self._e_N_arr_2 * self._MPN + self._e_t_vct_arr_2

    def _get_state_variables(self, sctx, eps_app_eng, sigma_kk):

        #e_vct_arr = self._get_e_vct_arr(eps_app_eng)
        e_N_arr = self._get_e_N_arr_2(eps_app_eng)
        e_T_vct_arr = self._get_e_T_vct_arr_2(eps_app_eng)

        sctx_arr = zeros((28, 18))

        for i in range(0, self.n_mp):

            sctx_N = self.get_normal_Law(
                e_N_arr[i], sctx[i, :])

            sctx_tangential = self.get_tangential_Law(
                e_T_vct_arr[i, :], sctx[i, :], sigma_kk)

            sctx_arr[i, 0:6] = sctx_N
            sctx_arr[i, 6:17] = sctx_tangential

        return sctx_arr

    def _get_stress_tns(self, sctx, eps_app_eng, sigma_kk):

        sigma_N_n = self._get_state_variables(
            sctx, eps_app_eng, sigma_kk)[:, 5]

        # print 'sigma_N_n', sigma_N_n

        sigma_T_ni = self._get_state_variables(
            sctx, eps_app_eng, sigma_kk)[:, 14:17]

        # print 'sigma_T_ni', sigma_T_ni

        MPNN_nij = self._get__MPNN()
        MPTT_nijk = self._get__MPTT()
        MPW_n = self._MPW

        sigma_ij = einsum('n,nij,n -> ij', MPW_n, MPNN_nij, sigma_N_n) + \
            einsum('n,nkij,nk -> ij', MPW_n, MPTT_nijk, sigma_T_ni)

        # print 'sigma ', sigma_ij

        return sigma_ij

    def _get_eps_N_p_arr(self, sctx, eps_app_eng, sigma_kk):
        '''
        Returns a list of the plastic normal strain  for all microplanes.

        '''
        eps_N_p = self._get_state_variables(sctx, eps_app_eng, sigma_kk)[:, 4]

        return eps_N_p

    def _get_eps_T_pi_arr(self, sctx, eps_app_eng, sigma_kk):
        '''
        Returns a list of the sliding strain vector for all microplanes.

        '''

        eps_T_pi_vct_arr = self._get_state_variables(
            sctx, eps_app_eng, sigma_kk)[:, 11:14]

        return eps_T_pi_vct_arr

    def _get_beta_N_arr(self, sctx, eps_app_eng, sigma_kk):

        # Returns a list of the normal integrity factors for all microplanes.

        beta_N_arr = 1 - \
            self._get_state_variables(sctx, eps_app_eng, sigma_kk)[:, 0]

        return beta_N_arr

    def _get_beta_T_arr(self, sctx, eps_app_eng, sigma_kk):

        # Returns a list of the tangential integrity factors for all
        # microplanes.

        beta_T_arr = 1 - \
            self._get_state_variables(sctx, eps_app_eng, sigma_kk)[:, 6]

        return beta_T_arr

    def _get_beta_tns_2(self, sctx, eps_app_eng, sigma_kk):

        # Returns the 4th order damage tensor 'beta4' using
        #(cf. [Baz99], Eq.(63))

        delta = identity(3)
        beta_N_n = self._get_beta_N_arr(sctx, eps_app_eng, sigma_kk)
        beta_T_n = self._get_beta_T_arr(sctx, eps_app_eng, sigma_kk)

        beta_ijkl = einsum('n, n, ni, nj, nk, nl -> ijkl', self._MPW, beta_N_n, self._MPN, self._MPN, self._MPN, self._MPN) + \
            0.25 * (einsum('n, n, ni, nk, jl -> ijkl', self._MPW, beta_T_n, self._MPN, self._MPN, delta) +
                    einsum('n, n, ni, nl, jk -> ijkl', self._MPW, beta_T_n, self._MPN, self._MPN, delta) +
                    einsum('n, n, nj, nk, il -> ijkl', self._MPW, beta_T_n, self._MPN, self._MPN, delta) +
                    einsum('n, n, nj, nl, ik -> ijkl', self._MPW, beta_T_n, self._MPN, self._MPN, delta) -
                    4 * einsum('n, n, ni, nj, nk, nl -> ijkl', self._MPW, beta_T_n, self._MPN, self._MPN, self._MPN, self._MPN))

        # print 'beta_ijkl =', beta_ijkl
        return beta_ijkl

    def _get_phi_arr(self, sctx, eps_app_eng, sigma_kk):
        # Returns a list of the integrity factors for all microplanes.

        w_n = self._get_state_variables(sctx, eps_app_eng, sigma_kk)[:, 0]
        w_T = self._get_state_variables(sctx, eps_app_eng, sigma_kk)[:, 6]

        w = np.zeros(self.n_mp)

        for i in range(0, self.n_mp):
            w[i] = np.maximum(w_n[i], w_T[i])

        phi_arr = 1. - w

#         phi_arr = sqrt(1 -
# self._get_state_variables(sctx, eps_app_eng, sigma_kk)[:, 7])

        return phi_arr

    def _get_phi_mtx(self, sctx, eps_app_eng, sigma_kk):
        # Returns the 2nd order damage tensor 'phi_mtx'

        # scalar integrity factor for each microplane
        phi_arr = self._get_phi_arr(sctx, eps_app_eng, sigma_kk)
        # integration terms for each microplanes

        phi_ij = einsum('n,n,nij->ij', phi_arr, self._MPW, self._MPNN)

        return phi_ij

    def _get_beta_tns(self, sctx, eps_app_eng, sigma_kk):
        '''
        Returns the 4th order damage tensor 'beta4' using sum-type symmetrization
        (cf. [Jir99], Eq.(21))
        '''
        delta = identity(3)

        phi_mtx = self._get_phi_mtx(sctx, eps_app_eng, sigma_kk)

        # use numpy functionality (einsum) to evaluate [Jir99], Eq.(21)
        beta_ijkl = 0.25 * (einsum('ik,jl->ijkl', phi_mtx, delta) +
                            einsum('il,jk->ijkl', phi_mtx, delta) +
                            einsum('jk,il->ijkl', phi_mtx, delta) +
                            einsum('jl,ik->ijkl', phi_mtx, delta))

        # print 'beta_ijkl', beta_ijkl

        return beta_ijkl

    def _get_eps_pi_mtx(self, sctx, eps_app_eng, sigma_kk):

        # Integration of the (inelastic) strains for each microplane
        eps_N_P_n = self._get_eps_N_p_arr(sctx, eps_app_eng, sigma_kk)
        eps_T_pi_ni = self._get_eps_T_pi_arr(sctx, eps_app_eng, sigma_kk)
        delta = identity(3)

        # print eps_N_P_n.shape

        eps_pi_ij = einsum('n,n,ni,nj -> ij', self._MPW, eps_N_P_n, self._MPN, self._MPN) + \
            0.5 * (einsum('n,nr,ni,rj->ij', self._MPW, eps_T_pi_ni, self._MPN, delta) +
                   einsum('n,nr,nj,ri->ij', self._MPW, eps_T_pi_ni, self._MPN, delta))

        return eps_pi_ij

    #-------------------------------------------------------------------------
    # Evaluation - get the corrector and predictor
    #-------------------------------------------------------------------------

    def get_corr_pred(self, sctx, eps_app_eng, sigma_kk):
        '''
        Corrector predictor computation.
        @param eps_app_eng input variable - engineering strain
        '''
        # -----------------------------------------------------------------------------------------------
        # check if average strain is to be used for damage evaluation
        # -----------------------------------------------------------------------------------------------
        # if eps_avg != None:
        #    pass
        # else:
        #    eps_avg = eps_app_eng

        # -----------------------------------------------------------------------------------------------
        # for debugging purposes only: if elastic_debug is switched on, linear elastic material is used
        # -----------------------------------------------------------------------------------------------
        if self.elastic_debug:
            # NOTE: This must be copied otherwise self.D2_e gets modified when
            # essential boundary conditions are inserted
            D2_e = copy(self.D2_e)
            sig_eng = tensordot(D2_e, eps_app_eng, [[1], [0]])
            return sig_eng, D2_e

        # print 'sctx_correc', sctx.shape
        # -----------------------------------------------------------------------------------------------
        # update state variables
        # -----------------------------------------------------------------------------------------------
        # if sctx.update_state_on:
        #    eps_n = eps_avg - d_eps
        #    e_max = self._get_state_variables(sctx, eps_n)
        #    sctx.mats_state_array[:] = e_max

        #----------------------------------------------------------------------
        # if the regularization using the crack-band concept is on calculate the
        # effective element length in the direction of principle strains
        #----------------------------------------------------------------------
        # if self.regularization:
        #    h = self.get_regularizing_length(sctx, eps_app_eng)
        #    self.phi_fn.h = h

        #------------------------------------------------------------------
        # Damage tensor (4th order) using product- or sum-type symmetrization:
        #------------------------------------------------------------------
        beta_ijkl = self._get_beta_tns(sctx, eps_app_eng, sigma_kk)

        #------------------------------------------------------------------
        # Damaged stiffness tensor calculated based on the damage tensor beta4:
        #------------------------------------------------------------------
        # (cf. [Jir99] Eq.(7): C = beta * D_e * beta^T),
        # minor symmetry is tacitly assumed ! (i.e. beta_ijkl = beta_jilk)
        # D4_mdm = tensordot(
        # beta_ijkl, tensordot(self.D4_e, beta_ijkl, [[2, 3], [2, 3]]),
        # [[2, 3], [0, 1]])

        D4_mdm_ijmn = einsum(
            'ijkl,klsr,mnsr->ijmn', beta_ijkl, self.D4_e, beta_ijkl)

        # print 'Elastic_stiffness', self.D4_e
        # print 'secant_stiffness', D4_mdm_ijmn
        #------------------------------------------------------------------
        # Reduction of the fourth order tensor to a matrix assuming minor and major symmetry:
        #------------------------------------------------------------------
        #D2_mdm = self.map_tns4_to_tns2(D4_mdm_ijmn)

        # print'D2_mdm', D2_mdm.shape

        #----------------------------------------------------------------------
        # Return stresses (corrector) and damaged secant stiffness matrix (predictor)
        #----------------------------------------------------------------------
        eps_pi_ij = self._get_eps_pi_mtx(sctx, eps_app_eng, sigma_kk)
        eps_e_mtx = eps_app_eng - eps_pi_ij

        # print'eps_e_mtx', eps_e_mtx.shape

        #sig_eng = tensordot(D2_mdm, eps_e_mtx, [[1], [0]])
        sig_eng = einsum('ijmn,mn -> ij', D4_mdm_ijmn, eps_e_mtx)

        # print 'sigma ', sig_eng
        return sig_eng, D4_mdm_ijmn


class MATS3DMicroplaneDamage(MATSXDMicroplaneDamageFatigue, MATS3DEval):

    implements(IMATSEval)
    # number of spatial dimensions
    #
    n_dim = Constant(3)

    # number of components of engineering tensor representation
    #
    n_eng = Constant(6)

    #-------------------------------------------------------------------------
    # PolarDiscr related data
    #-------------------------------------------------------------------------
    #
    # number of microplanes - currently fixed for 3D
    #
    n_mp = Constant(28)

    # get the normal vectors of the microplanes
    #
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

    # get the weights of the microplanes
    #
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

    elasticity_tensors = Property(
        depends_on='E, nu, dimensionality, stress_state')

    @cached_property
    def _get_elasticity_tensors(self):
        '''
        Intialize the fourth order elasticity tensor for 3D or 2D plane strain or 2D plane stress
        '''
        # ----------------------------------------------------------------------------
        # Lame constants calculated from E and nu
        # ----------------------------------------------------------------------------
        E = self.E
        nu = self.nu

        # first Lame paramter
        la = E * nu / ((1 + nu) * (1 - 2 * nu))
        # second Lame parameter (shear modulus)
        mu = E / (2 + 2 * nu)

        # -----------------------------------------------------------------------------------------------------
        # Get the fourth order elasticity and compliance tensors for the 3D-case
        # -----------------------------------------------------------------------------------------------------

        # The following line correspond to the tensorial expression:
        # (using numpy functionality in order to avoid the loop):
        #
        # D4_e_3D = zeros((3,3,3,3),dtype=float)
        # C4_e_3D = zeros((3,3,3,3),dtype=float)
        # delta = identity(3)
        # for i in range(0,3):
        #     for j in range(0,3):
        #         for k in range(0,3):
        #             for l in range(0,3):
        #                 # elasticity tensor (cf. Jir/Baz Inelastic analysis of structures Eq.D25):
        #                 D4_e_3D[i,j,k,l] = la * delta[i,j] * delta[k,l] + \
        #                                    mu * ( delta[i,k] * delta[j,l] + delta[i,l] * delta[j,k] )
        #                 # elastic compliance tensor (cf. Simo, Computational Inelasticity, Eq.(2.7.16) AND (2.1.16)):
        #                 C4_e_3D[i,j,k,l] = (1+nu)/(E) * \
        #                                    ( delta[i,k] * delta[j,l] + delta[i,l]* delta[j,k] ) - \
        #                                    nu / E * delta[i,j] * delta[k,l]
        # NOTE: swapaxes returns a reference not a copy!
        # (the index notation always refers to the initial indexing (i=0,j=1,k=2,l=3))
        #delta = identity(3)
        #delta_ijkl = outer(delta, delta).reshape(3, 3, 3, 3)
        #delta_ikjl = delta_ijkl.swapaxes(1, 2)
        #delta_iljk = delta_ikjl.swapaxes(2, 3)
        #D4_e_3D = la * delta_ijkl + mu * (delta_ikjl + delta_iljk)
        # C4_e_3D = -nu / E * delta_ijkl + \
        #    (1 + nu) / (2 * E) * (delta_ikjl + delta_iljk)

        # construct the elasticity tensor (using Numpy - einsum function)
        delta = identity(3)
        D_ijkl = (einsum(',ij,kl->ijkl', la, delta, delta) +
                  einsum(',ik,jl->ijkl', mu, delta, delta) +
                  einsum(',il,jk->ijkl', mu, delta, delta))
        # -----------------------------------------------------------------------------------------------------
        # Get the fourth order elasticity and compliance tensors for the 3D-case
        # -----------------------------------------------------------------------------------------------------
        D2_e_3D = self.map_tns4_to_tns2(D_ijkl)

        return D_ijkl, D2_e_3D

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

    n = 100
    s_levels = linspace(0, -0.005, 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= -1
    #s_levels[0] = 0
    #s_levels.reshape(-1, 2)[:, 1] = -0.006
    s_history = s_levels.flatten()
    s_arr = hstack([linspace(s_history[i], s_history[i + 1], n)
                    for i in range(len(s_levels) - 1)])

    # print s_arr
    m = 1.0
    n = 0.0
    eps = array([array([[m * s_arr[i], n * s_arr[i], 0],
                        [0, -0.0 * s_arr[i], 0],
                        [0, 0, -0.0 * s_arr[i]]]) for i in range(0, len(s_arr))])

    m1 = MATSEvalMicroplaneFatigue()
    m2 = MATS3DMicroplaneDamage()
    sigma = zeros_like(eps)
    sigma_2 = zeros_like(eps)
    sigma_kk = zeros(len(s_arr) + 1)
    w = zeros((len(eps[:, 0, 0]), 28))
    eps_P_N = zeros((len(eps[:, 0, 0]), 28))
    eps_Pi_T = zeros((len(eps[:, 0, 0]), 28, 3))
    e = zeros((len(eps[:, 0, 0]), 28, 3))
    e_T = zeros((len(eps[:, 0, 0]), 28, 3))
    e_N = zeros((len(eps[:, 0, 0]), 28))
    sigma_N = zeros((len(eps[:, 0, 0]), 28))
    sigma_T = zeros((len(eps[:, 0, 0]), 28, 3))
    sctx = zeros((len(eps[:, 0, 0]) + 1, 28, 18))

    for i in range(0, len(eps[:, 0, 0])):

        sigma[i, :] = m2.get_corr_pred(sctx[i, :], eps[i, :], sigma_kk[i])[0]
        sigma_2[i, :] = m2._get_stress_tns(sctx[i, :], eps[i, :], sigma_kk[i])
        sigma_kk[i + 1] = trace(sigma[i, :])

        sctx[
            i + 1] = m2._get_state_variables(sctx[i, :], eps[i, :], sigma_kk[i])
        w[i, :] = sctx[i, :, 6]
        eps_P_N[i, :] = sctx[i, :, 4]
        eps_Pi_T[i, :, :] = sctx[i, :, 10:13]
        sigma_N[i, :] = sctx[i, :, 5]
        sigma_T[i, :, :] = sctx[i, :, 14:17]
        # print eps_P_N

        e[i, :] = m2._get_e_vct_arr(eps[i, :])
        e_T[i, :] = m2._get_e_T_vct_arr_2(eps[i, :])
        e_N[i, :] = m2._get_e_N_arr(e[i, :])

    plt.subplot(221)
    plt.plot(eps[:, 0, 0], sigma[:, 0, 0],
             linewidth=1, label='sigma_11_jirasek')
    plt.plot(eps[:, 0, 0], sigma_2[:, 0, 0],
             linewidth=1, label='sigma_11_carol')
    #plt.plot(eps[:, 0, 0], sigma[:, 1, 1], linewidth=1, label='sigma_22')
    #plt.plot(eps[:, 0, 0], sigma_2[:, 0, 1], linewidth=1, label='sigma_12')
    #plt.plot(eps[:, 0, 0], sigma[:, 2, 2], linewidth=1, label='sigma_33')
    plt.xlabel('Strain')
    plt.ylabel('Stress(MPa)')
    plt.legend()

    plt.subplot(222)
    for i in range(0, 28):
        plt.plot(
            eps[:, 0, 0], w[:, i], linewidth=0.5, label='Damage', alpha=0.5)

        plt.xlabel('Strain')
        plt.ylabel('Damage')
        plt.legend()

    plt.subplot(223)
    for i in range(0, 28):
        #plt.plot(eps[:, 0, 0], eps_P_N[:, i], linewidth=1, label='plastic_strain')\
        plt.plot(
            eps[:, 0, 0], sigma_N[:, i], linewidth=1, label='plastic_strain')

        plt.xlabel('Strain')
        plt.ylabel('normal stress')
        # plt.legend()

    plt.subplot(224)
    for i in range(0, 28):
        #plt.plot(eps[:, 0, 0], eps_Pi_T[:, i, 1], linewidth=1, label='Normal_strain')\
        plt.plot(eps[:, 0, 0], sigma_T[:, i, 1],
                 linewidth=1, label='Normal_strain')
        plt.xlabel('Strain')
        plt.ylabel('tangential stress')

    plt.show()
