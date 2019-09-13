'''
Created on 29.03.2017

@author: abaktheer

Microplane Fatigue model

(compression plasticiy (CP) + Tensile Damage (TD) + Cumulative Damage Sliding (CSD))

Using Jirasek homogenization approach [1999]
'''

from numpy import \
    array, zeros, trace, \
    einsum, zeros_like,\
    identity, sign, linspace, hstack, maximum,\
    sqrt, linalg
from scipy.interpolate import interp1d
from scipy.linalg import \
    eigh
from traits.api import \
    Constant, implements,\
    Float, HasTraits, \
    Property, cached_property
from traitsui.api import \
    View,  Include

#from afxres import AFX_ID_PREVIEW_NEXT
#from ibvpy.mats.mats3D.mats3D_eval import MATS3DEval
# from ibvpy.mats.mats_eval import \
#     IMATSEval
from mathkit.mfn import MFnLineArray
import matplotlib.pyplot as plt
import numpy as np


class MATSEvalMicroplaneFatigue(HasTraits):
    #--------------------------
    # material model parameters
    #--------------------------

    E = Float(35000.,
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
    gamma_T = Float(80000.,
                    label="Gamma",
                    desc=" Tangential Kinematic hardening modulus",
                    enter_set=True,
                    auto_set=False)

    K_T = Float(10000.0,
                label="K",
                desc="Tangential Isotropic harening",
                enter_set=True,
                auto_set=False)

    S = Float(0.000010,
              label="S",
              desc="Damage strength",
              enter_set=True,
              auto_set=False)

    r = Float(1.21,
              label="r",
              desc="Damage cumulation parameter",
              enter_set=True,
              auto_set=False)

    c = Float(1.85,
              label="c",
              desc="Damage cumulation parameter",
              enter_set=True,
              auto_set=False)

    tau_pi_bar = Float(5.0,
                       label="Tau_bar",
                       desc="Reversibility limit",
                       enter_set=True,
                       auto_set=False)

    a = Float(0.001,
              label="a",
              desc="Lateral pressure coefficient",
              enter_set=True,
              auto_set=False)

    #-------------------------------------------
    # Normal_Tension constitutive law parameters
    #-------------------------------------------
    Ad = Float(6000.0,
               label="a",
               desc="brittleness coefficient",
               enter_set=True,
               auto_set=False)

#     eps_f = Float(0.0001,
#                   label="a",
#                   desc="brittleness coefficient",
#                   enter_set=True,
#                   auto_set=False)

    eps_0 = Float(8.0e-5,
                  label="a",
                  desc="threshold strain",
                  enter_set=True,
                  auto_set=False)

    #-----------------------------------------------
    # Normal_Compression constitutive law parameters
    #-----------------------------------------------
    K_N = Float(5000.,
                label="K_N",
                desc=" Normal isotropic harening",
                enter_set=True,
                auto_set=False)

    gamma_N = Float(20000.,
                    label="gamma_N",
                    desc="Normal kinematic hardening",
                    enter_set=True,
                    auto_set=False)

    sigma_0 = Float(18.,
                    label="sigma_0",
                    desc="Yielding stress",
                    enter_set=True,
                    auto_set=False)

    #--------------------------------------------------------------
    # microplane constitutive law (normal behavior CP + TD) (full thermodynamic)
    #--------------------------------------------------------------
    def get_normal_Law(self, eps, sctx):

        E_N = self.E / (1.0 - 2.0 * self.nu)

        w_N = sctx[:, 0]
        z_N = sctx[:, 1]
        alpha_N = sctx[:, 2]
        r_N = sctx[:, 3]
        eps_N_p = sctx[:, 4]

        pos = eps > 1e-6
        H = 1.0 * pos

        sigma_n_trial = (1. - H * w_N) * E_N * (eps - eps_N_p)
        Z = self.K_N * r_N
        X = self.gamma_N * alpha_N

        h = self.sigma_0 + Z
        pos_iso = h > 1e-6
        f_trial = abs(sigma_n_trial - X) - h * pos_iso

        thres_1 = f_trial > 1e-6

        delta_lamda = f_trial / \
            (E_N + abs(self.K_N) + self.gamma_N) * thres_1
        eps_N_p = eps_N_p + delta_lamda * sign(sigma_n_trial - X)
        r_N = r_N + delta_lamda
        alpha_N = alpha_N + delta_lamda * sign(sigma_n_trial - X)

        def Z_N(z_N): return 1.0 / self.Ad * (-z_N) / (1.0 + z_N)
        Y_N = 0.5 * H * E_N * eps ** 2.0
        Y_0 = 0.5 * E_N * self.eps_0 ** 2.0
        f = Y_N - (Y_0 + Z_N(z_N))

        thres_2 = f > 1e-6

        def f_w(Y): return 1.0 - 1.0 / (1.0 + self.Ad * (Y - Y_0))
        w_N = f_w(Y_N) * thres_2
        z_N = - w_N * thres_2

        new_sctx = zeros((28, 5))

        new_sctx[:, 0] = w_N
        new_sctx[:, 1] = z_N
        new_sctx[:, 2] = alpha_N
        new_sctx[:, 3] = r_N
        new_sctx[:, 4] = eps_N_p
        return new_sctx

    #--------------------------------------------------------------
    # microplane constitutive law (normal behavior CP + TD)
    #--------------------------------------------------------------
    def get_normal_Law_2(self, eps, sctx):

        E_N = self.E / (1.0 - 2.0 * self.nu)

        w_N = sctx[:, 0]
        eps_max = sctx[:, 1]
        alpha_N = sctx[:, 2]
        r_N = sctx[:, 3]
        eps_N_p = sctx[:, 4]
        #eps_N_p_cum = sctx[:, 5]

        pos = eps > 1e-6
        H = 1.0 * pos

        sigma_n_trial = (1. - H * w_N) * E_N * (eps - eps_N_p)
        Z = self.K_N * r_N
        X = self.gamma_N * alpha_N

        h = self.sigma_0 + Z
        pos_iso = h > 1e-6
        f_trial = abs(sigma_n_trial - X) - h * pos_iso

        thres_1 = f_trial > 1e-6

        delta_lamda = f_trial / \
            (E_N + abs(self.K_N) + self.gamma_N) * thres_1
        eps_N_p = eps_N_p + delta_lamda * sign(sigma_n_trial - X)
        r_N = r_N + delta_lamda
        alpha_N = alpha_N + delta_lamda * sign(sigma_n_trial - X)

        idx = np.where(eps_max >= self.eps_0)
        #eps_N_p_cum += eps_N_p

        w_N[idx] = (1. - sqrt(
            (self.eps_0 / eps_max[idx]) * np.exp(- (eps_max[idx] - self.eps_0) / (self.eps_f - self.eps_0))))

        eps_max = np.maximum(eps, eps_max)

        # print eps_max

#         def Z_N(z_N): return 1. / self.Ad * (-z_N) / (1. + z_N)
#         Y_N = 0.5 * H * E_N * eps ** 2.
#         Y_0 = 0.5 * E_N * self.eps_0 ** 2.
#         f = Y_N - (Y_0 + Z_N(z_N))
#
#         thres_2 = f > 1e-6
#
#         def f_w(Y): return 1. - 1. / (1. + self.Ad * (Y - Y_0))
#         w_N = f_w(Y_N) * thres_2
#         z_N = - w_N * thres_2

        new_sctx = zeros((28, 5))

        new_sctx[:, 0] = w_N
        new_sctx[:, 1] = eps_max
        new_sctx[:, 2] = alpha_N
        new_sctx[:, 3] = r_N
        new_sctx[:, 4] = eps_N_p
        #new_sctx[:, 5] = eps_N_p_cum
        return new_sctx

    #-------------------------------------------------------------------------
    # microplane constitutive law (Tangential CSD)-(Pressure sensitive cumulative damage)
    #-------------------------------------------------------------------------
    def get_tangential_Law(self, e_T, sctx, sigma_kk):

        E_T = self.E / (1. + self.nu)

        w_T = sctx[:, 5]
        z_T = sctx[:, 6]
        alpha_T = sctx[:, 7:10]
        eps_T_pi = sctx[:, 10:13]
        eps_T_pi_cum = sctx[:, 13]

        sig_pi_trial = E_T * (e_T - eps_T_pi)
        Z = self.K_T * z_T
        X = self.gamma_T * alpha_T
        norm_1 = sqrt(
            einsum('nj,nj -> n', (sig_pi_trial - X), (sig_pi_trial - X)))

        f = norm_1 - self.tau_pi_bar - \
            Z + self.a * sigma_kk / 3.0

        plas_1 = f > 1e-6
        elas_1 = f < 1e-6

        delta_lamda = f / \
            (E_T / (1.0 - w_T) + self.gamma_T + self.K_T) * plas_1

        norm_2 = 1.0 * elas_1 + sqrt(
            einsum('nj,nj -> n', (sig_pi_trial - X), (sig_pi_trial - X))) * plas_1

        eps_T_pi[:, 0] = eps_T_pi[:, 0] + plas_1 * delta_lamda * \
            ((sig_pi_trial[:, 0] - X[:, 0]) / (1.0 - w_T)) / norm_2
        eps_T_pi[:, 1] = eps_T_pi[:, 1] + plas_1 * delta_lamda * \
            ((sig_pi_trial[:, 1] - X[:, 1]) / (1.0 - w_T)) / norm_2
        eps_T_pi[:, 2] = eps_T_pi[:, 2] + plas_1 * delta_lamda * \
            ((sig_pi_trial[:, 2] - X[:, 2]) / (1.0 - w_T)) / norm_2

        Y = 0.5 * E_T * \
            einsum('nj,nj -> n', (e_T - eps_T_pi), (e_T - eps_T_pi))

        w_T += ((1 - w_T) ** self.c) * \
            (delta_lamda * (Y / self.S) ** self.r) * \
            (self.tau_pi_bar / (self.tau_pi_bar - self.a * sigma_kk / 3.0))

        alpha_T[:, 0] = alpha_T[:, 0] + plas_1 * delta_lamda *\
            (sig_pi_trial[:, 0] - X[:, 0]) / norm_2
        alpha_T[:, 1] = alpha_T[:, 1] + plas_1 * delta_lamda *\
            (sig_pi_trial[:, 1] - X[:, 1]) / norm_2
        alpha_T[:, 2] = alpha_T[:, 2] + plas_1 * delta_lamda *\
            (sig_pi_trial[:, 2] - X[:, 2]) / norm_2

        z_T = z_T + delta_lamda
        eps_T_pi_cum = eps_T_pi_cum + np.linalg.norm(eps_T_pi)

        new_sctx = zeros((28, 9))
        new_sctx[:, 0] = w_T
        new_sctx[:, 1] = z_T
        new_sctx[:, 2:5] = alpha_T
        new_sctx[:, 5:8] = eps_T_pi
        new_sctx[:, 8] = np.linalg.norm(eps_T_pi)
        return new_sctx


class MATSXDMicroplaneDamageFatigueJir(MATSEvalMicroplaneFatigue):

    '''
    Microplane Damage Fatigue Model.
    '''

    #-------------------------------------------------------------------------
    # Setup for computation within a supplied spatial context
    #-------------------------------------------------------------------------
    D4_e = Property

    def _get_D4_e(self):
        # Return the elasticity tensor
        return self.elasticity_tensors

    #-------------------------------------------------------------------------
    # MICROPLANE-Kinematic constraints
    #-------------------------------------------------------------------------

    # get the dyadic product of the microplane normals
    _MPNN = Property(depends_on='n_mp')

    @cached_property
    def _get__MPNN(self):
        # dyadic product of the microplane normals

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
        return einsum('nij,ij->n', self._MPNN, eps_eng)

    def _get_e_T_vct_arr_2(self, eps_eng):
        MPTT_ijr = self._get__MPTT()
        return einsum('nijr,ij->nr', MPTT_ijr, eps_eng)

    def _get_e_vct_arr_2(self, eps_eng):
        return self._e_N_arr_2 * self._MPN + self._e_t_vct_arr_2

    #--------------------------------------------------------
    # return the state variables (Damage , inelastic strains)
    #--------------------------------------------------------
    def _get_state_variables(self, sctx, eps_app_eng, sigma_kk):

        e_N_arr = self._get_e_N_arr_2(eps_app_eng)
        e_T_vct_arr = self._get_e_T_vct_arr_2(eps_app_eng)

        sctx_arr = zeros_like(sctx)

        sctx_N = self.get_normal_Law(e_N_arr, sctx)
        sctx_arr[:, 0:5] = sctx_N

        sctx_tangential = self.get_tangential_Law(e_T_vct_arr, sctx, sigma_kk)
        sctx_arr[:, 5:14] = sctx_tangential

        return sctx_arr

    #-----------------------------------------------------------------
    # Returns a list of the plastic normal strain  for all microplanes.
    #-----------------------------------------------------------------
    def _get_eps_N_p_arr(self, sctx, eps_app_eng, sigma_kk):

        eps_N_p = self._get_state_variables(sctx, eps_app_eng, sigma_kk)[:, 4]
        return eps_N_p

    #----------------------------------------------------------------
    # Returns a list of the sliding strain vector for all microplanes.
    #----------------------------------------------------------------
    def _get_eps_T_pi_arr(self, sctx, eps_app_eng, sigma_kk):

        eps_T_pi_vct_arr = self._get_state_variables(
            sctx, eps_app_eng, sigma_kk)[:, 10:13]

        return eps_T_pi_vct_arr

    #---------------------------------------------------------------------
    # Extra homogenization of damge tensor in case of two damage parameters
    #---------------------------------------------------------------------

    def _get_beta_N_arr(self, sctx, eps_app_eng, sigma_kk):

        # Returns a list of the normal integrity factors for all microplanes.

        beta_N_arr = sqrt(1 -
                          self._get_state_variables(sctx, eps_app_eng, sigma_kk)[:, 0])

        return beta_N_arr

    def _get_beta_T_arr(self, sctx, eps_app_eng, sigma_kk):

        # Returns a list of the tangential integrity factors for all
        # microplanes.

        beta_T_arr = sqrt(1 -
                          self._get_state_variables(sctx, eps_app_eng, sigma_kk)[:, 5])

        return beta_T_arr

    def _get_beta_tns(self, sctx, eps_app_eng, sigma_kk):

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

    #-------------------------------------------------------------
    # Returns a list of the integrity factors for all microplanes.
    #-------------------------------------------------------------

    def _get_phi_arr(self, sctx, eps_app_eng, sigma_kk):

        w_n = self._get_state_variables(sctx, eps_app_eng, sigma_kk)[:, 0]
        w_T = self._get_state_variables(sctx, eps_app_eng, sigma_kk)[:, 5]

        w = zeros(self.n_mp)

#         w = maximum(w_n, w_T)

        eig = np.linalg.eig(eps_app_eng)[0]

        ter_1 = np.sum(eig)

        if ter_1 > 0.0:
            w = w_n
        else:
            w = w_T

        phi_arr = sqrt(1.0 - w)

        return phi_arr

    #----------------------------------------------
    # Returns the 2nd order damage tensor 'phi_mtx'
    #----------------------------------------------
    def _get_phi_mtx(self, sctx, eps_app_eng, sigma_kk):

        # scalar integrity factor for each microplane
        phi_arr = self._get_phi_arr(sctx, eps_app_eng, sigma_kk)

        # integration terms for each microplanes
        phi_ij = einsum('n,n,nij->ij', phi_arr, self._MPW, self._MPNN)

        return phi_ij

    #----------------------------------------------------------------------
    # Returns the 4th order damage tensor 'beta4' using sum-type symmetrization
    #(cf. [Jir99], Eq.(21))
    #----------------------------------------------------------------------
    def _get_beta_tns_sum_type(self, sctx, eps_app_eng, sigma_kk):

        delta = identity(3)

        phi_mtx = self._get_phi_mtx(sctx, eps_app_eng, sigma_kk)

        # use numpy functionality (einsum) to evaluate [Jir99], Eq.(21)
        beta_ijkl = 0.25 * (einsum('ik,jl->ijkl', phi_mtx, delta) +
                            einsum('il,jk->ijkl', phi_mtx, delta) +
                            einsum('jk,il->ijkl', phi_mtx, delta) +
                            einsum('jl,ik->ijkl', phi_mtx, delta))

        return beta_ijkl

    #----------------------------------------------------------------------
    # Returns the 4th order damage tensor 'beta4' using product-type symmetrization
    #(cf. [Baz97], Eq.(87))
    #----------------------------------------------------------------------
    def _get_beta_tns_product_type(self, sctx, eps_app_eng, sigma_kk):

        delta = identity(3)

        phi_mtx = self._get_phi_mtx(sctx, eps_app_eng, sigma_kk)

        n_dim = 3
        phi_eig_value, phi_eig_mtx = eigh(phi_mtx)
        phi_eig_value_real = array([pe.real for pe in phi_eig_value])
        phi_pdc_mtx = zeros((n_dim, n_dim), dtype=float)
        for i in range(n_dim):
            phi_pdc_mtx[i, i] = phi_eig_value_real[i]
        # w_mtx = tensorial square root of the second order damage tensor:
        w_pdc_mtx = sqrt(phi_pdc_mtx)

        # transform the matrix w back to x-y-coordinates:
        w_mtx = einsum('ik,kl,lj -> ij', phi_eig_mtx, w_pdc_mtx, phi_eig_mtx)
        #w_mtx = dot(dot(phi_eig_mtx, w_pdc_mtx), transpose(phi_eig_mtx))

        beta_ijkl = 0.5 * \
            (einsum('ik,jl -> ijkl', w_mtx, w_mtx) +
             einsum('il,jk -> ijkl', w_mtx, w_mtx))

        return beta_ijkl

    #-----------------------------------------------------------
    # Integration of the (inelastic) strains for each microplane
    #-----------------------------------------------------------
    def _get_eps_p_mtx(self, sctx, eps_app_eng, sigma_kk):

        # plastic normal strains
        eps_N_P_n = self._get_eps_N_p_arr(sctx, eps_app_eng, sigma_kk)

        # sliding tangential strains
        eps_T_pi_ni = self._get_eps_T_pi_arr(sctx, eps_app_eng, sigma_kk)
        delta = identity(3)

        #eps_T_pi_ni = np.zeros_like(eps_T_pi_ni)

        # 2-nd order plastic (inelastic) tensor
        eps_p_ij = einsum('n,n,ni,nj -> ij', self._MPW, eps_N_P_n, self._MPN, self._MPN) + \
            0.5 * (einsum('n,nr,ni,rj->ij', self._MPW, eps_T_pi_ni, self._MPN, delta) +
                   einsum('n,nr,nj,ri->ij', self._MPW, eps_T_pi_ni, self._MPN, delta))

        # print 'eps_p', eps_p_ij

        return eps_p_ij

    #-------------------------------------------------------------------------
    # Evaluation - get the corrector and predictor
    #-------------------------------------------------------------------------

    def get_corr_pred(self, sctx, eps_app_eng, sigma_kk):

        # Corrector predictor computation.

        #------------------------------------------------------------------
        # Damage tensor (4th order) using product- or sum-type symmetrization:
        #------------------------------------------------------------------
        beta_ijkl = self._get_beta_tns_sum_type(
            sctx, eps_app_eng, sigma_kk)

#         beta_ijkl = self._get_beta_tns(
#             sctx, eps_app_eng, sigma_kk)

        #------------------------------------------------------------------
        # Damaged stiffness tensor calculated based on the damage tensor beta4:
        #------------------------------------------------------------------
        D4_mdm_ijmn = einsum(
            'ijkl,klsr,mnsr->ijmn', beta_ijkl, self.D4_e, beta_ijkl)

        #----------------------------------------------------------------------
        # Return stresses (corrector) and damaged secant stiffness matrix (predictor)
        #----------------------------------------------------------------------
        # plastic strain tensor
        eps_p_ij = self._get_eps_p_mtx(sctx, eps_app_eng, sigma_kk)

        # elastic strain tensor
        eps_e_mtx = eps_app_eng - eps_p_ij
        # print 'eps_e', eps_e_mtx

        # calculation of the stress tensor
        sig_eng = einsum('ijmn,mn -> ij', D4_mdm_ijmn, eps_e_mtx)

        return sig_eng, D4_mdm_ijmn


class MATS3DMicroplaneDamageJir(MATSXDMicroplaneDamageFatigueJir):

    # implements(IMATSEval)

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
    #==========================================================================
    # Check the model behavior at the single material point
    #==========================================================================

    model = MATS3DMicroplaneDamageJir()

    p = 1.0  # ratio of strain eps_11 (for bi-axial loading)
    m = 0.0  # ratio of strain eps_22 (for bi-axial loading)

    #------------------------------------
    # monotonic loading - comp
    #------------------------------------
    n = 500  # number of increments
    s_levels = linspace(0, -0.002, 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= 0
    s_history_1 = s_levels.flatten()
    s_arr_1 = hstack([linspace(s_history_1[i], s_history_1[i + 1], n)
                      for i in range(len(s_levels) - 1)])

    eps_1 = array([array([[p * s_arr_1[i], 0, 0],
                          [0, m * s_arr_1[i], 0],
                          [0, 0, m * s_arr_1[i]]]) for i in range(0, len(s_arr_1))])

    #--------------------------------------
    # construct the arrays
    #--------------------------------------
    sigma_1 = zeros_like(eps_1)
    sigma_kk_1 = zeros(len(s_arr_1) + 1)
    w_1_N = zeros((len(eps_1[:, 0, 0]), 28))
    w_1_T = zeros((len(eps_1[:, 0, 0]), 28))
    eps_P_N_1 = zeros((len(eps_1[:, 0, 0]), 28))
    eps_Pi_T_1 = zeros((len(eps_1[:, 0, 0]), 28, 3))
    eps_Pi_T_norm_1 = zeros((len(eps_1[:, 0, 0]), 28))
    sctx_1 = zeros((len(eps_1[:, 0, 0]) + 1, 28, 14))

    for i in range(0, len(eps_1[:, 0, 0])):
        sigma_1[i, :] = model.get_corr_pred(
            sctx_1[i, :], eps_1[i, :], sigma_kk_1[i])[0]
        sigma_kk_1[i + 1] = trace(sigma_1[i, :])
        sctx_1[
            i + 1] = model._get_state_variables(sctx_1[i, :], eps_1[i, :], sigma_kk_1[i])

        w_1_N[i, :] = sctx_1[i, :, 0]
        w_1_T[i, :] = sctx_1[i, :, 5]
        eps_P_N_1[i, :] = sctx_1[i, :, 4]
        eps_Pi_T_1[i, :, :] = sctx_1[i, :, 10:13]
        eps_Pi_T_norm_1[i, :] = sctx_1[i, :, 13]

    #-------------------------------------
    # cyclic loading - comp
    #-------------------------------------
    n = 100
    s_history_2 = [-0, -0.001, -0.00030, -
                   0.0015, -0.00060, -0.0020, -0.00090,
                   -0.0027, -0.0012, -0.004, -0.0018, -0.0055,
                   -0.0025, -0.007, -0.004, -0.008, -0.0046, -.009,
                   -0.0052, -0.01, -0.0058, -0.012, -0.007, -0.015]

    # manualy-2 (S=0.9)
    s_history_2 = [-0, -0.00169, -0.00065, -
                   0.00182, -0.000712, -0.00205, -0.0008, -0.0024, -0.00097, -0.0031, ]

    # manualy-2 (S=0.85)
    s_history_2 = [-0, -0.0015, -0.00065,
                   -0.001546, -0.000671,
                   -0.001599, -0.0006906,
                   -0.001667, -0.000726,
                   -0.001749, -0.0007676,
                   -0.00185, -0.000818,
                   -0.00198, -0.0008825,
                   -0.0021835, -0.0009826,
                   -0.00262, -0.0011942,
                   -0.0033, ]

    #s_history_2 = [0, -0.01]

    n = 200  # number of increments
    s_levels = linspace(0, -0.002, 200)
    s_levels.reshape(-1, 2)[:, 1] = -0.002
    s_levels.reshape(-1, 2)[:, 0] = -0.001
    s_levels[0] = 0
    s_history_2 = s_levels.flatten()

    s_arr_2 = hstack([linspace(s_history_2[i], s_history_2[i + 1], n)
                      for i in range(len(s_history_2) - 1)])

    idx_1 = np.where(s_arr_2 == -0.002)
    idx_2 = np.where(s_arr_2 == -0.001)

    eps_2 = array([array([[p * s_arr_2[i], 0, 0],
                          [0,  -0 * s_arr_2[i], 0],
                          [0, 0, -0 * s_arr_2[i]]]) for i in range(0, len(s_arr_2))])

    #--------------------------------------
    # construct the arrays
    #--------------------------------------
    sigma_2 = zeros_like(eps_2)
    sigma_kk_2 = zeros(len(s_arr_2) + 1)
    w_2_N = zeros((len(eps_2[:, 0, 0]), 28))
    w_2_T = zeros((len(eps_2[:, 0, 0]), 28))
    eps_P_N_2 = zeros((len(eps_2[:, 0, 0]), 28))
    eps_Pi_T_2 = zeros((len(eps_2[:, 0, 0]), 28, 3))
    eps_Pi_T_norm_2 = zeros((len(eps_2[:, 0, 0]), 28))
    sctx_2 = zeros((len(eps_2[:, 0, 0]) + 1, 28, 14))

    for i in range(0, len(eps_2[:, 0, 0])):

        sigma_2[i, :] = model.get_corr_pred(
            sctx_2[i, :], eps_2[i, :], sigma_kk_2[i])[0]
        sigma_kk_2[i + 1] = trace(sigma_2[i, :])
        sctx_2[
            i + 1] = model._get_state_variables(sctx_2[i, :], eps_2[i, :], sigma_kk_2[i])

        w_2_N[i, :] = sctx_2[i, :, 0]
        w_2_T[i, :] = sctx_2[i, :, 5]
        eps_P_N_2[i, :] = sctx_2[i, :, 4]
        eps_Pi_T_2[i, :, :] = sctx_2[i, :, 10:13]
        eps_Pi_T_norm_2[i, :] = sctx_2[i, :, 13]

#     #------------------------------------
#     # monotonic loading - ten
#     #------------------------------------
#     n = 100  # number of increments
#     s_levels = linspace(0, 0.0013, 2)
#     s_levels.reshape(-1, 2)[:, 0] = -0.0005
#     s_levels[0] = 0
#     s_history_3 = s_levels.flatten()
#     s_arr_3 = hstack([linspace(s_history_3[i], s_history_3[i + 1], n)
#                       for i in range(len(s_levels) - 1)])
#
#     eps_3 = array([array([[p * s_arr_3[i], 0, 0],
#                           [0, 0 * s_arr_3[i], 0],
#                           [0, 0, 0]]) for i in range(0, len(s_arr_3))])
#
#     # print s_history_3
#
#     #--------------------------------------
#     # construct the arrays
#     #--------------------------------------
#     sigma_3 = zeros_like(eps_3)
#     sigma_kk_3 = zeros(len(s_arr_3) + 1)
#     w_3_N = zeros((len(eps_3[:, 0, 0]), 28))
#     w_3_T = zeros((len(eps_3[:, 0, 0]), 28))
#     eps_P_N_3 = zeros((len(eps_3[:, 0, 0]), 28))
#     eps_Pi_T_3 = zeros((len(eps_3[:, 0, 0]), 28, 3))
#     eps_Pi_T_norm_3 = zeros((len(eps_3[:, 0, 0]), 28))
#     sctx_3 = zeros((len(eps_3[:, 0, 0]) + 1, 28, 14))
#
#     for i in range(0, len(eps_3[:, 0, 0])):
#         sigma_3[i, :] = model.get_corr_pred(
#             sctx_3[i, :], eps_3[i, :], sigma_kk_3[i])[0]
#         sigma_kk_3[i + 1] = trace(sigma_3[i, :])
#         sctx_3[
#             i + 1] = model._get_state_variables(sctx_3[i, :], eps_3[i, :], sigma_kk_3[i])
#
#         w_3_N[i, :] = sctx_3[i, :, 0]
#         w_3_T[i, :] = sctx_3[i, :, 5]
#         eps_P_N_3[i, :] = sctx_3[i, :, 4]
#         eps_Pi_T_3[i, :, :] = sctx_3[i, :, 10:13]
#         eps_Pi_T_norm_3[i, :] = sctx_3[i, :, 13]

    '''====================================================
    plotting
    ===================================================='''

    eps_c, sig_c = np.loadtxt(
        r'E:\Publishing\ECCM_2018\exp_results\uniaxial_comp.txt')
    eps_t, sig_t = np.loadtxt(
        r'E:\Publishing\ECCM_2018\exp_results\uniaxial_tension.txt')

    eps_cc, sig_cc = np.loadtxt(
        r'E:\Publishing\ECCM_2018\exp_results\uniaxial_comp_cyclic.txt')

    #------------------------------------------------------
    # stress -strain (compression)
    #------------------------------------------------------
    plt.subplot(221)
    plt.plot(-eps_1[:, 0, 0], -sigma_1[:, 0, 0], color='k',
             linewidth=1, label='sigma_11_(monotonic-compression)')
    plt.plot(-eps_2[:, 0, 0], -sigma_2[:, 0, 0], color='k',
             linewidth=1, label='sigma_11_(cyclic-compression)')

    plt.plot(eps_c / 100, sig_c, 'o', markersize=8, label='exp')

    plt.title('$\sigma - \epsilon$')
    plt.xlabel('strain')
    plt.ylabel('stress(MPa)')
    plt.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    plt.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.legend()

    np.savetxt(r'E:\Publishing\ECCM_2018\num_results\save\eps_monotonic.txt',
               np.transpose(-eps_1[:, 0, 0]), delimiter=" ", fmt="%s")
    np.savetxt(r'E:\Publishing\ECCM_2018\num_results\save\eps_cyclic.txt',
               np.transpose(-eps_2[:, 0, 0]), delimiter=" ", fmt="%s")
    np.savetxt(r'E:\Publishing\ECCM_2018\num_results\save\sig_monotonic.txt',
               np.transpose(-sigma_1[:, 0, 0]), delimiter=" ", fmt="%s")
    np.savetxt(r'E:\Publishing\ECCM_2018\num_results\save\sig_cyclic.txt',
               np.transpose(-sigma_2[:, 0, 0]), delimiter=" ", fmt="%s")

#     #------------------------------------------------------
#     # stress -strain (tension)
#     #------------------------------------------------------
#     plt.subplot(222)
#
#     plt.plot(eps_3[:, 0, 0], sigma_3[:, 0, 0], color='k',
#              linewidth=1, label='sigma_11_(monotonic-tension)')
#
#     plt.plot(eps_t / 10000, sig_t, 'o', label='exp')
#
#     plt.title('$\sigma - \epsilon$')
#     plt.xlabel('strain')
#     plt.ylabel('stress(MPa)')
#     plt.axhline(y=0, color='k', linewidth=1, alpha=0.5)
#     plt.axvline(x=0, color='k', linewidth=1, alpha=0.5)
#     plt.legend()

#     #------------------------------------------------------
#     # normal damage at the microplanes (TD)
#     #------------------------------------------------------
#     plt.subplot(222)
#     for i in range(0, 28):
#         plt.plot(
#             eps_1[:, 0, 0], w_1_N[:, i], linewidth=1.0, label='cyclic', alpha=1)
#         plt.plot(
#             eps_2[:, 0, 0], w_2_N[:, i], linewidth=1.0, label='monotonic', alpha=1)
#         plt.plot(
#             eps_3[:, 0, 0], w_3_N[:, i], linewidth=1.0, label='monotonic', alpha=1)
#
#         plt.xlabel('strain')
#         plt.ylabel('damage')
#         plt.title(' normal damage for all microplanes')

    #---------------------------------------------------------
    # tangential damage at the microplanes (CSD)
    #---------------------------------------------------------
    plt.subplot(222)

    d_array_1 = -sigma_2[idx_1, 0, 0][0]
    d_array_2 = -sigma_2[idx_2, 0, 0][0]
    n_max_1 = len(idx_1[0])
    n_max_2 = len(idx_2[0])
    n_arr_1 = np.arange(1, n_max_1, 2)
    n_arr_2 = np.arange(1, n_max_2, 2)
    d_array_1 = d_array_1[n_arr_1]
    d_array_2 = d_array_2[n_arr_2]

    n_arr_1 = np.arange(1, len(n_arr_1) + 1, 1)
    n_arr_2 = np.arange(1, len(n_arr_2) + 1, 1)
    # print n_arr, d_array_1
    time_func_1 = interp1d(n_arr_1, d_array_1)
    time_func_2 = interp1d(n_arr_2, d_array_2)
    # print time_func
    plt.plot(
        n_arr_1, time_func_1(n_arr_1))
    plt.plot(
        n_arr_2, time_func_2(n_arr_2))

    plt.fill_between(n_arr_1,  time_func_1(n_arr_1),
                     time_func_2(n_arr_2), facecolor='gray', alpha=0.2)

    plt.xlabel('strain')
    plt.ylabel('damage')
    plt.title(' tangential damage for all microplanes')

    np.savetxt(r'E:\Publishing\ECCM_2018\num_results\save\n_max.txt',
               np.transpose(n_arr_1), delimiter=" ", fmt="%s")
    np.savetxt(r'E:\Publishing\ECCM_2018\num_results\save\n_min.txt',
               np.transpose(n_arr_2), delimiter=" ", fmt="%s")
    np.savetxt(r'E:\Publishing\ECCM_2018\num_results\save\sig_max.txt',
               np.transpose(time_func_1(n_arr_1)), delimiter=" ", fmt="%s")
    np.savetxt(r'E:\Publishing\ECCM_2018\num_results\save\sig_min.txt',
               np.transpose(time_func_2(n_arr_2)), delimiter=" ", fmt="%s")

    #---------------------------------------------------------
    # tangential damage at the microplanes (CSD)
    #---------------------------------------------------------
    plt.subplot(223)
    for i in range(0, 28):
        plt.plot(
            eps_1[:, 0, 0], w_1_T[:, i], linewidth=1.0, label='cyclic', alpha=1)
        plt.plot(
            eps_2[:, 0, 0], w_2_T[:, i], linewidth=1.0, label='monotonic', alpha=1)

        plt.xlabel('strain')
        plt.ylabel('damage')
        plt.title(' tangential damage for all microplanes')

    #-----------------------------------------------------------
    # damage with sliding strains at the microplanes (CSD)
    #-----------------------------------------------------------
    plt.subplot(224)
    for i in range(0, 28):

        #         plt.plot(eps_Pi_T_1[:, i, 0], w_1_T[:, i])
        #         plt.plot(eps_Pi_T_2[:, i, 0], w_2_T[:, i])
        #         plt.plot(eps_Pi_T_3[:, i, 0], w_3_T[:, i])

        print(eps_Pi_T_1[:, i].shape)

        #plt.plot(eps_Pi_T_norm_1[:, i], w_1_T[:, i])
        plt.plot(eps_Pi_T_norm_2[:, i], w_2_T[:, i])

        plt.xlabel('sliding strain')
        plt.ylabel('damage')
        #plt.xlim(-0.00025, 0.00025)
        plt.ylim(0, 1.0)

    plt.show()
