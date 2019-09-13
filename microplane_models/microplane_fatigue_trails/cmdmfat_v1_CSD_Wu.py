'''
Created on 24.03.2017

@author: abaktheer

New implementation of microplane model - (Wu)
'''

from ibvpy.mats.mats3D.mats3D_eval import MATS3DEval
from ibvpy.mats.mats_eval import \
    IMATSEval
from numpy import \
    array, zeros, dot, trace, \
    tensordot, einsum, zeros_like,\
    identity, linspace, hstack, \
    sqrt as arr_sqrt, copy

from numpy.linalg import norm

from traits.api import \
    Constant,  Property, cached_property, implements,\
    Bool, Callable, Enum, Float, HasTraits, \
    Int, Trait, on_trait_change, \
    Dict, Property, cached_property
from traitsui.api import \
    Item, View, Group, Spring, Include

import matplotlib.pyplot as plt


class MATSEvalMicroplaneFatigue(HasTraits):

    E = Float(34000,
              label="G",
              desc="Elastic modulus",
              enter_set=True,
              auto_set=False)

    poisson = Float(0.2,
                    label="G",
                    desc="poisson's ratio",
                    enter_set=True,
                    auto_set=False)

    gamma = Float(100.0,
                  label="Gamma",
                  desc="Kinematic hardening modulus",
                  enter_set=True,
                  auto_set=False)

    K = Float(0.0,
              label="K",
              desc="Isotropic harening",
              enter_set=True,
              auto_set=False)

    S = Float(0.000001,
              label="S",
              desc="Damage cumulation parameter",
              enter_set=True,
              auto_set=False)

    r = Float(1.0,
              label="r",
              desc="Damage cumulation parameter",
              enter_set=True,
              auto_set=False)

    c = Float(1.0,
              label="c",
              desc="Damage cumulation parameter",
              enter_set=True,
              auto_set=False)

    tau_pi_bar = Float(1,
                       label="Tau_pi_bar",
                       desc="Reversibility limit",
                       enter_set=True,
                       auto_set=False)

    a = Float(0.0,
              label="a",
              desc="Lateral pressure coefficient",
              enter_set=True,
              auto_set=False)

    def get_phi_epspi(self, eps, sctx, sigma_kk):

        G = self.E / (2 * (1 + self.poisson))
        w = sctx[0]
        z = sctx[1]
        alpha = sctx[2:5]
        xs_pi = sctx[5:9]

        sig_pi_trial = G * (eps - xs_pi)

        Z = self.K * z
        X = self.gamma * alpha

        f = norm(sig_pi_trial - X) - self.tau_pi_bar - \
            Z + self.a * sigma_kk / 3.

        if f > 1e-6:

            delta_lamda = f / \
                (G / (1 - w) + self.gamma + self.K)

            xs_pi = xs_pi + delta_lamda * \
                ((sig_pi_trial - X) / (1 - w)) / norm(sig_pi_trial - X)

            Y = 0.5 * G * dot((eps - xs_pi), (eps - xs_pi))

            w += ((1 - w) ** self.c) * \
                (delta_lamda * (Y / self.S) ** self.r)

            # print 'w', w

            X = X + self.gamma * delta_lamda * \
                (sig_pi_trial - X) / norm(sig_pi_trial - X)
            alpha = alpha + delta_lamda * \
                (sig_pi_trial - X) / norm(sig_pi_trial - X)
            z = z + delta_lamda

        new_sctx = zeros(8)
        new_sctx[0:2] = w, z
        new_sctx[2:5] = alpha
        new_sctx[5:9] = xs_pi

        return new_sctx


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

    #-------------------------------------------------------------------------
    # MICROPLANE-DISCRETIZATION RELATED METHOD
    #-------------------------------------------------------------------------
    # get the dyadic product of the microplane normals
    _MPNN = Property(depends_on='n_mp')

    @cached_property
    def _get__MPNN(self):
        # dyadic product of the microplane normals

        MPNN_nij = einsum('ni,nj->nij', self._MPN, self._MPN)
        return MPNN_nij

    # get Third order tangential tensor (operator) for each microplane

    _MPTT = Property(depends_on='n_mp')

    @cached_property
    def _get__MPTT(self):
        # Third order tangential tensor for each microplane
        delta = identity(3)
        MPTT_nijr = 0.5 * (einsum('ni,jr -> nijr', self._MPN, delta) +
                           einsum('nj,ir -> njir', self._MPN, delta) - 2.0 *
                           einsum('ni,nj,nr -> nijr', self._MPN, self._MPN, self._MPN))
        return MPTT_nijr

    def _get_P_vol(self):
        delta = identity(3)
        P_vol_ij = (1 / 3.0) * delta
        return P_vol_ij

    def _get_P_dev(self):
        delta = identity(3)
        P_dev_njkl = 0.5 * einsum('ni,ij,kl -> njkl', self._MPN, delta, delta)
        return P_dev_njkl

    def _get_I_vol_4(self):
        # The fourth order volumetric-identity tensor
        delta = identity(3)
        I_vol_ijkl = (1.0 / 3.0) * einsum('ij,kl -> ijkl', delta, delta)
        return I_vol_ijkl

    def _get_I_dev_4(self):
        # The fourth order deviatoric-identity tensor
        delta = identity(3)
        I_dev_ijkl = 0.5 * (einsum('ik,jl -> ijkl', delta, delta) +
                            einsum('il,jk -> ijkl', delta, delta)) \
            - (1 / 3.0) * einsum('ij,kl -> ijkl', delta, delta)
        return I_dev_ijkl

    def _get_PP_vol_4(self):
        # outer product of P_vol
        delta = identity(3)
        PP_vol_ijkl = (1 / 9.) * einsum('ij,kl -> ijkl', delta, delta)
        return PP_vol_ijkl

    def _get_PP_dev_4(self):
        # inner product of P_dev
        delta = identity(3)
        PP_dev_ijkl = 0.5 * (0.5 * (einsum('ni,nk,jl -> nijkl', self._MPN, self._MPN, delta) +
                                    einsum('ni,nl,jk -> nijkl', self._MPN, self._MPN, delta)) +
                             0.5 * (einsum('ik,nj,nl -> nijkl',  delta, self._MPN, self._MPN) +
                                    einsum('il,nj,nk -> nijkl',  delta, self._MPN, self._MPN))) -\
                            (1 / 3.) * (einsum('ni,nj,kl -> nijkl', self._MPN, self._MPN, delta) +
                                        einsum('ij,nk,nl -> nijkl', delta, self._MPN, self._MPN)) +\
                            (1 / 9.) * einsum('ij,kl -> ijkl', delta, delta)

        return PP_dev_ijkl

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

    def _get_state_variables(self, sctx, eps_app_eng, sigma_kk):

        e_vct_arr = self._get_e_vct_arr(eps_app_eng)
        # print e_vct_arr.shape

        sctx_arr = zeros((28, 8))
        for i in range(0, self.n_mp):
            sctx_i = self.get_phi_epspi(
                e_vct_arr[i, :], sctx[i, :], sigma_kk)

            sctx_arr[i, :] = sctx_i

        return sctx_arr

    def _get_phi_arr(self, sctx, eps_app_eng, sigma_kk):
        '''
        Returns a list of the integrity factors for all microplanes.
        '''
        phi_arr = 1. - \
            self._get_state_variables(sctx, eps_app_eng, sigma_kk)[:, 0]
        return phi_arr

    def _get_phi_mtx(self, sctx, eps_app_eng, sigma_kk):
        '''
        Returns the 2nd order damage tensor 'phi_mtx'
        '''
        # scalar integrity factor for each microplane
        phi_arr = self._get_phi_arr(sctx, eps_app_eng, sigma_kk)
        # integration terms for each microplanes

        phi_ij = einsum(
            'n,n,nij->ij', phi_arr, self._MPW, self._MPNN)

        return phi_ij

    def _get_d_scalar(self, sctx, eps_app_eng, sigma_kk):

        # scalar integrity factor for each microplane
        phi_arr = self._get_phi_arr(sctx, eps_app_eng, sigma_kk)

        d_arr = 1 - phi_arr

        d = (1.0 / 3.0) * einsum('n,n->',  d_arr, self._MPW)

        print d

        return d

    def _get_M_vol_tns(self, sctx, eps_app_eng, sigma_kk):

        d = self._get_d_scalar(sctx, eps_app_eng, sigma_kk)
        delta = identity(3)

        I_4th_ijkl = einsum('ik,jl -> ijkl', delta, delta)

        return (1 - d) * I_4th_ijkl

    def _get_M_dev_tns(self, phi_mtx):
        '''
        Returns the 4th order deviatoric damage tensor
        '''
        n_dim = self.n_dim
        delta = identity(3)

        # use numpy functionality (einsum) to evaluate [Jir99], Eq.(21)
        # M_dev_ijkl = 0.25 * (einsum('ik,jl->ijkl', phi_mtx, delta) +
        #                    einsum('il,jk->ijkl', phi_mtx, delta) +
        #                     einsum('jk,il->ijkl', phi_mtx, delta) +
        #                     einsum('jl,ik->ijkl', phi_mtx, delta))

        M_dev_ijkl = 0.5 * (0.5 * (einsum('ik,jl->ijkl', delta, phi_mtx) +
                                   einsum('il,jk->ijkl', delta, phi_mtx)) +
                            0.5 * (einsum('ik,jl->ijkl', phi_mtx, delta) +
                                   einsum('il,jk->ijkl', phi_mtx, delta)))

        # print 'M_dev_ijkl', M_dev_ijkl

        return M_dev_ijkl

    def _get_eps_pi_arr(self, sctx, eps_app_eng, sigma_kk):
        # Returns a list of the sliding strain vector for all microplanes.
        eps_pi_arr = self._get_state_variables(
            sctx, eps_app_eng, sigma_kk)[:, 5:9]

        return eps_pi_arr

    def _get_eps_pi_mtx(self, sctx, eps_app_eng, sigma_kk):

        # Vector integration of sliding (inelastic) strain for each microplane

        eps_pi_ni = self._get_eps_pi_arr(sctx, eps_app_eng, sigma_kk)

        eps_pi_ij = 0.5 * (einsum('n,ni,nj -> ij', self._MPW, eps_pi_ni, self._MPN) +
                           einsum('n,nj,ni -> ij', self._MPW, eps_pi_ni, self._MPN))

        return eps_pi_ij

    #-------------------------------------------------------------------------
    # Secant stiffness (irreducible decomposition based on ODFs)
    #-------------------------------------------------------------------------

    def _get_S_tns(self, sctx, eps_app_eng, sigma_kk):

        K0 = self.E / (1. - 2. * self.poisson)
        G0 = self.E / (1. + self.poisson)

        I_vol_ijkl = self._get_I_vol_4()
        I_dev_ijkl = self._get_I_dev_4()

        # print 'I_vol_ijkl', I_vol_ijkl
        # print 'I_dev_ijkl', I_dev_ijkl

        phi_mtx = self._get_phi_mtx(sctx, eps_app_eng, sigma_kk)

        M_vol_ijkl = self._get_M_vol_tns(sctx, eps_app_eng, sigma_kk)
        M_dev_ijkl = self._get_M_dev_tns(phi_mtx)

        S_ijkl = K0 * einsum('ijmn,mnrs,rskl -> ijkl', I_vol_ijkl, M_vol_ijkl, I_vol_ijkl ) \
            + G0 * einsum('ijmn,mnrs,rskl -> ijkl', I_dev_ijkl, M_dev_ijkl, I_dev_ijkl)\

        S_0_ijkl = K0 * I_vol_ijkl + G0 * I_dev_ijkl

        print 'S_0_ijkl', S_0_ijkl

        return S_ijkl

    #-------------------------------------------------------------------------
    # Evaluation - get the corrector and predictor
    #-------------------------------------------------------------------------

    def get_corr_pred(self, sctx, eps_app_eng, sigma_kk):
        '''
        Corrector predictor computation.
        @param eps_app_eng input variable - engineering strain
        '''

        # -----------------------------------------------------------------------------------------------
        # for debugging purposes only: if elastic_debug is switched on, linear elastic material is used
        # -----------------------------------------------------------------------------------------------
        if self.elastic_debug:
            # NOTE: This must be copied otherwise self.D2_e gets modified when
            # essential boundary conditions are inserted
            D2_e = copy(self.D2_e)
            sig_eng = tensordot(D2_e, eps_app_eng, [[1], [0]])
            return sig_eng, D2_e

        #----------------------------------------------------------------------
        # if the regularization using the crack-band concept is on calculate the
        # effective element length in the direction of principle strains
        #----------------------------------------------------------------------
        # if self.regularization:
        #    h = self.get_regularizing_length(sctx, eps_app_eng)
        #    self.phi_fn.h = h

        #----------------------------------------------------------------------
        # Return stresses (corrector) and damaged secant stiffness matrix (predictor)
        #----------------------------------------------------------------------
        eps_pi_ij = self._get_eps_pi_mtx(sctx, eps_app_eng, sigma_kk)
        eps_e_mtx = eps_app_eng - eps_pi_ij

        S_ijkl = self._get_S_tns(sctx, eps_app_eng, sigma_kk)

        print 'S_ijkl', S_ijkl

        sig_eng = einsum('ijkl,kl -> ij', S_ijkl, eps_e_mtx)

        return sig_eng


class MATS3DMicroplaneDamage(MATSXDMicroplaneDamageFatigue, MATS3DEval):

    implements(IMATSEval)

    # number of spatial dimensions
    n_dim = Constant(3)

    # number of components of engineering tensor representation
    n_eng = Constant(6)

    # number of microplanes - currently fixed for 3D
    n_mp = Constant(28)

    # get the normal vectors of the microplanes

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
    s_levels = linspace(0, 0.002, 2)
    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= -1
    #s_levels.reshape(-1, 2)[:, 1] *= -1
    s_history = s_levels.flatten()
    s_arr = hstack([linspace(s_history[i], s_history[i + 1], n)
                    for i in range(len(s_levels) - 1)])

    # print s_arr
    m = 0.0
    eps = array([array([[s_arr[i], m * s_arr[i], 0],
                        [m * s_arr[i], 0, 0],
                        [0, 0, 0]]) for i in range(0, len(s_arr))])

    # print eps.shape
    # print eps[i, 0, 0]
    # print len(eps[:, 0, 0])

    m1 = MATSEvalMicroplaneFatigue()
    m2 = MATS3DMicroplaneDamage()
    sigma = zeros_like(eps)
    sigma_kk = zeros(len(s_arr) + 1)
    w = zeros((len(eps[:, 0, 0]), 28))
    e_pi = zeros((len(eps[:, 0, 0]), 28))
    e = zeros((len(eps[:, 0, 0]), 28, 3))
    e_T = zeros((len(eps[:, 0, 0]), 28, 3))
    e_N = zeros((len(eps[:, 0, 0]), 28))
    sctx = zeros((len(eps[:, 0, 0]) + 1, 28, 8))
    d = zeros(len(eps[:, 0, 0]))

    for i in range(0, len(eps[:, 0, 0])):

        sigma[i, :] = m2.get_corr_pred(sctx[i, :], eps[i, :], sigma_kk[i])
        sigma_kk[i + 1] = trace(sigma[i, :])

        sctx[
            i + 1] = m2._get_state_variables(sctx[i, :], eps[i, :], sigma_kk[i])
        w[i, :] = sctx[i, :, 0]
        e_pi[i, :] = sctx[i, :, 5]

        d[i] = m2._get_d_scalar(sctx[i, :], eps[i, :], sigma_kk[i])
        e[i, :] = m2._get_e_vct_arr(eps[i, :])
        e_T[i, :] = m2._get_e_T_vct_arr(e[i, :])
        #e_N[i, :] = m2._get_e_N_arr(e[i, :])

    # print sigma_kk
    print e_T.shape

    plt.subplot(221)
    plt.plot(eps[:, 0, 0], sigma[:, 0, 0], linewidth=1, label='sigma_11')
    plt.plot(eps[:, 0, 0], sigma[:, 1, 1], linewidth=1, label='sigma_22')
    plt.plot(eps[:, 0, 0], sigma[:, 0, 1], linewidth=1, label='sigma_12')
    #plt.plot(eps[:, 0, 0], sigma[:, 2, 2], linewidth=1, label='sigma_33')
    plt.xlabel('Strain')
    plt.ylabel('Stress(MPa)')
    plt.legend()

    plt.subplot(222)
    for i in range(0, 28):
        plt.plot(
            eps[:, 0, 0], w[:, i], linewidth=1, label='Damage', alpha=1)

        plt.xlabel('Strain')
        plt.ylabel('Damage')
        plt.legend()

    plt.subplot(223)
    for i in range(0, 28):
        plt.plot(
            eps[:, 0, 0], e[:, i, 0], linewidth=1, label='Tangential_strain')

        plt.xlabel('Strain')
        plt.ylabel('Tangential_strain')
        # plt.legend()

    plt.subplot(224)

    plt.plot(eps[:, 0, 0], d, linewidth=1, label='macro damage')

    plt.xlabel('Strain')
    plt.ylabel('Normal_strain')

    plt.show()
