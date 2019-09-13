
'''
Created on 30.10.2018

@author: abaktheer
'''

from traits.api import HasTraits
from traits.api import \
    Float
import matplotlib.pyplot as plt
import numpy as np
import traits.api as tr


class MATS3DDesmorat(HasTraits):
    '''Damage - plasticity model by Desmorat.
    '''

    #-------------------------------------------------------------------------
    # Material parameters
    #-------------------------------------------------------------------------

    E_1 = tr.Float(20.0e+3,
                   label="E_1",
                   desc="first Young's Modulus",
                   auto_set=False,
                   input=True)
    E_2 = tr.Float(15.0e+3,
                   label="E_2",
                   desc="second Young's Modulus",
                   auto_set=False,
                   input=True)

    nu = tr.Float(0.2,
                  label='nu',
                  desc="Poisson ratio",
                  auto_set=False,
                  input=True)

    def _get_lame_1_params(self):
        la = self.E_1 * self.nu / ((1. + self.nu) * (1. - 2. * self.nu))
        # second Lame parameter (shear modulus)
        mu = self.E_1 / (2. + 2. * self.nu)
        return la, mu

    D_1_abef = tr.Property(tr.Array, depends_on='+input')

    @tr.cached_property
    def _get_D_1_abef(self):
        la = self._get_lame_1_params()[0]
        mu = self._get_lame_1_params()[1]
        delta = np.identity(3)
        D_1_abef = (np.einsum(',ij,kl->ijkl', la, delta, delta) +
                    np.einsum(',ik,jl->ijkl', mu, delta, delta) +
                    np.einsum(',il,jk->ijkl', mu, delta, delta))

        return D_1_abef

    def _get_lame_2_params(self):
        la = self.E_2 * self.nu / ((1. + self.nu) * (1. - 2. * self.nu))
        # second Lame parameter (shear modulus)
        mu = self.E_2 / (2. + 2. * self.nu)
        return la, mu

    D_2_abef = tr.Property(tr.Array, depends_on='+input')

    @tr.cached_property
    def _get_D_2_abef(self):
        la = self._get_lame_2_params()[0]
        mu = self._get_lame_2_params()[1]
        delta = np.identity(3)
        D_2_abef = (np.einsum(',ij,kl->ijkl', la, delta, delta) +
                    np.einsum(',ik,jl->ijkl', mu, delta, delta) +
                    np.einsum(',il,jk->ijkl', mu, delta, delta))

        return D_2_abef

    gamma = Float(0.0,
                  label="Gamma",
                  desc="kinematic hardening modulus",
                  MAT=True,
                  symbol=r'\gamma',
                  unit='MPa/mm',
                  enter_set=True,
                  auto_set=False)

    K = Float(0.0,
              label="K",
              desc="isotropic hardening modulus",
              MAT=True,
              symbol='K',
              unit='MPa/mm',
              enter_set=True,
              auto_set=False)

    S = Float(324.0e-6,
              label="S",
              desc="damage strength",
              MAT=True,
              symbol='S',
              unit='MPa/mm',
              enter_set=True,
              auto_set=False)

    tau_bar = Float(9.0,
                    label="Tau_0 ",
                    desc="yield stress",
                    symbol=r'\bar{\tau}',
                    unit='MPa',
                    MAT=True,
                    enter_set=True,
                    auto_set=False)

    def _get_state_variables(self, eps_Emab, eps_pi_Emab, alpha_Emab, z_Ema, omega_Ema):

        D_1_abef = self.D_1_abef
        D_2_abef = self.D_2_abef

        sigma_pi_Emab_trial = (
            np.einsum('ijkl,kl->ij', D_2_abef, eps_Emab - eps_pi_Emab))

        a = sigma_pi_Emab_trial - self.gamma * alpha_Emab

        norm_a = np.sqrt(np.einsum('ij,ij', a, a))

        if norm_a == 0:
            n = 0 * a
        else:
            n = a / norm_a

        f = np.sqrt(np.einsum('ij,ij', a, a)
                    ) - self.tau_bar - self.K * z_Ema

        plas_1 = f > 1e-6
        elas_1 = f < 1e-6

        delta_pi = f / \
            (self.E_2 + (self.K + self.gamma) * (1. - omega_Ema)) * plas_1

        b = 1.0 * elas_1 + norm_a * plas_1

#         delta_lamda = f / (np.einsum('...ij,...ijkl,...kl',
#                                      n, D_2_abef, n) + self.gamma * np.einsum('...ij,...ij', n, n) + self.K)
#
#         delta_pi = delta_lamda / (1.0 - omega_Ema)

        eps_pi_Emab = eps_pi_Emab + (a * delta_pi / b)

        eps_diff_Emab = eps_Emab - eps_pi_Emab

        Y_Ema = 0.5 * (np.einsum('ij,ijkl,kl',
                                 eps_Emab, D_1_abef, eps_Emab)) + \
            0.5 * (np.einsum('ij,ijkl,kl',
                             eps_diff_Emab, D_2_abef, eps_diff_Emab))

        omega_Ema = omega_Ema + (Y_Ema / self.S) * delta_pi

        if omega_Ema >= 0.99:
            omega_Ema = 0.99

        alpha_Emab = alpha_Emab + plas_1 * \
            (a * delta_pi / b) * (1.0 - omega_Ema)

        z_Ema = z_Ema + delta_pi * (1.0 - omega_Ema)

        return eps_pi_Emab, alpha_Emab, z_Ema, omega_Ema

    def get_corr_pred(self, eps_Emab_n,
                      sigma_Emab, eps_pi_Emab,
                      alpha_Emab, z_Ema, omega_Ema):
        r'''
        Corrector predictor computation.
        '''

        D_1_abef = self.D_1_abef
        D_2_abef = self.D_2_abef

        eps_pi_Emab, alpha_Emab, z_Ema, omega_Ema = self._get_state_variables(
            eps_Emab_n, eps_pi_Emab, alpha_Emab, z_Ema, omega_Ema)

        phi_Emn = 1.0 - omega_Ema

        sigma_Emab = phi_Emn * (np.einsum('ijkl,kl->ij',
                                          D_1_abef, eps_Emab_n) +
                                np.einsum('ijkl,kl->ij',
                                          D_2_abef, eps_Emab_n - eps_pi_Emab))

        D_abef = phi_Emn * np.einsum('ijkl->ijkl',
                                     D_1_abef + D_2_abef)

        return D_abef, sigma_Emab


if __name__ == '__main__':

    model = MATS3DDesmorat()
#------------------------------------------------------------------------------
# monotonic loading
#------------------------------------------------------------------------------
    n = 1000  # number of increments
    s_levels = np.linspace(0, 0.0028, 2)
    #s_levels[0] = 0
#     s_levels.reshape(-1, 2)[:, 0] = 0.0005
#     s_levels.reshape(-1, 2)[:, 1] = 0.002
    s_levels[0] = 0
    s_history_1 = s_levels.flatten()
    s_arr_1 = np.hstack([np.linspace(s_history_1[i], s_history_1[i + 1], n)
                         for i in range(len(s_history_1) - 1)])

    eps_1 = np.array([np.array([[s_arr_1[i], 0, 0],
                                [0, 0, 0],
                                [0, 0,  0]]) for i in range(0, len(s_arr_1))])

    #--------------------------------------
    # construct the arrays
    #--------------------------------------
    sigma_1 = np.zeros_like(eps_1)
    w_1 = np.zeros(len(eps_1[:, 0, 0]) + 1)
    z_1 = np.zeros(len(eps_1[:, 0, 0]) + 1)
    eps_pi_1 = np.zeros((len(eps_1[:, 0, 0]) + 1, 3, 3))
    alpha_1 = np.zeros((len(eps_1[:, 0, 0]) + 1, 3, 3))

    for i in range(0, len(eps_1[:, 0, 0])):

        sigma_1[i, :] = model.get_corr_pred(eps_1[i, :], sigma_1[i, :], eps_pi_1[i, :],
                                            alpha_1[i, :], z_1[i], w_1[i])[1]

        eps_pi_1[i + 1, :], alpha_1[i + 1, :], z_1[i + 1], w_1[i + 1] = model._get_state_variables(eps_1[i, :], eps_pi_1[i, :],
                                                                                                   alpha_1[i, :], z_1[i], w_1[i])

    print(sigma_1)
#------------------------------------------------------------------------------
# cyclic loading
#------------------------------------------------------------------------------
    n = 1000  # number of increments
    s_levels = np.linspace(0, 0.005, 2)
    #s_levels[0] = 0
#     s_levels.reshape(-1, 2)[:, 0] = 0.0005
#     s_levels.reshape(-1, 2)[:, 1] = 0.002

    s_levels[0] = 0
    s_history_2 = s_levels.flatten()
    s_history_2 = [0, 0.0020, 0.0006,
                   0.0028]

    s_arr_2 = np.hstack([np.linspace(s_history_2[i], s_history_2[i + 1], n)
                         for i in range(len(s_history_2) - 1)])

    eps_2 = np.array([np.array([[s_arr_2[i], 0, 0],
                                [0, 0, 0],
                                [0, 0,  0]]) for i in range(0, len(s_arr_2))])

    #--------------------------------------
    # construct the arrays
    #--------------------------------------
    sigma_2 = np.zeros_like(eps_2)
    w_2 = np.zeros(len(eps_2[:, 0, 0]) + 1)
    z_2 = np.zeros(len(eps_2[:, 0, 0]) + 1)
    eps_pi_2 = np.zeros((len(eps_2[:, 0, 0]) + 1, 3, 3))
    alpha_2 = np.zeros((len(eps_2[:, 0, 0]) + 1, 3, 3))

    for i in range(0, len(eps_2[:, 0, 0])):

        sigma_2[i, :] = model.get_corr_pred(eps_2[i, :], sigma_2[i, :], eps_pi_2[i, :],
                                            alpha_2[i, :], z_2[i], w_2[i])[1]

        eps_pi_2[i + 1, :], alpha_2[i + 1, :], z_2[i + 1], w_2[i + 1] = model._get_state_variables(eps_2[i, :], eps_pi_2[i, :],
                                                                                                   alpha_2[i, :], z_2[i], w_2[i])
    #------------------------------------------------------
    # stress -strain
    #------------------------------------------------------
    plt.subplot(221)
    plt.plot(eps_1[:, 0, 0], sigma_1[:, 0, 0], color='k',
             linewidth=1, label='sigma_11_(monotonic-compression)')
    plt.plot(eps_2[:, 0, 0], sigma_2[:, 0, 0], color='r',
             linewidth=1, label='sigma_11_(cyclic-compression)')

    plt.subplot(222)
    plt.plot(eps_1[:, 0, 0], w_1[:-1], color='k',
             linewidth=1, label='sigma_11_(monotonic-compression)')
    plt.plot(eps_2[:, 0, 0], w_2[:-1], color='r',
             linewidth=1, label='sigma_11_(cyclic-compression)')

    plt.show()
