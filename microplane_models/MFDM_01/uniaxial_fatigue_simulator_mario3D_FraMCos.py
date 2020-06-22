#!/usr/bin/env python
# coding: utf-8

# # Simulation of fatigue for uniaxial stress state

# Assume a uniaxial stress state with $\sigma_{11} = \bar{\sigma}(\theta)$ representing the loading function. All other components of the stress tensor are assumed zero
# \begin{align}
# \sigma_{ab} = 0; \forall a,b \in (0,1,2), a = b \neq 1
# \end{align}
#

# In[1]:


import os

import matplotlib

from apps.sandbox.mario.Framcos.Micro2Dplot import Micro2Dplot
from apps.sandbox.mario.Framcos.vmats3D_mpl_csd_eeq import MATS3DMplCSDEEQ
#import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DELTA = np.identity(3)

EPS = np.zeros((3, 3, 3), dtype='f')
EPS[(0, 1, 2), (1, 2, 0), (2, 0, 1)] = 1
EPS[(2, 1, 0), (1, 0, 2), (0, 2, 1)] = -1


DD = np.hstack([DELTA, np.zeros_like(DELTA)])
EEPS = np.hstack([np.zeros_like(EPS), EPS])

GAMMA = np.einsum(
    'ik,jk->kij', DD, DD
) + np.einsum(
    'ikj->kij', np.fabs(EEPS)
)


def get_eps_ab(eps_O): return np.einsum(
    'Oab,...O->...ab', GAMMA, eps_O
)[np.newaxis, ...]


GAMMA_inv = np.einsum(
    'aO,bO->Oab', DD, DD
) + 0.5 * np.einsum(
    'aOb->Oab', np.fabs(EEPS)
)


def get_sig_O(sig_ab): return np.einsum(
    'Oab,...ab->...O', GAMMA_inv, sig_ab
)[0, ...]


GG = np.einsum(
    'Oab,Pcd->OPabcd', GAMMA_inv, GAMMA_inv
)


def get_K_OP(D_abcd):
    return np.einsum(
        'OPabcd,abcd->OP', GG, D_abcd
    )

# The above operators provide the three mappings
# map the primary variable to from vector to field
# map the residuum field to evctor (assembly operator)
# map the gradient of the residuum field to system matrix


m = MATS3DMplCSDEEQ()
plot = Micro2Dplot()


def get_UF_t(F, n_t, load, factor_H1, factor_H2, factor_L):

    n_mp = 61
    omegaN = np.zeros((n_mp, ))
    z_N_Emn = np.zeros((n_mp, ))
    alpha_N_Emn = np.zeros((n_mp, ))
    r_N_Emn = np.zeros((n_mp, ))
    eps_N_p_Emn = np.zeros((n_mp, ))
    sigma_N_Emn = np.zeros((n_mp, ))
    Y_n = np.zeros((n_mp, ))
    R_n = np.zeros((n_mp, ))
    w_T_Emn = np.zeros((n_mp, ))
    z_T_Emn = np.zeros((n_mp, ))
    alpha_T_Emna = np.zeros((n_mp, 3))
    eps_T_pi_Emna = np.zeros((n_mp, 3))
    sigma_T_Emna = np.zeros((n_mp, 3))
    X_T = np.zeros((n_mp, 3))
    Y_T = np.zeros((n_mp, ))
    sctx = np.zeros((n_mp, 23))
    sctx = sctx[np.newaxis, :, :]
    D = np.zeros((3, 3, 3, 3))
    D = D[np.newaxis, :, :, :, :]

    # total number of DOFs
    n_O = 6
    # Global vectors
    F_ext = np.zeros((n_O,), np.float_)
    F_O = np.zeros((n_O,), np.float_)
    U_k_O = np.zeros((n_O,), dtype=np.float_)
    eps_aux = get_eps_ab(U_k_O)
    # Setup the system matrix with displacement constraints
    # Time stepping parameters
    t_aux, t_n1, t_max, t_step = 0, 0, len(F), 1 / n_t
    # Iteration parameters
    k_max, R_acc = 1000, 1e-3
    # Record solutions
    U_t_list, F_t_list = [np.copy(U_k_O)], [np.copy(F_O)]

    # Load increment loop
    while t_n1 <= t_max - 1:
        #print('t:', t_n1)
        F_ext[0] = F[t_n1]
        F_ext[1] = -1. * F[t_n1]
        F_ext[2] = 0. * F[t_n1]

        k = 0
        # Equilibrium iteration loop
        while k < k_max:
            # Transform the primary vector to field
            eps_ab = get_eps_ab(U_k_O).reshape(3, 3)
            # Stress and material stiffness

            sig_ab, D_abcd = m.get_corr_pred(
                eps_ab, 1, omegaN, z_N_Emn,
                alpha_N_Emn, r_N_Emn, eps_N_p_Emn, sigma_N_Emn,
                w_T_Emn, z_T_Emn, alpha_T_Emna, eps_T_pi_Emna, eps_aux
            )
            D_abcd = D_abcd.reshape(3, 3, 3, 3)
            # Internal force
            F_O = get_sig_O(sig_ab).reshape(6,)
            # Residuum
            R_O = F_ext - F_O
            # System matrix
            K_OP = get_K_OP(D_abcd)
            # Convergence criterion
            R_norm = np.linalg.norm(R_O)
            delta_U_O = np.linalg.solve(K_OP, R_O)
            U_k_O += delta_U_O
            if R_norm < R_acc:
                # Convergence reached
                break
            # Next iteration
            k += 1

        else:
            print('no convergence')

            break

        # Update states variables after convergence
        [omegaN, z_N_Emn, alpha_N_Emn, r_N_Emn, eps_N_p_Emn, sigma_N_Emn, Y_n, R_n, w_T_Emn, z_T_Emn,
            alpha_T_Emna, eps_T_pi_Emna, sigma_T_Emna, Y_T, X_T] = m._get_state_variables(
                eps_ab, 1, omegaN, z_N_Emn, alpha_N_Emn, r_N_Emn, eps_N_p_Emn, sigma_N_Emn,
                w_T_Emn, z_T_Emn, alpha_T_Emna, eps_T_pi_Emna, eps_aux)

        omegaN = omegaN.reshape(n_mp, )
        z_N_Emn = z_N_Emn.reshape(n_mp, )
        alpha_N_Emn = alpha_N_Emn.reshape(n_mp, )
        r_N_Emn = r_N_Emn.reshape(n_mp, )
        eps_N_p_Emn = eps_N_p_Emn.reshape(n_mp, )
        sigma_N_Emn = sigma_N_Emn.reshape(n_mp,)
        Y_n = Y_n.reshape(n_mp,)
        R_n = R_n.reshape(n_mp,)
        w_T_Emn = w_T_Emn.reshape(n_mp, )
        z_T_Emn = z_T_Emn.reshape(n_mp, )
        alpha_T_Emna = alpha_T_Emna.reshape(n_mp, 3)
        eps_T_pi_Emna = eps_T_pi_Emna.reshape(n_mp, 3)
        sigma_T_Emna = sigma_T_Emna.reshape(n_mp, 3)
        X_T = X_T.reshape(n_mp, 3)
        Y_T = Y_T.reshape(n_mp, )

        sctx_aux = np.concatenate((omegaN.reshape(n_mp, 1), z_N_Emn.reshape(n_mp, 1), alpha_N_Emn.reshape(n_mp, 1),
                                   r_N_Emn.reshape(n_mp, 1), eps_N_p_Emn.reshape(
                                       n_mp, 1), sigma_N_Emn.reshape(n_mp, 1), Y_n.reshape(n_mp, 1), R_n.reshape(n_mp, 1),
                                   w_T_Emn.reshape(n_mp, 1), z_T_Emn.reshape(n_mp, 1), alpha_T_Emna, eps_T_pi_Emna, sigma_T_Emna, Y_T.reshape(n_mp, 1), X_T), axis=1)

        sctx_aux = sctx_aux[np.newaxis, :, :]
        sctx = np.concatenate((sctx, sctx_aux))
        D_aux = D_abcd[np.newaxis, :, :, :, :]
        D = np.concatenate((D, D_aux))
        U_t_list.append(np.copy(U_k_O))
        F_t_list.append(F_O)
        eps_aux = get_eps_ab(U_k_O)

        t_n = t_n1
        t_n1 += 1

    U_t, F_t = np.array(U_t_list), np.array(F_t_list)
    return U_t, F_t, t_n1 / t_max, t_aux, D


# load = -91.14371574719883

# load = -118.16283056375052

load = 60.54301292467442

# load = -64.895457863834174

l_H1 = 0.5
cycles1 = 2618

l_H2 = 0.5
cycles2 = 30

factor_H1 = -1.5
factor_H2 = 0.8
factor_L = 0.1

max_load1 = load * factor_H1
max_load2 = load * 0.55
max_load3 = load * 0.60
max_load4 = load * 0.65
max_load5 = load * 0.70
max_load6 = load * 0.75
max_load7 = load * 0.80
max_load8 = load * 0.85
max_load9 = load * 0.90
max_load10 = load * 0.95
max_load11 = load * 1.0
min_load = load * factor_L


# n_cycles1 = 88
# n_cycles2 = 99912

n_cycles1 = 1
n_cycles2 = 10
n_cycles3 = 10
n_cycles4 = 10
n_cycles5 = 10
n_cycles6 = 10
n_cycles7 = 10
n_cycles8 = 10
n_cycles9 = 10
n_cycles10 = 10
n_cycles11 = 10


t_steps_cycle = 2000
monotonic = np.linspace(0, max_load1, t_steps_cycle)

first_load = np.concatenate((np.linspace(0, max_load1, t_steps_cycle), np.linspace(
    max_load1, min_load, t_steps_cycle)[1:]))
cycle1 = np.concatenate((np.linspace(min_load, max_load1, t_steps_cycle)[1:], np.linspace(max_load1, min_load, t_steps_cycle)[
                        1:]))
cycle1 = np.tile(cycle1, n_cycles1 - 1)

cycle2 = np.concatenate((np.linspace(min_load, max_load2, t_steps_cycle)[1:], np.linspace(max_load2, min_load, t_steps_cycle)[
                        1:]))
cycle2 = np.tile(cycle2, n_cycles2)

cycle3 = np.concatenate((np.linspace(min_load, max_load3, t_steps_cycle)[1:], np.linspace(max_load3, min_load, t_steps_cycle)[
                        1:]))
cycle3 = np.tile(cycle3, n_cycles3)

cycle4 = np.concatenate((np.linspace(min_load, max_load4, t_steps_cycle)[1:], np.linspace(max_load4, min_load, t_steps_cycle)[
                        1:]))
cycle4 = np.tile(cycle4, n_cycles4)

cycle5 = np.concatenate((np.linspace(min_load, max_load5, t_steps_cycle)[1:], np.linspace(max_load5, min_load, t_steps_cycle)[
                        1:]))
cycle5 = np.tile(cycle5, n_cycles5)

cycle6 = np.concatenate((np.linspace(min_load, max_load6, t_steps_cycle)[1:], np.linspace(max_load6, min_load, t_steps_cycle)[
                        1:]))
cycle6 = np.tile(cycle6, n_cycles6)

cycle7 = np.concatenate((np.linspace(min_load, max_load7, t_steps_cycle)[1:], np.linspace(max_load7, min_load, t_steps_cycle)[
                        1:]))
cycle7 = np.tile(cycle7, n_cycles7)

cycle8 = np.concatenate((np.linspace(min_load, max_load8, t_steps_cycle)[1:], np.linspace(max_load8, min_load, t_steps_cycle)[
                        1:]))
cycle8 = np.tile(cycle8, n_cycles8)

cycle9 = np.concatenate((np.linspace(min_load, max_load9, t_steps_cycle)[1:], np.linspace(max_load9, min_load, t_steps_cycle)[
                        1:]))
cycle9 = np.tile(cycle9, n_cycles9)

cycle10 = np.concatenate((np.linspace(min_load, max_load10, t_steps_cycle)[1:], np.linspace(max_load10, min_load, t_steps_cycle)[
    1:]))
cycle10 = np.tile(cycle10, n_cycles10)

cycle11 = np.concatenate((np.linspace(min_load, max_load11, t_steps_cycle)[1:], np.linspace(max_load11, min_load, t_steps_cycle)[
    1:]))
cycle11 = np.tile(cycle11, n_cycles11)

sin_load = np.concatenate((first_load, cycle1, cycle2, cycle3,
                           cycle4, cycle5, cycle6, cycle7, cycle8, cycle9, cycle10, cycle11))

sin_load = np.concatenate((first_load, cycle1))

sin_load = monotonic


t_steps = len(sin_load)

t = np.linspace(0, 1, len(sin_load))


U, F, cyc, number_cyc, D = get_UF_t(
    sin_load, t_steps, load, factor_H1, factor_H2, factor_L)

font = {'family': 'normal',
        'size': 18}

matplotlib.rc('font', **font)

f, (ax1) = plt.subplots(1, 1, figsize=(5, 4))

ax1.plot(t[0:], np.abs(sin_load / load)[0:], linewidth=2.5)
ax1.set_xlabel('pseudotime [-]', fontsize=25)
ax1.set_ylabel(r'$|S_{max}$| [-]', fontsize=25)
ax1.set_title('L-H')

print(np.max(np.abs(F[:, 0])), 'sigma1')
print(np.max(np.abs(F[:, 1])), 'sigma2')

f, (ax2) = plt.subplots(1, 1, figsize=(5, 4))

ax2.plot(np.abs(U[:, 0]), np.abs(F[:, 0]), linewidth=2.5)
ax2.set_xlabel(r'$|\varepsilon_{11}$|', fontsize=25)
ax2.set_ylabel(r'$|\sigma_{11}$| [-]', fontsize=25)
ax2.set_title(str((n_cycles1)) + ',' + str(cyc))
# ax2.set_ylim(0.00, 140)
# ax2.set_xlim(-0.0005, 0.00333)
plt.show()


f, (ax) = plt.subplots(1, 1, figsize=(5, 4))

# ax.plot(np.arange(len(U[(t_steps_cycle)::2 * (t_steps_cycle - 1), 0])) + 1,
# np.abs(U[(t_steps_cycle)::2 * (t_steps_cycle - 1), 0]), linewidth=2.5)

ax.plot(np.arange(len(U[2::2, 0])) + 1,
        np.abs(U[2::2, 0]), linewidth=2.5)

#ax.set_xlim(0, ((n_cycles1)) + 1)
ax.set_xlabel('number of cycles [N]', fontsize=25)
ax.set_ylabel(r'$|\varepsilon_{11}^{max}$|', fontsize=25)
plt.title('creep fatigue Smax = 0.85')
plt.show()

# print(U[:, 0])
# print(U[(t_steps_cycle)::2 * (t_steps_cycle - 1), 0])


f, (ax) = plt.subplots(1, 1, figsize=(5, 4))


ax.plot((np.arange(len(U[(t_steps_cycle)::2 * (t_steps_cycle - 1), 0])) + 1),
        np.abs(U[(t_steps_cycle)::2 * (t_steps_cycle - 1), 0]), linewidth=2.5)
ax.plot((np.arange(len(U[0::2 * (t_steps_cycle - 1), 0])) + 1),
        np.abs(U[1::2 * (t_steps_cycle - 1), 0]), linewidth=2.5)

# ax.plot((np.arange(len(U[(t_steps_cycle)::2 * (t_steps_cycle - 1), 0])) + 1) / len(U[(t_steps_cycle)::2 * (t_steps_cycle - 1), 0]),
#         np.abs(U[(t_steps_cycle)::2 * (t_steps_cycle - 1), 0]), linewidth=2.5)
# ax.plot((np.arange(len(U[0::2 * (t_steps_cycle - 1), 0])) + 1) / len(U[0::2 * (t_steps_cycle - 1), 0]),
#         np.abs(U[1::2 * (t_steps_cycle - 1), 0]), linewidth=2.5)

# ax.plot((np.arange(len(U[2::2, 0])) + 1) / len(U[1::2, 0]),
#         np.abs(U[2::2, 0]), linewidth=2.5)
# ax.plot((np.arange(len(U[1::2, 0])) + 1) / len(U[0::2, 0]),
#         np.abs(U[1::2, 0]), linewidth=2.5)


# X_axis1 = (np.arange(l_H1 * cycles1) + 1) / cycles1
# Y_axis1 = np.abs(U[2:np.int(2 * l_H1 * cycles1) + 2:2, 0])
#
# X_axis2 = (np.arange(len(U[2::2, 0]) -
#                      (l_H1 * cycles1)) + 1) / (cycles2) + l_H1
# Y_axis2 = np.abs(U[np.int(2 * l_H1 * cycles1) + 2::2, 0])
#
# ax.plot(X_axis1, Y_axis1, 'k', linewidth=2.5)
# ax.plot(X_axis2, Y_axis2, 'k', linewidth=2.5)
# ax.plot([X_axis1[-1], X_axis2[0]],
#         [Y_axis1[-1], Y_axis2[0]], 'k', linewidth=2.5)


# ax.plot(np.arange(cycles1) / l_H1 * cycles1,
#         np.abs(U[2::2, 0]), linewidth=2.5)


ax.set_ylim(0.00, 0.004)
ax.set_xlim(-2.00, 100)
#ax.set_xlim(0, ((n_cycles1 + n_cycles2)) + 1)
ax.set_xlabel('N/Nf', fontsize=25)
ax.set_ylabel('strain', fontsize=25)
plt.title('creep fatigue Smax = 0.85')
plt.show()
