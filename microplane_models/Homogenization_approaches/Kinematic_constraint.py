'''
Created on 07.06.2020
@author: Abdul
plotting tool for microplane models
'''

from math import \
    cos, sin, pi

from numpy import \
    array, zeros, trace, \
    einsum, zeros_like,\
    identity, linspace, hstack,\
    sqrt, ones
    
from traits.api import \
    Constant, \
    HasTraits, \
    Property, cached_property


import matplotlib.pyplot as plt

import numpy as np


class MATS2DMicroplane(HasTraits):

    #-----------------------------------------------
    # number of microplanes
    #-----------------------------------------------
    n_mp = Constant(360)

    #-----------------------------------------------
    # get the normal vectors of the microplanes
    #-----------------------------------------------
    _MPN = Property(depends_on='n_mp')

    @cached_property
    def _get__MPN(self):
        # microplane normals:

        alpha_list = linspace(0, 2 * pi, self.n_mp)

        MPN = array([[cos(alpha), sin(alpha)] for alpha in alpha_list])

        return MPN

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
        MPW = ones(self.n_mp) / self.n_mp * 2

        return MPW

    #-------------------------------------------------------------------------
    # Cached elasticity tensors
    #-------------------------------------------------------------------------

    elasticity_tensors = Property(
        depends_on='E, nu, dimensionality, stress_state')

    @cached_property
    def _get_elasticity_tensors(self):
        '''
        Intialize the fourth order elasticity tensor for 2D or 2D plane strain or 2D plane stress
        '''
        # ----------------------------------------------------------------------------
        # Lame constants calculated from E and nu
        # ----------------------------------------------------------------------------

        # first Lame paramter
        la = self.E * self.nu / ((1 + self.nu) * (1 - 2 * self.nu))
        # second Lame parameter (shear modulus)
        mu = self.E / (2 + 2 * self.nu)

        # -----------------------------------------------------------------------------------------------------
        # Get the fourth order elasticity and compliance tensors for the 2D-case
        # -----------------------------------------------------------------------------------------------------

        # construct the elasticity tensor (using Numpy - einsum function)
        delta = identity(2)
        D_ijkl = (einsum(',ij,kl->ijkl', la, delta, delta) +
                  einsum(',ik,jl->ijkl', mu, delta, delta) +
                  einsum(',il,jk->ijkl', mu, delta, delta))

        return D_ijkl


    #-------------------------------------------------------------------------
    # Setup for computation within a supplied spatial context
    #-------------------------------------------------------------------------
    D4_e = Property

    def _get_D4_e(self):
        # Return the elasticity tensor
        return self.elasticity_tensors

    #-------------------------------------------------------------------------
    '''Kinematic constraint N-T split '''
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
        delta = identity(2)
        MPTT_nijr = 0.5 * (einsum('ni,jr -> nijr', self._MPN, delta) +
                           einsum('nj,ir -> njir', self._MPN, delta) - 2 *
                           einsum('ni,nj,nr -> nijr', self._MPN, self._MPN, self._MPN))
        return MPTT_nijr

    def _get_e_N_arr(self, eps_eng):
        return einsum('nij,ij->n', self._MPNN, eps_eng)

    def _get_e_T_vct_arr(self, eps_eng):
        MPTT_ijr = self._get__MPTT()
        return einsum('nijr,ij->nr', MPTT_ijr, eps_eng)
    
    
    #-------------------------------------------------------------------------
    '''Kinematic constraint V-D-T split '''
    #-------------------------------------------------------------------------


    # get the dyadic product of the microplane normals
    _MPVV_1 = Property(depends_on='n_mp')

    @cached_property
    def _get__MPVV_1(self):
        # dyadic product of the microplane normals
        I_ij = identity(2)
        MPVV_ij = 1/3 * I_ij 
        return MPVV_ij

    
    # get the dyadic product of the microplane normals
    _MPDD_1 = Property(depends_on='n_mp')

    @cached_property
    def _get__MPDD_1(self):
        # dyadic product of the microplane normals
        I_ij = identity(2)
        MPVV_ij = 1/3 * I_ij 
        
        MPNN_nij = einsum('ni,nj->nij', self._MPN, self._MPN)
        
        MPDD_ij = einsum('nij->nij', MPNN_nij - MPVV_ij) 
        
        return MPDD_ij 

    # get the third order tangential tensor (operator) for each microplane
    _MPTT_1 = Property(depends_on='n_mp')

    @cached_property
    def _get__MPTT_1(self):
        # Third order tangential tensor for each microplane
        delta = identity(2)
        MPTT_nijr = 0.5 * (einsum('ni,jr -> nijr', self._MPN, delta) +
                           einsum('nj,ir -> njir', self._MPN, delta) - 2 *
                           einsum('ni,nj,nr -> nijr', self._MPN, self._MPN, self._MPN))
        return MPTT_nijr

    def _get_e_V1_arr(self, eps_eng):
        return einsum('ij,ij->', self._MPVV_1, eps_eng)

    def _get_e_D1_arr(self, eps_eng):
        return einsum('nij,ij->n', self._MPDD_1, eps_eng)
    
    
    def _get_e_T1_vct_arr(self, eps_eng):
        MPTT_ijr = self._get__MPTT_1()
        return einsum('nijr,ij->nr', MPTT_ijr, eps_eng)
    

    #-------------------------------------------------------------------------
    '''Kinematic constraint V-D split '''
    #-------------------------------------------------------------------------

    # get the dyadic product of the microplane normals
    _MPVV = Property(depends_on='n_mp')

    @cached_property
    def _get__MPVV(self):
        # dyadic product of the microplane normals
        I_ij = identity(2)
        MPVV_ij = 1/3 * I_ij 
        return MPVV_ij

    # get the third order tangential tensor (operator) for each microplane
    _MPDD = Property(depends_on='n_mp')

    @cached_property
    def _get__MPDD(self):
        # Third order tangential tensor for each microplane
        delta = identity(2)
        MPDD_nijr = 0.5 * (einsum('ni,jr -> nijr', self._MPN, delta) +
                           einsum('nj,ir -> njir', self._MPN, delta)) -\
                             1/3. * einsum('nk,ki,jr -> nijr', self._MPN, delta, delta)
        return MPDD_nijr

    def _get_e_V_arr(self, eps_eng):
        return einsum('ij,ij->', self._MPVV, eps_eng)

    def _get_e_D_vct_arr(self, eps_eng):
        MPDD_ijr = self._get__MPDD()
        return einsum('nijr,ij->nr', MPDD_ijr, eps_eng)
    
    


if __name__ == '__main__':
    #==========================================================================
    # Check the model behavior at the single material point
    #==========================================================================

    model = MATS2DMicroplane()

    n_mp = model.n_mp

    p = 1.0  # ratio of strain eps_11 (for bi-axial loading)
    m = -0.2  # ratio of strain eps_22 (for bi-axial loading)

    # monotonic loading - comp
    n = 25  # number of increments
    s_levels = linspace(0, 0.001, 2)
    #s_levels[0] = 0
    eps_max = 0.001

    s_levels.reshape(-1, 2)[:, 0] = eps_max
    s_levels[0] = 0
    #s_levels.reshape(-1, 2)[:, 1] = eps_min
    s_history_1 = s_levels.flatten()
    s_arr_1 = hstack([linspace(s_history_1[i], s_history_1[i + 1], n)
                      for i in range(len(s_levels) - 1)])

    eps_1 = array([array([[0 * s_arr_1[i],  p* s_arr_1[i]],
                          [ p* s_arr_1[i] , 0 * s_arr_1[i]]]) for i in range(0, len(s_arr_1))])

    idx = np.where(eps_1[:, 0, 0] == eps_max)
    print(idx[0])

    #--------------------------------------
    # construct the arrays
    #--------------------------------------

    sigma_kk_1 = zeros(len(s_arr_1) + 1)

    eps_N_1 = zeros((len(eps_1[:, 0]), n_mp))
    eps_V_1 = zeros((len(eps_1[:, 0]), n_mp))
    eps_D_1 = zeros((len(eps_1[:, 0]), n_mp))

    eps_TT_1 = zeros((len(eps_1[:, 0]), n_mp, 2))
    eps_DD_1 = zeros((len(eps_1[:, 0]), n_mp, 2))


    for i in range(0, len(eps_1[:, 0])):
        
        eps_N_1[i, :] = model._get_e_N_arr(eps_1[i, :])
        eps_TT_1[i, :, :] = model._get_e_T_vct_arr(eps_1[i, :])
        
        eps_V_1[i, :] = model._get_e_V1_arr(eps_1[i, :])
        eps_D_1[i, :] = model._get_e_D1_arr(eps_1[i, :])
        
        eps_V_1[i, :] = model._get_e_V_arr(eps_1[i, :])
        eps_DD_1[i, :, :] = model._get_e_D_vct_arr(eps_1[i, :])



#     '''====================================================
#     plotting
#     ===================================================='''


    rads = np.arange(0, (2 * np.pi), (2 * np.pi) / n_mp)

    '=================================================================='
    'Plotting in increments'
    '=================================================================='
    
    #===================
    # Normal strain
    #===================
    plt.subplot(111)

    for i in range(0, len(eps_1[:, 0])):
        plt.polar(rads, eps_N_1[i, :])
    plt.title(' Normal strain for all microplanes')
    plt.show()
    
    #===================
    # Tangential  strain
    #===================
    plt.subplot(111)
    for i in range(0, len(eps_1[:, 0])):

        norm = sqrt(
            einsum('...i,...i->... ', eps_TT_1[i, :], eps_TT_1[i, :]))

        plt.polar(rads, norm)
    plt.title('Tangential strain for all microplanes')
    plt.show()
    
    #===================
    # Volumetric strain
    #===================
    plt.subplot(111)

    for i in range(0, len(eps_1[:, 0])):
        plt.polar(rads, eps_V_1[i, :])
    plt.title(' Volumetric strain for all microplanes')
    plt.show()
    
    

    #===================
    # Deviatoric strain (V-D-T)
    #===================
    plt.subplot(111)

    for i in range(0, len(eps_1[:, 0])):
        plt.polar(rads, eps_D_1[i, :])
    plt.title(' Deviatoric (V-D-T) strain for all microplanes')
    plt.show()

    #===================
    # Deviatoric  strain
    #===================
    plt.subplot(111)
    for i in range(0, len(eps_1[:, 0])):

        norm = sqrt(
            einsum('...i,...i->... ', eps_DD_1[i, :], eps_DD_1[i, :]))

        plt.polar(rads, norm)
    plt.title('Deviatoric (V-D) strain for all microplanes')
    plt.show()

    plt.show()
