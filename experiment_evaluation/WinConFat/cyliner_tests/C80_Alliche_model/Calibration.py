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
        #print('D=', D_i)

        # Helmholtz free energy
        phi_i = 0.5 * lamda * (eps_1_i + 2.0 * eps_2_i)**2.0 + mu * ((eps_1_i)**2.0 + 2.0 * eps_2_i**2.0) + 2.0 * g * w_i * eps_2_i + alpha * \
            (2.0 * w_i * eps_1_i * eps_2_i + 4.0 * w_i *
             eps_2_i**2.0) + 4.0 * beta * w_i * eps_2_i**2.0

#         if w_i > 5.0:
#             print(' ----------> No Convergence any more')
#             print(i)
#             break
        
        if d_w  > 0.5:
            print(' ----------> No Convergence any more')
            print(i)
            break

#         if abs(eps_1_i) > 0.005:
#             print(' ----------> No Convergence any more')
#             print(i)
#             break

        eps_1_arr[i] = eps_1_i
        eps_2_arr[i] = eps_2_i
        w_arr[i] = w_i
        f_arr[i] = f_i
        D_arr[i] = D_i
        phi_arr[i] = phi_i

    return sigma_1_arr, eps_1_arr, eps_2_arr, w_arr, f_arr, D_arr, i, phi_arr


if __name__ == '__main__':

    m = 200  # number of increments in each cycle

    n1 = 100000
    n2 = 100000

    sigma_u = - 99
    
    
    
#------------------------------------------------------------------------------ 
# simulation 1 (S=0.85)
#------------------------------------------------------------------------------ 
    stress_level_1_max = 0.85 * sigma_u
    stress_level_1_min = 0.2 * sigma_u

    d_0 = np.zeros(1)
    d_1 = np.linspace(0, stress_level_1_max, n1 * 2)
    d_1.reshape(-1, 2)[:, 0] = stress_level_1_max
    d_1.reshape(-1, 2)[:, 1] = stress_level_1_min
    d_history_1 = d_1.flatten()
    d_arr_1 = np.hstack((d_0, d_history_1))
    sigma_arr_1 = np.hstack([np.linspace(d_arr_1[i], d_arr_1[i + 1], m,  endpoint=False)
                           for i in range(len(d_arr_1) - 1)])

    t_arr_1 = np.linspace(0, 1, len(sigma_arr_1))
    
    
#------------------------------------------------------------------------------ 
# simulation 1 (S=0.75)
#------------------------------------------------------------------------------ 
    stress_level_2_max = 0.75 * sigma_u
    stress_level_2_min = 0.2 * sigma_u

    d_0 = np.zeros(1)
    d_2 = np.linspace(0, stress_level_2_max, n1 * 2)
    d_2.reshape(-1, 2)[:, 0] = stress_level_2_max
    d_2.reshape(-1, 2)[:, 1] = stress_level_2_min
    d_history_2 = d_2.flatten()
    d_arr_2 = np.hstack((d_0, d_history_2))
    sigma_arr_2 = np.hstack([np.linspace(d_arr_2[i], d_arr_2[i + 1], m,  endpoint=False)
                           for i in range(len(d_arr_2) - 1)])

    t_arr_2 = np.linspace(0, 1, len(sigma_arr_2))

    #================================================================

#     # C120
#     sigma_arr, eps_1_arr, eps_2_arr, w_arr, f_arr, D_arr, inc, phi_arr = get_stress_strain(
#         sigma_arr, lamda=12500, mu=18750, alpha=2237.5, beta=-2216.5, g=-10.0,
#         C0=0.00, C1=0.0019, K=0.00485, n=10)


#     # C80 - alliche paper
#     sigma_arr, eps_1_arr, eps_2_arr, w_arr, f_arr,  D_arr, inc, phi_arr = get_stress_strain(
#         sigma_arr, lamda=10555.55, mu=15833.33, alpha=2237.5, beta=-2216.5, g=-9.788, C0=0.00, C1=0.002033, K=0.003345, n=10)
    
    
    
#------------------------------------------------------------------------------ 
# material parameters
#------------------------------------------------------------------------------ 

    lamda=10555.55
    mu=15833.33 
    alpha=2000 
    beta=-2216.5 
    g=-9.8 
    C0=0.00
    C1=0.002
    K=0.005
    n=18
    
    
    # C80 - calibration
    sigma_arr_1, eps_1_arr_1, eps_2_arr_1, w_arr_1, f_arr_1,  D_arr_1, inc_1, phi_arr_1 = get_stress_strain(
        sigma_arr_1, lamda, mu, alpha, beta, g, C0, C1, K, n)
    
    sigma_arr_2, eps_1_arr_2, eps_2_arr_2, w_arr_2, f_arr_2,  D_arr_2, inc_2, phi_arr_2 = get_stress_strain(
        sigma_arr_2, lamda, mu, alpha, beta, g, C0, C1, K, n)
    
    
    Nf_1 = inc_1/(m*2)
    Nf_2 = inc_2/(m*2)
    
    
    
    
    

    

    #------------------------------------------------------------------------------ 
    # Experimental data - C80
    #------------------------------------------------------------------------------ 
     
    #------------------------------------------------------------------------------ 
    # fatigue creep curves (H)
    #------------------------------------------------------------------------------ 
    #=============================================
    ''' H (S= 0.85) (CT_80-39) '''
    #=============================================
    eps_max_1 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\High\NPY\CT80-39_6322_Zykl_avg_WA_1_WA_2_WA_3_max.npy')    
    N_1 = len(eps_max_1)
    N_max_1 = np.arange(1, N_1 +1, 1)
    
    
    #=============================================
    ''' H (S= 0.85) (CT_80-40) '''
    #=============================================
    eps_max_2 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\High\NPY\CT80-40_13713_Zykl_avg_WA_1_WA_2_WA_3_max.npy')    
    N_2 = len(eps_max_2)
    N_max_2 = np.arange(1, N_2 +1, 1)   
    
    #=============================================
    ''' H (S= 0.85) (CT_80-42) '''
    #=============================================
    eps_max_3 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\High\NPY\CT80-42_3610_Zykl_avg_WA_1_WA_2_WA_3_min.npy')    
    N_3 = len(eps_max_3)
    N_max_3 = np.arange(1, N_3 +1, 1) 
    
    
    #=============================================
    ''' H (S= 0.85) (CT_80-43) '''
    #=============================================
    eps_max_4 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\High\NPY\CT80-43_10949_Zykl_avg_WA_1_WA_2_WA_3_max.npy')    
    N_4 = len(eps_max_4)
    N_max_4 = np.arange(1, N_4 +1, 1) 
    
    #=============================================
    ''' H (S= 0.85) (CT_80-48) '''
    #=============================================
    eps_max_5 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\High\NPY\CT80-48_818_Zykl_avg_WA_1_WA_2_WA_3_max.npy')    
    N_5 = len(eps_max_5)
    N_max_5 = np.arange(1, N_5 +1, 1)
    
    #=============================================
    ''' H (S= 0.85) (CT_80-52) '''
    #=============================================
    eps_max_6 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\High\NPY\CT80-52_647_Zykl_avg_WA_1_WA_2_WA_3_max.npy')    
    N_6 = len(eps_max_6)
    N_max_6 = np.arange(1, N_6 +1, 1)
    
    
    #------------------------------------------------------------------------------ 
    # fatigue creep curves (L)
    #------------------------------------------------------------------------------ 
    #=============================================
    ''' H (S= 0.85) (CT_80-44) '''
    #=============================================
    eps_max_11 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\Low\NPY\CT80-44_74190_Zykl_avg_WA_1_WA_2_WA_3_max.npy')    
    N_11 = len(eps_max_11)
    N_max_11 = np.arange(1, N_11 +1, 1)
    
    #=============================================
    ''' H (S= 0.85) (CT_80-45) '''
    #=============================================
    eps_max_22 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\Low\NPY\CT80-45_32437_Zykl_avg_WA_1_WA_2_WA_3_max.npy')    
    N_22 = len(eps_max_22)
    N_max_22 = np.arange(1, N_22 +1, 1)   
    
    #=============================================
    ''' H (S= 0.85) (CT_80-46) '''
    #=============================================
    eps_max_33 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\Low\NPY\CT80-46_416585_Zykl_avg_WA_1_WA_2_WA_3_max.npy')    
    N_33 = len(eps_max_33)
    N_max_33 = np.arange(1, N_33 +1, 1) 
    
    #=============================================
    ''' H (S= 0.85) (CT_80-47) '''
    #=============================================
    eps_max_44 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\Low\NPY\CT80-47_253449_Zykl_avg_WA_1_WA_2_WA_3_max.npy')    
    N_44 = len(eps_max_44)
    N_max_44 = np.arange(1, N_44 +1, 1) 
    
    #=============================================
    ''' H (S= 0.85) (CT_80-59) '''
    #=============================================
    eps_max_55 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\Low\NPY\CT80-59_142547_Zykl_avg_WA_1_WA_2_WA_3_min.npy')    
    N_55 = len(eps_max_55)
    N_max_55 = np.arange(1, N_55 +1, 1)
    
    #=============================================
    ''' H (S= 0.85) (CT_80-60) '''
    #=============================================
    eps_max_66 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\Low\NPY\CT80-60_538303_Zykl2_avg_WA_1_WA_2_WA_3_max.npy')    
    N_66 = len(eps_max_66)
    N_max_66 = np.arange(1, N_66 +1, 1)
    
    
    #=============================================
    ''' H (S= 0.85) (CT_80-61) '''
    #=============================================
    eps_max_77 = np.load( r'H:\IMB-work\WinConFat\AP 1.5\Rohdaten\C80 Charge 2\Low\NPY\CT80-61_2012546_Zykl_avg_WA_1_WA_2_WA_3_max.npy')    
    N_77 = len(eps_max_77)
    N_max_77= np.arange(1, N_77 +1, 1)
    
    


    plt.subplot(223)
    n = (n1 + n2)
    eps_1_max_1 = np.zeros(n)
    eps_1_min_1 = np.zeros(n)
    cycle = np.zeros(n)
    for i in range(0, n, 1):

        idx_1 = m + 2 * i * m
        idx_2 = 2 * i * m
        if idx_1 <= len(eps_1_arr_1[0:inc_1]):
            idx_1 = idx_1
        else:
            idx_1 = m + 2 * (i - 1.0) * m
            break

        if idx_2 <= len(eps_1_arr_1[0:inc_1]):
            idx_2 = idx_2
        else:
            idx_2 = 1 * (i - 1.0) * m
            break

        eps_1_max_1[i] = eps_1_arr_1[int(idx_1)]
        eps_1_min_1[i] = eps_1_arr_1[int(idx_2)]
        cycle[i] = i + 1

    plt.plot(cycle[0:i]/Nf_1, abs(eps_1_max_1[0:i]), 'r', linewidth=3, alpha=1)

    
         
    #------------------------------------------------------------------------------ 
    plt.subplot(223)
    plt.plot(N_max_1[1:]/N_1, abs((eps_max_1[1:])/300), "k")
    plt.plot(N_max_2[1:]/N_2, abs((eps_max_2[1:])/300), "k")
    plt.plot(N_max_3[1:]/N_3, abs((eps_max_3[1:])/300), "k")
    plt.plot(N_max_4[1:]/N_4, abs((eps_max_4[1:])/300), "k")
    plt.plot(N_max_5[1:]/N_5, abs((eps_max_5[1:])/300), "k")
    plt.plot(N_max_6[1:]/N_6, abs((eps_max_6[1:])/300), "k")

    plt.ylim(0.0015, 0.0045)
    plt.title('Fatigue creep curve normalized (H-L)')
    plt.xlabel('N/Nf')
    plt.ylabel('Displacement [mm]')
    
    
    
    plt.subplot(224)
    n = (n1 + n2)
    eps_1_max_2 = np.zeros(n)
    eps_1_min_2 = np.zeros(n)
    cycle = np.zeros(n)
    for i in range(0, n, 1):

        idx_1 = m + 2 * i * m
        idx_2 = 2 * i * m
        if idx_1 <= len(eps_1_arr_2[0:inc_2]):
            idx_1 = idx_1
        else:
            idx_1 = m + 2 * (i - 1.0) * m
            break

        if idx_2 <= len(eps_1_arr_2[0:inc_2]):
            idx_2 = idx_2
        else:
            idx_2 = 1 * (i - 1.0) * m
            break

        eps_1_max_2[i] = eps_1_arr_2[int(idx_1)]
        eps_1_min_2[i] = eps_1_arr_2[int(idx_2)]
        cycle[i] = i + 1

    plt.plot(cycle[0:i]/Nf_2, abs(eps_1_max_2[0:i]), 'r', linewidth=3, alpha=1)
    
    #------------------------------------------------------------------------------ 
    plt.subplot(224)
    plt.plot(N_max_11[1:]/N_11, abs(eps_max_11[1:]/300), "k")
    plt.plot(N_max_22[1:]/N_22, abs(eps_max_22[1:]/300), "k")
    plt.plot(N_max_33[1:]/N_33, abs(eps_max_33[1:]/300), "k")
    plt.plot(N_max_44[1:]/N_44, abs(eps_max_44[1:]/300), "k")
    plt.plot(N_max_55[1:]/N_55, abs(eps_max_55[1:]/300), "k")
    plt.plot(N_max_66[1:]/N_66, abs(eps_max_66[1:]/300), "k")
    plt.plot(N_max_77[1:]/N_77, abs(eps_max_77[1:]/300), "k")
    
    plt.ylim(0.0015, 0.0045)
    plt.title('Fatigue creep curve normalized (H-L)')
    plt.xlabel('N/Nf')
    plt.ylabel('Displacement [mm]')
    
    
    
        
    plt.subplot(222)
    
    #Exp-IMB
    # Points
    n_1 = np.array([6322, 13713, 3610, 10949, 818, 647, ])
    s_1 = np.array([0.85,  0.85, 0.85, 0.85, 0.85, 0.85 ])
    plt.plot(np.log10(n_1), s_1, 'ro', markersize=4, color='r')
    
    n_1 = np.array([74190, 416585, 253449, 142547, 538303, 2012546 ])
    s_1 = np.array([0.75,  0.75, 0.75, 0.75, 0.75, 0.75 ])
    plt.plot(np.log10(n_1), s_1, 'ro', markersize=4, color='r')
    
    # average
    n_1 = np.array([6010 , 572937])
    s_1 = np.array([0.85,  0.75])
    plt.plot(np.log10(n_1), s_1,  color='r',)
    
    
    # FIB model code 2010
    n_11 = np.array([32064472260,    2853787864,    22605643,    2011937,    179066,    15937,    1418,    126,    11,    1,])
    s_11 = np.array([0.5,    0.55,    0.65,    0.7,    0.75,    0.8,    0.85,    0.9,    0.95,    1])
    plt.plot(np.log10(n_11), s_11, color='k')
    
    # simulation
    n_sim = np.array([Nf_1, Nf_2])
    s_sim = np.array([0.85, 0.75])
    plt.plot(np.log10(n_sim), s_sim, color='g')
    
    plt.xlim(1, 7)
    plt.ylim(0.65, 0.9)




    plt.subplot(221)
    plt.plot(abs(eps_1_arr_1[0:inc_1]), abs(sigma_arr_1[0:inc_1]), 'k', linewidth=1, alpha=1.0)
    plt.title('loading history')

    plt.xlabel('Time')
    plt.ylabel('$\sigma_{1}$')




    plt.show()
