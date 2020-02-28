

import matplotlib.pyplot as plt
import numpy as np



def get_stress_strain_(sigma_1_arr, lamda, mu, alpha, beta, g, C0, C1, K, n):

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
        d_D = Y_norm   * d_w
        D_i += d_D


        # Helmholtz free energy
        phi_i = 0.5 * lamda * (eps_1_i + 2.0 * eps_2_i)**2.0 + mu * ((eps_1_i)**2.0 + 2.0 * eps_2_i**2.0) + 2.0 * g * w_i * eps_2_i + alpha * \
            (2.0 * w_i * eps_1_i * eps_2_i + 4.0 * w_i *
             eps_2_i**2.0) + 4.0 * beta * w_i * eps_2_i**2.0

        
        if d_w  > 0.5:
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

    sigma_u = - 100
    
    
    N_levels = 10
    s_max_last = 1.0
    s_max_first = 0.55
    N_cycles = 10
    s_min_all = 0.1

    sigma_arr_1 = np.zeros(1)

    for i in range(1, N_levels + 1, 1):
        #print(i)

        d_i = np.linspace(0, sigma_u * (s_max_first + (i - 1.0) * (s_max_last - s_max_first) / (N_levels - 1.0)),  N_cycles * 2 )

        d_i.reshape(-1, 2)[:, 0] = sigma_u * (s_max_first + (i - 1.0) * (s_max_last - s_max_first) / (N_levels - 1.0))
        d_i.reshape(-1, 2)[:, 1] = sigma_u* s_min_all

        d_history_i = d_i.flatten()

        
        sigma_arr_1 = np.hstack((sigma_arr_1, d_history_i))
        
        #print('d_arr_1', d_arr_1)
        
    sigma_arr_1 = np.hstack([np.linspace(sigma_arr_1[i], sigma_arr_1[i + 1], m,  endpoint=False)
                           for i in range(len(sigma_arr_1) - 1)])
        
    
        
        #sigma_arr_1 = d_arr_1
    
    #print(sigma_arr_1)


    t_arr_1 = np.linspace(0, 1, len(sigma_arr_1))
    
    print('sigma_arr_1', len(sigma_arr_1))
    
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
    sigma_arr_1, eps_1_arr_1, eps_2_arr_1, w_arr_1, f_arr_1,  D_arr_1, inc_1, phi_arr_1 = get_stress_strain_(
        sigma_arr_1, lamda, mu, alpha, beta, g, C0, C1, K, n)


    print('sigma_arr_1', len(sigma_arr_1))
    print('eps_1_arr_1', len(eps_1_arr_1[0:inc_1]))

    Nf= inc_1/(m*2)
    
    
    plt.subplot(232)
    plt.plot(-eps_1_arr_1, -sigma_arr_1, 'r', linewidth=3, alpha=1)
    
    
    
    plt.ylim(0.00, 120)
    plt.xlim(0.00, 0.004)
    


    plt.subplot(233)
    n = np.int(Nf)
    print('n', n)
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
        
        print(' idx_1 ',  idx_1 )
        
        print(' idx_2 ',  idx_2 )
        eps_1_max_1[i] = eps_1_arr_1[idx_1]
        eps_1_min_1[i] = eps_1_arr_1[int(idx_2)]
        cycle[i] = i + 1
        
        
    print('eps_1_arr_1', np.max(-eps_1_arr_1))    
    print(' eps_1_max_1 ',  eps_1_max_1)
    
    plt.plot(cycle[0:i], abs(eps_1_max_1[0:i]), 'r', linewidth=3, alpha=1)
    plt.plot(cycle[1:i], abs(eps_1_min_1[1:i]), 'r', linewidth=3, alpha=1)


    
    #plt.ylim(0.00, 0.004)
    #plt.xlim(0.00, 100)
    plt.title('Fatigue creep curve normalized (H-L)')
    plt.xlabel('N')
    plt.ylabel('strian')









    #------------------------------------------------------------------------------ 
    # Experimental data - C80
    #------------------------------------------------------------------------------ 
      
    #------------------------------------------------------------------------------ 
    # Load-displacement & fatigue creep curves (step-wise loading scenario)
    #------------------------------------------------------------------------------ 
    #=============================================
    ''' (CT_80-5) '''
    #=============================================
     
    eps_1 = np.load( r'H:\Publishing\Journal_papers\Journal_paper_loading_sequence_effect_02\results\Experimental_results\step_wise\CT_80-5-86Zyk_avg_WA_1_WA_2_WA_3.npy') 
    sig_1 = np.load( r'H:\Publishing\Journal_papers\Journal_paper_loading_sequence_effect_02\results\Experimental_results\step_wise\CT_80-5-86Zyk_Kraft.npy') 
     
    eps_max_1 = np.load( r'H:\Publishing\Journal_papers\Journal_paper_loading_sequence_effect_02\results\Experimental_results\step_wise\CT_80-5-86Zyk_avg_WA_1_WA_2_WA_3_max.npy') 
    eps_min_1 = np.load( r'H:\Publishing\Journal_papers\Journal_paper_loading_sequence_effect_02\results\Experimental_results\step_wise\CT_80-5-86Zyk_avg_WA_1_WA_2_WA_3_min.npy')   
    N_01 = len(eps_max_1)
    N_1 = np.arange(1, N_01 +1, 1)
     
     
     
    #=============================================
    ''' (CT_80-6) '''
    #=============================================
    eps_2 = np.load( r'H:\Publishing\Journal_papers\Journal_paper_loading_sequence_effect_02\results\Experimental_results\step_wise\CT_80-6-86Zyk_avg_WA_1_WA_2_WA_3.npy') 
    sig_2 = np.load( r'H:\Publishing\Journal_papers\Journal_paper_loading_sequence_effect_02\results\Experimental_results\step_wise\CT_80-6-86Zyk_Kraft.npy') 
     
    eps_max_2 = np.load( r'H:\Publishing\Journal_papers\Journal_paper_loading_sequence_effect_02\results\Experimental_results\step_wise\CT_80-6-86Zyk_avg_WA_1_WA_2_WA_3_max.npy') 
    eps_min_2 = np.load( r'H:\Publishing\Journal_papers\Journal_paper_loading_sequence_effect_02\results\Experimental_results\step_wise\CT_80-6-86Zyk_avg_WA_1_WA_2_WA_3_min.npy')   
    N_02 = len(eps_max_2)
    N_2 = np.arange(1, N_02 +1, 1)
     
     
    #=============================================
    ''' (CT_80-8) '''
    #=============================================
    eps_3 = np.load( r'H:\Publishing\Journal_papers\Journal_paper_loading_sequence_effect_02\results\Experimental_results\step_wise\CT_80-8-91Zyk_avg_WA_1_WA_2_WA_3.npy') 
    sig_3 = np.load( r'H:\Publishing\Journal_papers\Journal_paper_loading_sequence_effect_02\results\Experimental_results\step_wise\CT_80-8-91Zyk_Kraft.npy') 
     
    eps_max_3 = np.load( r'H:\Publishing\Journal_papers\Journal_paper_loading_sequence_effect_02\results\Experimental_results\step_wise\CT_80-8-91Zyk_avg_WA_1_WA_2_WA_3_max.npy') 
    eps_min_3 = np.load( r'H:\Publishing\Journal_papers\Journal_paper_loading_sequence_effect_02\results\Experimental_results\step_wise\CT_80-8-91Zyk_avg_WA_1_WA_2_WA_3_min.npy')   
    N_03 = len(eps_max_3)
    N_3 = np.arange(1, N_03 +1, 1)
     
     
     
     
    # plot stress-strain
    plt.subplot(231)
    plt.plot(eps_1/300,-1000* sig_1/(np.pi * 100**2 /4), "k")
    plt.plot(eps_2/300,-1000* sig_2/(np.pi * 100**2 /4), "r")
    plt.plot(eps_3/300,-1000* sig_3/(np.pi * 100**2 /4), "g")
     
    plt.ylim(0.00, 120)
    plt.xlim(0.00, 0.004)
    plt.title('Fatigue creep curve normalized (H-L)')
    plt.xlabel('strain')
    plt.ylabel('Stress')
     
     
     
     
    # plot ceep-fatigue
    plt.subplot(233)
    plt.plot(N_1[1:], abs((eps_max_1[1:])/300), "k")
    plt.plot(N_1[1:], abs((eps_min_1[1:])/300), "k")
     
    plt.plot(N_2[1:], abs((eps_max_2[1:])/300), "k")
    plt.plot(N_2[1:], abs((eps_min_2[1:])/300), "k")
     
    plt.plot(N_3[1:], abs((eps_max_3[1:])/300), "k")
    plt.plot(N_3[1:], abs((eps_min_3[1:])/300), "k")
     
     
    plt.ylim(0.00, 0.004)
    plt.title('Fatigue creep curve normalized (H-L)')
    plt.xlabel('N')
    plt.ylabel('strian')
    
    
    plt.show()