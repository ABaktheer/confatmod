'''
Created on 14.11.2016

@author: abaktheer

two dimensional Drucker-Prager
'''


import matplotlib.pyplot as plt
import numpy as np


def get_bond_slip(eps_N_arr, eps_T_arr, sig_0, E_T, E_N, m):

    # arrays to store the values
    sig_N_arr = np.zeros_like(eps_N_arr)
    sig_T_arr = np.zeros_like(eps_T_arr)
    sig_N_trial_arr = np.zeros_like(eps_N_arr)
    sig_T_trial_arr = np.zeros_like(eps_T_arr)
    eps_N_p_arr = np.zeros_like(eps_N_arr)
    eps_T_p_arr = np.zeros_like(eps_T_arr)


    # state variables
    eps_T_p_i = 0.0
    eps_N_p_i = 0.0
    delta_lamda = 0.0
    
    


    for i in range(1, len(eps_N_arr)):
        
        eps_N_i = eps_N_arr[i]
        eps_T_i = eps_T_arr[i]
        
        sig_N_i_trial = E_N * (eps_N_i - eps_N_p_i)
        sig_T_i_trial  = E_T* (eps_T_i - eps_T_p_i)
        
        # Threshold
        f_pi_i = np.fabs(sig_T_i_trial) - ( sig_0  - m * sig_N_i_trial)
        

        if f_pi_i > 1e-6:
            # Return mapping
            delta_lamda = f_pi_i / (E_T + E_N * m*m)
            
            if delta_lamda <= (np.fabs(sig_T_i_trial) / E_T) :
            # return to the smooth portion
        
                eps_T_p_i += delta_lamda * \
                np.sign(sig_T_i_trial)
                
                eps_N_p_i += delta_lamda * m 
                
                sig_N_i = E_N * (eps_N_i - eps_N_p_i)
            
                sig_T_i = E_T * (eps_T_i - eps_T_p_i) 

            else:
                delta_lamda = (m * sig_N_i_trial - sig_0) / (E_N * m**2)
                # return to  the apex
                
                eps_N_p_i += delta_lamda * m
                
                sig_T_i = 0
                sig_N_i = E_N * (eps_N_i - eps_N_p_i)
                
                  
        else:
        # elstic    
            sig_N_i = sig_N_i_trial 
            sig_T_i = sig_T_i_trial 
            

        sig_N_trial_arr[i] = sig_N_i_trial
        sig_T_trial_arr[i] = sig_T_i_trial
        sig_N_arr[i] = sig_N_i
        sig_T_arr[i] = sig_T_i

        eps_T_p_arr[i] = eps_T_p_i
        eps_N_p_arr[i] = eps_N_p_i

    return  sig_N_trial_arr, sig_T_trial_arr, sig_N_arr, sig_T_arr, eps_T_p_arr, eps_N_p_arr


if __name__ == '__main__':
    
    inc = 20

    s_history = np.array([0, 5.0])
    eps_N_arr = np.hstack([np.linspace(s_history[i], s_history[i + 1], inc)
                         for i in range(len(s_history) - 1)])
    
    t_N_arr= t_arr = np.linspace(0, 1, len(eps_N_arr))
    
    

    s_history = np.array([0, 5.0])
    eps_T_arr = np.hstack([np.linspace(s_history[i], s_history[i + 1], inc)
                         for i in range(len(s_history) - 1)])
    
    t_T_arr= t_arr = np.linspace(0, 1, len(eps_T_arr))


    
    sig_0=10.0
    E_N=50.
    E_T=50.
    m=0.1


    sig_N_trial_arr, sig_T_trial_arr, sig_N_arr, sig_T_arr, eps_T_p_arr, eps_N_p_arr = get_bond_slip(
        eps_N_arr, eps_T_arr, sig_0=sig_0, E_T=E_T, E_N=E_N, m=m)

    
    
    
    ax1 = plt.subplot(231)
    ax1.plot(t_N_arr, eps_N_arr, 'k', linewidth=2,
             label='normal')
    ax1.plot(t_T_arr, eps_T_arr, 'r', linewidth=2,
             label='tangential')
    ax1.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax1.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('loading history')
    plt.xlabel('time')
    plt.ylabel('strain')
    plt.legend(loc=4)
    
    
    

    ax1 = plt.subplot(232)
    ax1.plot(eps_N_arr, sig_N_arr, 'k', linewidth=2,
             label='normal')
    ax1.plot(eps_N_arr, sig_N_trial_arr, '--k', linewidth=2,
             label='normal')
    ax1.plot(eps_T_arr, sig_T_arr, 'r', linewidth=2,
             label='tangential')
    ax1.plot(eps_T_arr, sig_T_trial_arr, '--r', linewidth=2,
             label='tangential')

    
    ax1.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax1.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('stress-strain')
    plt.xlabel('strain')
    plt.ylabel('Stress(MPa)')
    plt.legend(loc=4)
    
    

    ax2 = plt.subplot(233)
    ax2.plot(eps_N_arr, eps_N_p_arr, 'k', linewidth=2,
             label='normal')
    ax2.plot(eps_T_arr, eps_T_p_arr, 'r', linewidth=2,
             label='tangential')

    ax2.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax2.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('plastic strain')
    plt.ylim(0, 1)
    plt.xlabel('strain')
    plt.ylabel('plastic strain')
    plt.legend(loc=4)
    
    
    
 

    sig_N = np.linspace(-50 , sig_0 / m, 1000) 
    
    f_old = (sig_0 - m * sig_N) 
    
    #f_t = (1.0 - np.heaviside((sig_N), 1) * ((sig_N)**2 / (sig_t2)**2  ))
    
    #f = (sig_0 - m * sig_N)**1  * f_t  #* f_c
    

     

     
    
    plt.subplot(236) 
    plt.plot(sig_N , f_old, 'k', linewidth=2.0, alpha=1.0)
    #plt.plot(sig_N , f, 'r', linewidth=2.0, alpha=1.0)
    
    
    plt.plot(sig_N_trial_arr[:], np.abs(sig_T_trial_arr)[:], 'go', alpha=1.0)
    plt.plot(sig_N_arr[:], np.abs(sig_T_arr)[:], 'yo', linewidth=2.0, alpha=1.0)
    #plt.plot(sig_N_trial_arr, sig_N_arr, 'g', alpha=1.0)   
    #plt.plot(np.abs(sig_T_trial_arr), np.abs(sig_T_arr), 'y', alpha=1.0)  
    
    #plt.plot(sig_N , 20*np.heaviside((sig_N - sig_t1), 1), '--r', linewidth=2.0, alpha=1.0)
    #plt.plot(sig_N , 20*np.heaviside((sig_c - sig_N),1), '--b', linewidth=2.0, alpha=1.0)
    #plt.ylim(0,50)
    #plt.plot(sig_N , -f, 'k', linewidth=2, alpha=1.0)
    
    a = np.array([np.abs(sig_T_trial_arr)[2] , sig_N_arr[2]])
    
    for i in np.arange(0, inc):
        
        plt.arrow(sig_N_trial_arr[i], np.abs(sig_T_trial_arr)[i] ,
                sig_N_arr[i] -  sig_N_trial_arr[i], 
                np.abs(sig_T_arr)[i]- np.abs(sig_T_trial_arr)[i], 
                head_width=0.2,head_length=0.5,color = 'r')
    
        
    #plt.xlim([-200,200]) 
    plt.ylim([0,120])
    
    #plt.axes().arrow(sig_N_trial_arr[2], np.abs(sig_T_trial_arr)[2] , a, a , head_width=0.2,head_length=0.4,color = 'r')
    
    
    plt.axhline(y=0, color='k', linewidth=0.5, alpha=0.5)
    plt.axvline(x=0, color='k', linewidth=0.5, alpha=0.5)
    
    plt.title('Threshold and return mapping')
    plt.xlabel('sigma_N')
    plt.ylabel('sigma_T')



    plt.show()

