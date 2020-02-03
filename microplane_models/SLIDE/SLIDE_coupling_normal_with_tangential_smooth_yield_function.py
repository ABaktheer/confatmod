'''
Created on 14.11.2016

@author: abaktheer
'''

#from scipy.optimize import newton

import matplotlib.pyplot as plt
import numpy as np


def get_bond_slip(eps_N_arr, eps_T_arr, sig_0, sig_t1, sig_t2, Rt, E_T, E_N,  S_N, c_N, S_T, c_T,  m, b):

    # arrays to store the values
    sig_N_arr = np.zeros_like(eps_N_arr)
    sig_T_arr = np.zeros_like(eps_N_arr)
    omega_N_arr = np.zeros_like(eps_N_arr)
    omega_T_arr = np.zeros_like(eps_T_arr)
    eps_T_pi_arr = np.zeros_like(eps_T_arr)
    eps_T_pi_cum_arr = np.zeros_like(eps_T_arr)


    # state variables
    eps_T_pi_i = 0.0
    omega_N_i = 0.0 
    omega_T_i = 0.0
    delta_lamda = 0.0
    eps_T_pi_cum_i = 0.0

    for i in range(1, len(eps_N_arr)):
        
        eps_N_i = eps_N_arr[i]
        eps_T_i = eps_T_arr[i]
        
        sig_N_i = (1.0 - omega_N_i) * E_N * eps_N_i
        sig_T_i = (1.0 - omega_T_i) * E_N * (eps_T_i - eps_T_pi_i)

        sig_N_i_eff = E_N * eps_N_i
        sig_T_i_eff = E_N * (eps_T_i - eps_T_pi_i)
        
        Y_N_i = 0.5 * E_N * eps_N_i ** 2
        Y_T_i = 0.5 * E_T * (eps_T_i -  eps_T_pi_i)**2

        # Threshold
            
        f_pi_i = np.fabs(sig_T_i_eff) - ( sig_0  - m * sig_N_i_eff)**1 *\
        (1.0 - np.heaviside((sig_N_i_eff - sig_t1), 1) * ((sig_N_i_eff - sig_t1)**2 / (sig_t2 - sig_t1)**2  ))
        #(1.0 - np.heaviside((sig_N_i_eff), 1) * ((sig_N_i_eff)**2 / (sig_t)**2  ))
        
        print('f_pi_i ', f_pi_i )
         

        if f_pi_i > 1e-6:
            # Return mapping
            delta_lamda = f_pi_i / (E_T / (1. - omega_T_i))
            
            #delta_lamda = (E_T * (eps_T_i - eps_T_arr[i-1]) * np.sign(sig_T_i_eff ) + m*E_N* (eps_N_i - eps_N_arr[i-1])) / (E_T / (1. - omega_T_i))
#             print(delta_lamda)
            # update all the state variables

            eps_T_pi_i += delta_lamda * \
                np.sign(sig_T_i_eff ) / (1 - omega_T_i)

            #Y_N_i = 0.5 * E_N * eps_N_i ** 2
            Y_T_i = 0.5 * E_T * (eps_T_i -  eps_T_pi_i)**2


            omega_T_i  += ((1 - omega_T_i) ** c_T) * delta_lamda * (b * Y_N_i/S_N + Y_T_i/S_T) 

            omega_N_i  += ((1 - omega_N_i) ** c_N) * delta_lamda * (Y_N_i/S_N + b* Y_T_i/S_T) * np.heaviside(sig_N_i_eff,1)

            sig_N_i = (1.0 - omega_N_i) * E_N * eps_N_i
            
            sig_T_i = (1.0 - omega_T_i) * E_N * (eps_T_i - eps_T_pi_i)

            eps_T_pi_cum_i +=  delta_lamda / (1 - omega_T_i)


        sig_N_arr[i] = sig_N_i
        sig_T_arr[i] = sig_T_i
        
        omega_N_arr[i] = omega_N_i
        omega_T_arr[i] = omega_T_i
        
        eps_T_pi_arr[i] = eps_T_pi_i

        eps_T_pi_cum_arr[i] = eps_T_pi_cum_i


    return  sig_N_arr, sig_T_arr, omega_N_arr, omega_T_arr, eps_T_pi_arr, eps_T_pi_cum_arr


if __name__ == '__main__':

    s_history = np.array([0, 0.005])
    eps_N_arr = np.hstack([np.linspace(s_history[i], s_history[i + 1], 1000)
                         for i in range(len(s_history) - 1)])
    
    t_N_arr= t_arr = np.linspace(0, 1, len(eps_N_arr))
    
    

    s_history = np.array([0, 0.01])
    eps_T_arr = np.hstack([np.linspace(s_history[i], s_history[i + 1], 1000)
                         for i in range(len(s_history) - 1)])
    
    t_T_arr= t_arr = np.linspace(0, 1, len(eps_T_arr))


    
    sig_0=2
    E_N=10000.
    E_T=5000.
    S_N=0.00001
    S_T = 0.005
    c_N= 1.8
    c_T =2.5
    K = 0.0
    gamma= 0.0
    
    sig_t1 = -0.0
    sig_t2 = 5.0
    
    Rt=1
    
    m=0.1
    b=0.2


    sig_N_arr, sig_T_arr, omega_N_arr, omega_T_arr, eps_T_pi_arr, eps_T_pi_cum_arr = get_bond_slip(
        eps_N_arr, eps_T_arr, sig_0=sig_0, sig_t1=sig_t1, sig_t2=sig_t2, Rt=Rt, E_T=E_T, E_N=E_N,  S_N=S_N, c_N=c_N, S_T=S_T, c_T=c_T,  m=m, b=b)

    
    
    ax1 = plt.subplot(231)
    ax1.plot(t_N_arr, eps_N_arr, 'k', linewidth=2,
             label='normal')
    ax1.plot(t_T_arr, eps_T_arr, 'r', linewidth=2,
             label='tangential')

    
    ax1.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax1.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('Bond_slip')
    plt.xlabel('Slip(mm)')
    plt.ylabel('Stress(MPa)')
    plt.legend(loc=4)
    
    
    

    ax1 = plt.subplot(232)
    ax1.plot(eps_N_arr, sig_N_arr, 'k', linewidth=2,
             label='normal')
    ax1.plot(eps_T_arr, sig_T_arr, 'r', linewidth=2,
             label='tangential')

    
    ax1.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax1.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('Bond_slip')
    plt.xlabel('Slip(mm)')
    plt.ylabel('Stress(MPa)')
    plt.legend(loc=4)

    ax2 = plt.subplot(233)

    ax2.plot(eps_N_arr, omega_N_arr, 'k', linewidth=2,
             label='normal')
    ax2.plot(eps_T_arr, omega_T_arr, 'r', linewidth=2,
             label='tangential')

    ax2.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax2.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('Damage evolution')
    plt.ylim(0, 1)
    plt.xlabel('Slip(mm)')
    plt.ylabel('Damage')
    plt.legend(loc=4)

    #'''
    plt.subplot(234)
 
    plt.plot(eps_T_arr, eps_T_pi_cum_arr, 'r', linewidth=2,
             label='tangential')
 
 
    plt.xlabel('eps_T_arr(mm)')
    plt.ylabel('cum_eps_T')

    plt.legend()
    
    
    
    ax1 = plt.subplot(235)
    ax1.plot(t_N_arr, omega_N_arr, 'k', linewidth=2,
             label='normal')
    ax1.plot(t_T_arr, omega_T_arr, 'r', linewidth=2,
             label='tangential')

    
    ax1.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax1.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('Bond_slip')
    plt.xlabel('time')
    plt.ylabel('damage')
    plt.legend(loc=4)

    plt.show()
