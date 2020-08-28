'''
Created on 14.11.2016

@author: abaktheer
'''

import matplotlib.pyplot as plt
import numpy as np


def get_bond_slip(eps_N_arr, eps_T_arr, sig_0, E_T, E_N,  K, gamma, S_N, c_N, S_T, c_T,  m, h):

    # arrays to store the values
    sig_N_arr = np.zeros_like(eps_N_arr)
    sig_T_arr = np.zeros_like(eps_T_arr)
    omega_N_arr = np.zeros_like(eps_N_arr)
    omega_T_arr = np.zeros_like(eps_T_arr)
    eps_T_pi_arr = np.zeros_like(eps_T_arr)
    eps_T_pi_cum_arr = np.zeros_like(eps_T_arr)
    alpha_arr = np.zeros_like(eps_T_arr)
    z_arr = np.zeros_like(eps_T_arr)

    # state variables
    eps_T_pi_i = 0
    omega_T_i = 0.0
    z_i = 0.0
    alpha_i = 0.0
    omega_N_i = 0.0
    delta_lamda = 0
    eps_T_pi_cum_i = 0

    for i in range(1, len(eps_N_arr)):

        eps_N_i = eps_N_arr[i]
        eps_T_i = eps_T_arr[i]
        
        sig_N_i = (1.0 - omega_N_i) * E_N * eps_N_i
        sig_T_i = (1.0 - omega_T_i) * E_T * (eps_T_i - eps_T_pi_i)

        sig_N_i_eff = E_N * eps_N_i
        sig_T_i_eff = E_T * (eps_T_i - eps_T_pi_i)
        
        Y_N_i = 0.5 * E_N * eps_N_i ** 2
        Y_T_i = 0.5 * E_T * (eps_T_i -  eps_T_pi_i)**2

        # Threshold
        f_pi_i = np.fabs(sig_T_i_eff - gamma*alpha_i) - \
            (sig_0 + K * z_i) + m * sig_N_i_eff   
             

        if f_pi_i > 1e-6:
            # Return mapping
            delta_lamda = f_pi_i / (E_T / (1. - omega_T_i))
            delta_lamda = (E_T * (eps_T_i - eps_T_arr[i-1]) * np.sign(sig_T_i_eff - gamma*alpha_i) + m * E_N * (eps_N_i - eps_N_arr[i-1])) / (E_T / (1. - omega_T_i)+ K +gamma)
            
#             delta_lamda = (E_T * (eps_T_i - eps_T_arr[i-1]) * np.sign(sig_T_i_eff - gamma*alpha_i) + m * (1- omega_N_i) *E_N * (eps_N_i - eps_N_arr[i-1])) /\
#              (E_T / (1. - omega_T_i)+ K +gamma+ ((1- omega_N_i)**c_N * ((Y_N_i/S_N + h* Y_T_i/S_T) * np.heaviside(sig_N_i,1) * E_N * eps_N_i)))

# 
#             delta_lamda = (E_T * (eps_T_i - eps_T_arr[i-1]) * np.sign(sig_T_i_eff - gamma*alpha_i) + m * (1- omega_N_i) *E_N * (eps_N_i - eps_N_arr[i-1])) /\
#              (E_T / (1. - omega_T_i)+ K +gamma+ ((1- omega_N_i)**c_N * ((Y_N_i/S_N + h* Y_T_i/S_T) * np.heaviside(sig_N_i,1) * E_N * eps_N_i)))
             

             
            # update all the state variables
            eps_T_pi_i += delta_lamda * \
                np.sign(sig_T_i_eff - gamma*alpha_i ) / (1 - omega_T_i)

            Y_T_i = 0.5 * E_T * (eps_T_i -  eps_T_pi_i)**2

                
            omega_T_i  += ((1 - omega_T_i) ** c_T) * delta_lamda * (h* Y_N_i/S_N + Y_T_i/S_T)

            omega_N_i  += ((1 - omega_N_i) ** c_N) * delta_lamda * (Y_N_i/S_N + h* Y_T_i/S_T) * np.heaviside(sig_N_i_eff,1)

            sig_N_i = (1.0 - omega_N_i) * E_N * eps_N_i
            
            sig_T_i = (1.0 - omega_T_i) * E_T * (eps_T_i - eps_T_pi_i)
            
            z_i += delta_lamda
            alpha_i += delta_lamda * np.sign(sig_T_i_eff - gamma*alpha_i )

            eps_T_pi_cum_i +=  delta_lamda / (1 - omega_T_i)
            
            print( 'delta_lamda', delta_lamda )
            print('range', (1.0 - omega_T_i) * sig_T_i_eff / E_T )
            
            
            f = np.fabs(sig_T_i/(1.0 - omega_T_i) - gamma*alpha_i) - \
            (sig_0 + K * z_i) + m * sig_N_i /(1.0 - omega_N_i) 
            
            
            
            print('f', f)

        sig_N_arr[i] = sig_N_i
        sig_T_arr[i] = sig_T_i
        omega_N_arr[i] = omega_N_i
        omega_T_arr[i] = omega_T_i
        eps_T_pi_arr[i] = eps_T_pi_i
        z_arr[i] = z_i
        alpha_arr[i] = alpha_i
        eps_T_pi_cum_arr[i] = eps_T_pi_cum_i


    return  sig_N_arr, sig_T_arr, omega_N_arr, omega_T_arr, eps_T_pi_arr, eps_T_pi_cum_arr, f


if __name__ == '__main__':
    # normal strain
#     s_levels_1 = np.linspace(0, 1.0, 2)
#     s_levels_1[0] = 0
#     s_levels_1.reshape(-1, 2)[:, 0] *= 0
#     s_levels_1.reshape(-1, 2)[:, 1] = 0
#     
#     s_history_1 = s_levels_1.flatten()
    
    #s_history = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.05])
    #s_history = np.array([0, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])*0.2
    #s_history = np.array([0, 0, 0, 0.125])
    #s_history = np.array([0, 0.0075, 0.0075, 0.0075])
    
    inc = 100

    s_history_1 = np.array([0, 0.01])

    eps_N_arr_1 = np.hstack([np.linspace(s_history_1[i], s_history_1[i + 1], inc)
                         for i in range(len(s_history_1) - 1)])
    
    
    
    s_history_2 = np.array([0, 0.01])

    eps_N_arr_2 = np.hstack([np.linspace(s_history_2[i], s_history_2[i + 1], inc)
                         for i in range(len(s_history_2) - 1)])
    
    
    
    s_history_3 = np.array([0, 0.01])

    eps_N_arr_3 = np.hstack([np.linspace(s_history_3[i], s_history_3[i + 1], inc)
                         for i in range(len(s_history_3) - 1)])
    
    t_N_arr= t_arr = np.linspace(0, 1, len(eps_N_arr_1))
    
    
    # tangential strain
#     s_levels_1 = np.linspace(0, 1.0, 10)
#     s_levels_1[0] = 0
#     s_levels_1.reshape(-1, 2)[:, 0] *= -1
#     #s_levels_1.reshape(-1, 2)[:, 1] = 0
#     s_history_1 = s_levels_1.flatten()

    #s_history = np.array([0,0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.1, 0.1, 0.1])
    #s_history = np.array([0,0.0, 0.0, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1])
    #s_history = np.array([0, 0.125, 0.125, 0.125])
    #s_history = np.array([0, 0, 0, 0.025])
    s_history_1 = np.array([0, 0.0])

    eps_T_arr_1 = np.hstack([np.linspace(s_history_1[i], s_history_1[i + 1], inc)
                         for i in range(len(s_history_1) - 1)])
    
    
    s_history_2 = np.array([0, 0.05])

    eps_T_arr_2 = np.hstack([np.linspace(s_history_2[i], s_history_2[i + 1], inc)
                         for i in range(len(s_history_2) - 1)])
    
    
    s_history_3 = np.array([0, 0.1])

    eps_T_arr_3 = np.hstack([np.linspace(s_history_3[i], s_history_3[i + 1], inc)
                         for i in range(len(s_history_3) - 1)])
    
    
    
    t_T_arr = np.linspace(0, 1, len(eps_T_arr_1))



    sig_0=10
    E_N=1000.
    E_T=1000.
    S_N=0.000005
    S_T = 0.005
    c_N=2.0
    c_T =2.5
    K = 0.0
    gamma= 0.0
    
    m=2.0
    h=0.2


    sig_N_arr_1, sig_T_arr_1, omega_N_arr_1, omega_T_arr_1, eps_T_pi_arr_1, eps_T_pi_cum_arr_1, f1 = get_bond_slip(
        eps_N_arr_1, eps_T_arr_1, sig_0=sig_0, E_T=E_T, E_N=E_N, K=K, gamma=gamma, S_N=S_N, c_N=c_N, S_T=S_T, c_T=c_T,  m=m, h=h)
    
    sig_N_arr_2, sig_T_arr_2, omega_N_arr_2, omega_T_arr_2, eps_T_pi_arr_2, eps_T_pi_cum_arr_2, f2 = get_bond_slip(
        eps_N_arr_2, eps_T_arr_2, sig_0=sig_0, E_T=E_T, E_N=E_N, K=K, gamma=gamma, S_N=S_N, c_N=c_N, S_T=S_T, c_T=c_T,  m=m, h=h)
    
    sig_N_arr_3, sig_T_arr_3, omega_N_arr_3, omega_T_arr_3, eps_T_pi_arr_3, eps_T_pi_cum_arr_3, f3 = get_bond_slip(
        eps_N_arr_3, eps_T_arr_3, sig_0=sig_0, E_T=E_T, E_N=E_N, K=K, gamma=gamma, S_N=S_N, c_N=c_N, S_T=S_T, c_T=c_T,  m=m, h=h)
    

    
    
    ax1 = plt.subplot(231)
    ax1.plot(t_N_arr, eps_N_arr_1, 'g', linewidth=2,
             label='normal')
    ax1.plot(t_T_arr, eps_T_arr_1, 'r', linewidth=2,
             label='tangential')
    
    ax1.plot(t_N_arr, eps_N_arr_2, '--g', linewidth=2,
             label='normal')
    ax1.plot(t_T_arr, eps_T_arr_2, '--r', linewidth=2,
             label='tangential')
    
    ax1.plot(t_N_arr, eps_N_arr_3, '-.g', linewidth=2,
             label='normal')
    ax1.plot(t_T_arr, eps_T_arr_3, '-.r', linewidth=2,
             label='tangential')
    
    
    
    
    plt.title('Loading history')
    plt.xlabel('Time')
    plt.ylabel('Slip(mm)')
    plt.legend(loc=4)
    
    
    

    ax1 = plt.subplot(232)
    ax1.plot(eps_N_arr_1, sig_N_arr_1, 'g', linewidth=2,
             label='normal')
#     ax1.plot(eps_T_arr_1, sig_T_arr_1, 'r', linewidth=2,
#              label='tangential')
    
    ax1.plot(eps_N_arr_2, sig_N_arr_2, '--g', linewidth=2,
             label='normal')
#     ax1.plot(eps_T_arr_1, sig_T_arr_2, '--r', linewidth=2,
#              label='tangential')

    ax1.plot(eps_N_arr_3, sig_N_arr_3, '-.g', linewidth=2,
             label='normal')
#     ax1.plot(eps_T_arr_1, sig_T_arr_3, '-.r', linewidth=2,
#              label='tangential')

    

    plt.title('stress_strain')
    plt.xlabel('Strain')
    plt.ylabel('Stress(MPa)')
    plt.legend(loc=4)


    ax2 = plt.subplot(233)
    ax2.plot(eps_N_arr_1, omega_N_arr_1, 'g', linewidth=2,
             label='normal')
#     ax2.plot(eps_T_arr_1, omega_T_arr_1, 'r', linewidth=2,
#              label='tangential')
    
    ax2.plot(eps_N_arr_2, omega_N_arr_2, '--g', linewidth=2,
             label='normal')
#     ax2.plot(eps_T_arr_1, omega_T_arr_2, '--r', linewidth=2,
#              label='tangential')
    
    ax2.plot(eps_N_arr_3, omega_N_arr_3, '-.g', linewidth=2,
             label='normal')
#     ax2.plot(eps_T_arr_1, omega_T_arr_3, '-.r', linewidth=2,
#              label='tangential')
    
    

    plt.title('Damage evolution')
    plt.ylim(0, 1)
    plt.xlabel('Strain')
    plt.ylabel('Damage')
    plt.legend(loc=4)

    #'''
    plt.subplot(234)

    plt.plot(t_T_arr, eps_T_pi_cum_arr_1, 'r', linewidth=2,
             label='tangential')
 
    plt.title('Cumulative sliding')
    plt.xlabel('Strain')
    plt.ylabel('cum_eps_T')

    plt.legend()
    
    
    
    ax1 = plt.subplot(235)
    ax1.plot(t_N_arr, omega_N_arr_1, 'g', linewidth=2,
             label='normal')
    ax1.plot(t_T_arr, omega_T_arr_1, 'r', linewidth=2,
             label='tangential')
    
    ax1.plot(t_N_arr, omega_N_arr_2, '--g', linewidth=2,
             label='normal')
    ax1.plot(t_T_arr, omega_T_arr_2, '--r', linewidth=2,
             label='tangential')
     
    ax1.plot(t_N_arr, omega_N_arr_3, '-.g', linewidth=2,
             label='normal')
    ax1.plot(t_T_arr, omega_T_arr_3, '-.r', linewidth=2,
             label='tangential')

    plt.title('Damage history')
    plt.xlabel('Time')
    plt.ylabel('Damage')
    plt.legend(loc=4)
    
    
    ax1 = plt.subplot(236)
    ax1.plot(t_N_arr, sig_N_arr_1, 'g', linewidth=2,
             label='normal')
    ax1.plot(t_T_arr, sig_T_arr_1, 'r', linewidth=2,
             label='tangential')

    plt.title('stress history')
    plt.xlabel('Time')
    plt.ylabel('Damage')
    plt.legend(loc=4)

    plt.show()
