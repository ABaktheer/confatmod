'''
Created on 14.11.2016

@author: abaktheer
'''

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt


def get_bond_slip(s_arr, tau_pi_bar, K, gamma, E_b, S, c, r, m, sigma_n):
    '''for plotting the bond slip fatigue - Initial version modified modified threshold with cumulation-2 implicit
    '''

    #=====================================================
    # arrays to store the values
    #=====================================================

    # nominal stress
    tau_arr = np.zeros_like(s_arr)
    # sliding stress
    tau_pi_arr = np.zeros_like(s_arr)
    # damage factor
    w_arr = np.zeros_like(s_arr)
    # sliding slip
    s_pi_arr = np.zeros_like(s_arr)
    # max sliding
    s_max = np.zeros_like(s_arr)
    # max stress
    tau_max = np.zeros_like(s_arr)
    # cumulative sliding
    s_pi_cum = np.zeros_like(s_arr)
    diss = np.zeros_like(s_arr)
    k_ = np.zeros_like(s_arr)
    f_ = np.zeros_like(s_arr)

    #=====================================================
    # material parameters
    #=====================================================
    # shear modulus [MPa]
    E_b = E_b
    # damage - brittleness [MPa^-1]
    K = K
    # Kinematic hardening modulus [MPa]
    gamma = gamma
    # constant in the sliding threshold function
    tau_pi_bar = tau_pi_bar
    # material parameters
    S = S
    c = c
    r = r
    m = m
    sigma_n = sigma_n

    #================================================
    # state variables
    #================================================
    tau_i = 0
    alpha_i = 0.
    s_pi_i = 0
    z_i = 0.
    w_i = 0.
    w_n = 0.0
    delta_lamda = 0.
    s_pi_cum_i = 0.
    diss_i = 0.
    k_i = 0.
    f_i = 0

    for i in range(1, len(s_arr)):
        print('increment', i)
        s_i = s_arr[i]
        s_max_i = np.fabs(s_i)

        tau_i = (1 - w_i) * E_b * (s_i - s_pi_i)

        tau_i_1 = E_b * (s_i - s_pi_i)

        Y_i = 0.5 * E_b * (s_i - s_pi_i) ** 2

        # Threshold
        f_pi_i = np.fabs(tau_i_1 - gamma * alpha_i) - \
            tau_pi_bar - K * z_i + m * sigma_n

        if f_pi_i > 1e-10:
            # Return mapping

            # update all the state variables

            def f_w_n(w_n):

                a1 = (tau_pi_bar / (tau_pi_bar - m * sigma_n))
                a2 = np.sign(tau_i_1 - gamma * alpha_i)
                a3 = (f_pi_i / (E_b / (1.0 - w_n) + K + gamma))

                return w_n - w_i - a3 * ((1.0 - w_n)**c) * a1 * ((0.5 * E_b *
                (s_i - s_pi_i - a3 * a2 / (1.0 - w_n))**2.0) / S)**r
                
                
            def f_w_n_2(w_n):

                a1 = (tau_pi_bar / (tau_pi_bar - m * sigma_n))
                a2 = np.sign(tau_i_1 - gamma * alpha_i)
                a3 = (f_pi_i / (E_b / (1.0 - w_n) + K + gamma))

                return   1.0/((-a2*a3 + (1.0 - w_n)*(-s_pi_i + s_i))/
                              (1.0 - w_n))*(-2.0*a1*a2*pow(a3, 2)*r*\
                                            pow(0.5*E_b*pow((-a2*a3 + (1.0 - w_n)*(-s_pi_i + s_i))/\
                                            (1.0 - w_n), 2.0)/S, r)*pow(1.0 - w_n, c + 1) -\
                                             a1*a3*c*pow((-a2*a3 + (1.0 - w_n)*(-s_pi_i + s_i))/\
                                                         (1.0 - w_n), 1.0)*pow(0.5*E_b*pow((-a2*a3 + (1.0 - w_n)*(-s_pi_i + s_i))/\
                                                                                           (1.0 - w_n), 2.0)/S, r)*pow(1.0 - w_n, c + 2) +\
                                                          pow((-a2*a3 + (1.0 - w_n)*(-s_pi_i + s_i))/(1.0 - w_n), 1.0)*pow(w_n - 1.0, 3))/pow(w_n - 1.0, 3)

                    

            #w_i = opt.newton(f_w_n,  0.0, fprime= f_w_n_2,  tol=1.0e-8, maxiter=20, rtol=0.0)
            
            k=0
            while  k < 500:
                if abs(f_w_n(w_n)) < 1e-10:
                    w_i = w_n
                    break 
                else:  
                    w_n = w_n - f_w_n(w_n)/f_w_n_2(w_n)
                    k += 1
            else:
                print('no convergence')
                break

            
            k_i = k
            

            delta_lamda = f_pi_i / (E_b / (1.0 - w_i) + gamma + K)

            s_pi_i +=  delta_lamda * \
                np.sign(tau_i_1 - gamma * alpha_i) / (1 - w_i)

            Y_i = 0.5 * E_b * (s_i - s_pi_i) ** 2.0


            tau_i = E_b * (1.0 - w_i) * (s_i - s_pi_i)

            alpha_i += delta_lamda * \
                np.sign(tau_i_1 - gamma * alpha_i)
                
            z_i += delta_lamda
            
            s_pi_cum_i = s_pi_cum_i + delta_lamda / (1.0 - w_i)
            
            f_i = np.fabs(tau_i/(1.0 - w_i) - gamma * alpha_i) - \
            tau_pi_bar - K * z_i + m * sigma_n
            
            
#             diss_i += (tau_pi_bar - m * sigma_n / 3.0) * delta_lamda + \
#                 Y_i * ((1.0 - w_i) ** c) * (delta_lamda * (Y_i / S) ** r)
                
                
            diss_i += tau_i* (s_pi_i - s_pi_arr[i-1]) - gamma*alpha_i *(delta_lamda * \
                np.sign(tau_i_1 - gamma * alpha_i)) - K * z_i * delta_lamda + Y_i * (w_i - w_arr[i-1])    
            
            print('f=', f_i)

        tau_max_i = np.fabs(tau_i)
        tau_arr[i] = tau_i
        w_arr[i] = w_i
        s_pi_arr[i] = s_pi_i
        s_max[i] = s_max_i
        tau_max[i] = tau_max_i
        s_pi_cum[i] = s_pi_cum_i
        diss[i] = diss_i
        k_[i] = k_i
        f_[i] = f_i

    return s_arr, tau_arr, w_arr, s_pi_arr, s_max, tau_max, s_pi_cum, diss, k_, f_


if __name__ == '__main__':
    s_levels = np.linspace(0, 8, 2)

    s_levels[0] = 0
    s_levels.reshape(-1, 2)[:, 0] *= -1
    # s_levels.reshape(-1, 2)[:, 1] = 2
    s_history = s_levels.flatten()

    # slip array as input
    s_arr_1 = np.hstack([np.linspace(s_history[i], s_history[i + 1], 10)
                         for i in range(len(s_levels) - 1)])

    s_arr_2 = np.hstack([np.linspace(s_history[i], s_history[i + 1], 20)
                         for i in range(len(s_history) - 1)])
 
    s_arr_3 = np.hstack([np.linspace(s_history[i], s_history[i + 1], 30)
                         for i in range(len(s_history) - 1)])
 
    s_arr_4 = np.hstack([np.linspace(s_history[i], s_history[i + 1], 40)
                         for i in range(len(s_history) - 1)])
     
    s_arr_5 = np.hstack([np.linspace(s_history[i], s_history[i + 1], 50)
                         for i in range(len(s_history) - 1)])
    
    
    tau_pi_bar=1
    K=0.1
    gamma=0.2
    E_b=1
    S=0.001
    c=1
    r=0.001
    m=0
    sigma_n=0


    s_arr_1, tau_arr_1, w_arr_1, s_pi_arr_1, s_max_1, tau_max_1, s_pi_cum_1, diss_1, k_1, f_1= get_bond_slip(
        s_arr_1, tau_pi_bar=tau_pi_bar, K=K, gamma=gamma, E_b=E_b, S=S, c=c, r=r, m=m, sigma_n=sigma_n)

    s_arr_2, tau_arr_2, w_arr_2, s_pi_arr_2, s_max_2, tau_max_2, s_pi_cum_2, diss_2, k_2, f_2 = get_bond_slip(
        s_arr_2, tau_pi_bar=tau_pi_bar, K=K, gamma=gamma, E_b=E_b, S=S, c=c, r=r, m=m, sigma_n=sigma_n)
 
    s_arr_3, tau_arr_3, w_arr_3, s_pi_arr_3, s_max_3, tau_max_3, s_pi_cum_3, diss_3, k_3, f_3 = get_bond_slip(
        s_arr_3, tau_pi_bar=tau_pi_bar, K=K, gamma=gamma, E_b=E_b, S=S, c=c, r=r, m=m, sigma_n=sigma_n)
 
    s_arr_4, tau_arr_4, w_arr_4, s_pi_arr_4, s_max_4, tau_max_4, s_pi_cum_4, diss_4, k_4, f_4 = get_bond_slip(
        s_arr_4, tau_pi_bar=tau_pi_bar, K=K, gamma=gamma, E_b=E_b, S=S, c=c, r=r, m=m, sigma_n=sigma_n)
     
    s_arr_5, tau_arr_5, w_arr_5, s_pi_arr_5, s_max_5, tau_max_5, s_pi_cum_5, diss_5, k_5, f_5 = get_bond_slip(
        s_arr_5, tau_pi_bar=tau_pi_bar, K=K, gamma=gamma, E_b=E_b, S=S, c=c, r=r, m=m, sigma_n=sigma_n)
     
    
        
#     s_arr_5, tau_arr_5, w_arr_5, s_pi_arr_5, s_max_5, tau_max_5, s_pi_cum_5, diss_5 = get_bond_slip(
#         s_arr_5, tau_pi_bar=5, K=0, gamma=100, E_b=1000, S=.05, c=1, r=1.0, m=1.7, sigma_n=0)
    
    

    ax1 = plt.subplot(231)

    ax1.plot(s_arr_1, tau_arr_1, 'b', linewidth=2,
             label=' n_steps = 10 ')
    ax1.plot(s_arr_2, tau_arr_2, 'r', linewidth=2,
             label=' n_steps = 20 ')
    ax1.plot(s_arr_3, tau_arr_3, 'g', linewidth=2,
             label=' n_steps = 30 ')
    ax1.plot(s_arr_4, tau_arr_4, 'k', linewidth=2,
             label=' n_steps = 40 ')
    ax1.plot(s_arr_5, tau_arr_5, 'y', linewidth=2,
             label=' n_steps = 50 ')

    ax1.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax1.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('Bond_slip')
    plt.xlabel('Slip(mm)')
    plt.ylabel('Stress(MPa)')
    plt.legend(loc=4)

    ax2 = plt.subplot(233)

    ax2.plot(s_arr_1, w_arr_1, 'b', linewidth=2,
             label=' n_steps = 10 ')
    ax2.plot(s_arr_2, w_arr_2, 'r', linewidth=2,
             label=' n_steps = 20 ')
    ax2.plot(s_arr_3, w_arr_3, 'g', linewidth=2,
             label=' n_steps = 30 ')
    ax2.plot(s_arr_4, w_arr_4, 'k', linewidth=2,
             label=' n_steps = 40 ')
    ax2.plot(s_arr_5, w_arr_5, 'y', linewidth=2,
             label=' n_steps = 50 ')
    
    ax2.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax2.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    plt.title('Damage evolution')
    plt.ylim(0, 1)
    plt.xlabel('Slip(mm)')
    plt.ylabel('Damage')
    plt.legend(loc=4)
    
    
    plt.subplot(232)
    plt.plot(s_arr_1, k_1, 'b', linewidth=2,
             label=' n_steps = 10 ')
    plt.plot(s_arr_2, k_2, 'r', linewidth=2,
             label=' n_steps = 20 ')
    plt.plot(s_arr_3, k_3, 'g', linewidth=2,
             label=' n_steps = 30 ')
    plt.plot(s_arr_4, k_4, 'k', linewidth=2,
             label=' n_steps = 40 ')
    plt.plot(s_arr_5, k_5, 'y', linewidth=2,
             label=' n_steps = 50 ')

    plt.xlabel('Slip(mm)')
    plt.ylabel('number of iterations')
    #plt.ylim(0, 1)
    plt.legend()
    
    

    #'''
    plt.subplot(234)
    plt.plot(s_arr_1, s_pi_cum_1, 'b', linewidth=2,
             label=' n_steps = 10 ')
    plt.plot(s_arr_2, s_pi_cum_2, 'r', linewidth=2,
             label=' n_steps = 20 ')
    plt.plot(s_arr_3, s_pi_cum_3, 'g', linewidth=2,
             label=' n_steps = 30 ')
    plt.plot(s_arr_4, s_pi_cum_4,  'k', linewidth=2,
             label=' n_steps = 40 ')
    plt.plot(s_arr_5, s_pi_cum_5, 'y', linewidth=2,
             label=' n_steps = 50 ')

    plt.xlabel('Slip(mm)')
    plt.ylabel('Cumulative sliding(mm)')
    #plt.ylim(0, 1)
    plt.legend()
    
    
    #'''
    plt.subplot(235)
    plt.plot(s_arr_1, diss_1, 'b', linewidth=2,
             label=' n_steps = 10 ')
    plt.plot(s_arr_2, diss_2, 'r', linewidth=2,
             label=' n_steps = 20 ')
    plt.plot(s_arr_3, diss_3, 'g', linewidth=2,
             label=' n_steps = 30 ')
    plt.plot(s_arr_4, diss_4,  'k', linewidth=2,
             label=' n_steps = 40 ')
    plt.plot(s_arr_5, diss_5, 'y', linewidth=2,
             label=' n_steps = 50 ')

    plt.xlabel('Slip(mm)')
    plt.ylabel('Energy dissipation')
    #plt.ylim(0, 1)
    plt.legend()
    
    
    plt.subplot(236)
    plt.plot(s_arr_1, f_1, 'b', linewidth=2,
             label=' n_steps = 10 ')
    plt.plot(s_arr_2, f_2, 'r', linewidth=2,
             label=' n_steps = 20 ')
    plt.plot(s_arr_3, f_3, 'g', linewidth=2,
             label=' n_steps = 30 ')
    plt.plot(s_arr_4, f_4,  'k', linewidth=2,
             label=' n_steps = 40 ')
    plt.plot(s_arr_5, f_5, 'y', linewidth=2,
             label=' n_steps = 50 ')

    plt.xlabel('Slip(mm)')
    plt.ylabel('f (threshold) [MPa]')
    plt.ylim(20e-16, -20e-16)
    plt.legend()




    plt.show()
