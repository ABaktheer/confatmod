'''
Created on 14.11.2016

@author: abaktheer

iso error map of the SLIDE model
Version: two dimensional SLIDE (with Desmorat short-cut return mapping)
'''

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import  root





def get_corr_pred_(eps_N_arr, eps_T_arr, sig_0, sig_t, E_T, E_N, m, S_N, S_T, c_N, c_T, K_T, gamma_T):

    # arrays to store the values
    sig_N_arr= np.zeros_like(eps_N_arr)
    sig_T_arr= np.zeros_like(eps_T_arr)
    sig_N_trial_arr = np.zeros_like(eps_N_arr)
    sig_T_trial_arr = np.zeros_like(eps_T_arr)
    
    eps_N_p_arr = np.zeros_like(eps_N_arr)
    eps_T_p_arr = np.zeros_like(eps_T_arr)
    
    omega_N_arr = np.zeros_like(eps_N_arr)
    omega_T_arr = np.zeros_like(eps_T_arr)
    
    z_T_arr = np.zeros_like(eps_T_arr)
    alpha_T_arr = np.zeros_like(eps_T_arr)
    
    f_arr = np.zeros_like(eps_T_arr)

    # state variables
    eps_N_p_i = 0.0
    omega_N_i = 0.0
    eps_T_p_i = 0.0
    omega_T_i = 0.0
    z_T_i = 0.0 
    alpha_T_i = 0.0
    
    delta_lamda = 0.0
    
    for i in range(1, len(eps_N_arr)):
        
        eps_N_i = eps_N_arr[i]
        eps_T_i = eps_T_arr[i]
        
        sig_N_i_eff_trial = E_N * (eps_N_i - eps_N_p_i)
        sig_T_i_eff_trial  = E_T* (eps_T_i - eps_T_p_i)
        
        
        f_trial = np.fabs(sig_T_i_eff_trial - gamma_T * alpha_T_i)  -\
         (sig_0  + K_T * z_T_i - m * sig_N_i_eff_trial) * (1.0 - np.heaviside((sig_N_i_eff_trial), 1) * ((sig_N_i_eff_trial)**2 / (sig_t)**2 ))
         
        if f_trial > 1e-8: 
            
            def f(vars):
                
                delta_lamda ,   sig_N_eff  = vars

                N =  (m  + np.heaviside((sig_N_eff), 1) * ((-m * sig_N_eff **2)/ (sig_t)**2    +\
                                                                       (sig_0  + K_T * (z_T_i + delta_lamda) - m * sig_N_eff) * (2 * sig_N_eff / (sig_t)**2 ) )) /(1. - omega_N_i)
            
                f1 = sig_0 + K_T * (z_T_i + delta_lamda) - m * sig_N_eff
                
                ft =  1.0 - np.heaviside(sig_N_eff, 1) * (sig_N_eff**2 / (sig_t)**2)
                
#                sign = np.sign(sig_T_i_eff_trial - gamma_T * alpha_T_i)
#                 Y_N = 0.5 * E_N * (eps_N_i - (eps_N_p_i + delta_lamda * N))**2.0
#                 Y_T = 0.5 * E_T * (eps_T_i - (eps_T_p_i + delta_lamda * sign / (1.0 - omega_T_i)))**2.0
                
                
                
                
                f_lamda  = np.fabs(sig_T_i_eff_trial - gamma_T * alpha_T_i) - delta_lamda * ( E_T/(1. - omega_T_i) + gamma_T)  - f1 * ft
                f_N = sig_N_eff - sig_N_i_eff_trial +  delta_lamda * E_N * N
#                 f_omega_N = omega_N - omega_N_i  - delta_lamda*(1 -omega_N)**c_N * (Y_N / S_N + b * Y_T /S_T) #* np.heaviside(sig_N_eff, 1)
#                 f_omega_T = omega_T - omega_T_i  - delta_lamda*(1 -omega_T)**c_T * (Y_T / S_T + b * Y_N /S_N)
                
                

                return [f_lamda, f_N]
            
            x0 = np.array([0, 0])
            
            sol   =  root(f, x0 = x0,  method='lm')
            delta_lamda,  sig_N_eff= sol.x
            

            eps_N_p_i +=  delta_lamda* (m  - np.heaviside((sig_N_eff), 1) * ((m * sig_N_eff **2)/ (sig_t)**2 - (sig_0  - m * sig_N_eff) * (2 * sig_N_eff / (sig_t)**2 ) ))
            eps_T_p_i += delta_lamda * np.sign(sig_T_i_eff_trial - gamma_T * alpha_T_i) / (1.0 - omega_T_i)
            

            
            z_T_i += delta_lamda
            alpha_T_i += delta_lamda * np.sign(sig_T_i_eff_trial - gamma_T * alpha_T_i)
            
            Y_N = 0.5 * E_N * (eps_N_i - eps_N_p_i)**2.0
            Y_T = 0.5 * E_T * (eps_T_i - eps_T_p_i)**2.0
            
            omega_N_i +=  delta_lamda*(1 - omega_N_i)**c_N * (Y_N / S_N + b * Y_T /S_T) * np.heaviside(eps_N_i, 1)
            omega_T_i +=  delta_lamda*(1 - omega_T_i)**c_T * (Y_T / S_T + b * Y_N /S_N)
            
            sig_N_i = (1.0 - omega_N_i) *   sig_N_eff
            
            sig_T_i = (1.0 - omega_T_i) * E_T * (eps_T_i - eps_T_p_i) 
            
            
            f = np.fabs(sig_T_i /(1 - omega_T_i)  - gamma_T * alpha_T_i)  -\
            (sig_0  + K_T * z_T_i - m * sig_N_i /(1 - omega_N_i)) * (1.0 - np.heaviside(( sig_N_i /(1 - omega_N_i)), 1) * (( sig_N_i /(1 - omega_N_i))**2 / (sig_t)**2 ))   
                       
            
        else: 
                
            sig_N_i = (1 - omega_N_i) * sig_N_i_eff_trial 
            
            sig_T_i = (1 - omega_T_i) * sig_T_i_eff_trial 
            

            f  = f_trial
            

        sig_N_trial_arr[i] = sig_N_i_eff_trial 
        sig_T_trial_arr[i] = sig_T_i_eff_trial
        
        
        sig_N_arr[i] = sig_N_i 
        sig_T_arr[i] = sig_T_i 

        
        eps_T_p_arr[i] = eps_T_p_i
        eps_N_p_arr[i] = eps_N_p_i

        omega_N_arr[i] = omega_N_i
        omega_T_arr[i] = omega_T_i
        
        z_T_arr[i] = z_T_i
        alpha_T_arr[i] = alpha_T_i
        
        f_arr[i] = f
        
    return  sig_N_trial_arr, sig_T_trial_arr, sig_N_arr, sig_T_arr, eps_T_p_arr, eps_N_p_arr, omega_N_arr, omega_T_arr, z_T_arr, alpha_T_arr, f_arr

#------------------------------------------------------------------------------ 

inc = 1000

s_history = np.array([0, 0.0001])
eps_N_arr = np.hstack([np.linspace(s_history[i], s_history[i + 1], inc)
                     for i in range(len(s_history) - 1)])

t_N_arr= t_arr = np.linspace(0, 1, len(eps_N_arr))



s_history = np.array([0, 0.0002])
eps_T_arr = np.hstack([np.linspace(s_history[i], s_history[i + 1], inc)
                     for i in range(len(s_history) - 1)])

t_T_arr= t_arr = np.linspace(0, 1, len(eps_T_arr))



sig_0=10.0
E_N=50000.
E_T=20000.

sig_t = 5.0
S_N = 0.000000025
S_T = 0.000001
c_N = 1.2
c_T = 2.2
K_T = 0
gamma_T = 0

Rt=1

m=0.2
b=0.2
    
sig_N_trial_arr, sig_T_trial_arr, sig_N_arr, sig_T_arr, eps_T_p_arr, eps_N_p_arr, omega_N_arr, omega_T_arr, z_T_arr, alpha_T_arr, f_arr = get_corr_pred_(
        eps_N_arr, eps_T_arr, sig_0, sig_t, E_T, E_N, m, S_N, S_T, c_N, c_T, K_T, gamma_T)   

state_point= np.array([eps_N_arr[-1], eps_N_p_arr[-1], omega_N_arr[-1],  eps_T_arr[-1], eps_T_p_arr[-1], omega_T_arr[-1], z_T_arr[-1], alpha_T_arr[-1] ])

print('omega_N', omega_N_arr[-1])
print('omega_T', omega_T_arr[-1])

#------------------------------------------------------------------------------ 

def get_corr_pred_exact(state, d_eps_N, d_eps_T, sig_0, sig_t, E_T, E_N, m, S_N, S_T, c_N, c_T, K_T, gamma_T, inc):
    
    # state variables
    eps_N_i = state[0]
    eps_N_p_i = state[1]
    omega_N_i = state[2]
    
    eps_T_i = state[3]
    eps_T_p_i = state[4]
    omega_T_i = state[5]
    z_T_i = state[6] 
    alpha_T_i = state[7]
    
    
    
    eps_N = eps_N_i + d_eps_N
    eps_T = eps_T_i + d_eps_T
    
    eps_N_arr = np.linspace(eps_N_i, eps_N, inc)
    eps_T_arr = np.linspace(eps_T_i, eps_T, inc)
    sig_N_arr= np.zeros_like(eps_N_arr)
    sig_T_arr= np.zeros_like(eps_T_arr)
    sig_N_trial_arr = np.zeros_like(eps_N_arr)
    sig_T_trial_arr = np.zeros_like(eps_T_arr)
    eps_N_p_arr = np.zeros_like(eps_N_arr)
    eps_T_p_arr = np.zeros_like(eps_T_arr)
    omega_N_arr = np.zeros_like(eps_N_arr)
    omega_T_arr = np.zeros_like(eps_T_arr)
    
    z_T_arr = np.zeros_like(eps_T_arr)
    alpha_T_arr = np.zeros_like(eps_T_arr)

    delta_lamda = 0.0
    
    for i in range(1, len(eps_N_arr)):
    
        eps_N_i = eps_N_arr[i]
        eps_T_i = eps_T_arr[i]
        
        sig_N_i_eff_trial = E_N * (eps_N_i - eps_N_p_i)
        sig_T_i_eff_trial  = E_T* (eps_T_i - eps_T_p_i)
        
        f_trial = np.fabs(sig_T_i_eff_trial - gamma_T * alpha_T_i)  -\
         (sig_0  + K_T * z_T_i - m * sig_N_i_eff_trial) * (1.0 - np.heaviside((sig_N_i_eff_trial), 1) * ((sig_N_i_eff_trial)**2 / (sig_t)**2 ))
         
        if f_trial > 1e-8: 
            
            def f(vars):
                
                delta_lamda ,   sig_N_eff  = vars
    
                N =  (m  + np.heaviside((sig_N_eff), 1) * ((-m * sig_N_eff **2)/ (sig_t)**2    +\
                                                                       (sig_0  + K_T * (z_T_i + delta_lamda) - m * sig_N_eff) * (2 * sig_N_eff / (sig_t)**2 ) )) /(1. - omega_N_i)
            
                f1 = sig_0 + K_T * (z_T_i + delta_lamda) - m * sig_N_eff
                
                ft =  1.0 - np.heaviside(sig_N_eff, 1) * (sig_N_eff**2 / (sig_t)**2)
    
                f_lamda  = np.fabs(sig_T_i_eff_trial - gamma_T * alpha_T_i) - delta_lamda * ( E_T/(1. - omega_T_i) + gamma_T)  - f1 * ft
                f_N = sig_N_eff - sig_N_i_eff_trial +  delta_lamda * E_N * N
    
                return [f_lamda, f_N]
            
            x0 = np.array([0, 0])
            
            sol   =  root(f, x0 = x0,  method='lm')
            delta_lamda,  sig_N_eff= sol.x
            
            eps_N_p_i +=  delta_lamda* (m  - np.heaviside((sig_N_eff), 1) * ((m * sig_N_eff **2)/ (sig_t)**2 - (sig_0  - m * sig_N_eff) * (2 * sig_N_eff / (sig_t)**2 ) ))
            eps_T_p_i += delta_lamda * np.sign(sig_T_i_eff_trial - gamma_T * alpha_T_i) / (1.0 - omega_T_i)
            z_T_i += delta_lamda
            alpha_T_i += delta_lamda * np.sign(sig_T_i_eff_trial - gamma_T * alpha_T_i)
            
            Y_N = 0.5 * E_N * (eps_N_i - eps_N_p_i)**2.0
            Y_T = 0.5 * E_T * (eps_T_i - eps_T_p_i)**2.0
            
            omega_N_i +=  delta_lamda*(1 - omega_N_i)**c_N * (Y_N / S_N + b * Y_T /S_T) * np.heaviside(eps_N_i, 1)
            omega_T_i +=  delta_lamda*(1 - omega_T_i)**c_T * (Y_T / S_T + b * Y_N /S_N)
            
            sig_N_i = (1.0 - omega_N_i) *   sig_N_eff
            
            sig_T_i = (1.0 - omega_T_i) * E_T * (eps_T_i - eps_T_p_i) 
            
            f = np.fabs(sig_T_i /(1 - omega_T_i)  - gamma_T * alpha_T_i)  -\
            (sig_0  + K_T * z_T_i - m * sig_N_i /(1 - omega_N_i)) * (1.0 - np.heaviside(( sig_N_i /(1 - omega_N_i)), 1) * (( sig_N_i /(1 - omega_N_i))**2 / (sig_t)**2 ))   
                       
        else:   
            sig_N_i = (1 - omega_N_i) * sig_N_i_eff_trial 
            sig_T_i = (1 - omega_T_i) * sig_T_i_eff_trial 
            f  = f_trial
            
        sig_N_trial_arr[i] = sig_N_i_eff_trial 
        sig_T_trial_arr[i] = sig_T_i_eff_trial
        sig_N_arr[i] = sig_N_i 
        sig_T_arr[i] = sig_T_i 

    return  sig_N_arr[-1], sig_T_arr[-1]


'''
example
'''

d_eps_N = 0.0002
d_eps_T = 0.0004

sig_N_ret = get_corr_pred_exact(state = state_point, d_eps_N=d_eps_N, d_eps_T= d_eps_T ,
                           sig_0=sig_0, sig_t=sig_t, E_T = E_T, E_N = E_N, m=m, S_N=S_N, S_T=S_T, c_N=c_N, c_T=c_T, K_T=K_T, gamma_T= gamma_T, inc=2)[0]

sig_T_ret = get_corr_pred_exact(state = state_point, d_eps_N=d_eps_N, d_eps_T = d_eps_T ,
                           sig_0=sig_0, sig_t=sig_t, E_T = E_T, E_N = E_N, m=m, S_N=S_N, S_T=S_T, c_N=c_N, c_T=c_T, K_T=K_T, gamma_T= gamma_T, inc=2)[1]  
                           
                           
sig_N_exact = get_corr_pred_exact(state = state_point, d_eps_N=d_eps_N, d_eps_T= d_eps_T ,
                           sig_0=sig_0, sig_t=sig_t, E_T = E_T, E_N = E_N, m=m, S_N=S_N, S_T=S_T, c_N=c_N, c_T=c_T, K_T=K_T, gamma_T= gamma_T, inc=100)[0]

sig_T_exact = get_corr_pred_exact(state = state_point, d_eps_N=d_eps_N, d_eps_T = d_eps_T ,
                           sig_0=sig_0, sig_t=sig_t, E_T = E_T, E_N = E_N, m=m, S_N=S_N, S_T=S_T, c_N=c_N, c_T=c_T, K_T=K_T, gamma_T= gamma_T, inc=100)[1]  
                           
print('sig_N_ret=', sig_N_ret)                           
print('sig_T_ret=', sig_T_ret)                                                      
print('sig_N_exact=', sig_N_exact)                           
print('sig_T_exact=', sig_T_exact) 



def get_isoerror_map(state, d_eps_N_arr, d_eps_T_arr, sig_0, sig_t, E_T, E_N, m, S_N, S_T, c_N, c_T, K_T, gamma_T):
    
    sig_N_ret = np.zeros(( len(d_eps_N_arr), len(d_eps_T_arr) ))
    sig_T_ret = np.zeros(( len(d_eps_N_arr), len(d_eps_T_arr) ))
    
    sig_N_exact = np.zeros(( len(d_eps_N_arr), len(d_eps_T_arr) ))
    sig_T_exact = np.zeros(( len(d_eps_N_arr), len(d_eps_T_arr) ))
    
    
    
    for i in range(0,  len(d_eps_N_arr)):
        
        for j in range(0, len(d_eps_T_arr)):
            
            sig_N_ret[i,j] = get_corr_pred_exact(state, d_eps_N_arr[i], d_eps_T_arr[j], sig_0, sig_t, E_T, E_N, m, S_N, S_T, c_N, c_T, K_T, gamma_T, 2)[0]
            sig_T_ret[i,j]= get_corr_pred_exact(state, d_eps_N_arr[i], d_eps_T_arr[j], sig_0, sig_t, E_T, E_N, m, S_N, S_T, c_N, c_T, K_T, gamma_T, 2)[1]
            
            sig_N_exact[i,j] = get_corr_pred_exact(state, d_eps_N_arr[i], d_eps_T_arr[j], sig_0, sig_t, E_T, E_N, m, S_N, S_T, c_N, c_T, K_T, gamma_T, 100)[0]
            sig_T_exact[i,j] = get_corr_pred_exact(state, d_eps_N_arr[i], d_eps_T_arr[j], sig_0, sig_t, E_T, E_N, m, S_N, S_T, c_N, c_T, K_T, gamma_T, 100)[1]
            
            #print(sig_N_exact)
            #print(sig_N_ret)
            
            error = np.sqrt((sig_N_ret - sig_N_exact)**2 + (sig_T_ret - sig_T_exact)**2) / np.sqrt((sig_N_exact)**2 + (sig_T_exact)**2) *100


    return error


d_eps_N_arr = np.linspace(0.0, d_eps_N, 10)
d_eps_T_arr = np.linspace(0.0, d_eps_T, 10)



error = get_isoerror_map(state = state_point, d_eps_N_arr=d_eps_N_arr, d_eps_T_arr = d_eps_T_arr ,
                          sig_0=sig_0, sig_t=sig_t, E_T = E_T, E_N = E_N, m=m, S_N=S_N, S_T=S_T, c_N=c_N, c_T=c_T, K_T=K_T, gamma_T= gamma_T)
  
print(error)
  
  
  
#delta = 0.025
x = d_eps_N_arr
y = d_eps_T_arr
X, Y = np.meshgrid(x, y)
Z = error
  
print(X.shape)
print(Z.shape)
fig, ax = plt.subplots()
CS = ax.contour(X, Y, Z, 20,  cmap='inferno')
ax.clabel(CS, inline=2, fontsize=10)


# CS = ax.contourf(X, Y, Z, 10, cmap='RdGy')
# ax.clabel(CS)


#fig.colorbar()
  
  
plt.show()
