'''
Created on 29.01.2020

@author: abaktheer
'''



import matplotlib.pyplot as plt
import numpy as np
  
m = np.array([0.05, ])
f = np.zeros_like(m)
tau_bar = 2.0
sig_c1 = -20.0
sig_c2 = -40.0
sig_t1 = -0.0
sig_t2 =  5.0
Rc = 1.0
Rt = 1.0


  
for i in np.arange(0, len(m)):  

    sig_N = np.linspace(sig_c2 , tau_bar / m[i], 1000) 
    
    f_old = (tau_bar - m[i] * sig_N) 
    
    #f_t = (1.0 - np.heaviside((sig_N - sig_t1), 1) * ((sig_N - sig_t1)**2 / (sig_t2 - sig_t1)**2  ))
    
    f_t = (1.0 - np.heaviside((sig_N), 1) * ((sig_N)**2 / (sig_t2)**2  ))
    
    #f_c = (1.0 - np.heaviside((sig_c - sig_N),1) *((sig_N - sig_c)**2 / (Rc*(tau_bar - m[i] * sig_c))**2)  )
    
    f_c = (1.0 - np.heaviside((sig_c1 - sig_N),1) * ((sig_N - sig_c1)**2 / (sig_c2 - sig_c1)**2  ))
    
    f = (tau_bar - m[i] * sig_N)**2  * f_t  * f_c
    

    
    plt.subplot(221) 
    plt.plot(sig_N , f_old, 'k', linewidth=2.0, alpha=1.0)
    plt.plot(sig_N , f, 'r', linewidth=2.0, alpha=1.0)

    plt.ylim(0,150)
    plt.axhline(y=0, color='k', linewidth=0.5, alpha=0.5)
    plt.axvline(x=0, color='k', linewidth=0.5, alpha=0.5)
    
    E_N = 1000
    E_T = 1000
    eps_N  = 0.1 # np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]) *0.1
    eps_T  = 0.2 # np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]) *0.1
    
    omega_N = 0.0
    omega_T = 0.0
    c_N = 2
    c_T = 2
    S_N = 0.01
    S_T = 0.01
    b= 0.1
    
    Y_T = 0.5 * E_T * (eps_T -  0)**2
    Y_N = 0.5 * E_N * (eps_N)**2
    
    
    Y_N = 0.5 * eps_N * sig_N 
    Y_T = 0.5 * eps_T * (f_old - 0)
    
    phi =   f_old + ((1 - omega_T) ** c_T) * (b * Y_N * Y_T /S_N + Y_T**2 /S_T) + ((1 - omega_N) ** c_N) * (Y_N * Y_T/S_N + b* Y_T**2 /S_T) * np.heaviside(sig_N,1)  
  
    plt.subplot(222) 
    plt.plot(sig_N , f_old, 'k', linewidth=2.0, alpha=1.0)
    plt.plot(sig_N , phi, 'r', linewidth=2.0, alpha=1.0)

    #plt.ylim(0,40)
    plt.axhline(y=0, color='k', linewidth=0.5, alpha=0.5)
    plt.axvline(x=0, color='k', linewidth=0.5, alpha=0.5)

plt.show()
