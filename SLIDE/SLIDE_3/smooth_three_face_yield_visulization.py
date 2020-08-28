'''
Created on 13.01.2020

@author: abaktheer
'''



 
import matplotlib.pyplot as plt
import numpy as np
  
m = np.array([0.1 , 0.2, 0.3, 0.4, 0.5]) *2 
f = np.zeros_like(m)
tau_bar = 10.0
sig_c1 = -20.0
sig_c2 = -40.0
sig_t =  5.0
#Rc = 1.0
#Rt = 1.0


  
for i in np.arange(0, len(m)):  

    sig_N = np.linspace(sig_c2 , tau_bar/m[i] , 100000) 
    
    f_old = (tau_bar  - m[i] * sig_N) 
    
    #f_t = (1.0 - np.heaviside((sig_N - sig_t1), 1) * ((sig_N - sig_t1)**2 / (sig_t2 - sig_t1)**2  ))
    
    f_t = np.sqrt(1.0 - np.heaviside((sig_N), 1) * ((sig_N)**2 / (sig_t)**2))
    
    #f_c = (1.0 - np.heaviside((sig_c - sig_N),1) *((sig_N - sig_c)**2 / (Rc*(tau_bar - m[i] * sig_c))**2)  )
    
    f_c = np.sqrt((1.0 - np.heaviside((sig_c1 - sig_N),1) * ((sig_N - sig_c1)**2 / (sig_c2 - sig_c1)**2  )))
    
    f = (tau_bar - m[i] * sig_N) * f_t  #* f_c
    

     

     
    
    plt.subplot(111) 
    plt.plot(sig_N , f_old, 'k', linewidth=2.0, alpha=1.0)
    plt.plot(sig_N , f, 'r', linewidth=2.0, alpha=1.0)
    plt.fill_between(sig_N, f, y2 = 0)

    #plt.plot(sig_N , 20*np.heaviside((sig_N - sig_t1), 1), '--r', linewidth=2.0, alpha=1.0)
    #plt.plot(sig_N , 20*np.heaviside((sig_c - sig_N),1), '--b', linewidth=2.0, alpha=1.0)
    #plt.ylim(0,150)
    #plt.plot(sig_N , -f, 'k', linewidth=2, alpha=1.0)
    plt.axhline(y=0, color='k', linewidth=0.5, alpha=0.5)
    plt.axvline(x=0, color='k', linewidth=0.5, alpha=0.5)


plt.show()
 


 

