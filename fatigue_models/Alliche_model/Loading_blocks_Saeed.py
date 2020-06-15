
import matplotlib.pyplot as plt
import numpy as np


class alliche_model:

#-----------------------------------------------------------------------
# material parameters     
#-----------------------------------------------------------------------
 
    lamda = 12500
    mu=18750
    alpha=2237.5
    beta=-2216.5
    g=-10.0
    C0=0.00
    C1=0.0019
    K=0.00485
    n=10
        
    def __init__(self, m, n1, n2, sigma_u, H_max_level, L_max_level, min_level, rep_num):


        self.m = m
        self.n1 = n1
        self.n2 = n2
        self.sigma_u = sigma_u
        self.H_max_level = H_max_level
        self.L_max_level = L_max_level
        self.min_level = min_level
        self.rep_num = rep_num
        
    def get_sigma_arr_2block(self):
            
        n  = self.rep_num * ( self.n1 + self.n2)
       
        stress_level_1_max = self.H_max_level * self.sigma_u
        stress_level_2_max = self.L_max_level * self.sigma_u
    
        stress_level_1_min = self.min_level * self.sigma_u
        stress_level_2_min = self.min_level * self.sigma_u
    
        d_0 = np.zeros(1)
    
        d_1 = np.linspace(0, stress_level_1_max, self.n1 * 2)
        d_1.reshape(-1, 2)[:, 0] = stress_level_1_max
        d_1.reshape(-1, 2)[:, 1] = stress_level_1_min
        d_history_1 = d_1.flatten()
        sig_1_arr = np.hstack([np.linspace(d_history_1[i], d_history_1[i + 1], self.m, dtype=np.float_)
                               for i in range(len(d_1) - 1)])
        
        d_2 = np.linspace(0, stress_level_2_max, self.n2 * 2)
        d_2.reshape(-1, 2)[:, 0] = stress_level_2_max
        d_2.reshape(-1, 2)[:, 1] = stress_level_2_min
        d_history_2 = d_2.flatten()
        sig_2_arr = np.hstack([np.linspace(d_history_2[i], d_history_2[i + 1], self.m, dtype=np.float_)
                               for i in range(len(d_2) - 1)])
    
    
        sig_0_1_arr = np.linspace(
            d_0[-1], d_history_1[0], self.m, dtype=np.float_)
        
        sig_1_2_arr = np.linspace(
            d_history_1[-1], d_history_2[0], self.m, dtype=np.float_)
    
        sig_1_1_arr = np.linspace(
            sig_1_arr[-1], sig_1_arr[0], self.m, dtype=np.float_)
        
        sig_2_2_arr = np.linspace(
            sig_2_arr[-1], sig_2_arr[0], self.m, dtype=np.float_)
        
       
        sigma_data = np.hstack(
            self.rep_num * (sig_1_arr[1:-1], sig_1_2_arr, sig_2_arr[1:-1]))
        
#         print('sigma_data', np.shape(sigma_data))

        sigma_1_arr = np.hstack(
            (d_0, sig_0_1_arr, sigma_data))
        
#         
#         sigma_1_arr = np.hstack(
#             (d_0, sig_0_1_arr,sig_1_arr[1:-1], sig_1_2_arr, sig_2_arr[1:-1]))
#         print('sigma_1_arr shape', np.shape(sigma_1_arr))
        sigma_arr = sigma_1_arr

        t_arr = np.linspace(0, 1, len(sigma_arr))

        return sigma_arr, t_arr   
    
     
    def get_sigma_arr_2step(self):
            
        n = self.n1 + self.n2
       
        stress_level_1_max = self.H_max_level * self.sigma_u
        stress_level_2_max = self.L_max_level * self.sigma_u
    
        stress_level_1_min = self.min_level * self.sigma_u
        stress_level_2_min = self.min_level * self.sigma_u
    
        d_0 = np.zeros(1)
    
        d_1 = np.linspace(0, stress_level_1_max, self.n1 * 2)
        d_1.reshape(-1, 2)[:, 0] = stress_level_1_max
        d_1.reshape(-1, 2)[:, 1] = stress_level_1_min
        d_history_1 = d_1.flatten()
        sig_1_arr = np.hstack([np.linspace(d_history_1[i], d_history_1[i + 1], self.m, dtype=np.float_)
                               for i in range(len(d_1) - 1)])
        
        d_2 = np.linspace(0, stress_level_2_max, self.n2 * 2)
        d_2.reshape(-1, 2)[:, 0] = stress_level_2_max
        d_2.reshape(-1, 2)[:, 1] = stress_level_2_min
        d_history_2 = d_2.flatten()
        sig_2_arr = np.hstack([np.linspace(d_history_2[i], d_history_2[i + 1], self.m, dtype=np.float_)
                               for i in range(len(d_2) - 1)])
    
    
        sig_0_1_arr = np.linspace(
            d_0[-1], d_history_1[0], self.m, dtype=np.float_)
        
        sig_1_2_arr = np.linspace(
            d_history_1[-1], d_history_2[0], self.m, dtype=np.float_)
    
        sig_1_1_arr = np.linspace(
            sig_1_arr[-1], sig_1_arr[0], self.m, dtype=np.float_)
        
        sig_2_2_arr = np.linspace(
            sig_2_arr[-1], sig_2_arr[0], self.m, dtype=np.float_)
        
    
        sigma_1_arr = np.hstack(
            (d_0, sig_0_1_arr,sig_1_arr[1:-1], sig_1_2_arr, sig_2_arr[1:-1]))
        
#         print('sigma_1_arr shape 2stepi', np.shape(sigma_1_arr))

        sigma_arr = sigma_1_arr
        t_arr = np.linspace(0, 1, len(sigma_arr))

        return sigma_arr, t_arr 
        
        
    def get_stress_strain(self,sigma_arr):
        #-----------------------------------------------------------------------
        # arrays to store the values
        #-----------------------------------------------------------------------
        # normal strain
        eps_1_arr = np.zeros_like(sigma_arr, dtype=np.float_)
        # lateral strain
        eps_2_arr = np.zeros_like(sigma_arr, dtype=np.float_)
        # damaself.ge factor
        w_arr = np.zeros_like(sigma_arr, dtype=np.float_)
        f_arr = np.zeros_like(sigma_arr, dtype=np.float_)
        D_arr = np.zeros_like(sigma_arr, dtype=np.float_)
        phi_arr = np.zeros_like(sigma_arr, dtype=np.float_)
    

        #-----------------------------------------------------------------------
        # state variables
        #-----------------------------------------------------------------------

        eps_1_i = 0.0
        eps_2_i = 0.0
        w_i = 0.0
        D_i = 0.0
    
        for i in range(1, len(sigma_arr)):
            cycle = round(i/(self.m*2), 0)
    #         print(cycle)
            sigma_1_i = sigma_arr[i]
    
            eps_2_i = -1.0 * ((self.lamda + self.alpha * w_i) * sigma_1_i + self.g * w_i * (self.lamda + 2.0 * self.mu)) / \
                ((self.lamda + 2.0 * self.mu) * (2.0 * (self.lamda + self.mu) + 4.0 *
                                       w_i * (self.alpha + self.beta)) - 2.0 * (self.lamda + self.alpha * w_i) ** 2)
    
            eps_1_i = sigma_1_i / \
                (self.lamda + 2.0 * self.mu) - 2.0 * eps_2_i *  \
                (self.lamda + self.alpha * w_i) / (self.lamda + 2.0 * self.mu)
    
            f_i = abs(self.g) * eps_2_i - (self.C0 + 2 * self.C1 * w_i)
    
            kappa_i = (self.lamda + 2.0 * self.mu) * (2.0 * (self.lamda + self.mu) +
                                            4.0 * w_i * (self.alpha + self.beta) -
                                            self.alpha * (self.g / (2.0 * self.C1)) *
                                            (2.0 * eps_2_i + eps_1_i) -
                                            (self.g**2.0 / (2.0 * self.C1))) - 2.0 * (self.lamda + self.alpha * w_i)**2
    
            d_sigma_1 = sigma_arr[i] - sigma_arr[i - 1]
            m = -1.0 * ((self.lamda + self.alpha * w_i) / kappa_i) * d_sigma_1
    
            # loadinself.g staself.ge (evolve of the fatiself.gue damaself.ge based on (Mariself.go.85)
            # model)
            if m > 0:
                d_w = m * abs(self.g) / (2.0 * self.C1) * (f_i / self.K)**self.n
            else:  # unloadinself.g staself.ge (no fatiself.gue damaself.ge)
                d_w = 0
    
            w_i = w_i + d_w
    
            # Enerself.gy release rate
            Y_norm = np.sqrt((-self.g * eps_1_i - self.alpha * (eps_1_i + 2. * eps_2_i) * eps_1_i
                              - 2. * self.beta * (eps_1_i**2.0))**2.0 + 2.0 * (-self.g * eps_2_i - self.alpha * (eps_1_i + 2. * eps_2_i) * eps_1_i
                                                                          - 2 * self.beta * (eps_1_i**2.0))**2.0)
            d_D = Y_norm  
            D_i += d_D
    
            # Helmholtz free enerself.gy
            phi_i = 0.5 * self.lamda * (eps_1_i + 2.0 * eps_2_i)**2.0 + self.mu * ((eps_1_i)**2.0 + 2.0 * eps_2_i**2.0) + 2.0 * self.g * w_i * eps_2_i + self.alpha * \
                (2.0 * w_i * eps_1_i * eps_2_i + 4.0 * w_i *
                 eps_2_i**2.0) + 4.0 * self.beta * w_i * eps_2_i**2.0
    
    
            if w_i > 5.0:
                print(' ----------> No Convergence any more')
                print(i)
                break
    
            if abs(eps_1_i) > 0.005:
                print(' ----------> No Convergence any more')
                print(i)
                break
    
            eps_1_arr[i] = eps_1_i
            eps_2_arr[i] = eps_2_i
            w_arr[i] = w_i
            f_arr[i] = f_i
            D_arr[i] = D_i
            phi_arr[i] = phi_i
            
            eps_1_i = eps_1_arr[i]
            eps_2_i = eps_2_arr[i] 
            w_i = w_arr[i]
            D_i = D_arr[i]
            
        return eps_1_i, eps_2_i, w_i, D_i, eps_1_arr, eps_2_arr, w_arr, f_arr, D_arr, i, cycle, phi_arr
    
    
    def plot(self,t_arr,sigma_arr,w_arr,D_arr,inc,eps_1_arr):
        b = 1
        #-----------------------------------------------------------------------
        # plot 1
        #-----------------------------------------------------------------------
        plt.subplot(231)
        plt.plot(t_arr[0:inc], abs(sigma_arr[0:inc]), 'k', color='blue', linewidth=0.1, alpha=1.0)
        plt.title('loading history')
     
        plt.xlabel('Time')
        plt.ylabel('$\sigma_{1}$')
     
        #-----------------------------------------------------------------------
        # plot 2
        #-----------------------------------------------------------------------
        plt.subplot(232)
        plt.plot(abs(eps_1_arr[0:inc]), abs(sigma_arr[0:inc]), 'k', color='gray', linewidth=0.1, alpha=1.0)
        plt.title('$ \epsilon_{11} - \sigma_{11}$')
        plt.xlabel('$\epsilon_{11}$')
        plt.ylabel('$\sigma_{11}$[MPa]')
     
        #-----------------------------------------------------------------------
        # plot 3
        #-----------------------------------------------------------------------
        plt.subplot(233)
        plt.plot(abs(sigma_arr[0:inc]), w_arr[0:inc], 'k', color='black', linewidth=0.1, alpha=1.0)
        plt.title('$ \epsilon_{11} - \sigma_{11}$')
        plt.xlabel('$\epsilon_{11}$')
        plt.ylabel('$\sigma_{11}$[MPa]')
     
        PlotCalibration_n = self.rep_num * (self.n1 + self.n2)
        #-----------------------------------------------------------------------
        # plot 4
        #-----------------------------------------------------------------------
        plt.subplot(234)
     
        eps_1_max = np.zeros(PlotCalibration_n)
        eps_1_min = np.zeros(PlotCalibration_n)
        cycle = np.zeros(PlotCalibration_n)
        for i in range(0, PlotCalibration_n, 1):
            idx_1 = self.m + 2 * i * self.m - 1
            idx_2 = 2 * i * self.m
            if idx_1 <= len(eps_1_arr[0:inc]):
                idx_1 = idx_1
            else:
                idx_1 = self.m + 2 * (i - 1.0) * self.m - 1
                break
     
            if idx_2 <= len(eps_1_arr[0:inc]):
                idx_2 = idx_2
            else:
                idx_2 = 1 * (i - 1.0) * self.m
                break
     
            eps_1_max[i] = eps_1_arr[int(idx_1)]
            eps_1_min[i] = eps_1_arr[int(idx_2)]
            cycle[i] = i + 1
            
        plt.plot(cycle[0:i], abs(eps_1_max[0:i]), 'k', color='pink', linewidth=1, alpha=1)
              
        plt.xlabel('number of cycles')
        plt.ylabel('max $\epsilon_{11}$')
        
     
        #-----------------------------------------------------------------------
        # plot 5
        #-----------------------------------------------------------------------
        plt.subplot(235)
        PlotCalibration_n = self.rep_num * (self.n1 + self.n2)
        w = np.zeros(PlotCalibration_n)
        cycle = np.zeros(PlotCalibration_n)
        for i in range(0, PlotCalibration_n, 1):
            idx = self.m + 2 * i * self.m - 1
            if idx <= len(w_arr[0:inc]):
                idx = idx
            else:
                idx = self.m + 2 * (i - 1.0) * self.m - 1
                break
     
            w[i] = w_arr[idx]
            cycle[i] = i + 1
     
        plt.plot(cycle[0:i], w[0:i], 'k',color='green', linewidth=1, alpha=1)
        plt.xlabel('number of cycles')
        plt.ylabel('Damage') 
     
        #-----------------------------------------------------------------------
        # plot 6
        #-----------------------------------------------------------------------
        plt.subplot(236)
        D = np.zeros(PlotCalibration_n)
        cycle = np.zeros(PlotCalibration_n)
        for i in range(0, PlotCalibration_n, 1):
            idx = self.m + 2 * i * self.m - 1
            if idx <= len(D_arr[0:inc]):
                idx = idx
            else:
                idx = self.m + 2 * (i - 1.0) * self.m - 1
                break
     
            D[i] = D_arr[idx]
            cycle[i] = i + 1
     
        plt.plot(cycle[0:i], D[0:i], 'k', color='black', linewidth=1, alpha=1)
        plt.xlabel('number of cycles')
        plt.ylabel('$Dissipation$')
         
        fig = plt.gcf()
    
        manager = plt.get_current_fig_manager()
        manager.window.showMaximized()
        plt.grid(True)
        plt.tight_layout()
        plt.tight_layout()
        plt.show()
         
        
            
Study = alliche_model(m = 50, n1 = 50, n2=20 , 
                      sigma_u= -120, H_max_level = 0.80, 
                      L_max_level = 0.75, min_level = 0.20, rep_num = 50)

sigma_arr,t_arr = Study.get_sigma_arr_2step()
sigma_arr,t_arr = Study.get_sigma_arr_2block()
# print(sigma_arr)

eps_1_i, eps_2_i, w_i, D_i, eps_1_arr, eps_2_arr, w_arr, f_arr, D_arr, inc, cycle, phi_arr = Study.get_stress_strain(sigma_arr)
  
Study.plot(t_arr,sigma_arr,w_arr,D_arr,inc,eps_1_arr)       
