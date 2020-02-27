import matplotlib.pyplot as plt
import numpy as np

#=========================================================================
# num- summary - equal n
#=========================================================================
ax1 = plt.subplot(111)

Test = np.arange(1,30,1)

tau_u_80 = np.array([31.48, 26.66, 25.21,30.92, 29.13, 30.26, 35.17, 29.77, 24.13, # C80_ds16_lb2.5
                    23.53, 28.84, 27.04, 23.87 ,                                  # C80_ds25_lb2.5
                    24.83, 21.95, 22.87, 22.79, 26.67, 23.00 ])                    # C80_ds16_lb                            


tau_u_120 = np.array([  33.20, 36.88, 33.75, 38.5, 36.57, 36.88,                        # C120_ds16_lb2.5
                    32.89, 32.18, 33.16, 38.34])                                     # C120_ds25_lb2.5


t_u_all = np.hstack((tau_u_80/84.17, tau_u_120/107.04 ))

F_av_80 = np.zeros_like(tau_u_80 )
F_av_120 = np.zeros_like(tau_u_120)

for i in range(0, len(F_av_80)) :
    F_av_80[i] = np.average(t_u_all)
    
for i in range(0, len(F_av_120)) :
    F_av_120[i] = np.average(t_u_all)    
                 

ax1.plot(Test[0:19], tau_u_80/84.17 , 'ro', markersize=7, color='k',)
ax1.plot(Test[19:30], tau_u_120/107.04 , 'ro', markersize=7, color='r',)
ax1.plot(Test[0:19], F_av_80, color='k',)
ax1.plot(Test[19:30], F_av_120, color='r',)

print(F_av_80)
ax1.plot(Test[0], 0.45 , 'ro', markersize=7, color='y',)
ax1.plot(Test[0], 2.5*np.sqrt(84.17)/84.17 , 'ro', markersize=7, color='y',)
ax1.plot(Test[19], 2.5*np.sqrt(107.04)/107.04 , 'ro', markersize=7, color='y',)

#ax1.plot(Test, F_av_1, color='k',)
#ax1.plot(Test, F_av_2, color='r',)
#ax1.plot(Test, F_av_3, color='g',)
#ax1.plot(Test, F_av_4, '--k',)
#ax1.plot(Test, F_av_5, '--r',)



ax1.set_xlabel('n')
ax1.set_ylabel('eps_max')
ax1.set_ylim(0.0, 1)
ax1.legend(loc=2)

plt.show()