import matplotlib.pyplot as plt
import numpy as np

#=========================================================================
# num- summary - equal n
#=========================================================================
ax1 = plt.subplot(111)

Test = np.arange(1,28,1)

F_crack = np.array([84.02, 90.99, 76.65, 82.74, 87.09, 95.34, 95.87, 99.25, 99.87, # C80_ds16_lb2.5
                    95.44, 72.24, 62.08, 73.94  ,                                 # C80_ds25_lb2.5
                    72.18, 90.69, 87.22, 80.06, 68.91, 86.36  ,                   # C80_ds16_lb5
                    98.20, 97.1, 81.5, 91.7   ,                                    # C120_ds16_lb2.5
                    88.13, 95.46, 98.1, 97.2])                                       # C120_ds25_lb2.5

print(np.sum(F_crack ))

F_av_1 = np.zeros_like(F_crack)
F_av_2 = np.zeros_like(F_crack)
F_av_3 = np.zeros_like(F_crack)
F_av_4 = np.zeros_like(F_crack)
F_av_5 = np.zeros_like(F_crack)
F_av_C80 = np.zeros_like(F_crack)
F_av_C120 = np.zeros_like(F_crack)


for i in range(0, len(F_av_1)) :
    F_av_1[i] = 90.20
    F_av_2[i] = 75.93
    F_av_3[i] = 80.91 
    F_av_4[i] = 92.13 
    F_av_5[i] = 94.72
    F_av_C80[i] = np.average(F_crack[0:19])
    F_av_C120[i] = np.average(F_crack[19:28])

print( len(F_av_1))     
print(F_av_C80)               
print(F_av_C120)

ax1.plot(Test, F_crack, 'ro', markersize=7, color='g',)
ax1.plot(Test, F_av_C80, color='k')
ax1.plot(Test, F_av_C120, color='r')

#ax1.plot(Test, F_av_1, color='k',)
#ax1.plot(Test, F_av_2, color='r',)
#ax1.plot(Test, F_av_3, color='g',)
#ax1.plot(Test, F_av_4, '--k',)
#ax1.plot(Test, F_av_5, '--r',)



ax1.set_xlabel('n')
ax1.set_ylabel('eps_max')
ax1.set_ylim(0.0, 110)
ax1.legend(loc=2)

plt.show()