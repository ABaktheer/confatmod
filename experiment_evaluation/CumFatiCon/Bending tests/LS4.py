'''
Created on 27.03.2019

@author: abaktheer
'''


import matplotlib.pyplot as plt
import numpy as np

#===================================
"LS2"
#===================================

"#C40"


u_1 = np.load( r'H:\DFG_ComFatiCon\rohedaten\NPY\B-L60-05-0000-430_ W10TK.npy')
F_1= np.load( r'H:\DFG_ComFatiCon\rohedaten\NPY\B-L60-05-0000-430_ Kraft.npy')



eps_max_1= np.load( r'H:\DFG_ComFatiCon\rohedaten\NPY\B-L60-05-0000-430_ W10TK_max.npy')
N_1 =  len(eps_max_1)
N_max_1 = np.arange(1, N_1 +1, 1)




plt.subplot(221)

plt.plot((abs(u_1))/1.14, abs(F_1), "k")

plt.xlim(-0.01, 0.2)
plt.ylim(-0.05, 18)
plt.title('Smax = 0.85')
plt.xlabel('crack opening (w) [mm]')
plt.ylabel('force [kN]')



plt.subplot(222)
plt.plot(N_max_1, (abs(eps_max_1))/1.14, "k")


plt.ylim(0.0, 0.2)
#plt.xlim(-100, 750000)
plt.title('Smax = 0.85')
plt.xlabel('number of cycles ')
plt.ylabel('max (w) [mm]')

#========================================================================

u_2 = np.load( r'H:\DFG_ComFatiCon\rohedaten\NPY\B-L60-15-4427_avg_Riss_1 mm_Riss_2 mm.npy')
F_2= np.load( r'H:\DFG_ComFatiCon\rohedaten\NPY\B-L60-15-4427_Kraft kN.npy')



eps_max_2= np.load( r'H:\DFG_ComFatiCon\rohedaten\NPY\B-L60-15-4427_avg_Riss_1 mm_Riss_2 mm_max.npy')
N_2 =  len(eps_max_2)
N_max_2 = np.arange(1, N_2 +1, 1)



plt.subplot(223)

plt.plot((abs(u_2)), abs(F_2), "k")

plt.xlim(-0.01, 0.2)
plt.ylim(-0.05, 18)
plt.title('Smax = 0.70')
plt.xlabel('crack opening (w) [mm]')
plt.ylabel('force [kN]')

plt.subplot(224)
plt.plot(N_max_2, (abs(eps_max_2))/1.14, "k")


plt.ylim(0.0, 0.2)
#plt.xlim(-100, 750000)
plt.title('Smax = 0.70')
plt.xlabel('number of cycles ')
plt.ylabel('max (w) [mm]')








plt.show()