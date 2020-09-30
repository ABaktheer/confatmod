'''
Created on 27.03.2019

@author: abaktheer
'''


import matplotlib.pyplot as plt
import numpy as np

#===================================
"LS1"
#===================================

"#B-L60"

u_1 = np.load( r'H:\DFG_ComFatiCon\rohedaten\NPY\B-L60-01_Wa_Riss mm.npy')
F_1 = np.load( r'H:\DFG_ComFatiCon\rohedaten\NPY\B-L60-01_Kraft kN.npy')


u_2 = np.load( r'H:\DFG_ComFatiCon\rohedaten\NPY\B-L60-02_Wa_Riss mm.npy')
F_2 = np.load( r'H:\DFG_ComFatiCon\rohedaten\NPY\B-L60-02_Kraft kN.npy')


u_3 = np.load( r'H:\DFG_ComFatiCon\rohedaten\NPY\B-L60-03_Wa_Riss mm.npy')
F_3 = np.load( r'H:\DFG_ComFatiCon\rohedaten\NPY\B-L60-03_Kraft kN.npy')

plt.subplot(111)
plt.plot((abs(u_1))/2.5, abs(F_1), "k")
plt.plot((abs(u_2))/2.5, abs(F_2), "r")
plt.plot((abs(u_3)), abs(F_3), "g")

plt.xlim(-0.01, 0.5)
plt.ylim(-0.05, 20)
plt.xlabel('eps')
plt.ylabel('sigma')

#========================================================================



plt.show()