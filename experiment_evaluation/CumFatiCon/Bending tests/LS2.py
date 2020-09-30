'''
Created on 27.03.2019

@author: abaktheer
'''


import matplotlib.pyplot as plt
import numpy as np

#===================================
"LS2"
#===================================

"#B-L60"


u_4 = np.load( r'H:\DFG_ComFatiCon\rohedaten\NPY\B-L60-04_Wa_Riss mm.npy')
F_4 = np.load( r'H:\DFG_ComFatiCon\rohedaten\NPY\B-L60-04_Kraft kN.npy')


plt.subplot(111)

plt.plot((abs(u_4)), abs(F_4), "r")
#plt.plot((abs(u_3)), abs(F_3), "k")

plt.xlim(-0.01, 0.5)
plt.ylim(-0.05, 20)
plt.xlabel('eps')
plt.ylabel('sigma')

plt.show()