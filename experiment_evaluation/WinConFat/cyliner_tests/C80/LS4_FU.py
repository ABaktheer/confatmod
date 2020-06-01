
''''
Created on 20.10.2018

@author: Abedulgader Baktheer
'''


import matplotlib.pyplot as plt
import numpy as np






u_1 = np.load( r'D:\Work\WinConFat\C80 Charge 2\High\NPY\CT80-42_3610_Zykl_avg_WA_1_WA_2_WA_3.npy')

F_1 = np.load( r'D:\Work\WinConFat\C80 Charge 2\High\NPY\CT80-42_3610_Zykl_Kraft.npy')




plt.subplot(111)
plt.plot(abs(u_1)/300., abs(F_1) *1000 /(np.pi * 100**2/4), 'k')

plt.xlim(0, 0.0032)
plt.ylim(0, 90)
plt.title('stress - strain (H)')
plt.xlabel('Displacement [mm]')
plt.ylabel('Force [kN]')





plt.show()
