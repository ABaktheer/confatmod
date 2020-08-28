
import matplotlib.pyplot as plt
import numpy as np



"===================="
"P-tau"
"===================="

p_ex = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Lateral pressure_Zhang\data\tp_p.txt')
t_ex = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Lateral pressure_Zhang\data\tp_t.txt')


p_num= np.array([0, 0.05, 0.2,   0.4, 0.8 ])
t_num= np.array([0.68, 0.75, 1.78, 3.12, 5.89 ])

plt.subplot(2,3,2)
plt.plot(p_ex, t_ex *np.sqrt(40) , "ko", linewidth=2)
plt.plot(p_num, t_num, "r", linewidth=2)

plt.xlim(-0.1 , 1.1 )
#plt.ylim(-0.1, 1.2)



"===================="
"P-s0"
"===================="

p_ex = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Lateral pressure_Zhang\data\ps0_p.txt')
s0_ex = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Lateral pressure_Zhang\data\ps0_s0.txt')

p_num= np.array([0, 0.05, 0.2, 0.4, 0.8 ])
s0_num= np.array([0.126,  0.009, 0.022, 0.037, 0.073 ])


plt.subplot(2,3,3)
plt.plot(p_ex, s0_ex, "ko", linewidth=2)
plt.plot(p_num, s0_num, "r", linewidth=2)

plt.xlim(-0.1 , 1.1 )
plt.ylim(0, 0.25)







plt.show()






