
import matplotlib.pyplot as plt
import numpy as np



"===================="
"P-tau"
"===================="

p_ex = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Lateral pressure_Zhang\data\tp_p.txt')
t_ex = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Lateral pressure_Zhang\data\tp_t.txt')


p_num= np.array([0, 0.05,  0.2, 0.4, 0.8 ])
t_num= np.array([0.68, 0.75, 1.78, 3.12, 5.89 ])

p_num= np.array([0, 0.1,  0.2, 0.3, 0.4, 0.5, 0.6, 0.8 ])
t_num= [0.69255112, 1.07610105, 1.7499946,  2.41739373, 3.12625251, 3.79241707, 4.47294589, 5.81963928]

t_num= [0.58669774, 0.99714549, 1.58717435, 2.2004008,  2.80300464, 3.39078156, 3.96829484, 5.19438878]

plt.subplot(2,3,2)
plt.plot(p_ex, t_ex *np.sqrt(40) , "ko", linewidth=2)
plt.plot(p_num, t_num, "-ro", linewidth=2)

plt.xlim(-0.1 , 1.1 )
#plt.ylim(-0.1, 1.2)



"===================="
"P-s0"
"===================="

p_ex = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Lateral pressure_Zhang\data\ps0_p.txt')
s0_ex = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Lateral pressure_Zhang\data\ps0_s0.txt')

p_num= np.array([0, 0.05, 0.2, 0.4, 0.8 ])
s0_num= np.array([0.126,  0.009, 0.022, 0.037, 0.073 ])

p_num= np.array([0, 0.1,  0.2, 0.3, 0.4, 0.5, 0.6, 0.8 ])
s0_num= np.array([0.12985972, 0.01382766, 0.02224449, 0.03066132, 0.03907816, 0.04749499, 0.05591182, 0.07274549])
s0_num= np.array([0.12234469, 0.01683367, 0.02645291, 0.03667335, 0.04689379, 0.05651303, 0.06673347, 0.08657315])

plt.subplot(2,3,3)
plt.plot(p_ex, s0_ex, "ko", linewidth=2)
plt.plot(p_num, s0_num, "-ro", linewidth=2)

plt.xlim(-0.1 , 1.1 )
plt.ylim(0, 0.25)







plt.show()






