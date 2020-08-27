
import matplotlib.pyplot as plt
import numpy as np



"===================="
"Bond-slip response"
"===================="

# monotonic
s_mon = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\monotonic_s.txt')
t_mon = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\monotonic_t.txt')

# S=0.5
s_05 = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\ts_05_2_s.txt')
t_05 = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\ts_05_2_t.txt')


# S=2.0
s_2 = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\ts_2_1_s.txt')
t_2 = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\ts_2_1_t.txt')


# S=4.0
s_4 = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\ts_4_1_s.txt')
t_4 = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\ts_4_1_t.txt')


# S=8.0
s_8 = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\ts_8_1_s.txt')
t_8 = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\ts_8_1_t.txt')





plt.subplot(4,4,1)
plt.plot(s_mon, t_mon, "--k")
plt.plot(s_05 -0.25 , t_05, "k")
plt.xlim(-9 , 9 )
plt.ylim(-2, 2)


plt.subplot(4,4,5)
plt.plot(s_mon, t_mon, "--k")
plt.plot(s_2, t_2, "k")
plt.xlim(-9 , 9 )
plt.ylim(-2, 2)


plt.subplot(4,4,9)
plt.plot(s_mon, t_mon, "--k")
plt.plot(s_4, t_4, "k")
plt.xlim(-9 , 9 )
plt.ylim(-2, 2)


plt.subplot(4,4,13)
plt.plot(s_mon, t_mon, "--k")
plt.plot(s_8, t_8, "k")
plt.xlim(-9 , 9 )
plt.ylim(-2, 2)




"===================="
"Bond-slip numerical"
"===================="





"===================="
"stress over cycles"
"===================="

N = np.array([1, 2, 3, 4, 5])

# S=0.5
t_05_1 = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\tN_05_1.txt')
t_05_2 = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\tN_05_2.txt')

# S=2
t_2_1 = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\tN_2_1.txt')
t_2_2 = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\tN_2_2.txt')

# S=4
t_4_1 = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\tN_4_1.txt')
t_4_2 = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\tN_4_2.txt')

# S=8
t_8_1 = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\tN_8_1.txt')
t_8_2 = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\tN_8_2.txt')



plt.subplot(4,4,3)
plt.plot(N , t_05_1, "k")
plt.plot(N , t_05_2, "k")
plt.xlim(0 , 6 )
plt.ylim(0, 1.5)

plt.subplot(4,4,7)
plt.plot(N , t_2_1, "k")
plt.plot(N , t_2_2, "k")
plt.xlim(0 , 6 )
plt.ylim(0, 1.5)


plt.subplot(4,4,11)
plt.plot(N , t_4_1, "k")
plt.plot(N , t_4_2, "k")
plt.xlim(0 , 6 )
plt.ylim(0, 1.5)

plt.subplot(4,4,15)
plt.plot(N , t_8_1, "k")
plt.plot(N , t_8_2, "k")
plt.xlim(0 , 6 )
plt.ylim(0, 1.5)


"===================="
"Energy dissipation"
"===================="


# S=0.5
DE_05_1 = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\DE_05_1.txt')
DE_05_2 = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\DE_05_2.txt')


# S=2
DE_2_1 = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\DE_2_1.txt')
DE_2_2 = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\DE_2_2.txt')


# S=4
DE_4_1 = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\DE_4_1.txt')
DE_4_2 = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\DE_4_2.txt')


# S=8
DE_8_1 = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\DE_8_1.txt')
DE_8_2 = np.loadtxt( r'D:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\DE_8_2.txt')




plt.subplot(4,4,4)
plt.plot(N , DE_05_1, "k")
plt.plot(N , DE_05_2, "k")
plt.xlim(0 , 6 )
plt.ylim(0, 30)

plt.subplot(4,4,8)
plt.plot(N , DE_2_1, "k")
plt.plot(N , DE_2_2, "k")
plt.xlim(0 , 6 )
plt.ylim(0, 30)


plt.subplot(4,4,12)
plt.plot(N , DE_4_1, "k")
plt.plot(N , DE_4_2, "k")
plt.xlim(0 , 6 )
plt.ylim(0, 30)

plt.subplot(4,4,16)
plt.plot(N , DE_8_1, "k")
plt.plot(N , DE_8_2, "k")
plt.xlim(0 , 6 )
plt.ylim(0, 30)








plt.show()






