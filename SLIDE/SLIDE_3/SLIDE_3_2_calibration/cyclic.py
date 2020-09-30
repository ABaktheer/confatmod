
import matplotlib.pyplot as plt
import numpy as np



"===================="
"Bond-slip response"
"===================="

# monotonic
s_mon = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\monotonic_s.txt')
t_mon = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\monotonic_t.txt')

# S=0.5
s_05_1 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\ts_05_1_s.txt')
t_05_1 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\ts_05_1_t.txt')

s_05_2 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\ts_05_2_s.txt')
t_05_2 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\ts_05_2_t.txt')


# S=2.0
s_2_1 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\ts_2_1_s.txt')
t_2_1 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\ts_2_1_t.txt')

s_2_2 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\ts_2_2_s.txt')
t_2_2 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\ts_2_2_t.txt')


# S=4.0
s_4 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\ts_4_1_s.txt')
t_4 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\ts_4_1_t.txt')


# S=8.0
s_8_1 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\ts_8_1_s.txt')
t_8_1 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\ts_8_1_t.txt')

s_8_2 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\ts_8_2_s.txt')
t_8_2 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\ts_8_2_t.txt')



plt.subplot(4,4,1)
plt.plot(s_mon, t_mon, "--k", linewidth=2)
plt.plot(s_05_1, t_05_1, "k")
plt.plot(s_05_2, t_05_2, "gray")
plt.xlim(-2 , 2 )
plt.ylim(-2, 2.5)


plt.subplot(4,4,5)
plt.plot(s_mon, t_mon, "--k")
plt.plot(s_2_1, t_2_1, "k")
plt.plot(s_2_2, t_2_2, "gray")
plt.xlim(-4 , 4 )
plt.ylim(-2, 2)


plt.subplot(4,4,9)
plt.plot(s_mon, t_mon, "--k")
plt.plot(s_4, t_4, "k")
plt.xlim(-6 , 6 )
plt.ylim(-2, 2)


plt.subplot(4,4,13)
plt.plot(s_mon, t_mon, "--k")
plt.plot(s_8_1, t_8_1, "k")
plt.plot(s_8_2, t_8_2, "gray")
plt.xlim(-10 , 10 )
plt.ylim(-2, 2)




"===================="
"Bond-slip numerical"
"===================="

# monotonic
s_mon = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\model\monotonic_s.txt')
t_mon = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\model\monotonic_t.txt')

# S=0.5
s_05 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\model\ts_05_s.txt')
t_05 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\model\ts_05_t.txt')


# S=2
s_2 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\model\ts_2_s.txt')
t_2 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\model\ts_2_t.txt')
 
 
# S=4
s_4 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\model\ts_4_s.txt')
t_4 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\model\ts_4_t.txt')
 
 
# S=8
s_8 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\model\ts_8_s.txt')
t_8 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\model\ts_8_t.txt')


plt.subplot(4,4,2)
plt.plot(s_mon, t_mon, "--r", linewidth=2)
plt.plot(s_05 , t_05, "r")
plt.xlim(-2 , 2 )
plt.ylim(-2, 2.5)


plt.subplot(4,4,6)
plt.plot(s_mon, t_mon, "--r")
plt.plot(s_2 , t_2, "r")
plt.xlim(-4 , 4 )
plt.ylim(-2, 2)


plt.subplot(4,4,10)
plt.plot(s_mon, t_mon, "--r")
plt.plot(s_4 , t_4, "r")
plt.xlim(-6 , 6 )
plt.ylim(-2, 2)


plt.subplot(4,4,14)
plt.plot(s_mon, t_mon, "--r")
plt.plot(s_8 , t_8, "r")
plt.xlim(-10 , 10 )
plt.ylim(-2, 2)


"===================="
"stress over cycles"
"===================="

N = np.array([1, 2, 3, 4, 5])

# S=0.5 EXP
t_05_1 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\tN_05_1.txt')
t_05_2 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\tN_05_2.txt')

# S=0.5 Mod
tau_05 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\model\ts_05_t.txt')
t_05= np.array([tau_05[5001], tau_05[8335], tau_05[11667], tau_05[15001], tau_05[18335] ])

# S=2
t_2_1 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\tN_2_1.txt')
t_2_2 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\tN_2_2.txt')

# S=2 Mod
tau_2 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\model\ts_2_t.txt')
t_2= np.array([tau_2[5001], tau_2[8335], tau_2[11667], tau_2[15001], tau_2[18335] ])

# S=4
t_4_1 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\tN_4_1.txt')
t_4_2 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\tN_4_2.txt')

# S=4 Mod
tau_4 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\model\ts_4_t.txt')
t_4= np.array([tau_4[5001], tau_4[8335], tau_4[11667], tau_4[15001], tau_4[18335] ])

# S=8
t_8_1 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\tN_8_1.txt')
t_8_2 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\tN_8_2.txt')

# S=8 Mod
tau_8 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\model\ts_8_t.txt')
t_8= np.array([tau_8[5001], tau_8[8335], tau_8[11667], tau_8[15001], tau_8[18335] ])

plt.subplot(4,4,3)
plt.plot(N , t_05_1, "k")
plt.plot(N , t_05_2, "k")
plt.fill_between(N, t_05_1, t_05_2 , facecolor='black', alpha=0.2)
plt.plot(N , t_05, "r", linewidth=2)
plt.xlim(0 , 6 )
plt.ylim(0, 1.2)

plt.subplot(4,4,7)
plt.plot(N , t_2_1, "k")
plt.plot(N , t_2_2, "k")
plt.fill_between(N, t_2_1, t_2_2 , facecolor='black', alpha=0.2)
plt.plot(N , t_2, "r", linewidth=2)
plt.xlim(0 , 6 )
plt.ylim(0, 1.2)


plt.subplot(4,4,11)
plt.plot(N , t_4_1, "k")
plt.plot(N , t_4_2, "k")
plt.fill_between(N, t_4_1, t_4_2 , facecolor='black', alpha=0.2)
plt.plot(N , t_4, "r", linewidth=2)
plt.xlim(0 , 6 )
plt.ylim(0, 1.2)

plt.subplot(4,4,15)
plt.plot(N , t_8_1, "k")
plt.plot(N , t_8_2, "k")
plt.fill_between(N, t_8_1, t_8_2 , facecolor='black', alpha=0.2)
plt.plot(N , t_8, "r", linewidth=2)
plt.xlim(0 , 6 )
plt.ylim(0, 1.2)


"===================="
"Energy dissipation"
"===================="


# S=0.5 EXP
DE_05_1 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\DE_05_1.txt')
DE_05_2 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\DE_05_2.txt')

# S=0.5 Mod
t_05 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\model\t_05.txt')
DE_05 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\model\DE_05.txt')
DE_05= np.array([DE_05[5001], DE_05[8335], DE_05[11667], DE_05[15001], DE_05[18335] ])

# S=2
DE_2_1 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\DE_2_1.txt')
DE_2_2 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\DE_2_2.txt')

# S=2 Mod
DE_2 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\model\DE_2.txt')
DE_2 = np.array([DE_2[5001], DE_2[8335], DE_2[11667], DE_2[15001], DE_2[18335] ])


# S=4
DE_4_1 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\DE_4_1.txt')
DE_4_2 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\DE_4_2.txt')

# S=4 Mod
DE_4 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\model\DE_4.txt')
DE_4 = np.array([DE_4[5001], DE_4[8335], DE_4[11667], DE_4[15001], DE_4[18335] ])



# S=8
DE_8_1 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\DE_8_1.txt')
DE_8_2 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\data\DE_8_2.txt')


# S=8 Mod
DE_8 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Cyclic_Verndermonde\model\DE_8.txt')
DE_8 = np.array([DE_8[5001], DE_8[8335], DE_8[11667], DE_8[15001], DE_8[18335] ])


plt.subplot(4,4,4)
plt.plot(N , DE_05_1, "k")
plt.plot(N , DE_05_2, "k")
plt.fill_between(N, DE_05_1, DE_05_2 , facecolor='black', alpha=0.2)
plt.plot(N , DE_05, "r", linewidth=2)
plt.xlim(0 , 6 )
plt.ylim(0, 20)

plt.subplot(4,4,8)
plt.plot(N , DE_2_1, "k")
plt.plot(N , DE_2_2, "k")
plt.fill_between(N, DE_2_1, DE_2_2 , facecolor='black', alpha=0.2)
plt.plot(N , DE_2, "r", linewidth=2)
plt.xlim(0 , 6 )
plt.ylim(0, 20)


plt.subplot(4,4,12)
plt.plot(N , DE_4_1, "k")
plt.plot(N , DE_4_2, "k")
plt.fill_between(N, DE_4_1, DE_4_2 , facecolor='black', alpha=0.2)
plt.plot(N , DE_4, "r", linewidth=2)
plt.xlim(0 , 6 )
plt.ylim(0, 20)

plt.subplot(4,4,16)
plt.plot(N , DE_8_1, "k")
plt.plot(N , DE_8_2, "k")
plt.fill_between(N, DE_8_1, DE_8_2 , facecolor='black', alpha=0.2)
plt.plot(N , DE_8, "r", linewidth=2)
plt.xlim(0 , 6 )
plt.ylim(0, 30)








plt.show()






