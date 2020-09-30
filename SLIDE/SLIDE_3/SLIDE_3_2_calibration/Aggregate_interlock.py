
import matplotlib.pyplot as plt
import numpy as np



"===================="
"tau-s response exp."
"===================="

# w=0.013
s_013 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Agreggate interlock_Pauly\data\ts_w013_s.txt')
t_013 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Agreggate interlock_Pauly\data\ts_w013_t.txt')

#w=0.025
s_025 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Agreggate interlock_Pauly\data\ts_w025_s.txt')
t_025 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Agreggate interlock_Pauly\data\ts_w025_t.txt')

#w=0.051
s_051 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Agreggate interlock_Pauly\data\ts_w051_s.txt')
t_051 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Agreggate interlock_Pauly\data\ts_w051_t.txt')



"===================="
"tau-s response model."
"===================="
# w=0.013
s_013_ = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Agreggate interlock_Pauly\model\ts_w013_s.txt')
t_013_ = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Agreggate interlock_Pauly\model\ts_w013_t.txt')

#w=0.025
s_025_ = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Agreggate interlock_Pauly\model\ts_w025_s.txt')
t_025_ = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Agreggate interlock_Pauly\model\ts_w025_t.txt')

#w=0.051
s_051_ = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Agreggate interlock_Pauly\model\ts_w051_s.txt')
t_051_ = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Agreggate interlock_Pauly\model\ts_w051_t.txt')


plt.subplot(2,2,1)
plt.plot(s_013, t_013, "--k", linewidth=2)
plt.plot(s_025, t_025, "--k", linewidth=2)
plt.plot(s_051, t_051, "--k", linewidth=2)

plt.plot(s_013_, t_013_, "r", linewidth=2)
plt.plot(s_025_, t_025_, "r", linewidth=2)
plt.plot(s_051_, t_051_, "r", linewidth=2)

plt.xlim(-0.0 , 0.8 )
plt.ylim(-0.0 , 8 )





"===================="
"sigma-tau exp."
"===================="

# w=0.025
sig_025 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Agreggate interlock_Pauly\data\tau-sig_w025_sig.txt')
tau_025 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Agreggate interlock_Pauly\data\tau-sig_w025_tau.txt')

#w=0.051
sig_051 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Agreggate interlock_Pauly\data\tau-sig_w051_sig.txt')
tau_051 = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Agreggate interlock_Pauly\data\tau-sig_w051_tau.txt')


"===================="
"sigma-tau model."
"===================="

# w=0.025
sig_025_ = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Agreggate interlock_Pauly\model\tau-sig_w025_sig.txt')
tau_025_ = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Agreggate interlock_Pauly\model\ts_w025_t.txt')

#w=0.051
sig_051_ = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Agreggate interlock_Pauly\model\tau-sig_w051_sig.txt')
tau_051_ = np.loadtxt( r'H:\Heimarbeit\DFG_ComFatiCon\SLIDE\SLIDE_calibration\Agreggate interlock_Pauly\model\ts_w051_t.txt')



plt.subplot(2,2,2)
plt.plot(sig_025, tau_025, "--k", linewidth=2)
plt.plot(sig_051, tau_051, "--k", linewidth=2)

plt.plot(-sig_025_ , tau_025_, "r", linewidth=2)
plt.plot(-sig_051_ , tau_051_, "r", linewidth=2)



plt.xlim(-0.9 , 6 )
plt.ylim(-0.0 , 9 )


plt.show()






