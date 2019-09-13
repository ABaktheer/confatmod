
''''
Created on 20.10.2018

@author: Abedulgader Baktheer
'''


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import numpy as np
    from io import StringIO
    import os


#=========================================================================
''' C120 - LS4 - Constant'''
#=========================================================================
#
#=============================================
''' PSP 1000 - Test CT_120_1_20'''
#=============================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'CT')
path = os.path.join(path, 'C120')
path = os.path.join(path, 'CT_120_1_20_LS4L_6750Zyk')

N_max_1 = np.load(os.path.join(
    path, 'CT_120_1_20_LS4L_6750Zyk_Creep_n_load_max.npy'))
S_max_1 = np.load(os.path.join(
    path, 'CT_120_1_20_LS4L_6750Zyk_Creep_displacement_max2.npy'))
SS_max_1 = np.load(os.path.join(
    path, 'CT_120_1_20_LS4L_6750Zyk_Creep_displacement_max3.npy'))
SSS_max_1 = np.load(os.path.join(
    path, 'CT_120_1_20_LS4L_6750Zyk_Creep_displacement_max4.npy'))

u_1 = np.load(os.path.join(
    path, 'CT_120_1_20_LS4L_6750Zyk_Displacement2.npy'))
uu_1 = np.load(os.path.join(
    path, 'CT_120_1_20_LS4L_6750Zyk_Displacement3.npy'))
uuu_1 = np.load(os.path.join(
    path, 'CT_120_1_20_LS4L_6750Zyk_Displacement4.npy'))
F_1 = np.load(os.path.join(
    path, 'CT_120_1_20_LS4L_6750Zyk_Force.npy'))

#=============================================
''' PSP 1000 - Test CT_120_1_19'''
#=============================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'CT')
path = os.path.join(path, 'C120')
path = os.path.join(path, 'CT_120_1_21_LS4L_530Zyk')

N_max_2 = np.load(os.path.join(
    path, 'CT_120_1_21_LS4L_530Zyk_Creep_n_load_max.npy'))
S_max_2 = np.load(os.path.join(
    path, 'CT_120_1_21_LS4L_530Zyk_Creep_displacement_max2.npy'))
SS_max_2 = np.load(os.path.join(
    path, 'CT_120_1_21_LS4L_530Zyk_Creep_displacement_max3.npy'))
SSS_max_2 = np.load(os.path.join(
    path, 'CT_120_1_21_LS4L_530Zyk_Creep_displacement_max4.npy'))

u_2 = np.load(os.path.join(
    path, 'CT_120_1_21_LS4L_530Zyk_Displacement2.npy'))
uu_2 = np.load(os.path.join(
    path, 'CT_120_1_21_LS4L_530Zyk_Displacement3.npy'))
uuu_2 = np.load(os.path.join(
    path, 'CT_120_1_21_LS4L_530Zyk_Displacement4.npy'))
F_2 = np.load(os.path.join(
    path, 'CT_120_1_21_LS4L_530Zyk_Force.npy'))


#=============================================
''' PSP 1000 - Test CT_120_1_24'''
#=============================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'CT')
path = os.path.join(path, 'C120')
path = os.path.join(path, 'CT_120_1_24_13504Zyklen')

N_max_3 = np.load(os.path.join(
    path, 'CT_120_1_24_13504Zyklen_Creep_n_load_max.npy'))
S_max_3 = np.load(os.path.join(
    path, 'CT_120_1_24_13504Zyklen_Creep_displacement_max2.npy'))
SS_max_3 = np.load(os.path.join(
    path, 'CT_120_1_24_13504Zyklen_Creep_displacement_max3.npy'))
SSS_max_3 = np.load(os.path.join(
    path, 'CT_120_1_24_13504Zyklen_Creep_displacement_max4.npy'))

u_3 = np.load(os.path.join(
    path, 'CT_120_1_24_13504Zyklen_Displacement2.npy'))
uu_3 = np.load(os.path.join(
    path, 'CT_120_1_24_13504Zyklen_Displacement3.npy'))
uuu_3 = np.load(os.path.join(
    path, 'CT_120_1_24_13504Zyklen_Displacement4.npy'))
F_3 = np.load(os.path.join(
    path, 'CT_120_1_24_13504Zyklen_Force.npy'))


#=============================================
''' PSP 1000 - Test CT_120_1_25'''
#=============================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'CT')
path = os.path.join(path, 'C120')
path = os.path.join(path, 'CT_120_1_25_23256Zyklen')

N_max_4 = np.load(os.path.join(
    path, 'CT_120_1_25_23256Zyklen_Creep_n_load_max.npy'))
S_max_4 = np.load(os.path.join(
    path, 'CT_120_1_25_23256Zyklen_Creep_displacement_max2.npy'))
SS_max_4 = np.load(os.path.join(
    path, 'CT_120_1_25_23256Zyklen_Creep_displacement_max3.npy'))
SSS_max_4 = np.load(os.path.join(
    path, 'CT_120_1_25_23256Zyklen_Creep_displacement_max4.npy'))

u_4 = np.load(os.path.join(
    path, 'CT_120_1_25_23256Zyklen_Displacement2.npy'))
uu_4 = np.load(os.path.join(
    path, 'CT_120_1_25_23256Zyklen_Displacement3.npy'))
uuu_4 = np.load(os.path.join(
    path, 'CT_120_1_25_23256Zyklen_Displacement4.npy'))
F_4 = np.load(os.path.join(
    path, 'CT_120_1_25_23256Zyklen_Force.npy'))


#------------------------------------------------------------------------------
# Plot Force-displacement
#------------------------------------------------------------------------------
plt.subplot(241)
plt.plot(abs(u_1), abs(F_1), 'k')
plt.plot(abs(uu_1), abs(F_1), 'r')
plt.plot(abs(uuu_1), abs(F_1), 'g')
plt.xlim(0, 1.2)
#plt.ylim(1.6, 2.5)
plt.title('stress - strain (H)')
plt.xlabel('Displacement [mm]')
plt.ylabel('Force [kN]')

plt.subplot(242)
plt.plot(abs(u_2), abs(F_2), 'k')
plt.plot(abs(uu_2), abs(F_2), 'r')
plt.plot(abs(uuu_2), abs(F_2), 'g')
plt.xlim(0, 1.2)
#plt.ylim(1.6, 2.5)
plt.title('stress - strain (H)')
plt.xlabel('Displacement [mm]')
plt.ylabel('Force [kN]')


plt.subplot(243)
plt.plot(abs(u_3), abs(F_3), 'k')
plt.plot(abs(uu_3), abs(F_3), 'r')
plt.plot(abs(uuu_3), abs(F_3), 'g')
plt.xlim(0, 1.2)
#plt.ylim(1.6, 2.5)
plt.title('stress - strain (H)')
plt.xlabel('Displacement [mm]')
plt.ylabel('Force [kN]')


plt.subplot(244)
plt.plot(abs(u_4), abs(F_4), 'k')
plt.plot(abs(uu_4), abs(F_4), 'r')
plt.plot(abs(uuu_4), abs(F_4), 'g')
plt.xlim(0, 1.8)
#plt.ylim(1.6, 2.5)
plt.title('stress - strain (H)')
plt.xlabel('Displacement [mm]')
plt.ylabel('Force [kN]')


#------------------------------------------------------------------------------
# Plot fatigue-creep
#------------------------------------------------------------------------------
plt.subplot(245)
plt.plot(N_max_1, abs(S_max_1), 'k')
plt.plot(N_max_1, abs(SS_max_1), 'k')
plt.plot(N_max_1, abs(SSS_max_1), 'k')
plt.title('Fatigue creep curve (H)')
plt.xlabel('N')
plt.ylabel('Displacement [mm]')

plt.subplot(246)
plt.plot(N_max_2, abs(S_max_2), 'r')
plt.plot(N_max_2, abs(SS_max_2), 'r')
plt.plot(N_max_2, abs(SSS_max_2), 'r')
#plt.ylim(0.4, 1.3)
plt.title('Fatigue creep curve (H)')
plt.xlabel('N')
plt.ylabel('Displacement [mm]')

plt.subplot(247)
plt.plot(N_max_3, abs(S_max_3), 'r')
plt.plot(N_max_3, abs(SS_max_3), 'r')
plt.plot(N_max_3, abs(SSS_max_3), 'r')
#plt.ylim(0.4, 1.3)
plt.title('Fatigue creep curve (H)')
plt.xlabel('N')
plt.ylabel('Displacement [mm]')

plt.subplot(248)
plt.plot(N_max_4, abs(S_max_4), 'r')
plt.plot(N_max_4, abs(SS_max_4), 'r')
plt.plot(N_max_4, abs(SSS_max_4), 'r')
#plt.ylim(0.4, 1.3)
plt.title('Fatigue creep curve (H)')
plt.xlabel('N')
plt.ylabel('Displacement [mm]')


plt.show()
