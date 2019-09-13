
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
''' C120 - LS2'''
#=========================================================================
#
#=============================================
''' PSP 1000 - Test CT_120_1_16'''
#=============================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'CT')
path = os.path.join(path, 'C120')
path = os.path.join(path, 'CT_120_1_16_LS2')

N_max_1 = np.load(os.path.join(
    path, 'CT_120_1_16_LS2_Creep_n_load_max.npy'))
S_max_1 = np.load(os.path.join(
    path, 'CT_120_1_16_LS2_Creep_displacement_max2.npy'))
SS_max_1 = np.load(os.path.join(
    path, 'CT_120_1_16_LS2_Creep_displacement_max3.npy'))
SSS_max_1 = np.load(os.path.join(
    path, 'CT_120_1_16_LS2_Creep_displacement_max4.npy'))

u_1 = np.load(os.path.join(
    path, 'CT_120_1_16_LS2_Displacement2.npy'))
uu_1 = np.load(os.path.join(
    path, 'CT_120_1_16_LS2_Displacement3.npy'))
uuu_1 = np.load(os.path.join(
    path, 'CT_120_1_16_LS2_Displacement4.npy'))
F_1 = np.load(os.path.join(
    path, 'CT_120_1_16_LS2_Force.npy'))

#=============================================
''' PSP 1000 - Test CT_120_1_17'''
#=============================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'CT')
path = os.path.join(path, 'C120')
path = os.path.join(path, 'CT_120_1_17_LS2')

N_max_2 = np.load(os.path.join(
    path, 'CT_120_1_17_LS2_Creep_n_load_max.npy'))
S_max_2 = np.load(os.path.join(
    path, 'CT_120_1_17_LS2_Creep_displacement_max2.npy'))
SS_max_2 = np.load(os.path.join(
    path, 'CT_120_1_17_LS2_Creep_displacement_max3.npy'))
SSS_max_2 = np.load(os.path.join(
    path, 'CT_120_1_17_LS2_Creep_displacement_max4.npy'))

u_2 = np.load(os.path.join(
    path, 'CT_120_1_17_LS2_Displacement2.npy'))
uu_2 = np.load(os.path.join(
    path, 'CT_120_1_17_LS2_Displacement3.npy'))
uuu_2 = np.load(os.path.join(
    path, 'CT_120_1_17_LS2_Displacement4.npy'))
F_2 = np.load(os.path.join(
    path, 'CT_120_1_17_LS2_Force.npy'))


#------------------------------------------------------------------------------
# Plot Force-displacement
#------------------------------------------------------------------------------
plt.subplot(221)
plt.plot(abs(u_1), abs(F_1), 'k')
plt.plot(abs(uu_1), abs(F_1), 'r')
plt.plot(abs(uuu_1), abs(F_1), 'g')
plt.xlim(0, 1.2)
#plt.ylim(1.6, 2.5)
plt.title('stress - strain (H)')
plt.xlabel('Displacement [mm]')
plt.ylabel('Force [kN]')

plt.subplot(222)
plt.plot(abs(u_2), abs(F_2), 'k')
plt.plot(abs(uu_2), abs(F_2), 'r')
plt.plot(abs(uuu_2), abs(F_2), 'g')
plt.xlim(0, 1.2)
#plt.ylim(1.6, 2.5)
plt.title('stress - strain (H)')
plt.xlabel('Displacement [mm]')
plt.ylabel('Force [kN]')

#------------------------------------------------------------------------------
# Plot fatigue-creep
#------------------------------------------------------------------------------
plt.subplot(223)
plt.plot(N_max_1, abs(S_max_1), 'k')
plt.plot(N_max_1, abs(SS_max_1), 'k')
plt.plot(N_max_1, abs(SSS_max_1), 'k')
plt.title('Fatigue creep curve (H)')
plt.xlabel('N')
plt.ylabel('Displacement [mm]')

plt.subplot(224)
plt.plot(N_max_2, abs(S_max_2), 'r')
plt.plot(N_max_2, abs(SS_max_2), 'r')
plt.plot(N_max_2, abs(SSS_max_2), 'r')
#plt.ylim(0.4, 1.3)
plt.title('Fatigue creep curve (H)')
plt.xlabel('N')
plt.ylabel('Displacement [mm]')


plt.show()
