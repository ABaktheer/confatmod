
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
''' Beam-End-Test - C120 - LS2'''
#=========================================================================
#
#=============================================
'''BE_C120_16_04'''
#=============================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'BET')
path = os.path.join(path, 'C120')
path = os.path.join(path, 'BE_C120_16_04')

N_max_1 = np.load(os.path.join(
    path, 'BE_C120_16_04_Creep_n_load_max.npy'))
S_max_1 = np.load(os.path.join(
    path, 'BE_C120_16_04_Creep_displacement_max2.npy'))
SS_max_1 = np.load(os.path.join(
    path, 'BE_C120_16_04_Creep_displacement_max3.npy'))


u_1 = np.load(os.path.join(
    path, 'BE_C120_16_04_Displacement2.npy'))
uu_1 = np.load(os.path.join(
    path, 'BE_C120_16_04_Displacement3.npy'))

F_1 = np.load(os.path.join(
    path, 'BE_C120_16_04_Force.npy'))

#=============================================
'''BE_C120_16_05'''
#=============================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'BET')
path = os.path.join(path, 'C120')
path = os.path.join(path, 'BE_C120_16_05')

N_max_2 = np.load(os.path.join(
    path, 'BE_C120_16_05_Creep_n_load_max.npy'))
S_max_2 = np.load(os.path.join(
    path, 'BE_C120_16_05_Creep_displacement_max2.npy'))
SS_max_2 = np.load(os.path.join(
    path, 'BE_C120_16_05_Creep_displacement_max3.npy'))


u_2 = np.load(os.path.join(
    path, 'BE_C120_16_05_Displacement2.npy'))
uu_2 = np.load(os.path.join(
    path, 'BE_C120_16_05_Displacement3.npy'))

F_2 = np.load(os.path.join(
    path, 'BE_C120_16_05_Force.npy'))


#=============================================
'''BE_C120_16_06'''
#=============================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'BET')
path = os.path.join(path, 'C120')
path = os.path.join(path, 'BE_C120_16_06')

N_max_3 = np.load(os.path.join(
    path, 'BE_C120_16_06_Creep_n_load_max.npy'))
S_max_3 = np.load(os.path.join(
    path, 'BE_C120_16_06_Creep_displacement_max2.npy'))
SS_max_3 = np.load(os.path.join(
    path, 'BE_C120_16_06_Creep_displacement_max3.npy'))


u_3 = np.load(os.path.join(
    path, 'BE_C120_16_06_Displacement2.npy'))
uu_3 = np.load(os.path.join(
    path, 'BE_C120_16_06_Displacement3.npy'))

F_3 = np.load(os.path.join(
    path, 'BE_C120_16_06_Force.npy'))


#------------------------------------------------------------------------------
# Plot Force-displacement
#------------------------------------------------------------------------------
plt.subplot(231)
plt.plot(abs(u_1), abs(F_1), 'k')
plt.plot(abs(uu_1), abs(F_1), 'r')
plt.xlim(0, 1.2)
#plt.ylim(1.6, 2.5)
plt.title('stress - strain (H)')
plt.xlabel('Displacement [mm]')
plt.ylabel('Force [kN]')


plt.subplot(232)
plt.plot(abs(u_2), abs(F_2), 'k')
plt.plot(abs(uu_2), abs(F_2), 'r')
#plt.xlim(0, 1.2)
#plt.ylim(1.6, 2.5)
plt.title('stress - strain (H)')
plt.xlabel('Displacement [mm]')
plt.ylabel('Force [kN]')


plt.subplot(233)
plt.plot(abs(u_3), abs(F_3), 'k')
plt.plot(abs(uu_3), abs(F_3), 'r')
#plt.xlim(0, 1.2)
#plt.ylim(1.6, 2.5)
plt.title('stress - strain (H)')
plt.xlabel('Displacement [mm]')
plt.ylabel('Force [kN]')

#------------------------------------------------------------------------------
# Plot fatigue-creep
#------------------------------------------------------------------------------
plt.subplot(234)
plt.plot(N_max_1, abs(S_max_1), 'k')
plt.plot(N_max_1, abs(SS_max_1), 'r')
plt.title('Fatigue creep curve (H)')
plt.xlabel('N')
plt.ylabel('Displacement [mm]')


plt.subplot(235)
plt.plot(N_max_2, abs(S_max_2), 'k')
plt.plot(N_max_2, abs(SS_max_2), 'r')
#plt.ylim(0.4, 1.3)
plt.title('Fatigue creep curve (H)')
plt.xlabel('N')
plt.ylabel('Displacement [mm]')


plt.subplot(236)
plt.plot(N_max_3, abs(S_max_3), 'k')
plt.plot(N_max_3, abs(SS_max_3), 'r')
#plt.ylim(0.4, 1.3)
plt.title('Fatigue creep curve (H)')
plt.xlabel('N')
plt.ylabel('Displacement [mm]')


plt.show()
