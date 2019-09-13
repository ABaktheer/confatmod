
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
''' Beam-End-Test - C120 - LS1'''
#=========================================================================
#
#=============================================
'''BE_C120_16_01'''
#=============================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'BET')
path = os.path.join(path, 'C120')
path = os.path.join(path, 'BE_C120_16_01')

s_1 = np.load(os.path.join(
    path, 'BE_C120_16_01_Displacement1.npy'))
ss_1 = np.load(os.path.join(
    path, 'BE_C120_16_01_Displacement2.npy'))

F_1 = np.load(os.path.join(
    path, 'BE_C120_16_01_Force.npy'))


#=============================================
'''BE_C120_16_02'''
#=============================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'BET')
path = os.path.join(path, 'C120')
path = os.path.join(path, 'BE_C120_16_02')

s_2 = np.load(os.path.join(
    path, 'BE_C120_16_02_Displacement1.npy'))
ss_2 = np.load(os.path.join(
    path, 'BE_C120_16_02_Displacement2.npy'))

F_2 = np.load(os.path.join(
    path, 'BE_C120_16_02_Force.npy'))


#=============================================
'''BE_C120_16_03'''
#=============================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'BET')
path = os.path.join(path, 'C120')
path = os.path.join(path, 'BE_C120_16_03')

s_3 = np.load(os.path.join(
    path, 'BE_C120_16_03_Displacement1.npy'))
ss_3 = np.load(os.path.join(
    path, 'BE_C120_16_03_Displacement2.npy'))

F_3 = np.load(os.path.join(
    path, 'BE_C120_16_03_Force.npy'))


#=============================================
'''BE_C120_16_07'''
#=============================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'BET')
path = os.path.join(path, 'C120')
path = os.path.join(path, 'BE_C120_16_07')

s_4 = np.load(os.path.join(
    path, 'BE_C120_16_07_Displacement2.npy'))
ss_4 = np.load(os.path.join(
    path, 'BE_C120_16_07_Displacement3.npy'))

F_4 = np.load(os.path.join(
    path, 'BE_C120_16_07_Force.npy'))


#=============================================
'''BE_C120_16_09'''
#=============================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'BET')
path = os.path.join(path, 'C120')
path = os.path.join(path, 'BE_C120_16_09')

s_5 = np.load(os.path.join(
    path, 'BE_C120_16_09_Displacement2.npy'))
ss_5 = np.load(os.path.join(
    path, 'BE_C120_16_09_Displacement3.npy'))

F_5 = np.load(os.path.join(
    path, 'BE_C120_16_09_Force.npy'))


#------------------------------------------------------------------------------
# Plot Force-displacement
#------------------------------------------------------------------------------
plt.subplot(231)
plt.plot(abs(s_1), abs(F_1), 'k')
plt.plot(abs(ss_1), abs(F_1), 'k')
#plt.xlim(0, 1.2)
#plt.ylim(1.6, 2.5)
plt.title('stress - strain (H)')
plt.xlabel('Displacement [mm]')
plt.ylabel('Force [kN]')

plt.subplot(232)
plt.plot(abs(s_2), abs(F_2), 'k')
plt.plot(abs(ss_2), abs(F_2), 'k')
#plt.xlim(0, 1.2)
#plt.ylim(1.6, 2.5)
plt.title('stress - strain (H)')
plt.xlabel('Displacement [mm]')
plt.ylabel('Force [kN]')

plt.subplot(233)
plt.plot(abs(s_3), abs(F_3), 'k')
plt.plot(abs(ss_3), abs(F_3), 'k')
#plt.xlim(0, 1.2)
#plt.ylim(1.6, 2.5)
plt.title('stress - strain (H)')
plt.xlabel('Displacement [mm]')
plt.ylabel('Force [kN]')

plt.subplot(234)
plt.plot(abs(s_4), abs(F_4), 'k')
plt.plot(abs(ss_4), abs(F_4), 'k')
#plt.xlim(0, 1.2)
#plt.ylim(1.6, 2.5)
plt.title('stress - strain (H)')
plt.xlabel('Displacement [mm]')
plt.ylabel('Force [kN]')


plt.subplot(235)
plt.plot(abs(s_5), abs(F_5), 'k')
plt.plot(abs(ss_5), abs(F_5), 'k')
#plt.xlim(0, 1.2)
#plt.ylim(1.6, 2.5)
plt.title('stress - strain (H)')
plt.xlabel('Displacement [mm]')
plt.ylabel('Force [kN]')

plt.show()
