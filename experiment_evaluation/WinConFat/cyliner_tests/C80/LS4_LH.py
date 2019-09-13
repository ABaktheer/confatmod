'''
Created on 11.03.2019

@author: abaktheer
'''

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from io import StringIO
    import os

#=============================================
''' L-H (S= 0.85 -0.75) (CT_80-18) '''
#=============================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'CT')
path = os.path.join(path, 'C80')
path = os.path.join(path, 'CT_80-18-0036719Zyk_g')

N_max_1 = np.load(os.path.join(
    path, 'CT_80-18-0036719Zyk_g_Creep_n_load_max.npy')) * 36719
S_max_1 = np.load(os.path.join(
    path, 'CT_80-18-0036719Zyk_g_Creep_displacement_max2.npy'))
SS_max_1 = np.load(os.path.join(
    path, 'CT_80-18-0036719Zyk_g_Creep_displacement_max3.npy'))
SSS_max_1 = np.load(os.path.join(
    path, 'CT_80-18-0036719Zyk_g_Creep_displacement_max4.npy'))

u_1 = np.load(os.path.join(
    path, 'CT_80-18-0036719Zyk_g_Displacement1.npy'))
# F_1 = np.load(os.path.join(
#     path, 'CT_80-18-0036719Zyk_g_Force.npy'))

F_max_1 = np.load(os.path.join(
    path, 'CT_80-18-0036719Zyk_g_Creep_n_load_max.npy'))

mpl.rcParams['agg.path.chunksize'] = 10000

#==========================================
''' Constant High (S=0.85) (CT_80-19) '''
#==========================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'CT')
path = os.path.join(path, 'C80')
path = os.path.join(path, 'CT_80-19-0033983Zyk_g')

N_max_2 = np.load(os.path.join(
    path, 'CT_80-19-0033983Zyk_g_Creep_n_load_max.npy')) * 33983
S_max_2 = np.load(os.path.join(
    path, 'CT_80-19-0033983Zyk_g_Creep_displacement_max2.npy'))
SS_max_2 = np.load(os.path.join(
    path, 'CT_80-19-0033983Zyk_g_Creep_displacement_max3.npy'))
SSS_max_2 = np.load(os.path.join(
    path, 'CT_80-19-0033983Zyk_g_Creep_displacement_max4.npy'))

u_2 = np.load(os.path.join(
    path, 'CT_80-19-0033983Zyk_g_Displacement1.npy'))
# F_2 = np.load(os.path.join(
#     path, 'CT_80-19-0033983Zyk_g_Force.npy'))

mpl.rcParams['agg.path.chunksize'] = 10000

#==========================================
''' Constant High (S=0.85) (CT_80-20) '''
#==========================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'CT')
path = os.path.join(path, 'C80')
path = os.path.join(path, 'CT_80-20-0033013Zyk_g')

N_max_3 = np.load(os.path.join(
    path, 'CT_80-20-0033013Zyk_g_Creep_n_load_max.npy')) * 33013
S_max_3 = np.load(os.path.join(
    path, 'CT_80-20-0033013Zyk_g_Creep_displacement_max2.npy'))
SS_max_3 = np.load(os.path.join(
    path, 'CT_80-20-0033013Zyk_g_Creep_displacement_max3.npy'))
SSS_max_3 = np.load(os.path.join(
    path, 'CT_80-20-0033013Zyk_g_Creep_displacement_max4.npy'))

u_3 = np.load(os.path.join(
    path, 'CT_80-20-0033013Zyk_g_Displacement1.npy'))
# F_3 = np.load(os.path.join(
#     path, 'CT_80-20-0033013Zyk_g_Force.npy'))

mpl.rcParams['agg.path.chunksize'] = 10000


#==========================================
''' Constant High (S=0.85) (CT_80-24) '''
#==========================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'CT')
path = os.path.join(path, 'C80')
path = os.path.join(path, 'CT_80-24-0-65705Zyk_g')

N_max_4 = np.load(os.path.join(
    path, 'CT_80-24-0-65705Zyk_g_Creep_n_load_max.npy')) * 65705
S_max_4 = np.load(os.path.join(
    path, 'CT_80-24-0-65705Zyk_g_Creep_displacement_max2.npy'))
SS_max_4 = np.load(os.path.join(
    path, 'CT_80-24-0-65705Zyk_g_Creep_displacement_max3.npy'))
SSS_max_4 = np.load(os.path.join(
    path, 'CT_80-24-0-65705Zyk_g_Creep_displacement_max4.npy'))

u_4 = np.load(os.path.join(
    path, 'CT_80-24-0-65705Zyk_g_Displacement1.npy'))
# F_4 = np.load(os.path.join(
#     path, 'CT_80-24-0-65705Zyk_g_Force.npy'))

mpl.rcParams['agg.path.chunksize'] = 10000


#==========================================
''' Constant High (S=0.85) (CT_80-25) '''
#==========================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'CT')
path = os.path.join(path, 'C80')
path = os.path.join(path, 'CT_80-25-0-65097yk_g')

N_max_5 = np.load(os.path.join(
    path, 'CT_80-25-0-65097yk_g_Creep_n_load_max.npy')) * 65097
S_max_5 = np.load(os.path.join(
    path, 'CT_80-25-0-65097yk_g_Creep_displacement_max2.npy'))
SS_max_5 = np.load(os.path.join(
    path, 'CT_80-25-0-65097yk_g_Creep_displacement_max3.npy'))
SSS_max_5 = np.load(os.path.join(
    path, 'CT_80-25-0-65097yk_g_Creep_displacement_max4.npy'))

u_5 = np.load(os.path.join(
    path, 'CT_80-25-0-65097yk_g_Displacement1.npy'))
# F_5 = np.load(os.path.join(
#     path, 'CT_80-25-0-65097yk_g_Force.npy'))


#==============================================================
# Normalized fatigue creep curves
#==============================================================


n_01 = 32000
n_02 = 32000
n_03 = 32000
n_04 = 64000
n_05 = 64000

n_hf = 1443.0
n_lf = 160817.0

# normalized accumulative damage
N_01 = np.hstack(
    (N_max_1[0:n_01 + 1] / n_lf, n_01 / n_lf + (N_max_1[n_01 + 1:] - N_max_1[n_01]) / n_hf))
N_02 = np.hstack(
    (N_max_2[0:n_02 + 1] / n_lf, n_02 / n_lf + (N_max_2[n_02 + 1:] - N_max_2[n_02]) / n_hf))
N_03 = np.hstack(
    (N_max_3[0:n_03 + 1] / n_lf, n_03 / n_lf + (N_max_3[n_03 + 1:] - N_max_3[n_03]) / n_hf))
N_04 = np.hstack(
    (N_max_4[0:n_04 + 1] / n_lf, n_04 / n_lf + (N_max_4[n_04 + 1:] - N_max_4[n_04]) / n_hf))
N_05 = np.hstack(
    (N_max_5[0:n_05 + 1] / n_lf, n_05 / n_lf + (N_max_5[n_05 + 1:] - N_max_5[n_05]) / n_hf))


#========================================
# Plotting
#========================================

# plt.subplot(121)
# plt.plot(abs(S_max_1), abs(F_max_1), 'k')
#
# #plt.xlim(0.0, 3.0)
# plt.xlabel('Displacement [mm]')
# plt.ylabel('Force [kN]')


plt.subplot(221)
plt.plot(N_max_1, abs(S_max_1), 'k')
plt.plot(N_max_1, abs(SS_max_1), 'k')
plt.plot(N_max_1, abs(SSS_max_1), 'k')
plt.plot(N_max_2, abs(S_max_2), 'g')
plt.plot(N_max_2, abs(SS_max_2), 'g')
plt.plot(N_max_2, abs(SSS_max_2), 'g')
plt.plot(N_max_3, abs(S_max_3), 'b')
plt.plot(N_max_3, abs(SS_max_3), 'b')
plt.plot(N_max_3, abs(SSS_max_3), 'b')
plt.plot(N_max_4, abs(S_max_4), 'm')
plt.plot(N_max_4, abs(SS_max_4), 'm')
plt.plot(N_max_4, abs(SSS_max_4), 'm')
plt.plot(N_max_5, abs(S_max_5), 'r')
plt.plot(N_max_5, abs(SS_max_5), 'r')
plt.plot(N_max_5, abs(SSS_max_5), 'r')
plt.ylim(0.5, 1.5)
plt.title('Fatigue creep curve (L-H)')
plt.xlabel('N')
plt.ylabel('Displacement [mm]')

plt.subplot(222)
plt.plot(N_01[1:], abs(S_max_1[1:]), "k")
plt.plot(N_01[1:], abs(SS_max_1[1:]), "k")
plt.plot(N_01[1:], abs(SSS_max_1[1:]), "k")
plt.plot(N_02[1:], abs(S_max_2[1:]), "g")
plt.plot(N_02[1:], abs(SS_max_2[1:]), "g")
plt.plot(N_02[1:], abs(SSS_max_2[1:]), "g")
plt.plot(N_03[1:], abs(S_max_3[1:]), "b")
plt.plot(N_03[1:], abs(SS_max_3[1:]), "b")
plt.plot(N_03[1:], abs(SSS_max_3[1:]), "b")
plt.plot(N_04[1:], abs(S_max_4[1:]), "m")
plt.plot(N_04[1:], abs(SS_max_4[1:]), "m")
plt.plot(N_04[1:], abs(SSS_max_4[1:]), "m")
plt.plot(N_05[1:], abs(S_max_5[1:]), "r")
plt.plot(N_05[1:], abs(SS_max_5[1:]), "r")
plt.plot(N_05[1:], abs(SSS_max_5[1:]), "r")
plt.ylim(0.5, 1.5)
plt.title('Fatigue creep curve normalized (L-H)')
plt.xlabel('N/Nf')
plt.ylabel('Displacement [mm]')


plt.show()
