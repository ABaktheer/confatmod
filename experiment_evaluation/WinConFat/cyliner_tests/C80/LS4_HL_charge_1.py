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
''' H-L (S= 0.85 -0.75) (CT_80-30) '''
#=============================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'CT')
path = os.path.join(path, 'C80')
path = os.path.join(path, 'CT_80-30-0-0196959zykl_g')

N_max_1 = np.load(os.path.join(
    path, 'CT_80-30-0-0196959zykl_g_Creep_n_load_max.npy')) * 196959
S_max_1 = np.load(os.path.join(
    path, 'CT_80-30-0-0196959zykl_g_Creep_displacement_max2.npy'))
SS_max_1 = np.load(os.path.join(
    path, 'CT_80-30-0-0196959zykl_g_Creep_displacement_max3.npy'))
SSS_max_1 = np.load(os.path.join(
    path, 'CT_80-30-0-0196959zykl_g_Creep_displacement_max4.npy'))

u_1 = np.load(os.path.join(
    path, 'CT_80-30-0-0196959zykl_g_Displacement1.npy'))
F_1 = np.load(os.path.join(
    path, 'CT_80-30-0-0196959zykl_g_Force.npy'))


#==========================================
''' Constant High (S=0.85) (CT_80-31) '''
#==========================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'CT')
path = os.path.join(path, 'C80')
path = os.path.join(path, 'CT_80-31-0-0093991zykl_g')

N_max_2 = np.load(os.path.join(
    path, 'CT_80-31-0-0093991zykl_g_Creep_n_load_max.npy')) * 93991
S_max_2 = np.load(os.path.join(
    path, 'CT_80-31-0-0093991zykl_g_Creep_displacement_max2.npy'))
SS_max_2 = np.load(os.path.join(
    path, 'CT_80-31-0-0093991zykl_g_Creep_displacement_max3.npy'))
SSS_max_2 = np.load(os.path.join(
    path, 'CT_80-31-0-0093991zykl_g_Creep_displacement_max4.npy'))

u_2 = np.load(os.path.join(
    path, 'CT_80-31-0-0093991zykl_g_Displacement1.npy'))
# F_2 = np.load(os.path.join(
#     path, 'CT_80-31-0-0093991zykl_g_Force.npy'))


#==========================================
''' Constant High (S=0.85) (CT_80-32) '''
#==========================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'CT')
path = os.path.join(path, 'C80')
path = os.path.join(path, 'CT_80-32-0-00289376zykl_g')

N_max_3 = np.load(os.path.join(
    path, 'CT_80-32-0-00289376zykl_g_Creep_n_load_max.npy')) * 289376

S_max_3 = np.load(os.path.join(
    path, 'CT_80-32-0-00289376zykl_g_Creep_displacement_max2.npy'))
SS_max_3 = np.load(os.path.join(
    path, 'CT_80-32-0-00289376zykl_g_Creep_displacement_max3.npy'))
SSS_max_3 = np.load(os.path.join(
    path, 'CT_80-32-0-00289376zykl_g_Creep_displacement_max4.npy'))

u_3 = np.load(os.path.join(
    path, 'CT_80-32-0-00289376zykl_g_Displacement1.npy'))
F_3 = np.load(os.path.join(
    path, 'CT_80-32-0-00289376zykl_g_Force.npy'))


#==========================================
''' Constant High (S=0.85) (CT_80-27) '''
#==========================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'CT')
path = os.path.join(path, 'C80')
path = os.path.join(path, 'CT_80-27-0-194030zyk_g')

N_max_4 = np.load(os.path.join(
    path, 'CT_80-27-0-194030zyk_g_Creep_n_load_max.npy')) * 194030
S_max_4 = np.load(os.path.join(
    path, 'CT_80-27-0-194030zyk_g_Creep_displacement_max2.npy'))
SS_max_4 = np.load(os.path.join(
    path, 'CT_80-27-0-194030zyk_g_Creep_displacement_max3.npy'))
SSS_max_4 = np.load(os.path.join(
    path, 'CT_80-27-0-194030zyk_g_Creep_displacement_max4.npy'))

u_4 = np.load(os.path.join(
    path, 'CT_80-27-0-194030zyk_g_Displacement1.npy'))
# F_4 = np.load(os.path.join(
#     path, 'CT_80-27-0-194030zyk_g_Force.npy'))


#==============================================================
# Normalized fatigue creep curves
#==============================================================


n_01 = 1000
n_02 = 1000
n_03 = 1000
n_04 = 3000

n_hf = 1443.0
n_lf = 160817.0

# normalized accumulative damage
N_01 = np.hstack(
    (N_max_1[0:n_01 + 1] / n_hf, n_01 / n_hf + (N_max_1[n_01 + 1:] - N_max_1[n_01]) / n_lf))
N_02 = np.hstack(
    (N_max_3[0:n_02 + 1] / n_hf, n_02 / n_hf + (N_max_2[n_02 + 1:] - N_max_2[n_02]) / n_lf))
N_03 = np.hstack(
    (N_max_3[0:n_03 + 1] / n_hf, n_03 / n_hf + (N_max_3[n_03 + 1:] - N_max_3[n_03]) / n_lf))
N_04 = np.hstack(
    (N_max_4[0:n_04 + 1] / n_hf, n_04 / n_hf + (N_max_4[n_04 + 1:] - N_max_4[n_04]) / n_lf))

#========================================
# Plotting
#========================================

# plt.subplot(221)
# #plt.plot(abs(u_1), abs(F_1), 'k')
# #plt.plot(abs(u_2), abs(F_2), 'r')
# #plt.plot(abs(u_3), abs(F_3), 'g')
# #plt.plot(abs(u_4), abs(F_4), 'b')
# plt.xlim(0.0, 2.5)
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
plt.plot(N_max_4, abs(S_max_4), 'r')
plt.plot(N_max_4, abs(SS_max_4), 'r')
plt.plot(N_max_4, abs(SSS_max_4), 'r')
plt.ylim(0.5, 1.5)
plt.title('Fatigue creep curve (H-L)')
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
plt.plot(N_04[1:], abs(S_max_4[1:]), "r")
plt.plot(N_04[1:], abs(SS_max_4[1:]), "r")
plt.plot(N_04[1:], abs(SSS_max_4[1:]), "r")
plt.ylim(0.5, 1.5)
plt.title('Fatigue creep curve normalized (H-L)')
plt.xlabel('N/Nf')
plt.ylabel('Displacement [mm]')


plt.show()
