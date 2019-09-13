
''''
Created on 20.10.2018

@author: Abedulgader Baktheer
'''

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from io import StringIO
    import os


#=========================================================================
''' loading level 0.85 '''
#=========================================================================
#=============================================
''' Constant High (S_max=0.85) (CT_80-9) '''
#=============================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'CT')
path = os.path.join(path, 'C80')
path = os.path.join(path, 'CT_80-9-0001969Zyk_g')

N_max_1 = np.load(os.path.join(
    path, 'CT_80-9-0001969Zyk_g_Creep_n_load_max.npy')) * 1969
S_max_1 = np.load(os.path.join(
    path, 'CT_80-9-0001969Zyk_g_Creep_displacement_max2.npy'))
SS_max_1 = np.load(os.path.join(
    path, 'CT_80-9-0001969Zyk_g_Creep_displacement_max3.npy'))
SSS_max_1 = np.load(os.path.join(
    path, 'CT_80-9-0001969Zyk_g_Creep_displacement_max4.npy'))

u_1 = np.load(os.path.join(
    path, 'CT_80-9-0001969Zyk_g_Displacement2.npy'))
uu_1 = np.load(os.path.join(
    path, 'CT_80-9-0001969Zyk_g_Displacement3.npy'))
uuu_1 = np.load(os.path.join(
    path, 'CT_80-9-0001969Zyk_g_Displacement4.npy'))
F_1 = np.load(os.path.join(
    path, 'CT_80-9-0001969Zyk_g_Force.npy'))


#==========================================
''' Constant High (S_max=0.85) (CT_80-28) '''
#==========================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'CT')
path = os.path.join(path, 'C80')
path = os.path.join(path, 'CT_80-28-0-000045zykl_g')

N_max_2 = np.load(os.path.join(
    path, 'CT_80-28-0-000045zykl_g_Creep_n_load_max.npy')) * 45
S_max_2 = np.load(os.path.join(
    path, 'CT_80-28-0-000045zykl_g_Creep_displacement_max2.npy'))
SS_max_2 = np.load(os.path.join(
    path, 'CT_80-28-0-000045zykl_g_Creep_displacement_max3.npy'))
SSS_max_2 = np.load(os.path.join(
    path, 'CT_80-28-0-000045zykl_g_Creep_displacement_max4.npy'))
u_2 = np.load(os.path.join(
    path, 'CT_80-28-0-000045zykl_g_Displacement1.npy'))
F_2 = np.load(os.path.join(
    path, 'CT_80-28-0-000045zykl_g_Force.npy'))

#==========================================
''' Constant High (S_max=0.85) (CT_80-29) '''
#==========================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'CT')
path = os.path.join(path, 'C80')
path = os.path.join(path, 'CT_80-29-0-000916zykl_g')

N_max_3 = np.load(os.path.join(
    path, 'CT_80-29-0-000916zykl_g_Creep_n_load_max.npy')) * 916
S_max_3 = np.load(os.path.join(
    path, 'CT_80-29-0-000916zykl_g_Creep_displacement_max2.npy'))
SS_max_3 = np.load(os.path.join(
    path, 'CT_80-29-0-000916zykl_g_Creep_displacement_max3.npy'))
SSS_max_3 = np.load(os.path.join(
    path, 'CT_80-29-0-000916zykl_g_Creep_displacement_max4.npy'))

u_3 = np.load(os.path.join(
    path, 'CT_80-29-0-000916zykl_g_Displacement2.npy'))
uu_3 = np.load(os.path.join(
    path, 'CT_80-29-0-000916zykl_g_Displacement3.npy'))
uuu_3 = np.load(os.path.join(
    path, 'CT_80-29-0-000916zykl_g_Displacement4.npy'))
F_3 = np.load(os.path.join(
    path, 'CT_80-29-0-000916zykl_g_Force.npy'))


plt.subplot(223)
plt.plot(abs(u_1), abs(F_1), 'k')
plt.plot(abs(uu_1), abs(F_1), 'r')
plt.plot(abs(uuu_1), abs(F_1), 'g')
# plt.plot(abs(u_3), abs(F_3), 'k')
# plt.plot(abs(uu_3), abs(F_3), 'r')
# plt.plot(abs(uuu_3), abs(F_3), 'g')
plt.xlim(0, 1.2)
#plt.ylim(1.6, 2.5)
plt.title('stress - strain (H)')
plt.xlabel('Displacement [mm]')
plt.ylabel('Force [kN]')


plt.subplot(224)
# plt.plot(abs(u_1), abs(F_1), 'k')
# plt.plot(abs(uu_1), abs(F_1), 'r')
# plt.plot(abs(uuu_1), abs(F_1), 'g')
plt.plot(abs(u_3), abs(F_3), 'k')
plt.plot(abs(uu_3), abs(F_3), 'r')
plt.plot(abs(uuu_3), abs(F_3), 'g')
plt.xlim(0, 1.2)
#plt.ylim(1.6, 2.5)
plt.title('stress - strain (H)')
plt.xlabel('Displacement [mm]')
plt.ylabel('Force [kN]')

plt.subplot(221)
plt.plot(N_max_1, abs(S_max_1), 'k')
plt.plot(N_max_1, abs(SS_max_1), 'k')
plt.plot(N_max_1, abs(SSS_max_1), 'k')
plt.plot(N_max_2, abs(S_max_2), 'r')
plt.plot(N_max_2, abs(SS_max_2), 'r')
plt.plot(N_max_2, abs(SSS_max_2), 'r')
plt.plot(N_max_3, abs(S_max_3), 'g')
plt.plot(N_max_3, abs(SS_max_3), 'g')
plt.plot(N_max_3, abs(SSS_max_3), 'g')
plt.ylim(0.4, 1.3)
plt.title('Fatigue creep curve (H)')
plt.xlabel('N')
plt.ylabel('Displacement [mm]')


# #=========================================================================
# ''' loading level 0.80 '''
# #=========================================================================
# #=============================================
# ''' Constant High (S_max=0.80) (CT_80-12) '''
# #=============================================
# home_dir = os.path.expanduser('~')
# path = os.path.join(home_dir, 'Data Processing')
# path = os.path.join(path, 'CT_80-12-0032173Zyk_g')
#
# N_max_11 = np.load(os.path.join(
#     path, 'CT_80-12-0032173Zyk_g_Creep_n_load_max.npy'))  # * 32173
# S_max_11 = np.load(os.path.join(
#     path, 'CT_80-12-0032173Zyk_g_Creep_displacement_max2.npy'))
# SS_max_11 = np.load(os.path.join(
#     path, 'CT_80-12-0032173Zyk_g_Creep_displacement_max3.npy'))
# SSS_max_11 = np.load(os.path.join(
#     path, 'CT_80-12-0032173Zyk_g_Creep_displacement_max4.npy'))
#
# u_11 = np.load(os.path.join(
#     path, 'CT_80-12-0032173Zyk_g_Displacement1.npy'))
# F_11 = np.load(os.path.join(
#     path, 'CT_80-12-0032173Zyk_g_Force.npy'))
#
# mpl.rcParams['agg.path.chunksize'] = 10000
#
# #==========================================
# ''' Constant High (S_max=0.80) (CT_80-13) '''
# #==========================================
# home_dir = os.path.expanduser('~')
# path = os.path.join(home_dir, 'Data Processing')
# path = os.path.join(path, 'CT_80-13-0000276Zyk_g')
#
# N_max_22 = np.load(os.path.join(
#     path, 'CT_80-13-0000276Zyk_g_Creep_n_load_max.npy'))  # * 276
# S_max_22 = np.load(os.path.join(
#     path, 'CT_80-13-0000276Zyk_g_Creep_displacement_max2.npy'))
# SS_max_22 = np.load(os.path.join(
#     path, 'CT_80-13-0000276Zyk_g_Creep_displacement_max3.npy'))
# SSS_max_22 = np.load(os.path.join(
#     path, 'CT_80-13-0000276Zyk_g_Creep_displacement_max4.npy'))
#
# u_22 = np.load(os.path.join(
#     path, 'CT_80-13-0000276Zyk_g_Displacement1.npy'))
# F_22 = np.load(os.path.join(
#     path, 'CT_80-13-0000276Zyk_g_Force.npy'))
#
# mpl.rcParams['agg.path.chunksize'] = 10000
#
#
# #==========================================
# ''' Constant High (S_max=0.80) (CT_80-14) '''
# #==========================================
# home_dir = os.path.expanduser('~')
# path = os.path.join(home_dir, 'Data Processing')
# path = os.path.join(path, 'CT_80-14-00011339Zyk_g')
#
# N_max_33 = np.load(os.path.join(
#     path, 'CT_80-14-00011339Zyk_g_Creep_n_load_max.npy'))  # * 11339
# S_max_33 = np.load(os.path.join(
#     path, 'CT_80-14-00011339Zyk_g_Creep_displacement_max2.npy'))
# SS_max_33 = np.load(os.path.join(
#     path, 'CT_80-14-00011339Zyk_g_Creep_displacement_max3.npy'))
# SSS_max_33 = np.load(os.path.join(
#     path, 'CT_80-14-00011339Zyk_g_Creep_displacement_max4.npy'))
#
# u_33 = np.load(os.path.join(
#     path, 'CT_80-14-00011339Zyk_g_Displacement1.npy'))
# F_33 = np.load(os.path.join(
#     path, 'CT_80-14-00011339Zyk_g_Force.npy'))
#
# mpl.rcParams['agg.path.chunksize'] = 10000


#=========================================================================
''' loading level 0.75 '''
#=========================================================================
#
#=============================================
''' Constant low (S_max=0.75) (CT_80-15) '''
#=============================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'CT')
path = os.path.join(path, 'C80')
path = os.path.join(path, 'CT_80-15-0163931Zyk_g')

N_max_111 = np.load(os.path.join(
    path, 'CT_80-15-0163931Zyk_g_Creep_n_load_max.npy')) * 163931
S_max_111 = np.load(os.path.join(
    path, 'CT_80-15-0163931Zyk_g_Creep_displacement_max2.npy'))
SS_max_111 = np.load(os.path.join(
    path, 'CT_80-15-0163931Zyk_g_Creep_displacement_max3.npy'))
SSS_max_111 = np.load(os.path.join(
    path, 'CT_80-15-0163931Zyk_g_Creep_displacement_max4.npy'))

u_111 = np.load(os.path.join(
    path, 'CT_80-15-0163931Zyk_g_Displacement1.npy'))
F_111 = np.load(os.path.join(
    path, 'CT_80-15-0163931Zyk_g_Force.npy'))

mpl.rcParams['agg.path.chunksize'] = 10000

#==========================================
''' Constant low (S_max=0.75) (CT_80-16) '''
#==========================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'CT')
path = os.path.join(path, 'C80')
path = os.path.join(path, 'CT_80-16-00084205Zyk_g')

N_max_222 = np.load(os.path.join(
    path, 'CT_80-16-00084205Zyk_g_Creep_n_load_max.npy')) * 84205
S_max_222 = np.load(os.path.join(
    path, 'CT_80-16-00084205Zyk_g_Creep_displacement_max2.npy'))
SS_max_222 = np.load(os.path.join(
    path, 'CT_80-16-00084205Zyk_g_Creep_displacement_max3.npy'))
SSS_max_222 = np.load(os.path.join(
    path, 'CT_80-16-00084205Zyk_g_Creep_displacement_max4.npy'))

u_222 = np.load(os.path.join(
    path, 'CT_80-16-00084205Zyk_g_Displacement1.npy'))
F_222 = np.load(os.path.join(
    path, 'CT_80-16-00084205Zyk_g_Force.npy'))

mpl.rcParams['agg.path.chunksize'] = 10000


#==========================================
''' Constant low (S_max=0.75) (CT_80-17) '''
#==========================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'CT')
path = os.path.join(path, 'C80')
path = os.path.join(path, 'CT_80-17-0234316Zyk_g')

N_max_333 = np.load(os.path.join(
    path, 'CT_80-17-0234316Zyk_g_Creep_n_load_max.npy')) * 234316
S_max_333 = np.load(os.path.join(
    path, 'CT_80-17-0234316Zyk_g_Creep_displacement_max2.npy'))
SS_max_333 = np.load(os.path.join(
    path, 'CT_80-17-0234316Zyk_g_Creep_displacement_max3.npy'))
SSS_max_333 = np.load(os.path.join(
    path, 'CT_80-17-0234316Zyk_g_Creep_displacement_max4.npy'))

u_333 = np.load(os.path.join(
    path, 'CT_80-17-0234316Zyk_g_Displacement1.npy'))
F_333 = np.load(os.path.join(
    path, 'CT_80-17-0234316Zyk_g_Force.npy'))


#==========================================
''' Constant low (S_max=0.75) (CT_80-26) '''
#==========================================
home_dir = os.path.expanduser('~')
path = os.path.join(home_dir, 'Data Processing')
path = os.path.join(path, 'CT')
path = os.path.join(path, 'C80')
path = os.path.join(path, 'CT_80-26-0-1085zyk_g')

N_max_444 = np.load(os.path.join(
    path, 'CT_80-26-0-1085zyk_g_Creep_n_load_max.npy')) * 1085
S_max_444 = np.load(os.path.join(
    path, 'CT_80-26-0-1085zyk_g_Creep_displacement_max2.npy'))
SS_max_444 = np.load(os.path.join(
    path, 'CT_80-26-0-1085zyk_g_Creep_displacement_max3.npy'))
SSS_max_444 = np.load(os.path.join(
    path, 'CT_80-26-0-1085zyk_g_Creep_displacement_max4.npy'))

u_444 = np.load(os.path.join(
    path, 'CT_80-26-0-1085zyk_g_Displacement1.npy'))
F_444 = np.load(os.path.join(
    path, 'CT_80-26-0-1085zyk_g_Force.npy'))

mpl.rcParams['agg.path.chunksize'] = 10000

# plt.subplot(223)
# plt.plot(abs(u_111), abs(F_111), 'k')
# plt.plot(abs(u_222), abs(F_222), 'r')
# plt.plot(abs(u_333), abs(F_333), 'g')
# plt.plot(abs(u_444), abs(F_444), 'b')
# plt.xlim(0, 2.5)
# #plt.ylim(1.6, 2.5)
# plt.title('stress - strain (L)')
# plt.xlabel('Displacement [mm]')
# plt.ylabel('Force [kN]')

plt.subplot(222)
plt.plot(N_max_111, abs(S_max_111), 'k')
plt.plot(N_max_111, abs(SS_max_111), 'k')
plt.plot(N_max_111, abs(SSS_max_111), 'k')
plt.plot(N_max_222, abs(S_max_222), 'r')
plt.plot(N_max_222, abs(SS_max_222), 'r')
plt.plot(N_max_222, abs(SSS_max_222), 'r')
plt.plot(N_max_333, abs(S_max_333), 'g')
plt.plot(N_max_333, abs(SS_max_333), 'g')
plt.plot(N_max_333, abs(SSS_max_333), 'g')
plt.plot(N_max_444, abs(S_max_444), 'b')
plt.plot(N_max_444, abs(SS_max_444), 'b')
plt.plot(N_max_444, abs(SSS_max_444), 'b')
plt.ylim(0.4, 1.3)
plt.title('Fatigue creep curve (L)')
plt.xlabel('N')
plt.ylabel('Displacement [mm]')


plt.show()
