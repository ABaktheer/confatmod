'''
Created on 13.01.2020

@author: abaktheer
'''



 
import matplotlib.pyplot as plt
import numpy as np
  
m = np.linspace(0, 1, 6)
f = np.zeros_like(m)
tau_bar = 2.5
  
for i in np.arange(0, len(m)):  
#tau = np.linspace(-10, 10, 100)
    print (i)
    sig = np.linspace(-30, tau_bar/m[i], 100)
    print (sig)
  
    f = tau_bar - m[i] * sig
    
    plt.subplot(111) 
    plt.plot(sig , f, 'k', linewidth=2.0, alpha=1.0)
    plt.plot(sig , -f, 'k', linewidth=2, alpha=1.0)
    plt.axhline(y=0, color='k', linewidth=0.5, alpha=0.5)
    plt.axvline(x=0, color='k', linewidth=0.5, alpha=0.5)


plt.show()
 
tau_bar = np.linspace(0, 10, 5)
f = np.zeros_like(m)
m = 1.0 
  
for i in np.arange(0, len(tau_bar)):  
#tau = np.linspace(-10, 10, 100)
    print (i)
    sig = np.linspace(-30, tau_bar[i]/m, 100)
    print (sig)
  
    f = tau_bar[i] - m * sig
    
    plt.subplot(111) 
    plt.plot(sig , f, 'k', linewidth=2.0, alpha=1.0)
    plt.plot(sig , -f, 'k', linewidth=2, alpha=1.0)
    plt.axhline(y=0, color='k', linewidth=0.5, alpha=0.5)
    plt.axvline(x=0, color='k', linewidth=0.5, alpha=0.5) 
  
plt.show()

 
 
#  
# from mpl_toolkits.mplot3d import Axes3D
#  
# import matplotlib.pyplot as plt
# from matplotlib import cm
# from matplotlib.ticker import LinearLocator, FormatStrFormatter
# import numpy as np
#  
#  
# fig = plt.figure()
# ax = fig.gca(projection='3d')
#  
# # Make data.
# m = 1.0
# tau_bar = 5. 
# tau = np.linspace(-10, 10, 100)
# sig = np.linspace(-5, tau_bar/m, 100)
# X, Y = np.meshgrid(sig, sig)
#  
# f =  np.abs(X) - tau_bar + m * Y 
# Z = f 
#  
# # Plot the surface.
# surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)
#  
# # Customize the z axis.
# #ax.set_zlim(-10.0, 0.0)
# ax.zaxis.set_major_locator(LinearLocator(10))
# ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#  
# # Add a color bar which maps values to colors.
# fig.colorbar(surf, shrink=0.5, aspect=5)
# 
# plt.show()


# import matplotlib.pyplot as plt
# import numpy as np
# 
# fig, ax = plt.subplots()
# x1 = np.arange(0,10)
# x2 = np.arange(1,11)
# y = np.arange(20,30)
# c = np.arange(5)+range(5)
# cmap = np.vectorize(lambda x: {0: 'red', 1: 'blue', 2: 'green', 3: 'cyan', 4: 'magenta', 5: 'yellow'}.get(x))
# 
# # marker plot with different color
# # for i in range(len(y)):
# #     ax.plot(x1[i], y[i], "x", color=cmap(c)[i])
# 
# # horizotal lines
# ax.hlines(y, x1, x2, color=cmap(c))
# 
# plt.show()
# 
