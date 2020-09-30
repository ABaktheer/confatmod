
import matplotlib.pyplot as plt
import numpy as np
import random as r
from scipy.stats import norm
import seaborn as sns
from scipy.stats import pearsonr 


#=========================
# experimental (monotonic)
#=========================
etas = np.loadtxt(
    r'H:\Heimarbeit\results_random\1500 samples_3\etas.txt')
DS = np.loadtxt(
    r'H:\Heimarbeit\results_random\1500 samples_3\DS.txt')
abs_DS = np.loadtxt(
    r'H:\Heimarbeit\results_random\1500 samples_3\abs_DS.txt')
Sm = np.loadtxt(
    r'H:\Heimarbeit\results_random\1500 samples_3\Sm.txt')


Smax_mean = np.loadtxt(r'H:\Heimarbeit\results_random\1500 samples_3\Smax_mean.txt')
eta_mean = np.loadtxt(r'H:\Heimarbeit\results_random\1500 samples_3\eta_mean.txt')

#================================================
# plot
#================================================

plt.subplot(331)
x = np.arange(1, len(etas) +1)
plt.plot(x[1:], etas[1:], 'ro', markersize=3, color='k', alpha=0.5)
plt.plot(np.array([0, len(etas)]), np.array([np.mean(etas[1:]), np.mean(etas[1:])]), color='r', linewidth=2)
plt.plot(np.array([0, len(etas)]), np.array([np.mean(etas[1:]) + np.std(etas[1:]), np.mean(etas[1:]) + np.std(etas[1:])]), color='g')
plt.plot(np.array([0, len(etas)]), np.array([np.mean(etas[1:]) - np.std(etas[1:]), np.mean(etas[1:]) - np.std(etas[1:])]), color='g')
plt.fill_between(np.array([0, len(etas)]), np.array([np.mean(etas[1:]) - np.std(etas[1:]), np.mean(etas[1:]) - np.std(etas[1:])]) ,
                  np.array([np.mean(etas[1:]) + np.std(etas[1:]), np.mean(etas[1:]) + np.std(etas[1:])]), facecolor='green', alpha=0.4)
                                                                                                     
plt.plot(np.array([0, len(etas)]), np.array([1, 1]), color='k', linewidth=2)

plt.xlim(0.0, len(etas) + 1)
plt.ylim(0.0, 2)
plt.xlabel('simulation')
plt.ylabel('$\eta$')
plt.title('Simulations')

#================================================
plt.subplot(332)
ax = sns.distplot(etas[1:],
              bins=50,
              color='red',
              hist_kws={"linewidth": 2,'alpha':0.5},
              norm_hist = True)


plt.plot( np.array([np.mean(etas[1:]), np.mean(etas[1:])]), np.array([0, 3]), color='r', linewidth=2)

plt.plot( np.array([1, 1]), np.array([0, 3]), color='k', linewidth=2)

ax.set(xlabel='Normal Distribution', ylabel='Frequency')
plt.title('Normal Distribution')

#================================================
plt.subplot(333)

#kwargs = {'cumulative': True, 'density': True}
#sns.distplot(etas[1:], hist_kws=kwargs, kde_kws=kwargs)
plt.hist(etas[1:], cumulative=True, density=True, bins=100,  alpha=0.3)


# evaluate the histogram
values, base = np.histogram(etas[1:], bins=100)
#evaluate the cumulative
cumulative = np.cumsum(values) /np.cumsum(values)[-1]
# plot the cumulative function
plt.plot(base[:-1], cumulative, c='red' , linewidth=2)
plt.plot( np.array([1, 1]), np.array([0, 1]), color='k', linewidth=2)
plt.title('Cumulative Distribution')
    
    
#================================================    
plt.subplot(334)
plt.plot(DS[1:], etas[1:], 'ro', markersize=3, color='k')
#plt.plot(np.array([0, len(etas)]), np.array([np.mean(etas[1:]), np.mean(etas[1:])]), color='r')
 
m, b = np. polyfit(DS[1:], etas[1:], 1) 
plt.plot(DS[1:], m*DS[1:] + b , color='r') 
 
 
plt.xlabel('Delta S')
plt.ylabel('$\eta$')
 
corr, _ = pearsonr(DS[1:], etas[1:])
print('Pearsons correlation DS: %.3f' % corr)




#================================================
plt.subplot(335)
#plt.plot(abs_DS[1:], etas[1:], 'ro', markersize=3, color='k')
plt.plot(Smax_mean[1:], etas[1:], 'ro', markersize=3, color='k', alpha=0.5)


m, b = np. polyfit(Smax_mean[1:], etas[1:], 1) 
plt.plot(Smax_mean[1:], m*Smax_mean[1:] + b , color='r') 

plt.xlabel('Smax_mean')
plt.ylabel('$\eta$')

corr, _ = pearsonr(Smax_mean[1:], etas[1:])
print('Pearsons correlation Smax : %.3f' % corr)

#================================================
plt.subplot(336)
plt.plot(Sm[1:], etas[1:], 'ro', markersize=3, color='k', alpha=0.5)
#plt.plot(np.array([0, len(etas)]), np.array([np.mean(etas[1:]), np.mean(etas[1:])]), color='r')

#ax = sns.regplot(Sm[1:], etas[1:], marker='o', color='black', scatter_kws={'s':3})

m, b = np. polyfit(Sm[1:], etas[1:], 1) 
plt.plot(Sm[1:], m*Sm[1:] + b , color='r') 

plt.xlabel('Sm')
plt.ylabel('$\eta$')


#================================================
plt.subplot(337)
plt.plot(eta_mean[1:], etas[1:], 'ro', markersize=3, color='k', alpha=0.5)
#plt.plot(np.array([0, len(etas)]), np.array([np.mean(etas[1:]), np.mean(etas[1:])]), color='r')


m, b = np. polyfit(eta_mean[1:], etas[1:], 1) 
plt.plot(eta_mean[1:], m*eta_mean[1:] + b , color='r') 

plt.xlabel('eta_mean')
plt.ylabel('$\eta$')


#================================================
plt.subplot(338)
plt.plot(eta_mean[1:], Smax_mean[1:], 'ro', markersize=3, color='k', alpha=0.5)
#plt.plot(np.array([0, len(etas)]), np.array([np.mean(etas[1:]), np.mean(etas[1:])]), color='r')


m, b = np. polyfit(eta_mean[1:], Smax_mean[1:], 1) 
plt.plot(eta_mean[1:], m*eta_mean[1:] + b , color='r') 

corr, _ = pearsonr(Smax_mean[1:],eta_mean[1:])
print('Pearsons correlation Smax_mean-eta_mean : %.3f' % corr)

plt.xlabel('eta_mean')
plt.ylabel('$Smax_mean$')


#================================================
plt.subplot(339)
plt.plot(Sm[1:], Smax_mean[1:], 'ro', markersize=3, color='k', alpha=0.5)
#plt.plot(np.array([0, len(etas)]), np.array([np.mean(etas[1:]), np.mean(etas[1:])]), color='r')


m, b = np. polyfit(Sm[1:], Smax_mean[1:], 1) 
plt.plot(Sm[1:], m*Sm[1:] + b , color='r') 

corr, _ = pearsonr(Sm[1:], Smax_mean[1:])
print('Pearsons correlation Sm-Smax_mean : %.3f' % corr)

plt.xlabel('Sm')
plt.ylabel('$Smax_mean$')

plt.show()
