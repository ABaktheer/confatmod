

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd
 
#====================================================================================== 
# Data (H-L) (experiments)
#======================================================================================
r = [0,1,2]
raw_data = {'greenBars': [5, 7, 4], 'orangeBars': [0, 0, 1],'blueBars': [0, 4, 3]}
df = pd.DataFrame(raw_data)
 
# From raw value to percentage
totals = [i+j+k for i,j,k in zip(df['greenBars'], df['orangeBars'], df['blueBars'])]
greenBars = [i / j * 100 for i,j in zip(df['greenBars'], totals)]
orangeBars = [i / j * 100 for i,j in zip(df['orangeBars'], totals)]
blueBars = [i / j * 100 for i,j in zip(df['blueBars'], totals)]
 
# plot
plt.subplot(221)
barWidth = 0.6
names = ('C30-Holmen','C40','C80')
# Create green Bars
plt.bar(r, greenBars, color='red', edgecolor='white', width=barWidth)
# Create orange Bars
plt.bar(r, orangeBars, bottom=greenBars, color='black', edgecolor='white', width=barWidth)
# Create blue Bars
plt.bar(r, blueBars, bottom=[i+j for i,j in zip(greenBars, orangeBars)], color='green', edgecolor='white', width=barWidth)
 
# Custom x axis
plt.xticks(r, names)
plt.ylim(0.0, 100)


#======================================================================================
# Data (L-H) (experiments)
#======================================================================================
r = [0,1,2]
raw_data = {'greenBars': [2, 2, 0], 'orangeBars': [0, 1, 0],'blueBars': [3, 7, 5]}
df = pd.DataFrame(raw_data)
 
# From raw value to percentage
totals = [i+j+k for i,j,k in zip(df['greenBars'], df['orangeBars'], df['blueBars'])]
greenBars = [i / j * 100 for i,j in zip(df['greenBars'], totals)]
orangeBars = [i / j * 100 for i,j in zip(df['orangeBars'], totals)]
blueBars = [i / j * 100 for i,j in zip(df['blueBars'], totals)]
 
# plot
plt.subplot(222)
barWidth = 0.6
names = ('C30-Holmen','C40','C80')
# Create green Bars
plt.bar(r, greenBars, color='red', edgecolor='white', width=barWidth)
# Create orange Bars
plt.bar(r, orangeBars, bottom=greenBars, color='black', edgecolor='white', width=barWidth)
# Create blue Bars
plt.bar(r, blueBars, bottom=[i+j for i,j in zip(greenBars, orangeBars)], color='green', edgecolor='white', width=barWidth)
 
# Custom x axis
plt.xticks(r, names)

plt.ylim(0.0, 100)


#====================================================================================== 
# Data (H-L) (rules)
#======================================================================================
r = [0,1,2]
raw_data = {'greenBars': [0, 1, 1], 'orangeBars': [0, 0, 0],'blueBars': [1, 0, 0]}
df = pd.DataFrame(raw_data)
 
# From raw value to percentage
totals = [i+j+k for i,j,k in zip(df['greenBars'], df['orangeBars'], df['blueBars'])]
greenBars = [i / j * 100 for i,j in zip(df['greenBars'], totals)]
orangeBars = [i / j * 100 for i,j in zip(df['orangeBars'], totals)]
blueBars = [i / j * 100 for i,j in zip(df['blueBars'], totals)]
 
# plot
plt.subplot(223)
barWidth = 0.6
names = ('Mayer','Shah','Baktheer')
# Create green Bars
plt.bar(r, greenBars, color='red', edgecolor='white', width=barWidth)
# Create orange Bars
plt.bar(r, orangeBars, bottom=greenBars, color='black', edgecolor='white', width=barWidth)
# Create blue Bars
plt.bar(r, blueBars, bottom=[i+j for i,j in zip(greenBars, orangeBars)], color='green', edgecolor='white', width=barWidth)
 
# Custom x axis
plt.xticks(r, names)
plt.ylim(0.0, 100)


#======================================================================================
# Data (L-H) (rules)
#======================================================================================
r = [0,1,2]
raw_data = {'greenBars': [1, 1, 0], 'orangeBars': [0, 0, 0],'blueBars': [0, 0, 1]}
df = pd.DataFrame(raw_data)
 
# From raw value to percentage
totals = [i+j+k for i,j,k in zip(df['greenBars'], df['orangeBars'], df['blueBars'])]
greenBars = [i / j * 100 for i,j in zip(df['greenBars'], totals)]
orangeBars = [i / j * 100 for i,j in zip(df['orangeBars'], totals)]
blueBars = [i / j * 100 for i,j in zip(df['blueBars'], totals)]
 
# plot
plt.subplot(224)
barWidth = 0.6
names = ('Mayer','Shah','Baktheer')
# Create green Bars
plt.bar(r, greenBars, color='red', edgecolor='white', width=barWidth)
# Create orange Bars
plt.bar(r, orangeBars, bottom=greenBars, color='black', edgecolor='white', width=barWidth)
# Create blue Bars
plt.bar(r, blueBars, bottom=[i+j for i,j in zip(greenBars, orangeBars)], color='green', edgecolor='white', width=barWidth)
 
# Custom x axis
plt.xticks(r, names)

plt.ylim(0.0, 100)

 
# Show graphic
plt.show()
