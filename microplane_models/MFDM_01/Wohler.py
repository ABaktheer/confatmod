'''
Created on 04.07.2019
@author: Abdul
plotting tool for microplane models
'''

import matplotlib
from plotly.validators.carpet.aaxis import LinewidthValidator
from sklearn.linear_model import LinearRegression

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd  # To read data


font = {'family': 'normal',
        'size': 18}

matplotlib.rc('font', **font)


s_min = 0.05
S = np.linspace(0.9, 0.6)

Y = (0.45 + 1.8 * s_min) / (1 + 1.8 * s_min - 0.3 * s_min**2)
cyc = (8 / (Y - 1)) * (S - 1)

x = np.array([1.78, 2.42, 3.1, 3.83, 4.614, 5.49]).reshape(-1, 1)
y = np.array([0.9, 0.85, 0.8, 0.75, 0.7, 0.65]).reshape(-1, 1)

x2 = np.array([3.78, 4.46, 4.32, 4.12, 4.22, 4.02, 5.25,
               4.79, 5.08, 4.25, 5.23, 5.12]).reshape(-1, 1)
y2 = np.array([0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.65,
               0.65, 0.65, 0.65, 0.65, 0.65]).reshape(-1, 1)

f, (ax) = plt.subplots(1, 1, figsize=(5, 4))
ax.plot(cyc, S, color='black', linewidth=3.5)

ax.scatter(x2, y2, color='black', linewidth=4.5)
ax.scatter(x, y, color='red', linewidth=4.5)


ax.set_xlim(1, 7)
ax.set_ylim(0.61, 0.94)

ax.set_xlabel('Log(N)', fontsize=25)
ax.set_ylabel('Smax', fontsize=25)
plt.title('Wohler Curve Smin=0.2')
plt.show()


s_min = 0.
S = np.linspace(0.97, 0.68)
N = np.linspace(0.5, 7)
# cyc = np.linspace(2, 7)
Y = (0.45 + 1.8 * s_min) / (1 + 1.8 * s_min - 0.3 * s_min**2)
cyc = (8 / (Y - 1)) * (S - 1)

f, (ax) = plt.subplots(1, 1, figsize=(5, 4))


# ax.plot(cyc, S, 'g', linewidth=2.5)

X = np.array([1.114, 2.05, 2.16, 2.24, 2.27, 2.39]).reshape(-1, 1)
Y = np.array([0.9, 0.85, 0.8, 0.75, 0.7, 0.65]).reshape(-1, 1)
ax.scatter(X, Y, color='green', linewidth=4.5)

linear_regressor = LinearRegression()  # create object for the class
linear_regressor.fit(X, Y)  # perform linear regression
Y_pred = linear_regressor.predict(X)  # make predictions


m = (Y_pred[2] - Y_pred[3]) / (X[2] - X[3])
Y_extend = m * (N - X[4]) + Y[4]

ax.plot(N, Y_extend, color='green', linewidth=4.5)

plt.show()
s_min = 0.1

# cyc = np.linspace(2, 7)
Y = (0.45 + 1.8 * s_min) / (1 + 1.8 * s_min - 0.3 * s_min**2)
cyc = (8 / (Y - 1)) * (S - 1)

f, (ax) = plt.subplots(1, 1, figsize=(5, 4))


# ax.plot(cyc, S, 'g', linewidth=2.5)

X = np.array([2.52, 3.39, 4.22, 4.95, 5.6]).reshape(-1, 1)
Y = np.array([0.9, 0.85, 0.8, 0.75, 0.7]).reshape(-1, 1)
ax.scatter(X, Y, color='red', linewidth=4.5)

# X = np.array([1.01, 2.02, 2.97, 3.96, 5.01]).reshape(-1, 1)
# Y = np.array([0.95, 0.9, 0.85, 0.8, 0.75]).reshape(-1, 1)
# ax.scatter(X, Y, color='blue', linewidth=4.5)

linear_regressor = LinearRegression()  # create object for the class
linear_regressor.fit(X, Y)  # perform linear regression
Y_pred = linear_regressor.predict(X)  # make predictions


m = (Y_pred[1] - Y_pred[0]) / (X[1] - X[0])
Y_extend = m * (N - X[2]) + Y[2]

ax.plot(N, Y_extend, color='red', linewidth=4.5)


ax.set_xlim(0.5, 7)
ax.set_ylim(0.68, 0.97)

# ax.set_xlabel('Log(N)', fontsize=25)
# ax.set_ylabel('Smax', fontsize=25)
# plt.title('Wohler Curve Smin=0.2')
plt.show()

f, (ax) = plt.subplots(1, 1, figsize=(5, 4))


# ax.plot(cyc, S, 'g', linewidth=2.5)
N = np.linspace(0.5, 8)

X = np.array([2.7686, 4.00599, 5.318]).reshape(-1, 1)
Y = np.array([0.9, 0.85, 0.8]).reshape(-1, 1)
ax.scatter(X, Y, color='green', linewidth=4.5)

# X = np.array([1.01, 2.02, 2.97, 3.96, 5.01]).reshape(-1, 1)
# Y = np.array([0.95, 0.9, 0.85, 0.8, 0.75]).reshape(-1, 1)
# ax.scatter(X, Y, color='blue', linewidth=4.5)

linear_regressor = LinearRegression()  # create object for the class
linear_regressor.fit(X, Y)  # perform linear regression
Y_pred = linear_regressor.predict(X)  # make predictions


m = (Y_pred[1] - Y_pred[0]) / (X[1] - X[0])
Y_extend = m * (N - X[2]) + Y[2]

ax.plot(N, Y_extend, color='red', linewidth=4.5)


ax.set_xlim(2, 8)
ax.set_ylim(0.5, 0.97)

# ax.set_xlabel('Log(N)', fontsize=25)
# ax.set_ylabel('Smax', fontsize=25)
# plt.title('Wohler Curve Smin=0.2')
plt.show()

s_min = 0.2

# cyc = np.linspace(2, 7)
# Y = (0.45 + 1.8 * s_min) / (1 + 1.8 * s_min - 0.3 * s_min**2)
# cyc = (8 / (Y - 1)) * (S - 1)

f, (ax2) = plt.subplots(1, 1, figsize=(5, 4))
# ax.plot(cyc, S, linewidth=3.5)


X = np.array([2.52, 3.39, 4.22, 4.95]).reshape(-1, 1)
Y = np.array([0.9, 0.85, 0.8, 0.75]).reshape(-1, 1)
ax.scatter(X, Y, color='red', linewidth=4.5)

linear_regressor = LinearRegression()  # create object for the class
linear_regressor.fit(X, Y)  # perform linear regression
Y_pred = linear_regressor.predict(X)  # make predictions


m = (Y_pred[1] - Y_pred[0]) / (X[1] - X[0])
Y_extend = m * (N - X[2]) + Y[2]

ax.plot(N, Y_extend, color='red', linewidth=4.5)

plt.plot()
plt.show()

s_min = 0.2

# cyc = np.linspace(2, 7)
# Y = (0.45 + 1.8 * s_min) / (1 + 1.8 * s_min - 0.3 * s_min**2)
# cyc = (8 / (Y - 1)) * (S - 1)

f, (ax) = plt.subplots(1, 1, figsize=(5, 4))
# ax.plot(cyc, S, linewidth=3.5)
N = np.linspace(0.5, 3)

X = np.array([0.69, 1.23, 2.13, 2.36]).reshape(-1, 1)
Y = np.array([0.95, 0.9, 0.85, 0.8]).reshape(-1, 1)
ax.scatter(X, Y, color='red', linewidth=4.5)

linear_regressor = LinearRegression()  # create object for the class
linear_regressor.fit(X, Y)  # perform linear regression
Y_pred = linear_regressor.predict(X)  # make predictions


m = (Y_pred[1] - Y_pred[0]) / (X[1] - X[0])
Y_extend = m * (N - 0.0) + 1.0

ax.plot(N, Y_extend, color='red', linewidth=4.5)

plt.plot()
plt.show()

# s_min = 0.3
#
# X = np.array([1.51, 2.49, 3.43, 4.42]).reshape(-1, 1)
# Y = np.array([0.95, 0.9, 0.85, 0.8]).reshape(-1, 1)
# ax.scatter(X, Y, color='orange', linewidth=4.5)
#
# linear_regressor = LinearRegression()  # create object for the class
# linear_regressor.fit(X, Y)  # perform linear regression
# Y_pred = linear_regressor.predict(X)  # make predictions
#
#
# m = (Y_pred[1] - Y_pred[0]) / (X[1] - X[0])
# Y_extend = m * (N - X[2]) + Y[2]
#
# ax.plot(N, Y_extend, color='orange', linewidth=4.5)
#
# s_min = 0.4
#
# X = np.array([2.49, 3.44, 4.43]).reshape(-1, 1)
# Y = np.array([0.9, 0.85, 0.8]).reshape(-1, 1)
# ax.scatter(X, Y, color='black', linewidth=4.5)
#
# linear_regressor = LinearRegression()  # create object for the class
# linear_regressor.fit(X, Y)  # perform linear regression
# Y_pred = linear_regressor.predict(X)  # make predictions
#
#
# m = (Y_pred[1] - Y_pred[0]) / (X[1] - X[0])
# Y_extend = m * (N - X[2]) + Y[2]
#
# ax.plot(N, Y_extend, color='black', linewidth=4.5)
#
#
# # ax.set_xlim(0, 7)
# #ax.set_ylim(0.68, 0.97)
#
# ax.set_xlabel('Log(N)', fontsize=25)
# ax.set_ylabel('Smax', fontsize=25)
# plt.title('Wohler Curve Smin=0.2')
# plt.show()
