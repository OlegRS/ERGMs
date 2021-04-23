#######################################
# This script demonstrates graphical  # 
# analysis for models of even orders. #
#######################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(23/2, 10.5/2))

degree = 4-1

x_artanh = [0, 0, 1, 1]
y_artanh = [-100, 0, 0, 100]

x = np.array([0, 0.35, 0.7, 0.9, 1])
y = np.array([-3.5, 0, 0, 0, .5])

#p = np.polyfit(x, y, degree)
p = np.array([9*4, -20*3, 15*2, -4.3])/6

axs[0].set_xlim([-.02, 1.015])
axs[0].set_ylim([-.8, .6])
title = r"$ t_1=$" + str(round(p[3]/2, 4)) + r"$,\  t_2=$" + str(round(p[2]/(2*2), 4)) + r"$,\  t_3=$" + str(round(p[1]/(2*3), 4)) + r"$,\  t_4=$" + str(round(p[0]/(2*4), 4))
axs[0].set_xlabel(r'$\bar L$', fontsize=15)
axs[0].set_ylabel(r'${\cal A}(\bar L),\ {\cal P}(\bar L)$', fontsize=16)
axs[0].axhline(0, color='black')
axs[0].axvline(0, color='black')
X = np.arange(0, 1, .001)

t = []
for i in range(0,4):
    t.append(p[4-i-1]/(i+1))

def P(L):
    s=0
    for i in range(degree+1):
        s += p[i]*L**(degree-i)
    return s

axs[0].plot(X, P(X), color='green', zorder=11, linewidth=2.5, label=r'${\cal P}(\bar L)$')
axs[0].plot(x_artanh, y_artanh, color='blue', zorder=10, linewidth=2.5, label=r'${\cal A}(\bar L)$' + ' when ' + r'$T\to 0$')

x = np.array([0, fsolve(P,.2)[0], fsolve(P,.6)[0], fsolve(P,.8)[0], 1])
y=P(x)
axs[0].plot(x[0:5:2], y[0:5:2], marker='^', linestyle='none', color='red', markersize=8, label=r'${\cal A}(\bar L)={\cal P}(\bar L)$' + ' stable solutions', zorder=11)
axs[0].plot(x[1:5:2], y[1:5:2], marker='v', linestyle='none', color='yellow', markersize=8, label=r'${\cal A}(\bar L)={\cal P}(\bar L)$' + ' unstable solutions', zorder=11)

axs[0].legend(loc='lower right', fontsize=14)
axs[0].set_title(title, fontsize=13, y=1.01)
axs[0].grid()
####################################################################################
degree = 6-1

x = np.array([0, 0.05, 0.33333333, 0.5, 0.83333333, 1])
y = np.array([-.7, 0, 0, 0, 0, .5])

p = np.polyfit(x, y, degree)

axs[1].set_xlim([-.015, 1.015])
axs[1].set_ylim([-.8, .6])
title = r"$ t_1=$" + str(round(p[5]/2, 1)) + r"$,\  t_2=$" + str(round(p[4]/(2*2), 1)) + r"$,\  t_3=$" + str(round(p[3]/(2*3), 1)) + r"$,\  t_4=$" + str(round(p[2]/(2*4), 1)) + r"$,\  t_5=$" + str(round(p[1]/(2*5), 1)) + r"$,\  t_6=$" + str(round(p[0]/(2*6), 1))
axs[1].set_xlabel(r'$\bar L$', fontsize=16)
axs[1].set_label(r'${\cal A}(\bar L),\ {\cal P}(\bar L)$')
axs[1].axhline(0, color='black')
axs[1].axvline(0, color='black')
X = np.arange(0, 1, .001)

t = []
for i in range(0,6):
    t.append(p[6-i-1]/(i+1))

def P(L):
    s=0
    for i in range(degree+1):
        s += p[i]*L**(degree-i)
    return s
axs[1].plot(X, P(X), color='green', zorder=11, linewidth=2.5, label=r'${\cal P}(\bar L)$')
axs[1].plot(x_artanh, y_artanh, color='blue', zorder=10, linewidth=2.5, label=r'${\cal A}(\bar L)$' + ' when ' + r'$T\to 0$')

x6 = 1
y6 = P(1)
x[0] = 0
y[0] = P(0)
x[5] = fsolve(P, 9)[0]
y[5] = 0
axs[1].plot(np.append(x,[x6])[0:7:2], np.append(y,[y6])[0:7:2], marker='^', linestyle='none', color='red', markersize=8, label=r'${\cal A}(\bar L)={\cal P}(\bar L)$' + ' stable solutions', zorder=11)
axs[1].plot(x[1:7:2], y[1:7:2], marker='v', linestyle='none', color='yellow', markersize=8, label=r'${\cal A}(\bar L)={\cal P}(\bar L)$' + ' unstable solutions', zorder=11)

axs[1].legend(loc='lower right', fontsize=14)
axs[1].set_title(title, fontsize=13, x=.472, y=1.01)
axs[1].grid()
plt.tight_layout()
plt.savefig('../figures/even_order_GA.png', dpi=300)
plt.show()
plt.close()
