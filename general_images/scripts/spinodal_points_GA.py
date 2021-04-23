###########################################################
# This script shows the appearance of spinodal solutions  #
# of the self-consistency equation in graphical analysis. #
###########################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(23/2, 10.5/2))

degree = 4-1

def A(L):
    return 1/2*np.log(L/(1-L))

x_artanh = np.arange(-1,1,.00001) #.00001
y_artanh = A(x_artanh)

x = np.array([0, 0.35, 0.7, 0.9, 1])
y = np.array([-3.5, 0, 0, 0, .5])

# p = np.polyfit(x, y, degree)
p = np.array([9*4, -20*3, 15*2, -5.0139352])

axs[0].set_xlim([-.02, 1.015])
axs[0].set_ylim([-5.2, 2.5])
title = r"$t_1=$" + str(round(t[0]/2,4)) + r"$,\ t_2=$" + str(round(t[1]/(2*2),4)) + r"$,\ t_3=$" + str(round(t[2]/(2*3),4)) + r"$,\ t_4=$" + str(round(t[3]/(2*4),4))
axs[0].set_xlabel(r'$\bar L$', fontsize=16)
axs[0].set_ylabel(r'${\cal A}(\bar L; T),\ {\cal P}(\bar L;\bf t)$', fontsize=16)
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
def SC(L):
    return P(L)-A(L)

axs[0].plot(X, P(X), color='green', zorder=11, linewidth=2.5, label=r'${\cal P}(\bar L;\bf t)$')
axs[0].plot(x_artanh, y_artanh, color='blue', zorder=10, linewidth=2.5, label=r'${\cal A}(\bar L; T)$')
L_stab = fsolve(SC, 0)[0]
axs[0].plot(L_stab, P(L_stab), marker='^', linestyle='none', color='red', markersize=8, label=r'${\cal A}={\cal P}$' + ' stable regular solution', zorder=11)
axs[0].plot(.32779475, P(.32779475), marker='D', linestyle='none', color='orange', markersize=8, label=r'${\cal A}={\cal P}$' + ' unstable singular solution', zorder=11)


axs[0].legend(loc='lower right', fontsize=14)
axs[0].set_title(title, fontsize=14)
axs[0].grid()
####################################################################################
x = np.array([0, 0.35, 0.7, 0.9, 1])
y = np.array([-3.5, 0, 0, 0, .5])

# p = np.polyfit(x, y, degree)
p = np.array([9*4, -20*3, 15*2, -3.76886804])

axs[1].set_xlim([-.02, 1.015])
axs[1].set_ylim([-5.1, 2.5])
title = r"$t_1=$" + str(round(t[0],4)) + r"$,\ t_2=$" + str(round(t[1],4)) + r"$,\ t_3=$" + str(round(t[2],4)) + r"$,\ t_4=$" + str(round(t[3],4))
axs[1].set_xlabel(r'$\bar L$', fontsize=16)
axs[1].axhline(0, color='black')
axs[1].axvline(0, color='black')
X = np.arange(0, 1, .001)

t = []
for i in range(0,4):
    t.append(p[4-i-1]/(i+1))

axs[1].plot(X, P(X), color='green', zorder=11, linewidth=2.5, label=r'${\cal P}(\bar L;\bf t)$')
axs[1].plot(x_artanh, y_artanh, color='blue', zorder=10, linewidth=2.5, label=r'${\cal A}(\bar L; T)$')

L_stab = np.array([fsolve(SC, 0)[0], fsolve(SC, .6)[0]])
axs[1].plot(L_stab, P(L_stab), marker='^', linestyle='none', color='red', markersize=8, label=r'${\cal A}={\cal P}$' + ' stable regular solution', zorder=11)
L_unstab = fsolve(SC, .1)[0]
axs[1].plot(L_unstab, P(L_unstab), marker='v', linestyle='none', color='yellow', markersize=8, label=r'${\cal A}={\cal P}$' + ' unstable regular solution', zorder=11)
axs[1].plot(0.9649008, P(0.9649008), marker='D', linestyle='none', color='orange', markersize=8, label=r'${\cal A}={\cal P}$' + ' unstable singular solution', zorder=11)

axs[1].legend(loc='lower right', fontsize=14)

axs[1].legend(loc='lower right', fontsize=14)
axs[1].set_title(title, fontsize=14)
axs[1].grid()
plt.tight_layout()
plt.savefig('../figures/spinodal_points_GA.png', dpi=300)
plt.show()
plt.close()
