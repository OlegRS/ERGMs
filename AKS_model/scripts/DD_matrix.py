#######################################################
# This script plots the degree distribution computed  #
# using matrix inversion against that of the ER model #
# of matching connectance.                            #
#######################################################

import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import matplotlib.pyplot as plt
import numpy as np
from sympy import exp, N, S
from sympy.matrices import Matrix
import scipy.special

fig, axs = plt.subplots(1, 3, figsize=(19./1.5, 6./1.5))

n = 5
L = .3
S = Matrix(np.zeros((n-1, n-1)))
for l in range(1, n):
    for k in range(1, n):
        S[l-1, k-1] = scipy.special.binom(k, l)
MF_S = Matrix(np.zeros(n-1))
for l in range(1, n):
    MF_S[l-1] = scipy.special.binom(n-1, l)*L**l

P = S.inv()*MF_S

P_full = list(P)
P_full.insert(0,1-sum(P))

k = np.arange(0, n)
axs[0].plot(k, P_full, marker='d', label="System solution", linewidth=2, zorder=2, color='blue', alpha=.7)

p_theor = scipy.special.binom(n-1, k)*L**k*(1-L)**(n-1-k)
axs[0].plot(k, p_theor, marker='*', label="ER model", linewidth=2, zorder=2, color='red', alpha=.7)

axs[0].legend(prop={'size': 12})

axs[0].set_xlabel(r'$k$', fontsize=17)
axs[0].set_ylabel(r'$\mathrm{p}(k)$', fontsize=17)
axs[0].set_xlim([0,n-1+.1])
axs[0].set_title(r'$\bar L=$' + str(L) + r'$;\ n=$' + str(n), fontsize=15)

axs[0].axhline(0, color='black', linewidth=.5)

axs[0].grid()
#########################################################
n = 10
L = .3
S = Matrix(np.zeros((n-1, n-1)))
for l in range(1, n):
    for k in range(1, n):
        S[l-1, k-1] = scipy.special.binom(k, l)
MF_S = Matrix(np.zeros(n-1))
for l in range(1, n):
    MF_S[l-1] = scipy.special.binom(n-1, l)*L**l

P = S.inv()*MF_S

P_full = list(P)
P_full.insert(0,1-sum(P))

k = np.arange(0, n)
axs[1].plot(k, P_full, marker='d', label="System solution", linewidth=2, zorder=2, color='blue', alpha=.7)
p_theor = scipy.special.binom(n-1, k)*L**k*(1-L)**(n-1-k)
axs[1].plot(k, p_theor, marker='*', label="ER model", linewidth=2, zorder=2, color='red', alpha=.7)
axs[1].legend(prop={'size': 12})
axs[1].set_xlabel(r'$k$', fontsize=17)
axs[1].set_xlim([0,n-1+.1])
axs[1].axhline(0, color='black', linewidth=.5)
axs[1].grid()

axs[1].set_title(r'$\bar L=$' + str(L) + r'$;\ n=$' + str(n), fontsize=15)
##################################################
n = 40
L = .3
S = Matrix(np.zeros((n-1, n-1)))
for l in range(1, n):
    for k in range(1, n):
        S[l-1, k-1] = scipy.special.binom(k, l)
MF_S = Matrix(np.zeros(n-1))
for l in range(1, n):
    MF_S[l-1] = scipy.special.binom(n-1, l)*L**l

P = S.inv()*MF_S

P_full = list(P)
P_full.insert(0,1-sum(P))

k = np.arange(0, n)
axs[2].plot(k, P_full, marker='d', label="System solution", linewidth=2, zorder=2, color='blue', alpha=.7)

p_theor = scipy.special.binom(n-1, k)*L**k*(1-L)**(n-1-k)
axs[2].plot(k, p_theor, marker='*', label="ER model", linewidth=2, zorder=2, color='red', alpha=.7)
axs[2].legend(prop={'size': 12})
axs[2].set_xlabel(r'$k$', fontsize=17)
axs[2].set_xlim([0,n-1+.1])
axs[2].axhline(0, color='black', linewidth=.5)
axs[2].grid()
axs[2].set_title(r'$\bar L=$' + str(L) + r'$;\ n=$' + str(n), fontsize=15)


plt.tight_layout()
plt.savefig('../figures/DD_matrix.png', dpi=300)
plt.show()
