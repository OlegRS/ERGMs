###############################################
# This script plots self-consistency equation #
# and exact solutions for MF4 model.          #
###############################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

t2_MF = 15
t3_MF = -20
t4_MF = 9

title = "$t_2=$" + str(t2_MF) + "$,\ t_3=$" + str(t3_MF) + "$,\ t_4=$" + str(t4_MF) + "$;\ n=\{2,5,10\}$"
fig, ax = plt.subplots(figsize=(10., 6.*1.1))

T=1
x_min = -5.1
x_max = -2.5

# Plotting self-consistency equation
delta_t1 = 0.05
delta_L = 0.000005
t1_range = np.arange(x_min, x_max+delta_t1, delta_t1)
L_range = np.arange(0, 1, delta_L)
L, T1 = np.meshgrid(L_range, t1_range)
SC = L - 1/2*(1 + np.tanh(T1 + 2*t2_MF*L + 3*t3_MF*L**2 + 4*t4_MF*L**3)) #Self-consistency
CS = plt.contour(T1, L, SC, [0], colors='cyan', zorder=1, linewidths=1.5)
CS.collections[0].set_label('Equation of state')

plt.axvline(0, color='black')
plt.axhline(0, color='black')

plt.ylim([-.01, 1.01])
plt.xlim([x_min, x_max])

plt.title(title, fontsize=22)
plt.xlabel(r'$t_1$', fontsize=22)
plt.ylabel(r'$\bar L$', fontsize=22)

# Finding ensemble averages in thermodynamic limit
n_nodes = 2
N_links_max = int(n_nodes*(n_nodes-1)/2)

x_min_ = x_min
x_max_ = x_max

T1_MF = np.arange(x_min_, x_max_, 0.0001)
L_MIN = []
L_EXPECTED = []
T1_MF_ = []
L = np.arange(0, N_links_max+1)/N_links_max
S = np.zeros(N_links_max+1)
for i in range(0, N_links_max+1): # Computing log of degeneracy factor
    for j in range(0, N_links_max-1):
        S[i] += np.log(N_links_max-j)
    for j in range(0,i):
        S[i] -= np.log((i-j))
    for j in range(0, N_links_max-i):
        S[i] -= np.log(N_links_max-i-j)
S /= N_links_max

for t1_MF in T1_MF:
    F = -2*(t1_MF*L + t2_MF*L**2 + t3_MF*L**3 + t4_MF*L**4) - T*S
    numerator, denominator = 0, 0
    for i in range(0, N_links_max+1):
        numerator += L[i]*np.exp(-N_links_max/T*F[i])
        denominator += np.exp(-N_links_max/T*F[i])
    L_EXPECTED.append(numerator/denominator)

plt.plot(T1_MF, L_EXPECTED, color='orange', label='Exact solution for n=2', zorder=3)
#######################################################################################
n_nodes = 5
N_links_max = int(n_nodes*(n_nodes-1)/2)

x_min_ = x_min
x_max_ = x_max

T1_MF = np.arange(x_min_, x_max_, 0.0001)
L_MIN = []
L_EXPECTED = []
T1_MF_ = []
L = np.arange(0, N_links_max+1)/N_links_max
S = np.zeros(N_links_max+1)
for i in range(0, N_links_max+1): # Computing log of degeneracy factor
    for j in range(0, N_links_max-1):
        S[i] += np.log(N_links_max-j)
    for j in range(0,i):
        S[i] -= np.log((i-j))
    for j in range(0, N_links_max-i):
        S[i] -= np.log(N_links_max-i-j)
S /= N_links_max

for t1_MF in T1_MF:
    F = -2*(t1_MF*L + t2_MF*L**2 + t3_MF*L**3 + t4_MF*L**4) - T*S
    numerator, denominator = 0, 0
    for i in range(0, N_links_max+1):
        numerator += L[i]*np.exp(-N_links_max/T*F[i])
        denominator += np.exp(-N_links_max/T*F[i])
    L_EXPECTED.append(numerator/denominator)

plt.plot(T1_MF, L_EXPECTED, color='green', label='Exact solution for n=5', zorder=4)
#############################################################################
n_nodes = 10
N_links_max = int(n_nodes*(n_nodes-1)/2)

x_min_ = x_min
x_max_ = x_max

T1_MF = np.arange(x_min_, x_max_, 0.0001)
L_MIN = []
L_EXPECTED = []
T1_MF_ = []
L = np.arange(0, N_links_max+1)/N_links_max
S = np.zeros(N_links_max+1)
for i in range(0, N_links_max+1): # Computing log of degeneracy factor
    for j in range(0, N_links_max-1):
        S[i] += np.log(N_links_max-j)
    for j in range(0,i):
        S[i] -= np.log((i-j))
    for j in range(0, N_links_max-i):
        S[i] -= np.log(N_links_max-i-j)
S /= N_links_max

for t1_MF in T1_MF:
    F = -2*(t1_MF*L + t2_MF*L**2 + t3_MF*L**3 + t4_MF*L**4) - T*S
    numerator, denominator = 0, 0
    for i in range(0, N_links_max+1):
        numerator += L[i]*np.exp(-N_links_max/T*F[i])
        denominator += np.exp(-N_links_max/T*F[i])
    L_EXPECTED.append(numerator/denominator)

plt.plot(T1_MF, L_EXPECTED, color='blue', label='Exact solution for n=10', zorder=5)
#############################################################################
x_min_ = x_min 
x_max_ = -2.5

T1_MF = np.arange(x_min_, x_max_, .0001)
L_MIN = []
L_EXPECTED = []
T1_MF_ = []
for t1_MF in T1_MF:
    def SC_equation(L):
        return L - 1/2*(1 + np.tanh(t1_MF + 2*t2_MF*L + 3*t3_MF*L**2 + 4*t4_MF*L**3))
    L_ = np.array([fsolve(SC_equation, -0.01)[0], fsolve(SC_equation, 0)[0], fsolve(SC_equation, .5)[0], fsolve(SC_equation, .999)[0], fsolve(SC_equation, 1.001)[0], fsolve(SC_equation, .99)[0], fsolve(SC_equation, 1)[0]])
    def F(L):
        return -2*(t1_MF*L + t2_MF*L**2 + t3_MF*L**3 + t4_MF*L**4) + L*np.log(L/(1-L)) + np.log(1-L)
    L0 = L_[~(np.triu(np.abs(L_[:,None] - L_) < 0.001,1)).any(0)] # Distinct solution of SC equation
    L_ = np.array([l for l in L0 if (l<1 and l>0)])        
    L_min = L_[list(F(L_)).index(min(F(L_)))]
    L_MIN.append(L_min)

plt.plot(T1_MF, L_MIN, color='red', label='Thermodynamic solution', zorder=2)

plt.grid()

plt.legend(loc="upper left", prop={'size': 16}).set_zorder(200)

plt.tight_layout()

plt.savefig('../figures/MF4_SC_with_jumps.png', dpi=300)

plt.show()
fig.clear()
