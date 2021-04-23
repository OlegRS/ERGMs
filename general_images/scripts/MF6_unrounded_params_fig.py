###############################################
# This script plots self-consistency equation #
# and exact solutions for MF6 model.          #
###############################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

degree = 6-1
x = np.array([0, 0.05, 0.33333333, 0.5, 0.83333333, 1])
y = np.array([-.7, 0, 0, 0, 0, .5])

plt.plot(-y, x, marker='*')
p = np.polyfit(x, y, 5)

plt.ylim([-.1, 1.1])
plt.xlim([-2.5, 2.5])
plt.axhline(0, color='black')
plt.axvline(0, color='black')
X = np.arange(0, 1, .001)

t= np.array([-7, 95.53, -384.21, 692.97, -580.21, 183.79])
print("t=", t)

Y = np.zeros(len(X))
for i in range(degree+1):
    Y += p[i]*X**(degree-i)
    
plt.plot(-Y, X)
plt.grid()

from scipy.optimize import fsolve
fig2 = plt.figure()

t2_MF = t[1] 
t3_MF = t[2] 
t4_MF = t[3] 
t5_MF = t[4] 
t6_MF = t[5]

title = "$t_2=$" + str(t2_MF) + "$,\ t_3=$" + str(t3_MF) + "$,\ t_4=$" + str(t4_MF) + "$,\ t_5=$" + str(t5_MF) + "$,\ t_6=$" + str(t6_MF)
fig, ax = plt.subplots(figsize=(10., 6.*1.1))

T=1
x_min = -12 
x_max = -2 

# Plotting self-consistency equation
delta_t1 = 0.05
delta_L = 0.000005
t1_range = np.arange(x_min, x_max+delta_t1, delta_t1)
L_range = np.arange(0, 1, delta_L)
L, T1 = np.meshgrid(L_range, t1_range)
SC = L - 1/2*(1 + np.tanh(T1 + 2*t2_MF*L + 3*t3_MF*L**2 + 4*t4_MF*L**3 + 5*t5_MF*L**4 + 6*t6_MF*L**5)) #Self-consistency
CS = plt.contour(T1, L, SC, [0], colors='cyan', zorder=2, linewidths=1.5)
CS.collections[0].set_label('Equation of state')

plt.axvline(0, color='black')
plt.axhline(0, color='black', alpha=.5)

plt.ylim([-.01, 1.01])
plt.xlim([x_min, x_max])

plt.title(title, fontsize=18.5)
plt.xlabel(r'$t_1$', fontsize=20)
plt.ylabel(r'$\bar L$', fontsize=20)

# Finding ensemble averages in thermodynamic limit
n_nodes = 2
N_links_max = int(n_nodes*(n_nodes-1)/2)

x_min_ = x_min
x_max_ = x_max

T1_MF = np.arange(x_min_, x_max_, 0.001) #0.0001
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
    F = -2*(t1_MF*L + t2_MF*L**2 + t3_MF*L**3 + t4_MF*L**4 + t5_MF*L**5 + t6_MF*L**6) - T*S
    numerator, denominator = 0, 0
    Z = 0
    for i in range(0, N_links_max+1):
        numerator += L[i]*np.exp(-N_links_max/T*F[i])
        denominator += np.exp(-N_links_max/T*F[i])
    L_EXPECTED.append(numerator/denominator)

plt.plot(T1_MF, L_EXPECTED, color='orange', label='Exact solution for n=2', zorder=3)

######################################################
n_nodes = 5
N_links_max = int(n_nodes*(n_nodes-1)/2)

x_min_ = x_min
x_max_ = x_max

T1_MF = np.arange(x_min_, x_max_, .001) #0.0001
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
    F = -2*(t1_MF*L + t2_MF*L**2 + t3_MF*L**3 + t4_MF*L**4 + t5_MF*L**5 + t6_MF*L**6) - T*S
    numerator, denominator = 0, 0
    Z = 0
    for i in range(0, N_links_max+1):
        numerator += L[i]*np.exp(-N_links_max/T*F[i])
        denominator += np.exp(-N_links_max/T*F[i])
    L_EXPECTED.append(numerator/denominator)

plt.plot(T1_MF, L_EXPECTED, color='green', label='Exact solution for n=5', zorder=4)
################################################################
n_nodes = 10
N_links_max = int(n_nodes*(n_nodes-1)/2)

x_min_ = x_min
x_max_ = x_max

T1_MF = np.arange(x_min_, x_max_, .001)#0.0001
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
    F = -2*(t1_MF*L + t2_MF*L**2 + t3_MF*L**3 + t4_MF*L**4 + t5_MF*L**5 + t6_MF*L**6) - T*S
    numerator, denominator = 0, 0
    Z = 0
    for i in range(0, N_links_max+1):
        numerator += L[i]*np.exp(-N_links_max/T*F[i])
        denominator += np.exp(-N_links_max/T*F[i])
    L_EXPECTED.append(numerator/denominator)
    
plt.plot(T1_MF, L_EXPECTED, color='blue', label='Exact solution for n=10', zorder=5)
##################################################################
x_min_ = x_min 
x_max_ = -3

T1_MF = np.arange(x_min_, x_max_, .0001)
L_MIN = []
L_EXPECTED = []
T1_MF_ = []
for t1_MF in T1_MF:
    def SC_equation(L):
        return L - 1/2*(1 + np.tanh(t1_MF + 2*t2_MF*L + 3*t3_MF*L**2 + 4*t4_MF*L**3 + 5*t5_MF*L**4 + 6*t6_MF*L**5))
    L_ = np.array([fsolve(SC_equation, 0)[0], fsolve(SC_equation, .8)[0], fsolve(SC_equation, .3)[0], fsolve(SC_equation, .999)[0], fsolve(SC_equation, 1.001)[0], fsolve(SC_equation, .99)[0], fsolve(SC_equation, 1)[0]])
    def F(L):
        return -2*(t1_MF*L + t2_MF*L**2 + t3_MF*L**3 + t4_MF*L**4 + t5_MF*L**5 + t6_MF*L**6) + L*np.log(L/(1-L)) + np.log(1-L)
    L0 = L_[~(np.triu(np.abs(L_[:,None] - L_) < 0.001,1)).any(0)] # Distinct solution of SC equation
    L_ = np.array([l for l in L0 if (l<1 and l>0)])        
    L_min = L_[list(F(L_)).index(min(F(L_)))]
    L_MIN.append(L_min)

plt.plot(T1_MF, L_MIN, color='red', label='Thermodynamic solution', zorder=2)

plt.grid()

plt.legend(loc=[.58, .05], prop={'size': 16}).set_zorder(200)

plt.tight_layout()

plt.savefig('../figures/MF6_SC_with_n_2_5_10_unrounded_params.png', dpi=300)

plt.show()

fig.clear()
plt.close()
