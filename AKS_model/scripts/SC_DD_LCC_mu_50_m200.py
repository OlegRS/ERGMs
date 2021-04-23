import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as ss
from scipy.optimize import fsolve
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

fig, axs = plt.subplots(2, 3, figsize=(19./1.5, 6./1.5*2))

mu = 50
lam = .04
n_nodes = 200

file_names = ["../data/AKS_model/tau1_scanning/AKS___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_200__N_ITERS_PER_LINK_10__tau1_start_-0.600000__tau1_fin_-3.300000__tau1_step_-0.050000__mu_50__lambda_0.040000__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv",
              "../data/AKS_model/tau1_scanning/AKS___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_200__N_ITERS_PER_LINK_10__tau1_start_-3.300000__tau1_fin_-0.600000__tau1_step_0.050000__mu_50__lambda_0.040000__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv"
]

t1_PS = []
L_avrg_PS = []
L_RMS_PS = []
for file_name in file_names:
    L_vs_t1 = np.genfromtxt(file_name, delimiter=',')

    # Finding averages
    L_avrg = []
    L_RMS = []
    t1 = []
    i = 0
    while i < L_vs_t1.shape[0]:
        prev_t1 = L_vs_t1[i][0]
        L = []
        while i < L_vs_t1.shape[0] and prev_t1 == L_vs_t1[i][0]:
            L.append(L_vs_t1[i][1])
            i+=1
        if np.sqrt(np.array(L).var()) < .2: # Filtering out transition points
            L_avrg.append(np.array(L).mean())
            L_RMS.append(np.sqrt(np.array(L).var()))
            t1.append(prev_t1)
    t1_PS.extend(t1)
    L_avrg_PS.extend(L_avrg)
    L_RMS_PS.extend(L_RMS)
axs[0,0].errorbar(t1_PS, L_avrg_PS, xerr=None, yerr=L_RMS_PS, color="blue", marker='o', markersize=4, ls='none', label="AKS simulation", alpha=.5, elinewidth=2, capsize=2, zorder=51)

##### PLOTTING SELF-CONSISTENCY EQUATION #####
x_min, x_max = -3.1, -.9
delta_tau1 = .05
delta_L = .001
tau1_range = np.arange(x_min, x_max+.5, delta_tau1)
L_range = np.arange(0, 1+2*delta_L, delta_L)
L, TAU1 = np.meshgrid(L_range, tau1_range)
SC = L - 1/2*(1 + np.tanh(TAU1 + mu*lam*(1 - (1-L/(n_nodes*lam))**(n_nodes-2)))) #Self-consistency
CS = axs[0,0].contour(TAU1, L, SC, [0], colors='cyan', zorder=50, linewidths=2.5, alpha=1)
CS.collections[0].set_label('SC equation')

Ld = .3
rho_d = 1/2*np.log(Ld/(1-Ld)) - mu*lam*(1-(1-Ld/(n_nodes*lam))**(n_nodes-2))
print("rho_d=", rho_d)
axs[0,0].plot(rho_d, Ld, marker='+', color='green', markeredgewidth=1.2, markersize=25, zorder=150)

title = r"$\tilde\mu=50$" + r"$;\ \tilde\lambda=0.04$" + r"$;\ n=200$"
axs[0,0].set_title(title, fontsize=13.5)

axs[0,0].set_ylabel(r'$\bar L$', fontsize=18)
# axs[0,0].set_xlabel(r'$\rho$', fontsize=18)
axs[0,0].grid()

axs[0,0].set_xlim([x_min, x_max])
axs[0,0].set_ylim([-0.015, 1.015])
axs[0,0].axhline(0, color='black')


##### Exact thermodynamic solution #####
x_min_, x_max_ = -2.9, -2.8
RHO = np.arange(x_min_, x_max_, 0.0001)
L_MIN = []
for rho in RHO:
    def SC_equation(L):
        return L - 1/2*(1 + np.tanh(rho + mu*lam*(1 - (1-L/(n_nodes*lam))**(n_nodes-2))))
    L_ = np.array([fsolve(SC_equation, .00001)[0], fsolve(SC_equation, .2)[0]])
    def F(L):
        return -2*rho*L - 2*mu*lam**2*(n_nodes/(n_nodes-1)*(1-L/(n_nodes*lam))**(n_nodes-1)+L/lam-n_nodes/(n_nodes-1)) + L*np.log(L/(1-L)) + np.log(1-L)
    L0 = L_[~(np.triu(np.abs(L_[:,None] - L_) < 0.001,1)).any(0)] # Distinct solution of SC equation
    L_ = np.array([l for l in L0 if (l<1 and l>0)])        
    L_min = L_[list(F(L_)).index(min(F(L_)))]
    L_MIN.append(L_min)
axs[0,0].plot(RHO, L_MIN, color="red", label='Phase transition', zorder=10)

axs[0,0].legend(loc='upper left', prop={'size': 12})

##### ZOOM ######
axins = zoomed_inset_axes(axs[0,0], 3, bbox_to_anchor=(.55, .07, .6, .48), bbox_transform=axs[0,0].transAxes, loc="lower left")

CS = axins.contour(TAU1, L, SC, [0], colors='cyan', zorder=100, linewidths=2.5, alpha=.7)
axins.plot(RHO, L_MIN, color="red", label='Phase transition', zorder=10)
axins.plot(t1_, L_ensemble, color='g', label='Transition formula', zorder=5, linewidth=2)
axins.errorbar(t1_PS, L_avrg_PS, xerr=None, yerr=L_RMS_PS, color="blue", marker='o', markersize=4, ls='none', label="AKS simulation", alpha=.5, elinewidth=2, capsize=2, zorder=51)
axins.grid()

x1, x2, y1, y2 = -3.001, -2.7, 0, .14 # specify the limits
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.axhline(0, color='black')

from mpl_toolkits.axes_grid1.inset_locator import mark_inset
mark_inset(axs[0,0], axins, loc1=3, loc2=4, fc="none", ec='0.5', alpha=.5)


file_name = "../data/AKS_model/degree_distributions/AKS__N_NODES_200__N_ITERS_PER_LINK_100__tau1_-2.422615__mu_50__lambda_0.040000__N_AVRGNG_10000__INITIALIZE_RANDOMLY_0.csv"

mu = -2
lam = 1
tau1 = 0.093069

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_PS = []
for line in lines:
    deg_seq_PS.append(int(line.split('\n')[0]))
N_lines = np.array(lines).shape[0]
L_AKS = np.average(deg_seq_PS)
print(L_AKS)


# Theoretical
L_AKS /= (n_nodes-1)
sigma = np.sqrt((n_nodes-1)*L_AKS*(1-L_AKS))
x_min, x_max = L_AKS*(n_nodes-1) - 5*sigma, L_AKS*(n_nodes-1) + 5*sigma
degrees = np.arange(x_min, x_max)

p = np.sqrt(1/(2*np.pi*(n_nodes-1)*L_AKS*(1-L_AKS)))*np.exp(-(degrees-L_AKS*(n_nodes-1))**2/(2*(n_nodes-1)*L_AKS*(1-L_AKS)))
axs[0,1].plot(degrees, p, color='r', label='MF theor')
axs[0,1].hist(deg_seq_PS, bins = degrees, color='b', density=True, label='AKS sim')
axs[0,1].legend(loc='upper left', prop={'size': 11})
axs[0,1].grid()

title = r"$\rho=-2.423$" + r"$,\ \tilde\mu=50$" + r"$;\ \tilde\lambda=0.04$" + r"$;\ n=200$"
axs[0,1].set_title(title, fontsize=13.5)

axs[0,1].set_ylabel(r'$\langle \mathrm{p}(k|\mathbf{A}) \rangle$', fontsize=18)
# axs[0,1].set_xlabel(r'$k$', fontsize=18)
axs[0,1].set_xlim([30, 90])


file_name = "../data/AKS_model/clustering_distributions/AKS__N_NODES_200__N_ITERS_PER_LINK_100__tau1_-2.422615__mu_50__lambda_0.040000__N_AVRGNG_10000__INITIALIZE_RANDOMLY_0.csv"

mu = 200
lam = .01
degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_PS = []
for line in lines:
    deg_seq_PS.append(float(line.split('\n')[0]))
    
N_lines = np.array(lines).shape[0]
print(np.average(deg_seq_PS))

L_sol = .3
cc_theor = L_sol*(1-(1-L_sol)**(n_nodes-1)-L_sol*(n_nodes-1)*(1-L_sol)**(n_nodes-2))
print("C_theor=", cc_theor)

axs[0,2].hist(deg_seq_PS, 100, color='b', density=True, label="AKS sim")
axs[0,2].grid()
axs[0,2].axvline(cc_theor, color='r', alpha=.65, label="MF theor")
axs[0,2].legend(fontsize=11, loc="upper left")
axs[0,2].set_xlim([.25,.35])
axs[0,2].set_ylabel(r'$\langle \mathrm{p}(c|\mathbf{A}) \rangle$', fontsize=18)
title = r"$\rho=-2.423$" + r"$,\ \tilde\mu=50$" + r"$;\ \tilde\lambda=0.04$" + r"$;\ n=200$"
axs[0,2].set_title(title, fontsize=13.5)


mu = -200
lam = .04
n_nodes = 200

file_names = ["../data/AKS_model/tau1_scanning/AKS___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_200__N_ITERS_PER_LINK_10__tau1_start_1.500000__tau1_fin_11.500000__tau1_step_0.100000__mu_-200__lambda_0.040000__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv"]

t1_PS = []
L_avrg_PS = []
L_RMS_PS = []
for file_name in file_names:
    L_vs_t1 = np.genfromtxt(file_name, delimiter=',')

    # Finding averages
    L_avrg = []
    L_RMS = []
    t1 = []
    i = 0
    while i < L_vs_t1.shape[0]:
        prev_t1 = L_vs_t1[i][0]
        L = []
        while i < L_vs_t1.shape[0] and prev_t1 == L_vs_t1[i][0]:
            L.append(L_vs_t1[i][1])
            i+=1
        if np.sqrt(np.array(L).var()) < .2: # Filtering out transition points
            L_avrg.append(np.array(L).mean())
            L_RMS.append(np.sqrt(np.array(L).var()))
            t1.append(prev_t1)
    t1_PS.extend(t1)
    L_avrg_PS.extend(L_avrg)
    L_RMS_PS.extend(L_RMS)
axs[1,0].errorbar(t1_PS, L_avrg_PS, xerr=None, yerr=L_RMS_PS, color="blue", marker='o', markersize=4, ls='none', label="AKS simulation", alpha=.5, elinewidth=2, capsize=2, zorder=51)


##### PLOTTING SELF-CONSISTENCY EQUATION #####
x_min, x_max = 2, 11
delta_tau1 = .05
delta_L = .001
tau1_range = np.arange(x_min, x_max+.5, delta_tau1)
L_range = np.arange(0, 1+2*delta_L, delta_L)
L, TAU1 = np.meshgrid(L_range, tau1_range)
SC = L - 1/2*(1 + np.tanh(TAU1 + mu*lam*(1 - (1-L/(n_nodes*lam))**(n_nodes-2)))) #Self-consistency
CS = axs[1,0].contour(TAU1, L, SC, [0], colors='cyan', zorder=50, linewidths=2.5, alpha=1)
CS.collections[0].set_label('SC equation')

Ld = .3
rho_d = 1/2*np.log(Ld/(1-Ld)) - mu*lam*(1-(1-Ld/(n_nodes*lam))**(n_nodes-2))
print("rho_d=", rho_d)
axs[1,0].plot(rho_d, Ld, marker='+', color='green', markeredgewidth=1.2, markersize=25, zorder=150)

title = r"$\tilde\mu=-200$" + r"$;\ \tilde\lambda=0.04$" + r"$;\ n=200$"
axs[1,0].set_title(title, fontsize=13)

axs[1,0].set_ylabel(r'$\bar L$', fontsize=18)
axs[1,0].set_xlabel(r'$\rho$', fontsize=18)
axs[1,0].grid()

axs[1,0].legend(loc='upper left', prop={'size': 12})

axs[1,0].set_xlim([x_min, x_max])
axs[1,0].set_ylim([-0.015, 1.015])
axs[1,0].axhline(0, color='black')


file_name = "../data/AKS_model/degree_distributions/AKS__N_NODES_200__N_ITERS_PER_LINK_100__tau1_7.572216__mu_-200__lambda_0.040000__N_AVRGNG_10000__INITIALIZE_RANDOMLY_0.csv"

mu = -2
lam = 1
tau1 = 0.093069

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_PS = []
for line in lines:
    deg_seq_PS.append(int(line.split('\n')[0]))
N_lines = np.array(lines).shape[0]
L_AKS = np.average(deg_seq_PS)
print(L_AKS)


# Theoretical
L_AKS /= (n_nodes-1)
sigma = np.sqrt((n_nodes-1)*L_AKS*(1-L_AKS))
x_min, x_max = L_AKS*(n_nodes-1) - 5*sigma, L_AKS*(n_nodes-1) + 5*sigma
degrees = np.arange(x_min, x_max)

p = np.sqrt(1/(2*np.pi*(n_nodes-1)*L_AKS*(1-L_AKS)))*np.exp(-(degrees-L_AKS*(n_nodes-1))**2/(2*(n_nodes-1)*L_AKS*(1-L_AKS)))
axs[1,1].plot(degrees, p, color='r', label='MF theor')
axs[1,1].hist(deg_seq_PS, bins = degrees, color='b', density=True, label='AKS sim')
axs[1,1].legend(loc='upper left', prop={'size': 11})
axs[1,1].grid()

title = r"$\rho=7.572$" + r"$,\ \tilde\mu=-200$" + r"$;\ \tilde\lambda=0.04$" + r"$;\ n=200$"
axs[1,1].set_title(title, fontsize=13)

axs[1,1].set_ylabel(r'$\langle \mathrm{p}(k|\mathbf{A}) \rangle$', fontsize=18)
axs[1,1].set_xlabel(r'$k$', fontsize=18)
axs[1,1].set_xlim([30, 90])

file_name = "../data/AKS_model/clustering_distributions/AKS__N_NODES_200__N_ITERS_PER_LINK_100__tau1_7.572216__mu_-200__lambda_0.040000__N_AVRGNG_10000__INITIALIZE_RANDOMLY_0.csv"

mu = 200
lam = .01
degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_PS = []
for line in lines:
    deg_seq_PS.append(float(line.split('\n')[0]))
    
N_lines = np.array(lines).shape[0]
print(np.average(deg_seq_PS))

L_sol = .3
cc_theor = L_sol*(1-(1-L_sol)**(n_nodes-1)-L_sol*(n_nodes-1)*(1-L_sol)**(n_nodes-2))
print("C_theor=", cc_theor)

axs[1,2].hist(deg_seq_PS, 100, color='b', density=True, label="AKS sim")
axs[1,2].grid()
axs[1,2].axvline(cc_theor, color='r', alpha=.65, label="MF theor")
axs[1,2].legend(fontsize=11, loc="upper left")
axs[1,2].set_xlabel(r'$c$', fontsize=18)
axs[1,2].set_xlim([.25,.35])
axs[1,2].set_ylabel(r'$\langle \mathrm{p}(c|\mathbf{A}) \rangle$', fontsize=18)
title = r"$\rho=7.572$" + r"$,\ \tilde\mu=-200$" + r"$;\ \tilde\lambda=0.04$" + r"$;\ n=200$"
axs[1,2].set_title(title, fontsize=13)


plt.tight_layout()

plt.savefig('../figures/SC_DD_LCC_mu_50_m200.png', dpi=300)
plt.show()

fig.clear()
