############################################################
# This script plots time averages of connectance for the   #
# AKS model together with solution of the self-consistency #
# equation and (0,1)-bistability transition formula.       #                                     #
############################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

fig, axs = plt.subplots(1, 2, figsize=(15., 6.))

n_nodes = 150
mu = 1
lam = 1
x_min = -3
x_max = 2
n_pairs = n_nodes*(n_nodes-1)/2

title = r"$\tilde\mu=$" + str(mu) + r"$;\ \tilde\lambda=$" + str(lam) + r"$;\ n=$" + str(n_nodes)

file_names = ["../data/AKS_model/tau1_scanning/AKS___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_150__N_ITERS_PER_LINK_10__tau1_start_-3.100000__tau1_fin_2__tau1_step_0.050000__mu_1__lambda_1__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv"]

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

axs[0].errorbar(t1_PS, L_avrg_PS, xerr=None, yerr=L_RMS_PS, color="blue", marker='o', markersize=4, ls='none', label="AKS simulation", alpha=.5, elinewidth=2, capsize=2, zorder=51)


# Plotting-consistency equation
delta_tau1 = .05
delta_L = .001
tau1_range = np.arange(x_min, x_max+.5, delta_tau1)
L_range = np.arange(0, 1+2*delta_L, delta_L)
L, TAU1 = np.meshgrid(L_range, tau1_range)
SC = L - 1/2*(1 + np.tanh(TAU1 + mu*lam*(1 - (1-L/(n_nodes*lam))**(n_nodes-2)))) #Self-consistency
CS = axs[0].contour(TAU1, L, SC, [0], colors='cyan', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('SC equation')

axs[0].axvline(0, color='black')
axs[0].axhline(0, color='black')

axs[0].set_ylim([-.01, 1.01])
axs[0].set_xlim([x_min, x_max])

axs[0].set_title(title, fontsize=22)
axs[0].set_xlabel(r'$\rho$', fontsize=22)
axs[0].set_ylabel(r'$\bar L$', fontsize=22)


axs[0].grid(zorder=0)
axs[0].legend(loc='upper left', prop={'size': 15})
# axs[0].legend(loc='upper right', prop={'size': 17})

###########################
n_nodes = 150
mu = 10
lam = 1
x_min = -5
x_max = -1.5
n_pairs = n_nodes*(n_nodes-1)/2

title = r"$\tilde\mu=$" + str(mu) + r"$;\ \tilde\lambda=$" + str(lam) + r"$;\ n=$" + str(n_nodes)

file_names = ["../data/AKS_model/tau1_scanning/AKS___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_150__N_ITERS_PER_LINK_10__tau1_start_-5__tau1_fin_-1.500000__tau1_step_0.010000__mu_10__lambda_1__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv",
              "../data/AKS_model/tau1_scanning/AKS___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_150__N_ITERS_PER_LINK_10__tau1_start_-1.500000__tau1_fin_-5__tau1_step_-0.010000__mu_10__lambda_1__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv"]

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

axs[1].errorbar(t1_PS, L_avrg_PS, xerr=None, yerr=L_RMS_PS, color="blue", marker='o', markersize=4, ls='none', label="AKS simulation", alpha=.5, elinewidth=2, capsize=2, zorder=51)


# Plotting-consistency equation
delta_tau1 = .05
delta_L = .001
tau1_range = np.arange(x_min, x_max+.5, delta_tau1)
L_range = np.arange(0, 1+2*delta_L, delta_L)
L, TAU1 = np.meshgrid(L_range, tau1_range)
SC = L - 1/2*(1 + np.tanh(TAU1 + mu*lam*(1 - (1-L/(n_nodes*lam))**(n_nodes-2)))) #Self-consistency
CS = axs[1].contour(TAU1, L, SC, [0], colors='cyan', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('SC equation')

axs[1].axvline(0, color='black')
axs[1].axhline(0, color='black')

axs[1].set_ylim([-.01, 1.01])
axs[1].set_xlim([x_min, x_max])

axs[1].set_title(title, fontsize=22)
axs[1].set_xlabel(r'$\rho$', fontsize=22)
# axs[1].set_ylabel(r'$\bar L$', fontsize=22)


t1_ = np.arange(x_min, x_max, .001)
L_ensemble = 1/2*(1 + np.tanh(t1_*n_nodes*(n_nodes-1)/2 + (n_nodes*lam)**2*mu/2*((1-1/(n_nodes*lam))**(n_nodes-1) + (n_nodes-1)/(n_nodes*lam) - 1)))
axs[1].plot(t1_, L_ensemble, color='g', label='Transition formula', zorder=5, linewidth=2)

##### Exact thermodynamic solution #####
x_min_, x_max_ = -4.5, -3 #35.6261
RHO = np.arange(x_min_, x_max_, 0.0001)
L_MIN = []
for rho in RHO:
    def SC_equation(L):
        return L - 1/2*(1 + np.tanh(rho + mu*lam*(1 - (1-L/(n_nodes*lam))**(n_nodes-2))))
    L_ = np.array([fsolve(SC_equation, .00001)[0], fsolve(SC_equation, .99)[0]])
    def F(L):
        return -2*rho*L - 2*mu*lam**2*(n_nodes/(n_nodes-1)*(1-L/(n_nodes*lam))**(n_nodes-1)+L/lam-n_nodes/(n_nodes-1)) + L*np.log(L/(1-L)) + np.log(1-L)
    L0 = L_[~(np.triu(np.abs(L_[:,None] - L_) < 0.001,1)).any(0)] # Distinct solution of SC equation
    L_ = np.array([l for l in L0 if (l<1 and l>0)])        
    L_min = L_[list(F(L_)).index(min(F(L_)))]
    L_MIN.append(L_min)

axs[1].plot(RHO, L_MIN, color="red", label='Phase transition', zorder=10)


axs[1].grid(zorder=0)
axs[1].legend(loc=[.575, .32], prop={'size': 15})

##### ZOOM ######
axins = zoomed_inset_axes(axs[1], 13, bbox_to_anchor=(.51, .63, .6, .48), bbox_transform=axs[1].transAxes, loc="lower left")

CS = axins.contour(TAU1, L, SC, [0], colors='cyan', zorder=100, linewidths=2.5, alpha=.7)
axins.plot(RHO, L_MIN, color="red", label='Phase transition', zorder=10)
axins.plot(t1_, L_ensemble, color='g', label='Transition formula', zorder=5, linewidth=2)
axins.errorbar(t1_PS, L_avrg_PS, xerr=None, yerr=L_RMS_PS, color="blue", marker='o', markersize=4, ls='none', label="AKS simulation", alpha=.5, elinewidth=2, capsize=2, zorder=51)
axins.grid()

# x1, x2, y1, y2 = -3.75, -3.55, .95, 1.005 # specify the limits
x1, x2, y1, y2 = -3.651-.05, -3.651+.05, .98, 1.001 # specify the limits
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
# plt.yticks(visible=False)
# plt.xticks(visible=False)

from mpl_toolkits.axes_grid1.inset_locator import mark_inset
mark_inset(axs[1], axins, loc1=2, loc2=3, fc="none", ec="0.5")


plt.tight_layout()

plt.savefig('../figures/time_avrg_connectance_merged.png', dpi=300)

plt.show()

fig.clear()
