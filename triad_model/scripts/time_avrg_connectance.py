###################################################################
# This script plots time averages of connectance for the          #
# corresponding MF and triad models together with the numerical   #
# solution of the self-consistency equation and (0,1)-bistability #
# transition formula.                                             #
###################################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

n_nodes = 1000
t2_MF = 2
t3_MF = 1
x_min = -5.5
x_max = -1
n_pairs = n_nodes*(n_nodes-1)/2


title = "$t_2=$" + str(t2_MF) + "$;\ t_3=$" + str(t3_MF) + "$;\ n=1000$"
fig, ax = plt.subplots(figsize=(10., 6.))

file_names = ["../data/triad_model/t1_scanning/TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_-1__t1_fin_-5.500000__t1_step_-0.010000__t2_MF_2__t3_MF_1__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_-5.500000__t1_fin_-1__t1_step_0.010000__t2_MF_2__t3_MF_1__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv"
]

t2 = 2 * 2*n_nodes/(n_nodes-1) # t2_MF*(CONVERSION_FACTOR)
t3 = 1 * 6*n_nodes**2/(n_nodes**2-3*n_nodes+2) # t3_MF*(CONVERSION_FACTOR)

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

plt.errorbar(t1_PS, L_avrg_PS, xerr=None, yerr=L_RMS_PS, color="blue", marker='o', markersize=4, ls='none', label="Triad simulation", alpha=.5, elinewidth=2, capsize=2, zorder=51)


file_names = ["../data/mean_field_model/t1_scanning/NEW_METROPOLIS_MF__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_-5.500000__t1_fin_-1__t1_step_0.010000__t2_2__t3_1__N_DYNAMIC_PAIRS_499500__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv",
              "../data/mean_field_model/t1_scanning/NEW_METROPOLIS_MF__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_-1__t1_fin_-5.500000__t1_step_-0.010000__t2_2__t3_1__N_DYNAMIC_PAIRS_499500__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv"
]

t1_MF = []
L_avrg_MF = []
L_RMS_MF = []
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
    t1_MF.extend(t1)
    L_avrg_MF.extend(L_avrg)
    L_RMS_MF.extend(L_RMS)
    
plt.errorbar(t1_MF, L_avrg_MF, xerr=None, yerr=L_RMS_MF, color="orange", marker='s', markersize=4, ls='none', label="MF3 simulation", alpha=1, elinewidth=2, capsize=2, zorder=50)


# Plotting-consistency equation
delta_t1 = 0.05
delta_L = 0.000005
t1_range = np.arange(x_min, x_max+.5, delta_t1)
L_range = np.arange(0, 1, delta_L)
L, T1 = np.meshgrid(L_range, t1_range)
SC = L - 1/2*(1 + np.tanh(T1 +t2*(n_nodes-1)/n_nodes*L + t3*(n_nodes**2-3*n_nodes+2)/(2*n_nodes**2)*L**2)) #Self-consistency
CS = plt.contour(T1, L, SC, [0], colors='cyan', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('SC equation')

plt.axvline(0, color='black')
plt.axhline(0, color='black')

plt.ylim([-.01, 1.01])
plt.xlim([x_min, x_max])

plt.title(title, fontsize=22)
plt.xlabel(r'$t_1$', fontsize=22)
plt.ylabel(r'$\bar L$', fontsize=22)


## THEORETICAL ##
t1_ = np.arange(x_min, x_max, .001)
L_ensemble = 1/2*(1 + np.tanh(n_pairs*(t1_ + 1/2*t2*(n_nodes-1)/n_nodes + t3*(n_nodes**2-3*n_nodes+2)/(6*n_nodes**2))))
plt.plot(t1_, L_ensemble, color='g', label='Transition formula', zorder=5, linewidth=2)


plt.grid(zorder=0)
plt.legend(loc=[.64, .63], prop={'size': 17})

##### ZOOM ######
axins = zoomed_inset_axes(ax, 7.9, bbox_to_anchor=(0, .03, .6, .48), bbox_transform=ax.transAxes, loc="lower left")

CS = axins.contour(T1, L, SC, [0], colors='cyan', zorder=100, linewidths=2.5, alpha=.7)
axins.plot(t1_MF, L_avrg_MF, 's',markersize=4, label="MF3 simulation", zorder=50, alpha=1, color="orange")
axins.plot(t1_PS, L_avrg_PS, 'h',markersize=4, label="Triad simulation", zorder=51, alpha=.5, color="blue")

x1, x2, y1, y2 = -5.06, -4.84, .92, 1 # specify the limits
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
plt.yticks(visible=False)
plt.xticks(visible=False)

from mpl_toolkits.axes_grid1.inset_locator import mark_inset
mark_inset(ax, axins, loc1=2, loc2=1, fc="none", ec="0.5")


axins = zoomed_inset_axes(ax, 6, bbox_to_anchor=(.895, .05, .6, .5), bbox_transform=ax.transAxes, loc="lower left")

CS = axins.contour(T1, L, SC, [0], colors='cyan', zorder=100, linewidths=2.5, alpha=.7)
axins.errorbar(t1_MF, L_avrg_MF, xerr=None, yerr=L_RMS_MF, color="orange", marker='s', markersize=4, ls='none', label="MF3 simulation", alpha=1, elinewidth=2, capsize=2, zorder=50)
axins.errorbar(t1_PS, L_avrg_PS, xerr=None, yerr=L_RMS_PS, color="blue", marker='o', markersize=4, ls='none', label="Triad simulation", alpha=.5, elinewidth=2, capsize=2, zorder=51)

x1, x2, y1, y2 = -1.58, -1.51, .07, .14 # specify the limits
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
plt.yticks(visible=False)
plt.xticks(visible=False)

mark_inset(ax, axins, loc1=2, loc2=3, fc="none", ec="0.5")

plt.tight_layout()

plt.savefig('../figures/time_avrg_connectance.png', dpi=300)

plt.show()

fig.clear()
