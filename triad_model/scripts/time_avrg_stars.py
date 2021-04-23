################################################################################
# This script plots time averages of 2-stars and triangles for large networks  #
# together  with the corresponding (0,1)-bistability transition formula.       #
################################################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

n_nodes = 1000
n_pairs = n_nodes*(n_nodes-1)/2
t2_MF = 2
t3_MF = 1
x_min = -5.5
x_max = -1

fig, axs = plt.subplots(1, 2, figsize=(10.*1.2, 3.*1.2))

file_names = ["../data/triad_model/t1_scanning/TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_-1__t1_fin_-5.500000__t1_step_-0.010000__t2_MF_2__t3_MF_1__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_-5.500000__t1_fin_-1__t1_step_0.010000__t2_MF_2__t3_MF_1__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv"
]

t1_PS = []
S2_avrg_PS = []
S2_RMS_PS = []
for file_name in file_names:
    S2_vs_t1 = np.genfromtxt(file_name, delimiter=',')

    # Finding averages
    S2_avrg = []
    S2_RMS = []
    t1 = []
    i = 0
    while i < S2_vs_t1.shape[0]:
        prev_t1 = S2_vs_t1[i][0]
        S2 = []
        while i < S2_vs_t1.shape[0] and prev_t1 == S2_vs_t1[i][0]:
            S2.append(S2_vs_t1[i][2]*2/((n_nodes-1)*(n_nodes-2)))
            i+=1
        if np.sqrt(np.array(S2).var()) < .2: # Filtering out transition points
            S2_avrg.append(np.array(S2).mean())
            S2_RMS.append(np.sqrt(np.array(S2).var()))
            t1.append(prev_t1)
    t1_PS.extend(t1)
    S2_avrg_PS.extend(S2_avrg)
    S2_RMS_PS.extend(S2_RMS)

axs[0].errorbar(t1_PS, S2_avrg_PS, xerr=None, yerr=S2_RMS_PS, color="blue", marker='o', markersize=4, ls='none', label="Triad simulation", alpha=.5, elinewidth=2, capsize=2, zorder=51)


# Plotting-consistency equation
delta_t1 = 0.05
delta_L = 0.000005
t1_range = np.arange(x_min, x_max+.5, delta_t1)
S2_range = np.arange(0, 1, delta_L)
L, T1 = np.meshgrid(S2_range, t1_range)
SC = L - 1/2*(1 + np.tanh(T1 +2*t2_MF*L + 3*t3_MF*L**2)) #Self-consistency
CS = axs[0].contour(T1, (L)**2, SC, [0], colors='cyan', zorder=100, linewidths=1.5)

CS.collections[0].set_label('MF solution')

axs[0].axvline(0, color='black')
axs[0].axhline(0, color='black')

axs[0].set_xlim([x_min, x_max])
axs[0].set_ylim([-.01, 1.01])


axs[0].set_xlabel(r'$t_1$', fontsize=22)
axs[0].set_ylabel(r'$\mathfrak{S}_2$', fontsize=22)


## THEORETICAL ##
t1_ = np.arange(x_min, x_max, .001)
L_ensemble = 1/2*(1 + np.tanh(n_pairs*(t1_ + t2_MF + t3_MF)))
axs[0].plot(t1_, L_ensemble**2, color='g', label='Transition formula', zorder=5, linewidth=2)


axs[0].grid(zorder=0)
axs[0].legend(loc=[.58, .59], prop={'size': 11})


###########################################################
###################### 3_STARS ############################
###########################################################
file_names = ["../data/triad_model/t1_scanning/TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_-1__t1_fin_-5.500000__t1_step_-0.010000__t2_MF_2__t3_MF_1__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_-5.500000__t1_fin_-1__t1_step_0.010000__t2_MF_2__t3_MF_1__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv"
]

t1_PS = []
S3_avrg_PS = []
S3_RMS_PS = []
for file_name in file_names:
    S3_vs_t1 = np.genfromtxt(file_name, delimiter=',')

    # Finding averages
    S3_avrg = []
    S3_RMS = []
    t1 = []
    i = 0
    while i < S3_vs_t1.shape[0]:
        prev_t1 = S3_vs_t1[i][0]
        S2 = []
        while i < S3_vs_t1.shape[0] and prev_t1 == S3_vs_t1[i][0]:
            S2.append(S3_vs_t1[i][3]*6/(n_nodes*(n_nodes-1)*(n_nodes-2)))
            i+=1
        if np.sqrt(np.array(S2).var()) < .2: # Filtering out transition points
            S3_avrg.append(np.array(S2).mean())
            S3_RMS.append(np.sqrt(np.array(S2).var()))
            t1.append(prev_t1)
    t1_PS.extend(t1)
    S3_avrg_PS.extend(S3_avrg)
    S3_RMS_PS.extend(S3_RMS)


axs[1].errorbar(t1_PS, S3_avrg_PS, xerr=None, yerr=S3_RMS_PS, color="blue", marker='o', markersize=4, ls='none', label="Triad simulation", alpha=.5, elinewidth=2, capsize=2, zorder=51)



# Plotting-consistency equation
delta_t1 = 0.05
delta_L = 0.000005
t1_range = np.arange(x_min, x_max+.5, delta_t1)
S3_range = np.arange(0, 1, delta_L)
L, T1 = np.meshgrid(S3_range, t1_range)
SC = L - 1/2*(1 + np.tanh(T1 +2*t2_MF*L + 3*t3_MF*L**2)) #Self-consistency
CS = axs[1].contour(T1, (L)**3, SC, [0], colors='cyan', zorder=100, linewidths=1.5)

CS.collections[0].set_label('MF solution')

axs[1].axvline(0, color='black')
axs[1].axhline(0, color='black')

axs[1].set_xlim([x_min, x_max])
axs[1].set_ylim([-.01, 1.01])

axs[1].set_xlabel(r'$t_1$', fontsize=22)
axs[1].set_ylabel(r'$\mathfrak{T}$', fontsize=22)


## THEORETICAL ##
t1_ = np.arange(x_min, x_max, .001)
L_ensemble = 1/2*(1 + np.tanh(n_pairs*(t1_ + t2_MF + t3_MF)))
axs[1].plot(t1_,  L_ensemble**3, color='g', label='Transition formula', zorder=5, linewidth=2)


axs[1].grid(zorder=0)
axs[1].legend(loc=[.58, .59], prop={'size': 11})

plt.tight_layout()
plt.savefig('../figures/time_avrg_stars.png', dpi=300)

plt.show()

fig.clear()
