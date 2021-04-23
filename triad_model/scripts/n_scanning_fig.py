#################################################################
# This script plots the dependencies of connectance and         #
# average local clustering coefficient on the number on nodes.  # 
#################################################################

import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt


fig, axs = plt.subplots(1, 2, figsize=(10.*1.4/1.1, 6./1.3))

t1_MF = 20
t2_MF = -75
t3_MF = 20

file_names = ["../data/triad_model/n_scanning/MATRIX_BASED__TRIAD__CONNECTANCE_2_STARS_and_LC__N_NODES_start_3__N_NODES_fin_1000__N_ITERS_PER_LINK_10__t1_20__t2_MF_-75__t3_MF_20__N_AVRGNG_10__INIT_RANDOM_0__N_DYNAMIC_PAIRS_1.csv"]

title = "$t_1=$" + str(t1_MF) + "$;\ t_2=$" + str(t2_MF) + "$;\ t_3=$" + str(t3_MF)

x_min = 0
x_max = 1000

for file_name in file_names:
    L_vs_t1 = np.genfromtxt(file_name, delimiter=',')

    # Finding averages
    L_avrg = []
    C_avrg = []
    L_RMS = []
    t1 = []
    t1_unaveraged = []
    L_unaveraged = []
    C_unaveraged = []
    i = 0
    while i < L_vs_t1.shape[0]:
        prev_t1 = L_vs_t1[i][0]
        L = []
        C = []
        count = 0
        while i < L_vs_t1.shape[0] and prev_t1 == L_vs_t1[i][0]:
            count += 1
            L.append(L_vs_t1[i][1])
            C.append(L_vs_t1[i][4])
            t1_unaveraged.append(L_vs_t1[i][0])
            L_unaveraged.append(L_vs_t1[i][1])
            C_unaveraged.append(L_vs_t1[i][4])
            i += 1
        L_avrg.append(np.array(L).mean())
        C_avrg.append(np.array(C).mean())
        L_RMS.append(np.sqrt(np.array(L).var()))
        t1.append(prev_t1)
    axs[0].plot(t1, L_avrg, 'h',markersize=4, label=r'$\bar L$')
    axs[0].plot(t1, C_avrg, 'h',markersize=4, label=r'$\bar c$')

axs[0].axvline(0, color='black')
axs[0].axhline(0, color='black')

axs[0].set_ylim([-.01, 1.01])
axs[0].set_xlim([x_min, x_max])
axs[0].legend(loc='best', prop={'size': 22})

axs[0].set_title(title, fontsize=20)
axs[0].set_xlabel(r'$n$', fontsize=22)
axs[0].set_ylabel(r'$\bar L,\ \bar c$', fontsize=22)
axs[0].grid()

######################################################################
t1_MF = 20
t2_MF = -75
t3_MF = 25

file_names = ["../data/triad_model/n_scanning/MATRIX_BASED__TRIAD__CONNECTANCE_2_STARS_and_LC__N_NODES_start_3__N_NODES_fin_1000__N_ITERS_PER_LINK_10__t1_20__t2_MF_-75__t3_MF_25__N_AVRGNG_10__INIT_RANDOM_0__N_DYNAMIC_PAIRS_1.csv"]

title = "$t_1=$" + str(t1_MF) + "$;\ t_2=$" + str(t2_MF) + "$;\ t_3=$" + str(t3_MF)

x_min = 0
x_max = 1000

for file_name in file_names:
    L_vs_t1 = np.genfromtxt(file_name, delimiter=',')

    # Finding averages
    L_avrg = []
    C_avrg = []
    L_RMS = []
    t1 = []
    t1_unaveraged = []
    L_unaveraged = []
    C_unaveraged = []
    i = 0
    while i < L_vs_t1.shape[0]:
        prev_t1 = L_vs_t1[i][0]
        L = []
        C = []
        count = 0
        while i < L_vs_t1.shape[0] and prev_t1 == L_vs_t1[i][0]:
            count += 1
            L.append(L_vs_t1[i][1])
            C.append(L_vs_t1[i][4])
            t1_unaveraged.append(L_vs_t1[i][0])
            L_unaveraged.append(L_vs_t1[i][1])
            C_unaveraged.append(L_vs_t1[i][4])
            i += 1
        L_avrg.append(np.array(L).mean())
        C_avrg.append(np.array(C).mean())
        L_RMS.append(np.sqrt(np.array(L).var()))
        t1.append(prev_t1)
    axs[1].plot(t1, L_avrg, 'h',markersize=4, label=r'$\bar L$')
    axs[1].plot(t1, C_avrg, 'h',markersize=4, label=r'$\bar c$')

axs[1].axvline(0, color='black')
axs[1].axhline(0, color='black')

axs[1].set_ylim([-.01, 1.01])
axs[1].set_xlim([x_min, x_max])
axs[1].legend(loc='best', prop={'size': 20})

axs[1].set_title(title, fontsize=20)
axs[1].set_xlabel(r'$n$', fontsize=20)
axs[1].grid()

plt.tight_layout()
plt.savefig('../figures/n_scanning.png', dpi=300)
plt.show()

fig.clear()
