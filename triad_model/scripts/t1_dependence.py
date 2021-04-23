##########################################################
# This script plots the dependencies of connectance and  #
# average local clustering coefficients on t1 for the    #
# graphs of different sizes.                             #
##########################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

n_nodes = 50
N_links_max = n_nodes*(n_nodes-1)/2

t2_MF = -75
t3_MF = 20

fig, axs = plt.subplots(2, 3, figsize=(10.*1.4, 6.*1.4))

file_names = ["../data/triad_model/t1_scanning/TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_50__N_ITERS_PER_LINK_10__t1_start_-1__t1_fin_101__t1_step_0.300000__t2_MF_-75__t3_MF_20__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv"]

title = "$t_2=$" + str(t2_MF) + "$,\ t_3=$" + str(t3_MF) + "$;\ n=$" + str(n_nodes)

x_min = 0
x_max = 100

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
    axs[0,0].plot(t1, L_avrg, 'h',markersize=4, label=r"$\bar L$")
    axs[0,0].plot(t1, C_avrg, 'h',markersize=4, label=r"$\bar c$")

## THEORETICAL ##
# Plotting self-consistency equation
delta_t1 = 0.05
delta_L = 0.00005
t1_range = np.arange(x_min, x_max + delta_t1, delta_t1)
L_range = np.arange(0, 1+delta_L*2, delta_L)
L, T1 = np.meshgrid(L_range, t1_range)
SC = L - 1/2*(1 + np.tanh(T1 + 2*t2_MF*L + 3*t3_MF*L**2)) #Self-consistency
CS = axs[0,0].contour(T1, L, SC, [0], colors='cyan', zorder=100, linewidths=1.5)

CS.collections[0].set_label('SC')

axs[0,0].axvline(0, color='black')
axs[0,0].axhline(0, color='black')

axs[0,0].set_ylim([-.01, 1.01])
axs[0,0].set_xlim([x_min, x_max])
axs[0,0].legend(loc='lower right', prop={'size': 16.5})

axs[0,0].set_title(title, fontsize=18)
axs[0,0].set_ylabel(r'$\bar L,\ \bar c$', fontsize=20)
axs[0,0].grid()

##########################################################
n_nodes = 100

file_names = ["../data/triad_model/t1_scanning/TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_100__N_ITERS_PER_LINK_10__t1_start_-1__t1_fin_101__t1_step_0.300000__t2_MF_-75__t3_MF_20__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv"]


title = "$t_2=$" + str(t2_MF) + "$,\ t_3=$" + str(t3_MF) + "$;\ n=$" + str(n_nodes)

x_min = 0
x_max = 100

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
    axs[0,1].plot(t1, L_avrg, 'h',markersize=4, label=r"$\bar L$")
    axs[0,1].plot(t1, C_avrg, 'h',markersize=4, label=r"$\bar c$")

## THEORETICAL ##
# Plotting self-consistency equation
delta_t1 = 0.05
delta_L = 0.00005
t1_range = np.arange(x_min, x_max + delta_t1, delta_t1)
L_range = np.arange(0, 1+delta_L*2, delta_L)
L, T1 = np.meshgrid(L_range, t1_range)
SC = L - 1/2*(1 + np.tanh(T1 + 2*t2_MF*L + 3*t3_MF*L**2)) #Self-consistency
CS = axs[0,1].contour(T1, L, SC, [0], colors='cyan', zorder=100, linewidths=1.5)

CS.collections[0].set_label('SC')

axs[0,1].axvline(0, color='black')
axs[0,1].axhline(0, color='black')

axs[0,1].set_ylim([-.01, 1.01])
axs[0,1].set_xlim([x_min, x_max])
axs[0,1].legend(loc='lower right', prop={'size': 16.5})

axs[0,1].set_title(title, fontsize=18)
axs[0,1].grid()

##########################################################
n_nodes = 200

file_names = ["../data/triad_model/t1_scanning/TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_200__N_ITERS_PER_LINK_10__t1_start_-1__t1_fin_101__t1_step_0.300000__t2_MF_-75__t3_MF_20__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv"]

title = "$t_2=$" + str(t2_MF) + "$,\ t_3=$" + str(t3_MF) + "$;\ n=$" + str(n_nodes)

x_min = 0
x_max = 100

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
    axs[0,2].plot(t1, L_avrg, 'h',markersize=4, label=r"$\bar L$")
    axs[0,2].plot(t1, C_avrg, 'h',markersize=4, label=r"$\bar c$")

## THEORETICAL ##
# Plotting self-consistency equation
delta_t1 = 0.05
delta_L = 0.00005
t1_range = np.arange(x_min, x_max + delta_t1, delta_t1)
L_range = np.arange(0, 1+delta_L*2, delta_L)
L, T1 = np.meshgrid(L_range, t1_range)
SC = L - 1/2*(1 + np.tanh(T1 + 2*t2_MF*L + 3*t3_MF*L**2)) #Self-consistency
CS = axs[0,2].contour(T1, L, SC, [0], colors='cyan', zorder=100, linewidths=1.5)

CS.collections[0].set_label('SC')

axs[0,2].axvline(0, color='black')
axs[0,2].axhline(0, color='black')

axs[0,2].set_ylim([-.01, 1.01])
axs[0,2].set_xlim([x_min, x_max])
axs[0,2].legend(loc='lower right', prop={'size': 16.5})

axs[0,2].set_title(title, fontsize=18)
axs[0,2].grid()

##########################################################
n_nodes = 400

file_names = ["../data/triad_model/t1_scanning/TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_400__N_ITERS_PER_LINK_10__t1_start_-1__t1_fin_101__t1_step_0.300000__t2_MF_-75__t3_MF_20__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv"]

title = "$t_2=$" + str(t2_MF) + "$,\ t_3=$" + str(t3_MF) + "$;\ n=$" + str(n_nodes)

x_min = 0
x_max = 100

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
    axs[1,0].plot(t1, L_avrg, 'h',markersize=4, label=r"$\bar L$")
    axs[1,0].plot(t1, C_avrg, 'h',markersize=4, label=r"$\bar c$")

## THEORETICAL ##
# Plotting self-consistency equation
delta_t1 = 0.05
delta_L = 0.00005
t1_range = np.arange(x_min, x_max + delta_t1, delta_t1)
L_range = np.arange(0, 1+delta_L*2, delta_L)
L, T1 = np.meshgrid(L_range, t1_range)
SC = L - 1/2*(1 + np.tanh(T1 + 2*t2_MF*L + 3*t3_MF*L**2)) #Self-consistency
CS = axs[1,0].contour(T1, L, SC, [0], colors='cyan', zorder=100, linewidths=1.5)

CS.collections[0].set_label('SC')

axs[1,0].axvline(0, color='black')
axs[1,0].axhline(0, color='black')

axs[1,0].set_ylim([-.01, 1.01])
axs[1,0].set_xlim([x_min, x_max])
axs[1,0].legend(loc='lower right', prop={'size': 16.5})

axs[1,0].set_title(title, fontsize=18)
axs[1,0].set_xlabel(r'$t_1$', fontsize=20)
axs[1,0].set_ylabel(r'$\bar L,\ \bar c$', fontsize=20)
axs[1,0].grid()

##########################################################
n_nodes = 800

file_names = ["../data/triad_model/t1_scanning/TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_800__N_ITERS_PER_LINK_10__t1_start_-1__t1_fin_101__t1_step_0.300000__t2_MF_-75__t3_MF_20__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv"]

title = "$t_2=$" + str(t2_MF) + "$,\ t_3=$" + str(t3_MF) + "$;\ n=$" + str(n_nodes)

x_min = 0
x_max = 100

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
    axs[1,1].plot(t1, L_avrg, 'h',markersize=4, label=r"$\bar L$")
    axs[1,1].plot(t1, C_avrg, 'h',markersize=4, label=r"$\bar c$")

## THEORETICAL ##
# Plotting self-consistency equation
delta_t1 = 0.05
delta_L = 0.00005
t1_range = np.arange(x_min, x_max + delta_t1, delta_t1)
L_range = np.arange(0, 1+delta_L*2, delta_L)
L, T1 = np.meshgrid(L_range, t1_range)
SC = L - 1/2*(1 + np.tanh(T1 + 2*t2_MF*L + 3*t3_MF*L**2)) #Self-consistency
CS = axs[1,1].contour(T1, L, SC, [0], colors='cyan', zorder=100, linewidths=1.5)

CS.collections[0].set_label('SC')

axs[1,1].axvline(0, color='black')
axs[1,1].axhline(0, color='black')

axs[1,1].set_ylim([-.01, 1.01])
axs[1,1].set_xlim([x_min, x_max])
axs[1,1].legend(loc='lower right', prop={'size': 16.5})

axs[1,1].set_title(title, fontsize=18)
axs[1,1].set_xlabel(r'$t_1$', fontsize=20)
axs[1,1].grid()

##########################################################
n_nodes = 1600

file_names = ["../data/triad_model/t1_scanning/TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_1600__N_ITERS_PER_LINK_10__t1_start_-1__t1_fin_101__t1_step_0.300000__t2_MF_-75__t3_MF_20__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv"]

title = "$t_2=$" + str(t2_MF) + "$,\ t_3=$" + str(t3_MF) + "$;\ n=$" + str(n_nodes)

x_min = 0
x_max = 100

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
    axs[1,2].plot(t1, L_avrg, 'h',markersize=4, label=r"$\bar L$")
    axs[1,2].plot(t1, C_avrg, 'h',markersize=4, label=r"$\bar c$")

## THEORETICAL ##
# Plotting self-consistency equation
delta_t1 = 0.05
delta_L = 0.00005
t1_range = np.arange(x_min, x_max + delta_t1, delta_t1)
L_range = np.arange(0, 1+delta_L*2, delta_L)
L, T1 = np.meshgrid(L_range, t1_range)
SC = L - 1/2*(1 + np.tanh(T1 + 2*t2_MF*L + 3*t3_MF*L**2)) #Self-consistency
CS = axs[1,2].contour(T1, L, SC, [0], colors='cyan', zorder=100, linewidths=1.5)

CS.collections[0].set_label('SC')

axs[1,2].axvline(0, color='black')
axs[1,2].axhline(0, color='black')

axs[1,2].plot(39.176351, .3, marker='+', markersize=40, markeredgewidth=1.5, color='green', zorder=500, alpha=.7)

axs[1,2].set_ylim([-.01, 1.01])
axs[1,2].set_xlim([x_min, x_max])
axs[1,2].legend(loc='lower right', prop={'size': 16.5})

axs[1,2].set_title(title, fontsize=18)
axs[1,2].set_xlabel(r'$t_1$', fontsize=20)
axs[1,2].grid()

##########################################################

plt.tight_layout()
plt.savefig('../figures/t1_scanning_joint.png', dpi=300)
plt.show()

fig.clear()
