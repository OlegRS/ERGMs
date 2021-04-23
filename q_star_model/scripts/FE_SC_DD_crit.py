################################################################################
# This script plots the free energy, t1-dependence and degree distributions    #
# of the corresponding MF3 and 3-star models at the critical point with L=0.3. #
################################################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

Ld = .3
fig, axs = plt.subplots(1, 3, figsize=(10.*1.4, 3.*1.5))

t2_MF = 1.8707482993197284
t3_MF = -0.7558578987150416
t1_MF = 1/2*np.log(Ld/(1-Ld)) - 2*t2_MF*Ld - 3*t3_MF*Ld*Ld


title = "$t_1=$" + str(round(t1_MF, 2)) + "$;\ t_2=$" + str(round(t2_MF,3)) + "$;\ t_3=$" + str(round(t3_MF,3))

L = np.arange(0,1, 1e-5)

F = -2*t1_MF*L - 2*t2_MF*L**2 - 2*t3_MF*L**3 + L*np.log(L/(1-L)) + np.log(1-L)
axs[0].plot(L, F, label="Specific MCFE", linewidth=2)
axs[0].plot(Ld, F[30000], marker='+', markersize=40, color='red', zorder=500)

axs[0].axvline(0, color='black')
axs[0].axhline(0, color='black')

#plt.ylim([-.03, 1.03])
axs[0].set_xlim([0, 1])
axs[0].set_ylim([-.15, .5])
axs[0].legend(loc='upper left', prop={'size': 13.5})

axs[0].set_title(title, fontsize=13)
axs[0].set_xlabel(r'$\bar L$', fontsize=20)
axs[0].set_ylabel(r'$\bar{\cal F}(\bar L;{\bf t})$', fontsize=20)

axs[0].grid()


#########################################################
################## t1 DEPENDENCE ########################
#########################################################
n_nodes = 1000

title = r"$n=10^3$" + "$;\ t_2=$" + str(round(t2_MF,3)) + "$;\ t_3=$" + str(round(t3_MF,3))

x_min, x_max =-2, .5
delta_t1 = 0.05
delta_L = 0.00001
t1_range = np.arange(x_min, x_max+delta_t1, delta_t1)
L_range = np.arange(0, 1+5*delta_L, delta_L)
L, T1 = np.meshgrid(L_range, t1_range)
SC = L - 1/2*(1 + np.tanh(T1 +2*t2_MF*L + 3*t3_MF*L**2)) #Self-consistency
CS = axs[1].contour(T1, L, SC, [0], colors='cyan', zorder=100, linewidths=1.5)
CS.collections[0].set_label('SC equation')

axs[1].set_xlabel(r'$t_1$', fontsize=20)
axs[1].set_ylabel(r'$\bar L$', fontsize=20)
axs[1].set_ylim([-.03, 1.03])
axs[1].set_xlim([x_min, x_max])
axs[1].axvline(0, color="black")
axs[1].axhline(0, color="black")
axs[1].plot(t1_MF, Ld, marker='+', markersize=40, color='red', zorder=500, alpha=.7)
axs[1].grid()

file_names = ["../data/p_star_model/t1_scanning/NEW_MATRIX_BASED_HAMILTONIAN__PS__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_-2__t1_fin_0.500000__t1_step_0.001000__t2_3.745242__t3_-4.548785__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv"]

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


axs[1].errorbar(t1_PS, L_avrg_PS, xerr=None, yerr=L_RMS_PS, color="blue", marker='o', markersize=4, ls='none', label="3-star", alpha=.5, elinewidth=2, capsize=2, zorder=51)

file_names = ["../data/mean_field_model/t1_scanning/NEW_METROPOLIS_MF__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_-2__t1_fin_0.500000__t1_step_0.001000__t2_3.745242__t3_-4.548785__N_DYNAMIC_PAIRS_499500__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv"]

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


axs[1].errorbar(t1_MF, L_avrg_MF, xerr=None, yerr=L_RMS_MF, color="orange", marker='s', markersize=4, ls='none', label="MF3", alpha=1, elinewidth=2, capsize=2, zorder=50)
axs[1].set_title(title, fontsize=13)
axs[1].legend(loc="lower right", prop={'size': 13.5})



#########################################################
############### DEGREE DISTRIBUTIONS ####################
#########################################################
n_nodes = 10000
t1_MF = 1/2*np.log(Ld/(1-Ld)) - 2*t2_MF*Ld - 3*t3_MF*Ld*Ld
title = r"$n=10^4$"
title += "$;\ t_1=$" + str(round(t1_MF, 2)) + "$;\ t_2=$" + str(round(t2_MF,3)) + "$;\ t_3=$" + str(round(t3_MF,3))

file_name = "../data/p_star_model/degree_distributions/MATRIX_BASED_PS__N_NODES_10000__N_ITERS_PER_LINK_10__t1_-1.342016__t2_3.741497__t3_-4.535147__N_AVRGNG_1000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_PS = []
for line in lines:
    deg_seq_PS.append(int(line.split('\n')[0]))
N_lines = np.array(lines).shape[0]
deg_seq_PS_last = deg_seq_PS[-n_nodes:]

L_PS = np.average(deg_seq_PS_last)
print(L_PS)


file_name = "../data/mean_field_model/degree_distributions/NEW_MF__N_NODES_10000__N_ITERS_PER_LINK_10__t1_-1.342016__t2_3.741497__t3_-4.535147__N_AVRGNG_1000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_MF = []
for lc in range(0, N_lines):
    deg_seq_MF.append(int(lines[lc].split('\n')[0]))

deg_seq_MF_last = deg_seq_MF[-n_nodes:]

L_MF = np.average(deg_seq_MF_last)
print(L_MF)


axs[2].set_title(title, fontsize=13)
# Theoretical
L_MF /= (n_nodes-1)
sigma = np.sqrt((n_nodes-1)*L_MF*(1-L_MF))
x_min, x_max = L_MF*(n_nodes-1) - 15*sigma, L_MF*(n_nodes-1) + 5*sigma
degrees = np.arange(np.rint(x_min), np.rint(x_max)+1)
p = []
p = np.sqrt(1/(2*np.pi*(n_nodes-1)*Ld*(1-Ld)))*np.exp(-(degrees-Ld*(n_nodes-1))**2/(2*(n_nodes-1)*Ld*(1-Ld)))
axs[2].plot(degrees, p, color='r', label="MF theory")
axs[2].hist(deg_seq_PS[-n_nodes:], bins = degrees, color='b', density=True, label="3-star")
axs[2].hist(deg_seq_MF[-n_nodes:], bins = degrees, color='orange', density=True, alpha=.65, label="MF3")
axs[2].set_xlabel(r"$k$", fontsize=20)
axs[2].set_ylabel(r"$\mathrm{p}(k|{\bf A})$", fontsize=20)


axs[2].grid()
axs[2].set_xlim([2450, 3190])
axs[2].legend(loc='best', prop={'size': 13.5})


plt.tight_layout()
plt.savefig("../figures/FE_SC_DD_crit.png", dpi=300)

plt.show()

