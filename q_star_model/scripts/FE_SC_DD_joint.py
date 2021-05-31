#########################################################
# This script plots what happens with t1-hysteresis and #
# degree distributions with decreasing temperature.     #
#########################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

Ld = .3
n_nodes = 10000

fig, axs = plt.subplots(3, 3, figsize=(10.*1.4, 6.*1.4))

#########################################################
################### Free Energy #########################
#########################################################
t2_MF = -.75
t3_MF = .45
t1_MF = 1/2*np.log(Ld/(1-Ld)) - 2*t2_MF*Ld - 3*t3_MF*Ld*Ld #2.7*10

title = "$t_1=$" + str(round(t1_MF, 4)) + "$;\ t_2=$" + str(t2_MF) + "$;\ t_3=$" + str(t3_MF)

L = np.arange(1e-5,1, 1e-5)

F = -2*t1_MF*L - 2*t2_MF*L**2 - 2*t3_MF*L**3 + L*np.log(L/(1-L)) + np.log(1-L)
axs[0,0].plot(L, F, label="SMFE")

axs[0,0].axvline(0, color='black')
axs[0,0].axhline(0, color='black')

axs[0,0].set_xlim([0, 1])
axs[0,0].legend(loc='upper left', prop={'size': 10})

axs[0,0].set_title(title, fontsize=12)
axs[0,0].set_ylabel(r'$\bar {\cal F}(\bar L ;{\bf t})$', fontsize=17)

axs[0,0].grid()

#########################################################
################## t1 DEPENDENCE ########################
#########################################################
title = r"$n=10^3$" + "$;\ t_2=$" + str(t2_MF) + "$;\ t_3=$" + str(t3_MF)

x_min, x_max =-2.5, 2.5
delta_t1 = 0.05
delta_L = 0.00001
t1_range = np.arange(x_min, x_max+delta_t1, delta_t1)
L_range = np.arange(0, 1+2*delta_L, delta_L)
L, T1 = np.meshgrid(L_range, t1_range)
SC = L - 1/2*(1 + np.tanh(T1 +2*t2_MF*L + 3*t3_MF*L**2)) #Self-consistency
CS = axs[0,1].contour(T1, L, SC, [0], colors='cyan', zorder=100, linewidths=1.5)
CS.collections[0].set_label('SC equation')

axs[0,1].set_ylabel(r'$\bar L$', fontsize=17)
axs[0,1].set_xlim([x_min, x_max])
axs[0,1].set_ylim([-.03, 1.03])
axs[0,1].axvline(0, color="black")
axs[0,1].axhline(0, color="black")
axs[0,1].plot(t1_MF, Ld, marker='+', markersize=50, color='green', zorder=500, alpha=.7)
axs[0,1].grid()

file_names = ["../data/p_star_model/t1_scanning/NEW_MATRIX_BASED_HAMILTONIAN__PS__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_2.500000__t1_fin_-2.500000__t1_step_-0.005000__t2_-1.501502__t3_2.708119__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv",
              "../data/p_star_model/t1_scanning/NEW_MATRIX_BASED_HAMILTONIAN__PS__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_-2.500000__t1_fin_2.500000__t1_step_0.005000__t2_-1.501502__t3_2.708119__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv"
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

axs[0,1].errorbar(t1_PS, L_avrg_PS, xerr=None, yerr=L_RMS_PS, color="blue", marker='o', markersize=4, ls='none', label="3-star simulation", alpha=.5, elinewidth=2, capsize=2, zorder=51)

file_names = ["../data/mean_field_model/t1_scanning/NEW_METROPOLIS_MF__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_2.500000__t1_fin_-2.500000__t1_step_-0.005000__t2_-1.501502__t3_2.708119__N_DYNAMIC_PAIRS_499500__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv",
              "../data/mean_field_model/t1_scanning/NEW_METROPOLIS_MF__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_-2.500000__t1_fin_2.500000__t1_step_0.005000__t2_-1.501502__t3_2.708119__N_DYNAMIC_PAIRS_499500__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv"
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

axs[0,1].errorbar(t1_MF, L_avrg_MF, xerr=None, yerr=L_RMS_MF, color="orange", marker='s', markersize=4, ls='none', label="MF3 simulation", alpha=1, elinewidth=2, capsize=2, zorder=50)

handles, labels = axs[0,1].get_legend_handles_labels()
axs[0,1].legend([handles[2], handles[1], handles[0]], [labels[2],labels[1],labels[0]], loc="upper left", prop={'size': 10})
# axs[0,1].legend(loc="upper left", prop={'size': 10})
axs[0,1].set_title(title, fontsize=12)

#########################################################
############### DEGREE DISTRIBUTIONS ####################
#########################################################
t1_MF = 1/2*np.log(Ld/(1-Ld)) - 2*t2_MF*Ld - 3*t3_MF*Ld*Ld

title = r"$n=10^4$"
title += "$;\ t_1=$" + str(round(t1_MF, 4)) + "$;\ t_2=$" + str(t2_MF) + "$;\ t_3=$" + str(t3_MF)
file_name = "../data/p_star_model/degree_distributions/MATRIX_BASED_PS__N_NODES_10000__N_ITERS_PER_LINK_10__t1_-0.752149__t2_1.500000__t3_-2.700000__N_AVRGNG_1000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_PS = []
for line in lines:
    deg_seq_PS.append(int(line.split('\n')[0]))
N_lines = np.array(lines).shape[0]
L_PS = np.average(deg_seq_PS)
print(L_PS)


file_name = "../data/mean_field_model/degree_distributions/NEW_MF__N_NODES_10000__N_ITERS_PER_LINK_10__t1_-0.752149__t2_1.500000__t3_-2.700000__N_AVRGNG_1000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_MF = []
for lc in range(0, N_lines):
    deg_seq_MF.append(int(lines[lc].split('\n')[0]))

L_MF = np.average(deg_seq_MF)
print(L_MF)


# Theoretical
n_nodes = 10000
L_MF /= (n_nodes-1)
sigma = np.sqrt((n_nodes-1)*L_MF*(1-L_MF))
x_min, x_max = L_MF*(n_nodes-1) - 3*sigma, L_MF*(n_nodes-1) + 3*sigma
degrees = np.arange(np.rint(x_min), np.rint(x_max)+1)
p = []
p = np.sqrt(1/(2*np.pi*(n_nodes-1)*L_MF*(1-L_MF)))*np.exp(-(degrees-L_MF*(n_nodes-1))**2/(2*(n_nodes-1)*L_MF*(1-L_MF)))
axs[0,2].plot(degrees, p, color='r', label="MF theor")
axs[0,2].hist(deg_seq_PS, bins = degrees, color='b', density=True, label="3-star sim")
axs[0,2].hist(deg_seq_MF, bins = degrees, color='orange', density=True, alpha=.65, label="MF3 sim")

axs[0,2].set_ylabel(r"$\langle \mathrm{p}(k|{\bf A})\rangle$", fontsize=17)

axs[0,2].grid()
axs[0,2].set_xlim([x_min, x_max])
axs[0,2].legend(loc=[.005,.66], prop={'size': 10})
axs[0,2].set_title(title, fontsize=12)

#########################################################
############### HIGHER COUPLINGS FE #####################
#########################################################

t2_MF = -7.5
t3_MF = 4.5
t1_MF = 1/2*np.log(Ld/(1-Ld)) - 2*t2_MF*Ld - 3*t3_MF*Ld*Ld

title = "$t_1=$" + str(round(t1_MF, 3)) + "$;\ t_2=$" + str(t2_MF) + "$;\ t_3=$" + str(t3_MF)

L = np.arange(1e-5,1, 1e-5)

F = -2*t1_MF*L - 2*t2_MF*L**2 - 2*t3_MF*L**3 + L*np.log(L/(1-L)) + np.log(1-L)
axs[1,0].plot(L, F, label="SMFE")

axs[1,0].axvline(0, color='black')
axs[1,0].axhline(0, color='black')

axs[1,0].set_xlim([0, 1])
axs[1,0].legend(loc='upper left', prop={'size': 10})

axs[1,0].set_title(title, fontsize=12)
axs[1,0].set_ylabel(r'$\bar{\cal F}(\bar L;{\bf t})$', fontsize=17)

axs[1,0].grid()

#########################################################
################## t1 DEPENDENCE ########################
#########################################################
title = r"$n=10^3$" + "$;\ t_2=$" + str(t2_MF) + "$;\ t_3=$" + str(t3_MF)

x_min, x_max =-1, 7
delta_t1 = 0.05
delta_L = 0.00001
t1_range = np.arange(x_min, x_max+delta_t1, delta_t1)
L_range = np.arange(0, 1+2*delta_L, delta_L)
L, T1 = np.meshgrid(L_range, t1_range)
SC = L - 1/2*(1 + np.tanh(T1 +2*t2_MF*L + 3*t3_MF*L**2)) #Self-consistency
CS = axs[1,1].contour(T1, L, SC, [0], colors='cyan', zorder=100, linewidths=1.5)
CS.collections[0].set_label('SC equation')

axs[1,1].set_ylabel(r'$\bar L$', fontsize=17)
axs[1,1].set_xlim([x_min, x_max])
axs[1,1].set_ylim([-.03, 1.03])
axs[1,1].axvline(0, color="black")
axs[1,1].axhline(0, color="black")
axs[1,1].plot(t1_MF, Ld, marker='+', markersize=50, color='green', zorder=500, alpha=.7)
axs[1,1].grid()


file_names = ["../data/p_star_model/t1_scanning/NEW_MATRIX_BASED_HAMILTONIAN__PS__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_-1__t1_fin_7__t1_step_0.005000__t2_-15.015015__t3_27.081189__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv",
              "../data/p_star_model/t1_scanning/NEW_MATRIX_BASED_HAMILTONIAN__PS__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_7__t1_fin_-1__t1_step_-0.005000__t2_-15.015015__t3_27.081189__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv"
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

ind = t1_PS.index(4.36) # Removing transition point
del t1_PS[ind]
del L_avrg_PS[ind]
del L_RMS_PS[ind]

axs[1,1].errorbar(t1_PS, L_avrg_PS, xerr=None, yerr=L_RMS_PS, color="blue", marker='o', markersize=4, ls='none', label="3-star simulation", alpha=.5, elinewidth=2, capsize=2, zorder=51)

file_names = ["../data/mean_field_model/t1_scanning/NEW_METROPOLIS_MF__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_-1__t1_fin_7__t1_step_0.005000__t2_-15.015015__t3_27.081189__N_DYNAMIC_PAIRS_499500__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv",
              "../data/mean_field_model/t1_scanning/NEW_METROPOLIS_MF__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_7__t1_fin_-1__t1_step_-0.005000__t2_-15.015015__t3_27.081189__N_DYNAMIC_PAIRS_499500__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv"
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

ind = t1_MF.index(4.36) # Removing transition point
del t1_MF[ind]
del L_avrg_MF[ind]
del L_RMS_MF[ind]

axs[1,1].errorbar(t1_MF, L_avrg_MF, xerr=None, yerr=L_RMS_MF, color="orange", marker='s', markersize=4, ls='none', label="MF3 simulation", alpha=1, elinewidth=2, capsize=2, zorder=50)
axs[1,1].set_title(title, fontsize=12)

##### Exact thermodynamic solution #####
x_min_, x_max_ = 3.8, 4.3
T1_MF = np.arange(x_min_, x_max_, 0.0001)
L_MIN = []
for t1_MF in T1_MF:
    def SC_equation(L):
        return L - 1/2*(1 + np.tanh(t1_MF + 2*t2_MF*L + 3*t3_MF*L**2))
    L_ = np.array([fsolve(SC_equation, .35)[0], fsolve(SC_equation, .99)[0]])
    def F(L):
        return -2*(t1_MF*L + t2_MF*L**2 + t3_MF*L**3) + L*np.log(L/(1-L)) + np.log(1-L)
    L0 = L_[~(np.triu(np.abs(L_[:,None] - L_) < 0.001,1)).any(0)] # Distinct solution of SC equation
    L = np.array([l for l in L0 if (l<1 and l>0)])        
    L_min = L_[list(F(L_)).index(min(F(L_)))]
    L_MIN.append(L_min)
    
axs[1,1].plot(T1_MF, L_MIN, color="red", label='Phase transition')
handles, labels = axs[1,1].get_legend_handles_labels()
axs[1,1].legend([handles[2], handles[3], handles[0], handles[1]], [labels[2],labels[3],labels[0],labels[1]], loc="upper left", prop={'size': 10})


#########################################################
############### DEGREE DISTRIBUTIONS ####################
#########################################################
t1_MF = 1/2*np.log(Ld/(1-Ld)) - 2*t2_MF*Ld - 3*t3_MF*Ld*Ld
title = r"$n=10^4$"
title += "$;\ t_1=$" + str(round(t1_MF, 3)) + "$;\ t_2=$" + str(t2_MF) + "$;\ t_3=$" + str(t3_MF)

file_name = "../data/p_star_model/degree_distributions/MATRIX_BASED_PS__N_NODES_10000__N_ITERS_PER_LINK_10__t1_-0.095149__t2_-1.500000__t3_2.700000__N_AVRGNG_1000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_PS = []
for line in lines:
    deg_seq_PS.append(int(line.split('\n')[0]))
N_lines = np.array(lines).shape[0]
L_PS = np.average(deg_seq_PS)
print(L_PS)

file_name = "../data/mean_field_model/degree_distributions/NEW_MF__N_NODES_10000__N_ITERS_PER_LINK_10__t1_-0.095149__t2_-1.500000__t3_2.700000__N_AVRGNG_1000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_MF = []
for lc in range(0, N_lines):
    deg_seq_MF.append(int(lines[lc].split('\n')[0]))

L_MF = np.average(deg_seq_MF)
print(L_MF)

axs[1,2].set_title(title, fontsize=12)
# Theoretical
n_nodes = 10000
L_MF /= (n_nodes-1)
sigma = np.sqrt((n_nodes-1)*L_MF*(1-L_MF))
x_min, x_max = L_MF*(n_nodes-1) - 3*sigma, L_MF*(n_nodes-1) + 3*sigma
degrees = np.arange(np.rint(x_min), np.rint(x_max)+1)
p = []
p = np.sqrt(1/(2*np.pi*(n_nodes-1)*L_MF*(1-L_MF)))*np.exp(-(degrees-L_MF*(n_nodes-1))**2/(2*(n_nodes-1)*L_MF*(1-L_MF)))
axs[1,2].plot(degrees, p, color='r', label="MF theor")
axs[1,2].hist(deg_seq_PS, bins = degrees, color='b', density=True, label="3-star sim")
axs[1,2].hist(deg_seq_MF, bins = degrees, color='orange', density=True, alpha=.65, label="MF3 sim")

axs[1,2].set_ylabel(r"$\langle \mathrm{p}(k|{\bf A})\rangle$", fontsize=17)

axs[1,2].grid()
axs[1,2].set_xlim([x_min, x_max])
axs[1,2].legend(loc=[.005,.66], prop={'size': 10})


#########################################################
############### HIGHER COUPLINGS FE #####################
#########################################################

t2_MF = -75
t3_MF = 45
t1_MF = 1/2*np.log(Ld/(1-Ld)) - 2*t2_MF*Ld - 3*t3_MF*Ld*Ld

title = "$t_1=$" + str(round(t1_MF, 3)) + "$;\ t_2=$" + str(t2_MF) + "$;\ t_3=$" + str(t3_MF)

L = np.arange(1e-5,1, 1e-5)

F = -2*t1_MF*L - 2*t2_MF*L**2 - 2*t3_MF*L**3 + L*np.log(L/(1-L)) + np.log(1-L)
axs[2,0].plot(L, F, label="SMFE")

axs[2,0].axvline(0, color='black')
axs[2,0].axhline(0, color='black')

axs[2,0].set_xlim([0, 1])
axs[2,0].legend(loc='upper left', prop={'size': 10})

axs[2,0].set_title(title, fontsize=12)
axs[2,0].set_xlabel(r'$\bar L$', fontsize=17)
axs[2,0].set_ylabel(r'$\bar{\cal F}(\bar L;{\bf t})$', fontsize=17)

axs[2,0].grid()


#########################################################
################## t1 DEPENDENCE ########################
#########################################################
title = r"$n=10^3$" + "$;\ t_2=$" + str(t2_MF) + "$;\ t_3=$" + str(t3_MF)

x_min, x_max =-1, 50
delta_t1 = 0.05
delta_L = 0.00001
t1_range = np.arange(x_min, x_max+delta_t1, delta_t1)
L_range = np.arange(0, 1+2*delta_L, delta_L)
L, T1 = np.meshgrid(L_range, t1_range)
SC = L - 1/2*(1 + np.tanh(T1 +2*t2_MF*L + 3*t3_MF*L**2)) #Self-consistency
CS = axs[2,1].contour(T1, L, SC, [0], colors='cyan', zorder=100, linewidths=1.5)
CS.collections[0].set_label('SC equation')

axs[2,1].set_xlabel(r'$t_1$', fontsize=17)
axs[2,1].set_ylabel(r'$\bar L$', fontsize=17)
axs[2,1].set_xlim([x_min, x_max])
axs[2,1].set_ylim([-.03, 1.03])
axs[2,1].axvline(0, color="black")
axs[2,1].axhline(0, color="black")
axs[2,1].plot(t1_MF, Ld, marker='+', markersize=50, color='green', zorder=500, alpha=.7)
axs[2,1].grid()

file_names = ["../data/p_star_model/t1_scanning/NEW_MATRIX_BASED_HAMILTONIAN__PS__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_-1__t1_fin_50__t1_step_0.050000__t2_-150.150150__t3_270.811894__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv",
              "../data/p_star_model/t1_scanning/NEW_MATRIX_BASED_HAMILTONIAN__PS__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_50__t1_fin_-1__t1_step_-0.050000__t2_-150.150150__t3_270.811894__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv"
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


axs[2,1].errorbar(t1_PS, L_avrg_PS, xerr=None, yerr=L_RMS_PS, color="blue", marker='o', markersize=4, ls='none', label="3-star simulation", alpha=.5, elinewidth=2, capsize=2, zorder=51)

file_names = ["../data/mean_field_model/t1_scanning/NEW_METROPOLIS_MF__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_50__t1_fin_-1__t1_step_-0.050000__t2_-150.150150__t3_270.811894__N_DYNAMIC_PAIRS_499500__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv",
              "../data/mean_field_model/t1_scanning/NEW_METROPOLIS_MF__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_-1__t1_fin_50__t1_step_0.050000__t2_-150.150150__t3_270.811894__N_DYNAMIC_PAIRS_499500__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv"
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

ind = t1_MF.index(18.2)  # Removing transition point
del t1_MF[ind]
del L_avrg_MF[ind]
del L_RMS_MF[ind]

axs[2,1].errorbar(t1_MF, L_avrg_MF, xerr=None, yerr=L_RMS_MF, color="orange", marker='s', markersize=4, ls='none', label="MF3 simulation", alpha=1, elinewidth=2, capsize=2, zorder=50)
axs[2,1].set_title(title, fontsize=12)

##### Exact thermodynamic solution #####
x_min_, x_max_ = 35, 38 #35.6261
T1_MF = np.arange(x_min_, x_max_, 0.0001)
L_MIN = []
for t1_MF in T1_MF:
    def SC_equation(L):
        return L - 1/2*(1 + np.tanh(t1_MF + 2*t2_MF*L + 3*t3_MF*L**2))
    L_ = np.array([fsolve(SC_equation, .3)[0], fsolve(SC_equation, .99)[0]])
    def F(L):
        return -2*(t1_MF*L + t2_MF*L**2 + t3_MF*L**3) + L*np.log(L/(1-L)) + np.log(1-L)
    L0 = L_[~(np.triu(np.abs(L_[:,None] - L_) < 0.001,1)).any(0)] # Distinct solution of SC equation
    L = np.array([l for l in L0 if (l<1 and l>0)])        
    L_min = L_[list(F(L_-.00001)).index(min(F(L_-.00001)))]
    L_MIN.append(L_min)

axs[2,1].plot(T1_MF, L_MIN, color="red", label='Phase transition')
handles, labels = axs[2,1].get_legend_handles_labels()
axs[2,1].legend([handles[2], handles[3], handles[0], handles[1]], [labels[2],labels[3],labels[0],labels[1]], loc=[.025,.44], prop={'size': 10})


#########################################################
############### DEGREE DISTRIBUTIONS ####################
#########################################################
t1_MF = 1/2*np.log(Ld/(1-Ld)) - 2*t2_MF*Ld - 3*t3_MF*Ld*Ld
title = r"$n=10^4$"
title += "$;\ t_1=$" + str(round(t1_MF, 3)) + "$;\ t_2=$" + str(t2_MF) + "$;\ t_3=$" + str(t3_MF)

file_name = "../data/p_star_model/degree_distributions/MATRIX_BASED_PS__N_NODES_10000__N_ITERS_PER_LINK_10__t1_32.426351__t2_-150__t3_270__N_AVRGNG_1000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_PS = []
for line in lines:
    deg_seq_PS.append(int(line.split('\n')[0]))
N_lines = np.array(lines).shape[0]
L_PS = np.average(deg_seq_PS)
print(L_PS)


file_name = "../data/mean_field_model/degree_distributions/NEW_MF__N_NODES_10000__N_ITERS_PER_LINK_10__t1_32.426351__t2_-150__t3_270__N_AVRGNG_1000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_MF = []
for lc in range(0, N_lines):
    deg_seq_MF.append(int(lines[lc].split('\n')[0]))

L_MF = np.average(deg_seq_MF)
print(L_MF)


axs[2,2].set_title(title, fontsize=12)
# Theoretical
n_nodes = 10000
L_MF /= (n_nodes-1)
sigma = np.sqrt((n_nodes-1)*L_MF*(1-L_MF))
x_min, x_max = L_MF*(n_nodes-1) - 3*sigma, L_MF*(n_nodes-1) + 3*sigma
degrees = np.arange(np.rint(x_min), np.rint(x_max)+1)
p = []
p = np.sqrt(1/(2*np.pi*(n_nodes-1)*L_MF*(1-L_MF)))*np.exp(-(degrees-L_MF*(n_nodes-1))**2/(2*(n_nodes-1)*L_MF*(1-L_MF)))
axs[2,2].plot(degrees, p, color='r', label="MF theor")
axs[2,2].hist(deg_seq_PS, bins = degrees, color='b', density=True, label="3-star sim")
axs[2,2].hist(deg_seq_MF, bins = degrees, color='orange', density=True, alpha=.65, label="MF3 sim")

axs[2,2].set_ylabel(r"$\langle \mathrm{p}(k|{\bf A})\rangle$", fontsize=17)

axs[2,2].grid()
axs[2,2].set_xlim([x_min, x_max])
axs[2,2].legend(loc=[.005,.66], prop={'size': 10})
axs[2,2].set_xlabel(r"$k$", fontsize=17)

plt.tight_layout()
plt.savefig("../figures/FE_SC_DD_joint.png", dpi=300)

plt.show()
