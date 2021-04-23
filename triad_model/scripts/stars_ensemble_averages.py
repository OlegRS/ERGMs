################################################################################
# This script plots ensemble averages of 2-star and triangle counts for small  #
# networks together with the corresponding numerical solutions of the Burgers' #
# hierarchy and (0,1)-bistability transition formula.                          #
################################################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

fig, axs = plt.subplots(1, 2, figsize=(10.*1.2, 3.*1.2))

#########################################
############ SIMULATIONS ################
#########################################
t2_MF = 2
t3_MF = 1
x_min = -3.3
x_max = -2.7

###################### 2_STARS ############################
n_nodes=5
n_pairs = n_nodes*(n_nodes-1)/2

file_names = [ "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_5__N_ITERS_PER_LINK_100__t1_start_-3.500000__t1_fin_-2.500000__t1_step_0.010000__t2_MF_2__t3_MF_1__N_AVRGNG_50000__ENSEMBLE_AVRG_0.csv"
]

t2 = t2_MF * 2*n_nodes/(n_nodes-2) # t2_MF*(CONVERSION_FACTOR)
t3 = t3_MF * 6*n_nodes**2/(n_nodes**2-5*n_nodes+6) # t3_MF*(CONVERSION_FACTOR)

t1_PS = []
S2_avrg_PS = []
S2_RMS_PS = []
S2_err_PS_squared = []

for file_name in file_names:
    S2_vs_t1 = np.genfromtxt(file_name, delimiter=',')
    # Finding averages
    S2_avrg = []
    S2_RMS = []
    t1 = []
    i = 0
    N_0, N_1 = [], [] # Numbers of data points close to L=0 and L=1
    S2_var_0, S2_var_1 = [], [] # Fluctuations around L=0 and L=1
    while i < S2_vs_t1.shape[0]:
        prev_t1 = S2_vs_t1[i][0]
        S2 = []
        S2_0, S2_1 = [], [] # Data points close to L=0 and L=1
        n_0, n_1 = 0, 0
        while i < S2_vs_t1.shape[0] and prev_t1 == S2_vs_t1[i][0]:
            S2.append(S2_vs_t1[i][2]*2/((n_nodes-1)*(n_nodes-2)))
            if(S2[-1] < 0.5):
                S2_0.append(S2[-1])
                n_0+=1
            else:
                S2_1.append(S2[-1])
                n_1+=1
            i+=1
        N_0.append(n_0)
        N_1.append(n_1)
        S2_avrg.append(np.array(S2).mean())
        S2_RMS.append(np.sqrt(np.array(S2).var()))
        if(S2_0!=[]):
            S2_var_0.append(np.array(S2_0).var())
        else:
            S2_var_0.append(0)
        if(S2_1!=[]):
            S2_var_1.append(np.array(S2_1).var())
        else:
            S2_var_1.append(0)
        t1.append(prev_t1)
    t1_PS.extend(t1)
    S2_avrg_PS.extend(S2_avrg)
    S2_RMS_PS.extend(S2_RMS)
    S2_err_PS_squared.extend((np.array(S2_var_0)*np.array(N_0) + np.array(S2_var_1)*np.array(N_1))/(np.array(N_0)+np.array(N_1)))

axs[0].errorbar(t1_PS, S2_avrg_PS, xerr=None, yerr=[np.minimum(np.sqrt(S2_err_PS_squared+(np.array(S2_avrg_PS)*(1-np.array(S2_avrg_PS)))**2/50000), np.array(S2_avrg_PS)), np.minimum(np.sqrt(S2_err_PS_squared+(np.array(S2_avrg_PS)*(1-np.array(S2_avrg_PS)))**2/50000), 1-np.array(S2_avrg_PS))], color="blue", marker='o', markersize=4, ls='none', label="Triad sim, n=5", alpha=.5, elinewidth=1, capsize=2, zorder=50)

## THEORETICAL ##
t1_ = np.arange(x_min, x_max, .001)
S2_ensemble = 1/2*(1 + np.tanh(n_pairs*(t1_ + t2_MF + t3_MF)))
axs[0].plot(t1_, S2_ensemble, color='g', label='Transition formula', zorder=5, linewidth=2)

n_nodes=10
n_pairs = n_nodes*(n_nodes-1)/2

file_names = ["../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-3.200000__t1_fin_-3.120000__t1_step_0.004000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_45__N_AVRGNG_10000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-3.120000__t1_fin_-3.040000__t1_step_0.004000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_45__N_AVRGNG_10000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-3.040000__t1_fin_-2.960000__t1_step_0.004000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_45__N_AVRGNG_10000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-2.960000__t1_fin_-2.880000__t1_step_0.004000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_45__N_AVRGNG_10000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-2.880000__t1_fin_-2.796000__t1_step_0.004000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_45__N_AVRGNG_10000__ENSEMBLE_AVRG_0.csv"
]

t2 = t2_MF * 2*n_nodes/(n_nodes-2) # t2_MF*(CONVERSION_FACTOR)
t3 = t3_MF * 6*n_nodes**2/(n_nodes**2-5*n_nodes+6) # t3_MF*(CONVERSION_FACTOR)

t1_PS = []
S2_avrg_PS = []
S2_RMS_PS = []
S2_err_PS_squared = []

for file_name in file_names:
    S2_vs_t1 = np.genfromtxt(file_name, delimiter=',')
    # Finding averages
    S2_avrg = []
    S2_RMS = []
    t1 = []
    i = 0
    N_0, N_1 = [], [] # Numbers of data points close to L=0 and L=1
    S2_var_0, S2_var_1 = [], [] # Fluctuations around L=0 and L=1
    while i < S2_vs_t1.shape[0]:
        prev_t1 = S2_vs_t1[i][0]
        S2 = []
        S2_0, S2_1 = [], [] # Data points close to L=0 and L=1
        n_0, n_1 = 0, 0
        while i < S2_vs_t1.shape[0] and prev_t1 == S2_vs_t1[i][0]:
            S2.append(S2_vs_t1[i][2]*2/((n_nodes-1)*(n_nodes-2)))
            if(S2[-1] < 0.5):
                S2_0.append(S2[-1])
                n_0+=1
            else:
                S2_1.append(S2[-1])
                n_1+=1
            i+=1
        N_0.append(n_0)
        N_1.append(n_1)
        S2_avrg.append(np.array(S2).mean())
        S2_RMS.append(np.sqrt(np.array(S2).var()))
        if(S2_0!=[]):
            S2_var_0.append(np.array(S2_0).var())
        else:
            S2_var_0.append(0)
        if(S2_1!=[]):
            S2_var_1.append(np.array(S2_1).var())
        else:
            S2_var_1.append(0)
        t1.append(prev_t1)
    t1_PS.extend(t1)
    S2_avrg_PS.extend(S2_avrg)
    S2_RMS_PS.extend(S2_RMS)
    S2_err_PS_squared.extend((np.array(S2_var_0)*np.array(N_0) + np.array(S2_var_1)*np.array(N_1))/(np.array(N_0)+np.array(N_1)))

axs[0].errorbar(t1_PS, S2_avrg_PS, xerr=None, yerr=[np.minimum(np.sqrt(S2_err_PS_squared+(np.array(S2_avrg_PS)*(1-np.array(S2_avrg_PS)))**2/1000), np.array(S2_avrg_PS)), np.minimum(np.sqrt(S2_err_PS_squared+(np.array(S2_avrg_PS)*(1-np.array(S2_avrg_PS)))**2/1000), 1-np.array(S2_avrg_PS))], color="royalblue", marker='o', markersize=4, ls='none', label="Triad sim, n=10", alpha=.5, elinewidth=1, capsize=2, zorder=50)

## THEORETICAL ##
S2_ensemble = 1/2*(1 + np.tanh(n_pairs*(t1_ + t2_MF + t3_MF)))
axs[0].plot(t1_, S2_ensemble, color='g', zorder=5, linewidth=2)


n_nodes=15
n_pairs = n_nodes*(n_nodes-1)/2

file_names = ["../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-3.040000__t1_fin_-3.032000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-3.032000__t1_fin_-3.024000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-3.024000__t1_fin_-3.016000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-3.016000__t1_fin_-3.008000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-3.008000__t1_fin_-3.000000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-3.000000__t1_fin_-2.992000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-2.992000__t1_fin_-2.984000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-2.984000__t1_fin_-2.976000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-2.976000__t1_fin_-2.968000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-2.968000__t1_fin_-2.959000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__ENSEMBLE_AVRG_0.csv"
]

t2 = t2_MF * 2*n_nodes/(n_nodes-2) # t2_MF*(CONVERSION_FACTOR)
t3 = t3_MF * 6*n_nodes**2/(n_nodes**2-5*n_nodes+6) # t3_MF*(CONVERSION_FACTOR)

t1_PS = []
S2_avrg_PS = []
S2_RMS_PS = []
S2_err_PS_squared = []

for file_name in file_names:
    S2_vs_t1 = np.genfromtxt(file_name, delimiter=',')
    # Finding averages
    S2_avrg = []
    S2_RMS = []
    t1 = []
    i = 0
    N_0, N_1 = [], [] # Numbers of data points close to L=0 and L=1
    S2_var_0, S2_var_1 = [], [] # Fluctuations around L=0 and L=1
    while i < S2_vs_t1.shape[0]:
        prev_t1 = S2_vs_t1[i][0]
        S2 = []
        S2_0, S2_1 = [], [] # Data points close to L=0 and L=1
        n_0, n_1 = 0, 0
        while i < S2_vs_t1.shape[0] and prev_t1 == S2_vs_t1[i][0]:
            S2.append(S2_vs_t1[i][2]*2/((n_nodes-1)*(n_nodes-2)))
            if(S2[-1] < 0.5):
                S2_0.append(S2[-1])
                n_0+=1
            else:
                S2_1.append(S2[-1])
                n_1+=1
            i+=1
        N_0.append(n_0)
        N_1.append(n_1)
        S2_avrg.append(np.array(S2).mean())
        S2_RMS.append(np.sqrt(np.array(S2).var()))
        if(S2_0!=[]):
            S2_var_0.append(np.array(S2_0).var())
        else:
            S2_var_0.append(0)
        if(S2_1!=[]):
            S2_var_1.append(np.array(S2_1).var())
        else:
            S2_var_1.append(0)
        t1.append(prev_t1)
    t1_PS.extend(t1)
    S2_avrg_PS.extend(S2_avrg)
    S2_RMS_PS.extend(S2_RMS)
    S2_err_PS_squared.extend((np.array(S2_var_0)*np.array(N_0) + np.array(S2_var_1)*np.array(N_1))/(np.array(N_0)+np.array(N_1)))

axs[0].errorbar(t1_PS, S2_avrg_PS, xerr=None, yerr=[np.minimum(np.sqrt(S2_err_PS_squared+(np.array(S2_avrg_PS)*(1-np.array(S2_avrg_PS)))**2/1000), np.array(S2_avrg_PS)), np.minimum(np.sqrt(S2_err_PS_squared+(np.array(S2_avrg_PS)*(1-np.array(S2_avrg_PS)))**2/1000), 1-np.array(S2_avrg_PS))], color="navy", marker='o', markersize=4, ls='none', label="Triad sim, n=15", alpha=.5, elinewidth=1, capsize=2, zorder=50)

## THEORETICAL ##
S2_ensemble = 1/2*(1 + np.tanh(n_pairs*(t1_ + t2_MF + t3_MF)))
axs[0].plot(t1_, S2_ensemble, color='g', zorder=5, linewidth=2)


# Plotting-consistency equation
delta_t1 = 0.05
delta_L = 0.000005
t1_range = np.arange(x_min, x_max+.5, delta_t1)
S2_range = np.arange(0, 1, delta_L)
L, T1 = np.meshgrid(S2_range, t1_range)
SC = L - 1/2*(1 + np.tanh(T1 + 2*t2_MF*L + 3*t3_MF*L**2)) #Self-consistency
CS = axs[0].contour(T1, L**2, SC, [0], colors='cyan', zorder=100, linewidths=1.5)

CS.collections[0].set_label('MF solution')

axs[0].axvline(0, color='black')
axs[0].axhline(0, color='black')

axs[0].set_xlim([x_min, x_max])
axs[0].set_ylim([-.01, 1.01])

axs[0].set_xlabel(r'$t_1$', fontsize=22)
axs[0].set_ylabel(r'$\mathfrak{S}_2$', fontsize=22)

axs[0].grid(zorder=0)
l = axs[0].legend(loc='upper left', prop={'size': 11})
l.set_zorder(110)


###################### 3_STARS ############################

n_nodes = 5
n_pairs = n_nodes*(n_nodes-1)/2

file_names = [ "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_5__N_ITERS_PER_LINK_100__t1_start_-3.500000__t1_fin_-2.500000__t1_step_0.010000__t2_MF_2__t3_MF_1__N_AVRGNG_50000__ENSEMBLE_AVRG_0.csv"]

t2 = t2_MF * 2*n_nodes/(n_nodes-2) # t2_MF*(CONVERSION_FACTOR)
t3 = t3_MF * 6*n_nodes**2/(n_nodes**2-5*n_nodes+6) # t3_MF*(CONVERSION_FACTOR)

t1_PS = []
S3_avrg_PS = []
S3_RMS_PS = []
S3_err_PS_squared = []

for file_name in file_names:
    S3_vs_t1 = np.genfromtxt(file_name, delimiter=',')
    # Finding averages
    S3_avrg = []
    S3_RMS = []
    t1 = []
    i = 0
    N_0, N_1 = [], [] # Numbers of data points close to L=0 and L=1
    S3_var_0, S3_var_1 = [], [] # Fluctuations around L=0 and L=1
    while i < S3_vs_t1.shape[0]:
        prev_t1 = S3_vs_t1[i][0]
        S3 = []
        S3_0, S3_1 = [], [] # Data points close to L=0 and L=1
        n_0, n_1 = 0, 0
        while i < S3_vs_t1.shape[0] and prev_t1 == S3_vs_t1[i][0]:
            S3.append(S3_vs_t1[i][3]*6/(n_nodes*(n_nodes-1)*(n_nodes-2)))
            if(S3[-1] < 0.5):
                S3_0.append(S3[-1])
                n_0+=1
            else:
                S3_1.append(S3[-1])
                n_1+=1
            i+=1
        N_0.append(n_0)
        N_1.append(n_1)
        S3_avrg.append(np.array(S3).mean())
        S3_RMS.append(np.sqrt(np.array(S3).var()))
        if(S3_0!=[]):
            S3_var_0.append(np.array(S3_0).var())
        else:
            S3_var_0.append(0)
        if(S3_1!=[]):
            S3_var_1.append(np.array(S3_1).var())
        else:
            S3_var_1.append(0)
        t1.append(prev_t1)
    t1_PS.extend(t1)
    S3_avrg_PS.extend(S3_avrg)
    S3_RMS_PS.extend(S3_RMS)
    S3_err_PS_squared.extend((np.array(S3_var_0)*np.array(N_0) + np.array(S3_var_1)*np.array(N_1))/(np.array(N_0)+np.array(N_1)))

axs[1].errorbar(t1_PS, S3_avrg_PS, xerr=None, yerr=[np.minimum(np.sqrt(S3_err_PS_squared+(np.array(S3_avrg_PS)*(1-np.array(S3_avrg_PS)))**2/50000), np.array(S3_avrg_PS)), np.minimum(np.sqrt(S3_err_PS_squared+(np.array(S3_avrg_PS)*(1-np.array(S3_avrg_PS)))**2/50000), 1-np.array(S3_avrg_PS))], color="blue", marker='o', markersize=4, ls='none', label="Triad sim, n=5", alpha=.5, elinewidth=1, capsize=2, zorder=50)

## THEORETICAL ##
S3_ensemble = 1/2*(1 + np.tanh(n_pairs*(t1_ + t2_MF + t3_MF)))
axs[1].plot(t1_, S3_ensemble, color='g', label='Transition formula', zorder=5, linewidth=2)


n_nodes = 10
n_pairs = n_nodes*(n_nodes-1)/2

file_names = ["../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-3.200000__t1_fin_-3.120000__t1_step_0.004000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_45__N_AVRGNG_10000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-3.120000__t1_fin_-3.040000__t1_step_0.004000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_45__N_AVRGNG_10000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-3.040000__t1_fin_-2.960000__t1_step_0.004000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_45__N_AVRGNG_10000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-2.960000__t1_fin_-2.880000__t1_step_0.004000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_45__N_AVRGNG_10000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-2.880000__t1_fin_-2.796000__t1_step_0.004000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_45__N_AVRGNG_10000__ENSEMBLE_AVRG_0.csv"
]

t2 = t2_MF * 2*n_nodes/(n_nodes-2) # t2_MF*(CONVERSION_FACTOR)
t3 = t3_MF * 6*n_nodes**2/(n_nodes**2-5*n_nodes+6) # t3_MF*(CONVERSION_FACTOR)

t1_PS = []
S3_avrg_PS = []
S3_RMS_PS = []
S3_err_PS_squared = []

for file_name in file_names:
    S3_vs_t1 = np.genfromtxt(file_name, delimiter=',')
    # Finding averages
    S3_avrg = []
    S3_RMS = []
    t1 = []
    i = 0
    N_0, N_1 = [], [] # Numbers of data points close to L=0 and L=1
    S3_var_0, S3_var_1 = [], [] # Fluctuations around L=0 and L=1
    while i < S3_vs_t1.shape[0]:
        prev_t1 = S3_vs_t1[i][0]
        S3 = []
        S3_0, S3_1 = [], [] # Data points close to L=0 and L=1
        n_0, n_1 = 0, 0
        while i < S3_vs_t1.shape[0] and prev_t1 == S3_vs_t1[i][0]:
            S3.append(S3_vs_t1[i][3]*6/(n_nodes*(n_nodes-1)*(n_nodes-2)))
            if(S3[-1] < 0.5):
                S3_0.append(S3[-1])
                n_0+=1
            else:
                S3_1.append(S3[-1])
                n_1+=1
            i+=1
        N_0.append(n_0)
        N_1.append(n_1)
        S3_avrg.append(np.array(S3).mean())
        S3_RMS.append(np.sqrt(np.array(S3).var()))
        if(S3_0!=[]):
            S3_var_0.append(np.array(S3_0).var())
        else:
            S3_var_0.append(0)
        if(S3_1!=[]):
            S3_var_1.append(np.array(S3_1).var())
        else:
            S3_var_1.append(0)
        t1.append(prev_t1)
    t1_PS.extend(t1)
    S3_avrg_PS.extend(S3_avrg)
    S3_RMS_PS.extend(S3_RMS)
    S3_err_PS_squared.extend((np.array(S3_var_0)*np.array(N_0) + np.array(S3_var_1)*np.array(N_1))/(np.array(N_0)+np.array(N_1)))

axs[1].errorbar(t1_PS, S3_avrg_PS, xerr=None, yerr=[np.minimum(np.sqrt(S3_err_PS_squared+(np.array(S3_avrg_PS)*(1-np.array(S3_avrg_PS)))**2/1000), np.array(S3_avrg_PS)), np.minimum(np.sqrt(S3_err_PS_squared+(np.array(S3_avrg_PS)*(1-np.array(S3_avrg_PS)))**2/1000), 1-np.array(S3_avrg_PS))], color="royalblue", marker='o', markersize=4, ls='none', label="Triad sim, n=10", alpha=.5, elinewidth=1, capsize=2, zorder=50)

## THEORETICAL ##
S3_ensemble = 1/2*(1 + np.tanh(n_pairs*(t1_ + t2_MF + t3_MF)))
axs[1].plot(t1_, S3_ensemble, color='g', zorder=5, linewidth=2)


n_nodes = 15
n_pairs = n_nodes*(n_nodes-1)/2

file_names = ["../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-3.040000__t1_fin_-3.032000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-3.032000__t1_fin_-3.024000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-3.024000__t1_fin_-3.016000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-3.016000__t1_fin_-3.008000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-3.008000__t1_fin_-3.000000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-3.000000__t1_fin_-2.992000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-2.992000__t1_fin_-2.984000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-2.984000__t1_fin_-2.976000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-2.976000__t1_fin_-2.968000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__ENSEMBLE_AVRG_0.csv",
              "../data/triad_model/t1_scanning/NEW_TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-2.968000__t1_fin_-2.959000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__ENSEMBLE_AVRG_0.csv"
]

t2 = t2_MF * 2*n_nodes/(n_nodes-2) # t2_MF*(CONVERSION_FACTOR)
t3 = t3_MF * 6*n_nodes**2/(n_nodes**2-5*n_nodes+6) # t3_MF*(CONVERSION_FACTOR)

t1_PS = []
S3_avrg_PS = []
S3_RMS_PS = []
S3_err_PS_squared = []

for file_name in file_names:
    S3_vs_t1 = np.genfromtxt(file_name, delimiter=',')
    # Finding averages
    S3_avrg = []
    S3_RMS = []
    t1 = []
    i = 0
    N_0, N_1 = [], [] # Numbers of data points close to L=0 and L=1
    S3_var_0, S3_var_1 = [], [] # Fluctuations around L=0 and L=1
    while i < S3_vs_t1.shape[0]:
        prev_t1 = S3_vs_t1[i][0]
        S3 = []
        S3_0, S3_1 = [], [] # Data points close to L=0 and L=1
        n_0, n_1 = 0, 0
        while i < S3_vs_t1.shape[0] and prev_t1 == S3_vs_t1[i][0]:
            S3.append(S3_vs_t1[i][3]*6/(n_nodes*(n_nodes-1)*(n_nodes-2)))
            if(S3[-1] < 0.5):
                S3_0.append(S3[-1])
                n_0+=1
            else:
                S3_1.append(S3[-1])
                n_1+=1
            i+=1
        N_0.append(n_0)
        N_1.append(n_1)
        S3_avrg.append(np.array(S3).mean())
        S3_RMS.append(np.sqrt(np.array(S3).var()))
        if(S3_0!=[]):
            S3_var_0.append(np.array(S3_0).var())
        else:
            S3_var_0.append(0)
        if(S3_1!=[]):
            S3_var_1.append(np.array(S3_1).var())
        else:
            S3_var_1.append(0)
        t1.append(prev_t1)
    t1_PS.extend(t1)
    S3_avrg_PS.extend(S3_avrg)
    S3_RMS_PS.extend(S3_RMS)
    S3_err_PS_squared.extend((np.array(S3_var_0)*np.array(N_0) + np.array(S3_var_1)*np.array(N_1))/(np.array(N_0)+np.array(N_1)))

axs[1].errorbar(t1_PS, S3_avrg_PS, xerr=None, yerr=[np.minimum(np.sqrt(S3_err_PS_squared+(np.array(S3_avrg_PS)*(1-np.array(S3_avrg_PS)))**2/1000), np.array(S3_avrg_PS)), np.minimum(np.sqrt(S3_err_PS_squared+(np.array(S3_avrg_PS)*(1-np.array(S3_avrg_PS)))**2/1000), 1-np.array(S3_avrg_PS))], color="navy", marker='o', markersize=4, ls='none', label="Triad sim, n=15", alpha=.5, elinewidth=1, capsize=2, zorder=50)

## THEORETICAL ##
S3_ensemble = 1/2*(1 + np.tanh(n_pairs*(t1_ + t2_MF + t3_MF)))
axs[1].plot(t1_, S3_ensemble, color='g', zorder=5, linewidth=2)


# Plotting-consistency equation
delta_t1 = 0.05
delta_L = 0.000005
t1_range = np.arange(x_min, x_max+.5, delta_t1)
S3_range = np.arange(0, 1, delta_L)
L, T1 = np.meshgrid(S3_range, t1_range)
SC = L - 1/2*(1 + np.tanh(T1 + 2*t2_MF*L + 3*t3_MF*L**2)) #Self-consistency
CS = axs[1].contour(T1, L**3, SC, [0], colors='cyan', zorder=100, linewidths=1.5)

CS.collections[0].set_label('MF solution')

axs[1].axvline(0, color='black')
axs[1].axhline(0, color='black')

axs[1].set_xlim([x_min, x_max])
axs[1].set_ylim([-.01, 1.01])

axs[1].set_xlabel(r'$t_1$', fontsize=22)
axs[1].set_ylabel(r'$\mathfrak{T}$', fontsize=22)

axs[1].grid(zorder=0)
axs[1].legend(loc='upper left', prop={'size': 11})

plt.tight_layout()
plt.savefig('../figures/stars_ensemble_averages.png', dpi=300)

plt.show()
fig.clear()
