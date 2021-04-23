################################################################
# This script plots theoretical prediction for t1-dependence   #
# of connectance together with the corresponding results of    #
# ensemble averaging MC simulations for MF4 and 4-star models. #
################################################################

import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from scipy.optimize import fsolve

t2_MF = 15
t3_MF = -20
t4_MF = 9

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10*1.6, 3.2*1.9))

title_5 = "$t_2=$" + str(t2_MF) + "$;\ t_3=$" + str(t3_MF) + "$;\ t_4=$" + str(t4_MF) + "$;\ n=5$"
title_10 = "$t_2=$" + str(t2_MF) + "$;\ t_3=$" + str(t3_MF) + "$;\ t_4=$" + str(t4_MF) + "$;\ n=10$"

x_min = -6
x_max = -2

file_names = ["../data/p_star_model/t1_scanning/PS__CONNECTANCE_and_STARS__N_NODES_5__N_ITERS_PER_LINK_1000__t1_start_-8__t1_fin_-7__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_10__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv",
              "../data/p_star_model/t1_scanning/PS__CONNECTANCE_and_STARS__N_NODES_5__N_ITERS_PER_LINK_1000__t1_start_-6__t1_fin_-5.500000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_10__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv",
              "../data/p_star_model/t1_scanning/PS__CONNECTANCE_and_STARS__N_NODES_5__N_ITERS_PER_LINK_1000__t1_start_-5.500000__t1_fin_-2__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_10__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv",
              "../data/p_star_model/t1_scanning/_PS__CONNECTANCE_and_STARS__N_NODES_5__N_ITERS_PER_LINK_1000__t1_start_-7__t1_fin_-6__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_10__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv",
              "../data/p_star_model/t1_scanning/_PS__CONNECTANCE_and_STARS__N_NODES_5__N_ITERS_PER_LINK_1000__t1_start_-2__t1_fin_-1__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_10__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv"
]
n_nodes = 5
N_links_max = n_nodes*(n_nodes-1)/2

t1_PS_5 = []
L_avrg_PS_5 = []
L_RMS_PS_5 = []
L_err_PS_5_squared = []

for file_name in file_names:
    L_vs_t1 = np.genfromtxt(file_name, delimiter=',')
    # Finding averages
    L_avrg = []
    L_RMS = []
    t1 = []
    i = 0
    N_0, N_1 = [], [] # Numbers of data points close to L=0 and L=1
    L_var_0, L_var_1 = [], [] # Fluctuations around L=0 and L=1
    while i < L_vs_t1.shape[0]:
        prev_t1 = L_vs_t1[i][0]
        L = []
        L_0, L_1 = [], [] # Data points close to L=0 and L=1
        n_0, n_1 = 0, 0
        while i < L_vs_t1.shape[0] and prev_t1 == L_vs_t1[i][0]:
            L.append(L_vs_t1[i][1])
            if(L[-1] < 0.5):
                L_0.append(L[-1])
                n_0+=1
            else:
                L_1.append(L[-1])
                n_1+=1
            i+=1
        N_0.append(n_0)
        N_1.append(n_1)
        L_avrg.append(np.array(L).mean())
        L_RMS.append(np.sqrt(np.array(L).var()))
        if(L_0!=[]):
            L_var_0.append(np.array(L_0).var())
        else:
            L_var_0.append(0)
        if(L_1!=[]):
            L_var_1.append(np.array(L_1).var())
        else:
            L_var_1.append(0)
        t1.append(prev_t1)
    t1_PS_5.extend(t1)
    L_avrg_PS_5.extend(L_avrg)
    L_RMS_PS_5.extend(L_RMS)
    L_err_PS_5_squared.extend((np.array(L_var_0)*np.array(N_0) + np.array(L_var_1)*np.array(N_1))/(np.array(N_0)+np.array(N_1)))


file_names = ["../data/p_star_model/t1_scanning/_PS__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-6__t1_fin_-5__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_100__RANDOM_INITIALISATION_0.csv",
              "../data/p_star_model/t1_scanning/_PS__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-5__t1_fin_-4__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_100__RANDOM_INITIALISATION_0_.csv",
              "../data/p_star_model/t1_scanning/_PS__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-4__t1_fin_-3__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_100__RANDOM_INITIALISATION_0.csv",
              "../data/p_star_model/t1_scanning/_PS__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-3__t1_fin_-2__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_100__RANDOM_INITIALISATION_0.csv",
               "../data/p_star_model/t1_scanning/_PS__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-2__t1_fin_-1__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_100__RANDOM_INITIALISATION_0.csv",
               "../data/p_star_model/t1_scanning/_PS__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-7__t1_fin_-6__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_100__RANDOM_INITIALISATION_0.csv"
]

file_names += ["../data/p_star_model/t1_scanning/_PS__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_2000000__t1_start_-4.500000__t1_fin_-4.480000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_500__RANDOM_INITIALISATION_0.csv",
               "../data/p_star_model/t1_scanning/_PS__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_2000000__t1_start_-4.480000__t1_fin_-4.460000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_500__RANDOM_INITIALISATION_0.csv",
               "../data/p_star_model/t1_scanning/_PS__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_2000000__t1_start_-4.460000__t1_fin_-4.440000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_500__RANDOM_INITIALISATION_0.csv",
               "../data/p_star_model/t1_scanning/_PS__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_2000000__t1_start_-4.440000__t1_fin_-4.420000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_500__RANDOM_INITIALISATION_0.csv",
               "../data/p_star_model/t1_scanning/_PS__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_2000000__t1_start_-4.420000__t1_fin_-4.400000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_500__RANDOM_INITIALISATION_0.csv",
               "../data/p_star_model/t1_scanning/_PS__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_2000000__t1_start_-4.400000__t1_fin_-4.380000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_500__RANDOM_INITIALISATION_0.csv",
               "../data/p_star_model/t1_scanning/_PS__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_2000000__t1_start_-4.380000__t1_fin_-4.360000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_500__RANDOM_INITIALISATION_0.csv",
               "../data/p_star_model/t1_scanning/_PS__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_2000000__t1_start_-4.360000__t1_fin_-4.340000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_500__RANDOM_INITIALISATION_0.csv",
               "../data/p_star_model/t1_scanning/_PS__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_2000000__t1_start_-4.340000__t1_fin_-4.320000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_500__RANDOM_INITIALISATION_0.csv",
               "../data/p_star_model/t1_scanning/_PS__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_2000000__t1_start_-4.320000__t1_fin_-4.300000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_500__RANDOM_INITIALISATION_0.csv",
               "../data/p_star_model/t1_scanning/_PS__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_2000000__t1_start_-4.300000__t1_fin_-4.280000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_500__RANDOM_INITIALISATION_0.csv",
               "../data/p_star_model/t1_scanning/_PS__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_2000000__t1_start_-4.280000__t1_fin_-4.260000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_500__RANDOM_INITIALISATION_0.csv"
]

n_nodes = 10
N_links_max = n_nodes*(n_nodes-1)/2

t1_PS_10 = []
L_avrg_PS_10 = []
L_RMS_PS_10 = []
L_err_PS_10_squared = []

for file_name in file_names:
    L_vs_t1 = np.genfromtxt(file_name, delimiter=',')
    # Finding averages
    L_avrg = []
    L_RMS = []
    t1 = []
    i = 0
    N_0, N_1 = [], [] # Numbers of data points close to L=0 and L=1
    L_var_0, L_var_1 = [], [] # Fluctuations around L=0 and L=1
    while i < L_vs_t1.shape[0]:
        prev_t1 = L_vs_t1[i][0]
        L = []
        L_0, L_1 = [], [] # Data points close to L=0 and L=1
        n_0, n_1 = 0, 0
        while i < L_vs_t1.shape[0] and prev_t1 == L_vs_t1[i][0]:
            L.append(L_vs_t1[i][1])
            if(L_vs_t1[i][1] < 0.5):
                L_0.append(L_vs_t1[i][1])
                n_0+=1
            else:
                L_1.append(L_vs_t1[i][1])
                n_1+=1
            i+=1
        N_0.append(n_0)
        N_1.append(n_1)
        L_avrg.append(np.array(L).mean())
        L_RMS.append(np.sqrt(np.array(L).var()))
        if(L_0!=[]):
            L_var_0.append(np.array(L_0).var())
        else:
            L_var_0.append(0)
        if(L_1!=[]):
            L_var_1.append(np.array(L_1).var())
        else:
            L_var_1.append(0)
        t1.append(prev_t1)
    t1_PS_10.extend(t1)
    L_avrg_PS_10.extend(L_avrg)
    L_RMS_PS_10.extend(L_RMS)
    L_err_PS_10_squared.extend((np.array(L_var_0)*np.array(N_0) + np.array(L_var_1)*np.array(N_1))/(np.array(N_0)+np.array(N_1)))


file_names = ["../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_5__N_ITERS_PER_LINK_1000__t1_start_-7__t1_fin_-1__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_10__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv"]

n_nodes = 5
N_links_max = n_nodes*(n_nodes-1)/2

t1_MF_5 = []
L_avrg_MF_5 = []
L_RMS_MF_5 = []
L_err_MF_5_squared = []

for file_name in file_names:
    L_vs_t1 = np.genfromtxt(file_name, delimiter=',')
    # Finding averages
    L_avrg = []
    L_RMS = []
    t1 = []
    i = 0
    N_0, N_1 = [], [] # Numbers of data points close to L=0 and L=1
    L_var_0, L_var_1 = [], [] # Fluctuations around L=0 and L=1
    while i < L_vs_t1.shape[0]:
        prev_t1 = L_vs_t1[i][0]
        L = []
        L_0, L_1 = [], [] # Data points close to L=0 and L=1
        n_0, n_1 = 0, 0
        while i < L_vs_t1.shape[0] and prev_t1 == L_vs_t1[i][0]:
            L.append(L_vs_t1[i][1])
            if(L_vs_t1[i][1] < 0.5):
                L_0.append(L_vs_t1[i][1])
                n_0+=1
            else:
                L_1.append(L_vs_t1[i][1])
                n_1+=1
            i+=1
        N_0.append(n_0)
        N_1.append(n_1)
        L_avrg.append(np.array(L).mean())
        L_RMS.append(np.sqrt(np.array(L).var()))
        if(L_0!=[]):
            L_var_0.append(np.array(L_0).var())
        else:
            L_var_0.append(0)
        if(L_1!=[]):
            L_var_1.append(np.array(L_1).var())
        else:
            L_var_1.append(0)
        t1.append(prev_t1)
    t1_MF_5.extend(t1)
    L_avrg_MF_5.extend(L_avrg)
    L_RMS_MF_5.extend(L_RMS)
    L_err_MF_5_squared.extend((np.array(L_var_0)*np.array(N_0) + np.array(L_var_1)*np.array(N_1))/(np.array(N_0)+np.array(N_1)))


file_names = ["../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-6__t1_fin_-5__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_100__RANDOM_INITIALISATION_0.csv",
               "../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-5__t1_fin_-4__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_100__RANDOM_INITIALISATION_0.csv",
               "../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-4__t1_fin_-3__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_100__RANDOM_INITIALISATION_0.csv",
               "../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-3__t1_fin_-2__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_100__RANDOM_INITIALISATION_0.csv",
               "../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-2__t1_fin_-1__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_100__RANDOM_INITIALISATION_0.csv"
]
file_names += ["../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_2000000__t1_start_-4.450000__t1_fin_-4.430000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_500__RANDOM_INITIALISATION_0.csv",
             "../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_2000000__t1_start_-4.430000__t1_fin_-4.410000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_500__RANDOM_INITIALISATION_0.csv",
             "../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_2000000__t1_start_-4.410000__t1_fin_-4.390000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_500__RANDOM_INITIALISATION_0.csv",
             "../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_2000000__t1_start_-4.390000__t1_fin_-4.370000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_500__RANDOM_INITIALISATION_0.csv",
             "../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_2000000__t1_start_-4.370000__t1_fin_-4.350000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_500__RANDOM_INITIALISATION_0.csv",
             "../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_2000000__t1_start_-4.350000__t1_fin_-4.330000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_500__RANDOM_INITIALISATION_0.csv",
             "../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_2000000__t1_start_-4.330000__t1_fin_-4.310000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_500__RANDOM_INITIALISATION_0.csv",
             "../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_2000000__t1_start_-4.310000__t1_fin_-4.290000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_500__RANDOM_INITIALISATION_0.csv",
             "../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_2000000__t1_start_-4.290000__t1_fin_-4.270000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_500__RANDOM_INITIALISATION_0.csv",
             "../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_2000000__t1_start_-4.270000__t1_fin_-4.250000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_500__RANDOM_INITIALISATION_0.csv",
             "../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_2000000__t1_start_-4.250000__t1_fin_-4.230000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_500__RANDOM_INITIALISATION_0.csv",
             "../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_2000000__t1_start_-4.230000__t1_fin_-4.210000__t1_step_0.010000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_45__N_AVRGNG_500__RANDOM_INITIALISATION_0.csv"
]

n_nodes = 10
N_links_max = n_nodes*(n_nodes-1)/2

t1_MF_10 = []
L_avrg_MF_10 = []
L_RMS_MF_10 = []
L_err_MF_10_squared = []

for file_name in file_names:
    L_vs_t1 = np.genfromtxt(file_name, delimiter=',')
    # Finding averages
    L_avrg = []
    L_RMS = []
    t1 = []
    i = 0
    N_0, N_1 = [], [] # Numbers of data points close to L=0 and L=1
    L_var_0, L_var_1 = [], [] # Fluctuations around L=0 and L=1
    while i < L_vs_t1.shape[0]:
        prev_t1 = L_vs_t1[i][0]
        L = []
        L_0, L_1 = [], [] # Data points close to L=0 and L=1
        n_0, n_1 = 0, 0
        while i < L_vs_t1.shape[0] and prev_t1 == L_vs_t1[i][0]:
            L.append(L_vs_t1[i][1])
            if(L_vs_t1[i][1] < 0.5):
                L_0.append(L_vs_t1[i][1])
                n_0+=1
            else:
                L_1.append(L_vs_t1[i][1])
                n_1+=1
            i+=1
        N_0.append(n_0)
        N_1.append(n_1)
        L_avrg.append(np.array(L).mean())
        L_RMS.append(np.sqrt(np.array(L).var()))
        if(L_0!=[]):
            L_var_0.append(np.array(L_0).var())
        else:
            L_var_0.append(0)
        if(L_1!=[]):
            L_var_1.append(np.array(L_1).var())
        else:
            L_var_1.append(0)
        t1.append(prev_t1)
    t1_MF_10.extend(t1)
    L_avrg_MF_10.extend(L_avrg)
    L_RMS_MF_10.extend(L_RMS)
    L_err_MF_10_squared.extend((np.array(L_var_0)*np.array(N_0) + np.array(L_var_1)*np.array(N_1))/(np.array(N_0)+np.array(N_1)))


# Plotting self-consistency equation
delta_t1 = 0.05
delta_L = 0.000005
t1_range = np.arange(x_min, x_max+delta_t1, delta_t1)
L_range = np.arange(0, 1, delta_L)
L, T1 = np.meshgrid(L_range, t1_range)
SC = L - 1/2*(1 + np.tanh(T1 + 2*t2_MF*L + 3*t3_MF*L**2 + 4*t4_MF*L**3)) #Self-consistency
CS_5 = axs[0].contour(T1, L, SC, [0], colors='cyan', zorder=50, linewidths=1.5)
CS_10 = axs[1].contour(T1, L, SC, [0], colors='cyan', zorder=50, linewidths=1.5)
CS_5.collections[0].set_label('Self-consistency equation')
CS_10.collections[0].set_label('Self-consistency equation')

# plt.title(file_names[0])
axs[0].axvline(0, color='black')
axs[0].axhline(0, color='black')
axs[1].axvline(0, color='black')
axs[1].axhline(0, color='black')


axs[0].set_ylim([-.01, 1.01])
axs[0].set_xlim([x_min, x_max])
axs[1].set_ylim([-.01, 1.01])
axs[1].set_xlim([x_min, x_max])

axs[0].set_title(title_5, fontsize=22)
axs[1].set_title(title_10, fontsize=22)

axs[0].set_xlabel(r'$t_1$', fontsize=22)
axs[1].set_xlabel(r'$t_1$', fontsize=22)
axs[0].set_ylabel(r'$\bar L$', fontsize=22)

##### Finding ensemble averages #####
n_nodes = 5
N_links_max = int(n_nodes*(n_nodes-1)/2)

delta = 0.0001

T = 1
T1_MF = np.arange(x_min, x_max+delta, delta)
L_EXPECTED_5 = []
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
    F = -2*(t1_MF*L + t2_MF*L**2 + t3_MF*L**3 + t4_MF*L**4) - T*S
    numerator, denominator = 0, 0
    Z = 0
    for i in range(0, N_links_max+1):
        numerator += L[i]*np.exp(-N_links_max/T*F[i])
        denominator += np.exp(-N_links_max/T*F[i])
    L_EXPECTED_5.append(numerator/denominator)

n_nodes = 10
N_links_max = int(n_nodes*(n_nodes-1)/2)

T1_MF = np.arange(x_min, x_max+delta, delta)
L_EXPECTED_10 = []
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
    F = -2*(t1_MF*L + t2_MF*L**2 + t3_MF*L**3 + t4_MF*L**4) - T*S
    numerator, denominator = 0, 0
    Z = 0
    for i in range(0, N_links_max+1):
        numerator += L[i]*np.exp(-N_links_max/T*F[i])
        denominator += np.exp(-N_links_max/T*F[i])
    L_EXPECTED_10.append(numerator/denominator)


axs[0].plot(t1_PS_5, L_avrg_PS_5, 'h',markersize=4, label="4-star simulation, n=5", zorder=100, color='blue', alpha=.5)
axs[0].plot(t1_MF_5, L_avrg_MF_5, 's',markersize=4, label="MF4 simulation, n=5", zorder=100, color='orange', alpha=.5)
axs[1].plot(t1_PS_10, L_avrg_PS_10, 'h',markersize=4, label="4-star simulation, n=10", zorder=101, color='blue', alpha=.5)
axs[1].plot(t1_MF_10, L_avrg_MF_10, 's',markersize=4, label="MF4 simulation, n=10", zorder=101, color='orange', alpha=.5)
axs[0].plot(T1_MF, L_EXPECTED_5, color='red', label='Exact solution', zorder=150)
axs[1].plot(T1_MF, L_EXPECTED_10, color='red', label='Exact solution', zorder=150)

axs[0].legend(loc='upper left', prop={'size': 16}).set_zorder(200)
axs[1].legend(loc='upper left', prop={'size': 16}).set_zorder(200)

axs[0].grid()
axs[1].grid()

plt.tight_layout()

plt.savefig('../figures/5_10_nodes_MF4_4S_ensemble_avrg.png', dpi=300)

plt.show()

fig.clear()
