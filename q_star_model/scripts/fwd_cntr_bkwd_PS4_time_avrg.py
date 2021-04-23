######################################################
# This script plots the results of time averaging MC #
# simulations for the 4-star model together with the #
# corresponding theoretical predictions.             #
######################################################
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

fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(10*1.4, 3.2*1.4))

x_min = -5.1
x_max = -2.2

file_names = ["../data/p_star_model/t1_scanning/PS__CONNECTANCE_and_STARS__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_-5.500000__t1_fin_-2__t1_step_0.005000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_499500__N_AVRGNG_10__RANDOM_INITIALISATION_0.csv"]

t1_fwd = []
L_avrg_fwd = []
L_RMS_fwd = []
L_err_fwd_squared = []

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
    t1_fwd.extend(t1)
    L_avrg_fwd.extend(L_avrg)
    L_RMS_fwd.extend(L_RMS)
    L_err_fwd_squared.extend((np.array(L_var_0)*np.array(N_0) + np.array(L_var_1)*np.array(N_1))/(np.array(N_0)+np.array(N_1)))



file_names = ["../data/p_star_model/t1_scanning/PS__CONNECTANCE_and_STARS__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_-2__t1_fin_-5.500000__t1_step_-0.005000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_499500__N_AVRGNG_10__RANDOM_INITIALISATION_0.csv"]

t1_bkwd = []
L_avrg_bkwd = []
L_RMS_bkwd = []
L_err_bkwd_squared = []

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
    t1_bkwd.extend(t1)
    L_avrg_bkwd.extend(L_avrg)
    L_RMS_bkwd.extend(L_RMS)
    L_err_bkwd_squared.extend((np.array(L_var_0)*np.array(N_0) + np.array(L_var_1)*np.array(N_1))/(np.array(N_0)+np.array(N_1)))



file_names = ["../data/p_star_model/t1_scanning/PS__CONNECTANCE_and_STARS__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_-4__t1_fin_-5.500000__t1_step_-0.005000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_499500__N_AVRGNG_10__RANDOM_INITIALISATION_0.csv",
              "../data/p_star_model/t1_scanning/PS__CONNECTANCE_and_STARS__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_-4__t1_fin_-2__t1_step_0.005000__t2_MF_15__t3_MF_-20__t4_MF_9__N_DYNAMIC_PAIRS_499500__N_AVRGNG_10__RANDOM_INITIALISATION_0.csv"
]

t1_center = []
L_avrg_center = []
L_RMS_center = []
L_err_center_squared = []

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
    t1_center.extend(t1)
    L_avrg_center.extend(L_avrg)
    L_RMS_center.extend(L_RMS)
    L_err_center_squared.extend((np.array(L_var_0)*np.array(N_0) + np.array(L_var_1)*np.array(N_1))/(np.array(N_0)+np.array(N_1)))


# Plotting self-consistency equation
delta_t1 = 0.05
delta_L = 0.000005
t1_range = np.arange(x_min, x_max+delta_t1, delta_t1)
L_range = np.arange(0, 1, delta_L)
L, T1 = np.meshgrid(L_range, t1_range)
SC = L - 1/2*(1 + np.tanh(T1 + 2*t2_MF*L + 3*t3_MF*L**2 + 4*t4_MF*L**3)) #Self-consistency
CS_fwd = axs[0].contour(T1, L, SC, [0], colors='cyan', zorder=50, linewidths=1.5)
CS_bkwd = axs[1].contour(T1, L, SC, [0], colors='cyan', zorder=50, linewidths=1.5)
CS_cntr = axs[2].contour(T1, L, SC, [0], colors='cyan', zorder=50, linewidths=1.5)
CS_fwd.collections[0].set_label('Self-consistency equation')
CS_bkwd.collections[0].set_label('Self-consistency equation')
CS_cntr.collections[0].set_label('Self-consistency equation')

axs[0].axvline(0, color='black')
axs[0].axhline(0, color='black')
axs[1].axvline(0, color='black')
axs[1].axhline(0, color='black')
axs[2].axvline(0, color='black')
axs[2].axhline(0, color='black')


axs[0].set_ylim([-.01, 1.01])
axs[0].set_xlim([x_min, x_max])
axs[1].set_ylim([-.01, 1.01])
axs[1].set_xlim([x_min, x_max])
axs[2].set_ylim([-.01, 1.01])
axs[2].set_xlim([x_min, x_max])

axs[0].set_title(r'$\longrightarrow$', fontsize=22)
axs[1].set_title(r'$\longleftarrow\longrightarrow$', fontsize=22)
axs[2].set_title(r'$\longleftarrow$', fontsize=22)

axs[0].set_xlabel(r'$t_1$', fontsize=22)
axs[1].set_xlabel(r'$t_1$', fontsize=22)
axs[2].set_xlabel(r'$t_1$', fontsize=22)
axs[0].set_ylabel(r'$\bar L$', fontsize=22)

# Finding ensemble averages in thermodynamic limit
x_min_ = -5.1
x_max_ = -2.5

T1_MF = np.arange(x_min_, x_max_, 0.001)
L_MIN = []
T1_MF_ = []
for t1_MF in T1_MF:
    def func(L):
        return L - 1/2*(1 + np.tanh(t1_MF + 2*t2_MF*L + 3*t3_MF*L**2 + 4*t4_MF*L**3))
    L_ = np.array([fsolve(func, 0)[0], fsolve(func, .5)[0], fsolve(func, 1)[0]])
    F = -2*(t1_MF*L_ + t2_MF*L_**2 + t3_MF*L_**3 + t4_MF*L_**4)+ L_*np.log(L_/(1-L_)) + np.log(1-L_)
    L_min = L_[list(F).index(min(F))]
    L_MIN.append(L_min)

axs[0].plot(t1_fwd, L_avrg_fwd, 'h',markersize=4, label="4-star simulation, n=1000", zorder=100, color='blue', alpha=.3)
axs[1].plot(t1_center, L_avrg_center, 'h',markersize=4, label="4-star simulation, n=1000", zorder=100, color='blue', alpha=.3)
axs[2].plot(t1_bkwd, L_avrg_bkwd, 'h',markersize=4, label="4-star simulation, n=1000", zorder=101, color='blue', alpha=.3)


axs[0].plot(T1_MF, L_MIN, color='red', zorder=51, label='Thermodynamic solution')
axs[1].plot(T1_MF, L_MIN, color='red', zorder=51, label='Thermodynamic solution')
axs[2].plot(T1_MF, L_MIN, color='red', zorder=51, label='Thermodynamic solution')

axs[0].legend(loc=[.295, .21], prop={'size': 11}).set_zorder(200)
axs[1].legend(loc=[.295, .21], prop={'size': 11}).set_zorder(200)
axs[2].legend(loc=[.295, .21], prop={'size': 11}).set_zorder(200)

axs[0].grid()
axs[1].grid()
axs[2].grid()

plt.tight_layout()

plt.savefig('../figures/fwd_cntr_bkwd_PS4_time_avrg.png', dpi=300)

plt.show()

fig.clear()
