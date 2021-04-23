###############################################################
# This script plots the degree distributions together with    #           
# the distributions of local clustering coefficients for the  #
# corresponding 3-star and MF models with desired connectance #
# Ld controlled by appropriate setting of t1.                 #
###############################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as ss

fig, axs = plt.subplots(2, 2, figsize=(10.*1.02, 6.*1.02))


#############################################################################
############################ LCC DISTRIBUTIONS ##############################
#############################################################################

N_nodes = 1000

file_name = "../data/p_star_model/clustering_distributions/PS__N_NODES_1000__N_ITERS_PER_LINK_10__t1_MF_3284.576351__t2_MF_-7500__t3_MF_4500__N_AVRGNG_100000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_PS = []
for line in lines:
    deg_seq_PS.append(float(line.split('\n')[0]))
    
N_lines = np.array(lines).shape[0]
print(np.average(deg_seq_PS))

file_name = "../data/mean_field_model/clustering_distributions/MF__N_NODES_1000__N_ITERS_PER_LINK_10__t1_MF_3284.576351__t2_MF_-7500__t3_MF_4500__N_AVRGNG_100000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_MF = []
for lc in range(0, N_lines):
    deg_seq_MF.append(float(lines[lc].split('\n')[0]))

print(np.average(deg_seq_MF))


#### Finding MF theory prediction for <C> ####
from scipy.optimize import fsolve
Ld = .3
t2_MF = -7500
t3_MF = 4500
t1_MF = 1/2*np.log(Ld/(1-Ld)) - 2*t2_MF*Ld - 3*t3_MF*Ld*Ld

SC = lambda L : L - .5*(1+np.tanh(t1_MF + 2*t2_MF*L + 3*t3_MF*L**2))

L_initial_guess = 0.4
L_sol = fsolve(SC, L_initial_guess)[0]
cc_theor = L_sol*(1-(1-L_sol)**(N_nodes-1)-L_sol*(N_nodes-1)*(1-L_sol)**(N_nodes-2))
print("C_theor=", cc_theor)

axs[0,0].hist(deg_seq_PS, 100, color='b', density=True, label="3-star sim")
axs[0,0].hist(deg_seq_MF, 100, color='orange', density=True, alpha=.65, label="MF3 sim")
axs[0,0].grid()
title = "$t_1=$" + str(round(t1_MF, 4)) + "$,\ t_2=$" + str(t2_MF) + "$,\ t_3=$" + str(t3_MF) + r"$;\ n=10^3$"
axs[0,0].set_title(title, fontsize=13)
axs[0,0].axvline(cc_theor, color='r', alpha=.65, label="MF theor")
axs[0,0].legend(fontsize=12, loc="upper right")



N_nodes = 5000

file_name = "../data/p_star_model/clustering_distributions/PS__N_NODES_5000__N_ITERS_PER_LINK_10__t1_MF_3284.576351__t2_MF_-7500__t3_MF_4500__N_AVRGNG_1000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_PS = []
for line in lines:
    deg_seq_PS.append(float(line.split('\n')[0]))
    
N_lines = np.array(lines).shape[0]
print(np.average(deg_seq_PS))


file_name = "../data/mean_field_model/clustering_distributions/MF__N_NODES_5000__N_ITERS_PER_LINK_10__t1_MF_3284.576351__t2_MF_-7500__t3_MF_4500__N_AVRGNG_1000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_MF = []
for lc in range(0, N_lines):
    deg_seq_MF.append(float(lines[lc].split('\n')[0]))

print(np.average(deg_seq_MF))

axs[0,1].hist(deg_seq_PS, 100, color='b', density=True, label="3-star sim")
axs[0,1].hist(deg_seq_MF, 100, color='orange', density=True, alpha=.65, label="MF3 sim")
axs[0,1].grid()
title = "$t_1=$" + str(round(t1_MF, 4)) + "$,\ t_2=$" + str(t2_MF) + "$,\ t_3=$" + str(t3_MF) + r"$;\ n=5\cdot 10^3$"
axs[0,1].set_title(title, fontsize=13)
axs[0,1].axvline(cc_theor, color='r', alpha=.65, label="MF theor")
axs[0,1].legend(fontsize=12, loc="upper right")


#############################################################################
######################### DEGREE DISTRIBUTIONS ##############################
#############################################################################
N_nodes = 1000
file_name = "../data/p_star_model/degree_distributions/PS__N_NODES_1000__N_ITERS_PER_LINK_10__t1_MF_3284.576351__t2_MF_-7500__t3_MF_4500__N_AVRGNG_100000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_PS = []
for line in lines:
    deg_seq_PS.append(int(line.split('\n')[0]))
N_lines = np.array(lines).shape[0]
L_PS = np.average(deg_seq_PS)
print(L_PS)


file_name = "../data/mean_field_model/degree_distributions/MF__N_NODES_1000__N_ITERS_PER_LINK_10__t1_MF_3284.576351__t2_MF_-7500__t3_MF_4500__N_AVRGNG_100000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_MF = []
for lc in range(0, N_lines):
    deg_seq_MF.append(int(lines[lc].split('\n')[0]))

L_MF = np.average(deg_seq_MF)
print(L_MF)

# Theoretical
L_MF /= (N_nodes-1)
sigma = np.sqrt((N_nodes-1)*L_MF*(1-L_MF))
x_min, x_max = L_MF*(N_nodes-1) - 3.5*sigma, L_MF*(N_nodes-1) + 3.5*sigma
degrees = np.arange(np.rint(x_min), np.rint(x_max)+1)
p = []
p = np.sqrt(1/(2*np.pi*(N_nodes-1)*L_MF*(1-L_MF)))*np.exp(-(degrees-L_MF*(N_nodes-1))**2/(2*(N_nodes-1)*L_MF*(1-L_MF)))
axs[1,0].plot(degrees, p, color='r', label="MF theor")
axs[1,0].hist(deg_seq_PS, bins = degrees, color='b', density=True, label="3-star sim")
axs[1,0].hist(deg_seq_MF, bins = degrees, color='orange', density=True, alpha=.65, label="MF3 sim")
axs[1,0].legend(fontsize=12, loc="upper right")

axs[1,0].grid()
axs[1,0].set_xlim([x_min, x_max])

################
N_nodes = 5000
file_name = "../data/p_star_model/degree_distributions/PS__N_NODES_5000__N_ITERS_PER_LINK_10__t1_MF_3284.576351__t2_MF_-7500__t3_MF_4500__N_AVRGNG_1000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_PS = []
for line in lines:
    deg_seq_PS.append(int(line.split('\n')[0]))
N_lines = np.array(lines).shape[0]
L_PS = np.average(deg_seq_PS)
print(L_PS)


file_name = "../data/mean_field_model/degree_distributions/MF__N_NODES_5000__N_ITERS_PER_LINK_10__t1_MF_3284.576351__t2_MF_-7500__t3_MF_4500__N_AVRGNG_1000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_MF = []
for lc in range(0, N_lines):
    deg_seq_MF.append(int(lines[lc].split('\n')[0]))

L_MF = np.average(deg_seq_MF)
print(L_MF)

# Theoretical
L_MF /= (N_nodes-1)
sigma = np.sqrt((N_nodes-1)*L_MF*(1-L_MF))
x_min, x_max = L_MF*(N_nodes-1) - 3.5*sigma, L_MF*(N_nodes-1) + 3.5*sigma
degrees = np.arange(np.rint(x_min), np.rint(x_max)+1)
p = []
p = np.sqrt(1/(2*np.pi*(N_nodes-1)*L_MF*(1-L_MF)))*np.exp(-(degrees-L_MF*(N_nodes-1))**2/(2*(N_nodes-1)*L_MF*(1-L_MF)))
axs[1,1].plot(degrees, p, color='r', label="MF theor")
axs[1,1].hist(deg_seq_PS, bins = degrees, color='b', density=True, label="3-star sim")
axs[1,1].hist(deg_seq_MF, bins = degrees, color='orange', density=True, alpha=.65, label="MF3 sim")
axs[1,1].legend(fontsize=12, loc="upper right")

axs[1,1].grid()
axs[1,1].set_xlim([x_min, x_max])


axs[1,0].set_ylabel(r'$\langle \mathrm{p}(k|\mathbf{A}) \rangle$', fontsize=20)
axs[0,0].set_ylabel(r'$\langle \mathrm{p}(c|\mathbf{A}) \rangle$', fontsize=20)
axs[0,0].set_xlabel(r'$c$', fontsize=20)
axs[0,1].set_xlabel(r'$c$', fontsize=20)
axs[1,0].set_xlabel(r'$k$', fontsize=20)
axs[1,1].set_xlabel(r'$k$', fontsize=20)

plt.tight_layout()

plt.savefig('../figures/low_temp_LCC_and_degree_distr.png', dpi=300)
plt.show()

fig.clear()
