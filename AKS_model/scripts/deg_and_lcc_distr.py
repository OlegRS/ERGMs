##################################################
# This script plots degree and LCC distributions #
# together with the corresponding theoretical    #
# predictions for different AKS models.          #
##################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as ss

fig, axs = plt.subplots(2, 3, figsize=(19./1.5, 6./1.5*2))

N_nodes = 200

file_name = "../data/AKS_model/degree_distributions/AKS__N_NODES_200__N_ITERS_PER_LINK_100__tau1_51.000681__mu_-200__lambda_1__N_AVRGNG_10000__INITIALIZE_RANDOMLY_0.csv"

mu = -200
lam = 1
tau1 = 51.000681

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_PS = []
for line in lines:
    deg_seq_PS.append(int(line.split('\n')[0]))
N_lines = np.array(lines).shape[0]
L_AKS = np.average(deg_seq_PS)
print(L_AKS)

# Theoretical
L_AKS /= (N_nodes-1)
sigma = np.sqrt((N_nodes-1)*L_AKS*(1-L_AKS))
x_min, x_max = L_AKS*(N_nodes-1) - 5*sigma, L_AKS*(N_nodes-1) + 5*sigma
degrees = np.arange(np.rint(x_min), np.rint(x_max)+1)
p = []

p = np.sqrt(1/(2*np.pi*(N_nodes-1)*L_AKS*(1-L_AKS)))*np.exp(-(degrees-L_AKS*(N_nodes-1))**2/(2*(N_nodes-1)*L_AKS*(1-L_AKS)))
axs[0,0].plot(degrees, p, color='r', label='MF theor')
axs[0,0].hist(deg_seq_PS, bins = degrees, color='b', density=True, label='AKS sim')

axs[0,0].grid()

title = r"$\rho=$" + str(round(tau1,4)) + r"$,\ \tilde\mu=$" + str(mu) + r"$;\ \tilde\lambda=$" + str(lam) + r"$;\ n=$" + str(N_nodes)
axs[0,0].set_title(title, fontsize=14)
axs[0,0].legend(loc='upper left', prop={'size': 11})


file_name = "../data/AKS_model/degree_distributions/AKS__N_NODES_200__N_ITERS_PER_LINK_100__tau1_0.090594__mu_-2__lambda_1__N_AVRGNG_10000__INITIALIZE_RANDOMLY_0.csv"

mu = -2
lam = 1
tau1 = 0.090594

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_PS = []
for line in lines:
    deg_seq_PS.append(int(line.split('\n')[0]))
N_lines = np.array(lines).shape[0]
L_AKS = np.average(deg_seq_PS)
print(L_AKS)


# Theoretical
L_AKS /= (N_nodes-1)
sigma = np.sqrt((N_nodes-1)*L_AKS*(1-L_AKS))
x_min, x_max = L_AKS*(N_nodes-1) - 5*sigma, L_AKS*(N_nodes-1) + 5*sigma
degrees = np.arange(x_min, x_max)

p = np.sqrt(1/(2*np.pi*(N_nodes-1)*L_AKS*(1-L_AKS)))*np.exp(-(degrees-L_AKS*(N_nodes-1))**2/(2*(N_nodes-1)*L_AKS*(1-L_AKS)))
axs[0,1].plot(degrees, p, color='r', label='MF theor')
axs[0,1].hist(deg_seq_PS, bins = degrees, color='b', density=True, label='AKS sim')
axs[0,1].legend(loc='upper left', prop={'size': 11})
axs[0,1].grid()

title = r"$\rho=$" + str(round(tau1,4)) + r"$,\ \tilde\mu=$" + str(mu) + r"$;\ \tilde\lambda=$" + str(lam) + r"$;\ n=$" + str(N_nodes)
axs[0,1].set_title(title, fontsize=14)

file_name = "../data/AKS_model/degree_distributions/AKS__N_NODES_200__N_ITERS_PER_LINK_100__tau1_-1.195014__mu_3__lambda_1__N_AVRGNG_10000__INITIALIZE_RANDOMLY_0.csv"

mu = 3
lam = 1
tau1 = -1.195014

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_PS = []
for line in lines:
    deg_seq_PS.append(int(line.split('\n')[0]))
N_lines = np.array(lines).shape[0]
L_AKS = np.average(deg_seq_PS)
print(L_AKS)

# Theoretical
L_AKS /= (N_nodes-1)
sigma = np.sqrt((N_nodes-1)*L_AKS*(1-L_AKS))
x_min, x_max = L_AKS*(N_nodes-1) - 5*sigma, L_AKS*(N_nodes-1) + 5*sigma
degrees = np.arange(np.rint(x_min), np.rint(x_max)+1)
p = []

p = np.sqrt(1/(2*np.pi*(N_nodes-1)*L_AKS*(1-L_AKS)))*np.exp(-(degrees-L_AKS*(N_nodes-1))**2/(2*(N_nodes-1)*L_AKS*(1-L_AKS)))
axs[0,2].plot(degrees, p, color='r', label='MF theor')
axs[0,2].hist(deg_seq_PS, bins = degrees, color='b', density=True, label='AKS sim')
axs[0,2].legend(loc='upper left', prop={'size': 11})
axs[0,2].grid()

title = r"$\rho=$" + str(round(tau1,4)) + r"$,\ \tilde\mu=$" + str(mu) + r"$;\ \tilde\lambda=$" + str(lam) + r"$;\ n=$" + str(N_nodes)
axs[0,2].set_title(title, fontsize=14)

axs[0,0].set_ylabel(r'$\langle \mathrm{p}(k|\mathbf{A}) \rangle$', fontsize=20)
axs[0,0].set_xlabel(r'$k$', fontsize=20)
axs[0,1].set_xlabel(r'$k$', fontsize=20)
axs[0,2].set_xlabel(r'$k$', fontsize=20)


#############################################################################
############################ LCC DISTRIBUTIONS ##############################
#############################################################################

file_name = "../data/AKS_model/clustering_distributions/AKS__N_NODES_100__N_ITERS_PER_LINK_10__tau1_50.586903__mu_-200__lambda_1__N_AVRGNG_10000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_PS = []
for line in lines:
    deg_seq_PS.append(float(line.split('\n')[0]))
    
N_lines = np.array(lines).shape[0]
print(np.average(deg_seq_PS))

L_sol = .3
cc_theor = L_sol*(1-(1-L_sol)**(N_nodes-1)-L_sol*(N_nodes-1)*(1-L_sol)**(N_nodes-2))
print("C_theor=", cc_theor)

axs[1,0].hist(deg_seq_PS, 100, color='b', density=True, label="AKS sim")
axs[1,0].grid()
axs[1,0].axvline(cc_theor, color='r', alpha=.65, label="MF theor")
axs[1,0].legend(fontsize=11, loc="upper left")
axs[1,0].set_ylabel(r'$\langle \mathrm{p}(c|\mathbf{A}) \rangle$', fontsize=20)
axs[1,0].set_xlabel(r'$c$', fontsize=20)
axs[1,0].set_xlim([.15, .45])


file_name = "../data/AKS_model/clustering_distributions/AKS__N_NODES_100__N_ITERS_PER_LINK_10__tau1_0.086457__mu_-2__lambda_1__N_AVRGNG_10000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_PS = []
for line in lines:
    deg_seq_PS.append(float(line.split('\n')[0]))
    
N_lines = np.array(lines).shape[0]
print(np.average(deg_seq_PS))

L_sol = .3
cc_theor = L_sol*(1-(1-L_sol)**(N_nodes-1)-L_sol*(N_nodes-1)*(1-L_sol)**(N_nodes-2))
print("C_theor=", cc_theor)

axs[1,1].hist(deg_seq_PS, 100, color='b', density=True, label="AKS sim")
axs[1,1].grid()
axs[1,1].axvline(cc_theor, color='r', alpha=.65, label="MF theor")
axs[1,1].legend(fontsize=11, loc="upper left")
axs[1,1].set_xlabel(r'$c$', fontsize=20)
axs[1,1].set_xlim([.15, .45])

file_name = "../data/AKS_model/clustering_distributions/AKS__N_NODES_100__N_ITERS_PER_LINK_10__tau1_-1.188807__mu_3__lambda_1__N_AVRGNG_10000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_PS = []
for line in lines:
    deg_seq_PS.append(float(line.split('\n')[0]))
    
N_lines = np.array(lines).shape[0]
print(np.average(deg_seq_PS))

L_sol = .3
cc_theor = L_sol*(1-(1-L_sol)**(N_nodes-1)-L_sol*(N_nodes-1)*(1-L_sol)**(N_nodes-2))
print("C_theor=", cc_theor)

axs[1,2].hist(deg_seq_PS, 100, color='b', density=True, label="AKS sim")
axs[1,2].grid()
axs[1,2].axvline(cc_theor, color='r', alpha=.65, label="MF theor")
axs[1,2].legend(fontsize=11, loc="upper left")
axs[1,2].set_xlabel(r'$c$', fontsize=20)
axs[1,2].set_xlim([.15, .45])


plt.tight_layout()

plt.savefig('../figures/deg_and_lcc_distr.png', dpi=300)
plt.show()

fig.clear()
