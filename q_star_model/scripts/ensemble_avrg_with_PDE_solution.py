#########################################################################
# This script plots ensemble averages of connectance for small networks #
# together with the corresponding numerical solutions of the Burgers'   #
# hierarchy and the transition formula for deep (0,1)-bistability.      #
#########################################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

t2_MF = 2
t3_MF = 1

title = "$t_2=$" + str(t2_MF) + "$;\ t_3=$" + str(t3_MF) + "$;\ n=\{5,10,15\}$"
fig, ax = plt.subplots(figsize=(10., 6.*1.1))

#### Numerical solution of PDEs ####
n_nodes = 5
x_min = -20
x_max = 20
nx = 10000
dx = (x_max-x_min)/nx
t2 = 2
nt2 = 50000
dt2 = t2/nt2
t3 = 1
nt3 = 100000
dt3 = t3/nt3
x = np.arange(x_min, x_max, dx)
u5 = np.zeros(nx)

nu = 1/(n_nodes*(n_nodes-1))

# Setting the initial conditions
u5 = np.exp(2*x)/(1+np.exp(2*x))

### Defining auxiliary functions ###
def d1x(y):
    D1x = np.zeros(nx)
    D1x[1:nx-1]=(y[2:nx]-y[0:nx-2])/(2*dx)
    D1x[0]=0
    D1x[nx-1]=0
    return D1x


def d2x(y):
    D2x = np.zeros(nx)
    D2x[1:nx-1]=(y[2:nx]-2*y[1:nx-1]+y[0:nx-2])/dx**2
    D2x[0]=0
    D2x[nx-1]=0
    return D2x

def F1(u5): # Right hand side of the 1st equation of the hierarchy
    return d1x(u5**2 + nu*d1x(u5))
def F1_(u5): # Right hand side of the 1st equation with t->-t
    return d1x(-u5**2 + nu*d1x(u5))

def F2(u5): # Right hand side of the 2nd equation of the hierarchy
    return d1x(u5**3 + 3*nu*u5*d1x(u5) + nu**2*d2x(u5))
def F2_(u5): # Right hand side of the 2nd equation with t->-t
    return d1x(-u5**3 + 3*nu*u5*d1x(u5) - nu**2*d2x(u5))

def RK4(u5, F, dt):
    K1 = F(u5)
    K2 = F(u5+.5*dt*K1)
    K3 = F(u5+.5*dt*K2)
    K4 = F(u5+dt*K3)
    return u5+dt/6.*(K1 + 2*K2 + 2*K3 + K4)

### RK4 time stepping for 1st eqn of hierarchy ###
for j in range(0, nt2):
    u5 = RK4(u5, F1, dt2)

### RK4 time stepping for 2nd eqn of hierarchy ###
for j in range(0, nt3):
    u5 = RK4(u5, F2, dt3)


n_nodes = 10
x_min = -20
x_max = 20
nx = 10000
dx = (x_max-x_min)/nx
t2 = 2
nt2 = 50000
dt2 = t2/nt2
t3 = 1
nt3 = 100000
dt3 = t3/nt3
x = np.arange(x_min, x_max, dx)
u10 = np.zeros(nx)

nu = 1/(n_nodes*(n_nodes-1))

# Setting the initial conditions
u10 = np.exp(2*x)/(1+np.exp(2*x))

### Defining auxiliary functions ###
def d1x(y):
    D1x = np.zeros(nx)
    D1x[1:nx-1]=(y[2:nx]-y[0:nx-2])/(2*dx)
    D1x[0]=0
    D1x[nx-1]=0
    return D1x

def d2x(y):
    D2x = np.zeros(nx)
    D2x[1:nx-1]=(y[2:nx]-2*y[1:nx-1]+y[0:nx-2])/dx**2
    D2x[0]=0
    D2x[nx-1]=0
    return D2x

def F1(u10): # Right hand side of the 1st equation of the hierarchy
    return d1x(u10**2 + nu*d1x(u10))
def F1_(u10): # Right hand side of the 1st equation with t->-t
    return d1x(-u10**2 + nu*d1x(u10))

def F2(u10): # Right hand side of the 2nd equation of the hierarchy
    return d1x(u10**3 + 3*nu*u10*d1x(u10) + nu**2*d2x(u10))
def F2_(u10): # Right hand side of the 2nd equation with t->-t
    return d1x(-u10**3 + 3*nu*u10*d1x(u10) - nu**2*d2x(u10))

def RK4(u10, F, dt):
    K1 = F(u10)
    K2 = F(u10+.5*dt*K1)
    K3 = F(u10+.5*dt*K2)
    K4 = F(u10+dt*K3)
    return u10+dt/6.*(K1 + 2*K2 + 2*K3 + K4)


### RK4 time stepping for 1st eqn of hierarchy ###
for j in range(0, nt2):
    u10 = RK4(u10, F1, dt2)

### RK4 time stepping for 2nd eqn of hierarchy ###
for j in range(0, nt3):
    u10 = RK4(u10, F2, dt3)


n_nodes = 15
x_min = -20
x_max = 20
nx = 10000
dx = (x_max-x_min)/nx
t2 = 2
nt2 = 50000
dt2 = t2/nt2
t3 = 1
nt3 = 100000
dt3 = t3/nt3
x = np.arange(x_min, x_max, dx)
u15 = np.zeros(nx)

nu = 1/(n_nodes*(n_nodes-1))

# Setting the initial conditions
u15 = np.exp(2*x)/(1+np.exp(2*x))

### Defining auxiliary functions ###
def d1x(y):
    D1x = np.zeros(nx)
    D1x[1:nx-1]=(y[2:nx]-y[0:nx-2])/(2*dx)
    D1x[0]=0
    D1x[nx-1]=0
    return D1x


def d2x(y):
    D2x = np.zeros(nx)
    D2x[1:nx-1]=(y[2:nx]-2*y[1:nx-1]+y[0:nx-2])/dx**2
    D2x[0]=0
    D2x[nx-1]=0
    return D2x

def F1(u15): # Right hand side of the 1st equation of the hierarchy
    return d1x(u15**2 + nu*d1x(u15))
def F1_(u15): # Right hand side of the 1st equation with t->-t
    return d1x(-u15**2 + nu*d1x(u15))


def F2(u15): # Right hand side of the 2nd equation of the hierarchy
    return d1x(u15**3 + 3*nu*u15*d1x(u15) + nu**2*d2x(u15))
def F2_(u15): # Right hand side of the 2nd equation with t->-t
    return d1x(-u15**3 + 3*nu*u15*d1x(u15) - nu**2*d2x(u15))

    
def RK4(u15, F, dt):
    K1 = F(u15)
    K2 = F(u15+.5*dt*K1)
    K3 = F(u15+.5*dt*K2)
    K4 = F(u15+dt*K3)
    return u15+dt/6.*(K1 + 2*K2 + 2*K3 + K4)


### RK4 time stepping for 1st eqn of hierarchy ###
for j in range(0, nt2):
    u15 = RK4(u15, F1, dt2)

### RK4 time stepping for 2nd eqn of hierarchy ###
for j in range(0, nt3):
    u15 = RK4(u15, F2, dt3)


#############################################
#############################################
x_min = -3.5
x_max = -2.6



file_names = ["../data/p_star_model/t1_scanning/PS__CONNECTANCE_and_STARS__N_NODES_5__N_ITERS_PER_LINK_100__t1_start_-3.500000__t1_fin_-2.500000__t1_step_0.010000__t2_MF_2__t3_MF_1__N_AVRGNG_50000__RANDOM_INITIALISATION_0.csv"
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


##### TRANSITION FORMULA #####
t1_ = np.arange(x_min, x_max, .001)
L_ensemble_TF_5 = 1/2*(1 + np.tanh(N_links_max*(t1_ + t2_MF + t3_MF)))



file_names = ["../data/p_star_model/t1_scanning/PS__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-3.200000__t1_fin_-3.120000__t1_step_0.004000__t2_MF_2__t3_MF_1__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv",
              "../data/p_star_model/t1_scanning/PS__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-3.120000__t1_fin_-3.040000__t1_step_0.004000__t2_MF_2__t3_MF_1__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv",
              "../data/p_star_model/t1_scanning/PS__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-3.040000__t1_fin_-2.960000__t1_step_0.004000__t2_MF_2__t3_MF_1__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv",
              "../data/p_star_model/t1_scanning/PS__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-2.960000__t1_fin_-2.880000__t1_step_0.004000__t2_MF_2__t3_MF_1__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv",
              "../data/p_star_model/t1_scanning/PS__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-2.880000__t1_fin_-2.796000__t1_step_0.004000__t2_MF_2__t3_MF_1__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv"
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
n_nodes = 10
N_links_max = n_nodes*(n_nodes-1)/2

t1_PS_10 = []
L_avrg_PS_10 = []
L_RMS_PS_10 = []
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
        L_avrg.append(np.array(L).mean())
        L_RMS.append(np.sqrt(np.array(L).var()))
        t1.append(prev_t1)
    t1_PS_10.extend(t1)
    L_avrg_PS_10.extend(L_avrg)
    L_RMS_PS_10.extend(L_RMS)

##### TRANSITION FORMULA #####
L_ensemble_TF_10 = 1/2*(1 + np.tanh(N_links_max*(t1_ + t2_MF + t3_MF)))



file_names = ["../data/p_star_model/t1_scanning/PS__CONNECTANCE_and_STARS__N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-2.976000__t1_fin_-2.959000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv",
              "../data/p_star_model/t1_scanning/PS__CONNECTANCE_and_STARS__N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-2.992000__t1_fin_-2.976000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv",
              "../data/p_star_model/t1_scanning/PS__CONNECTANCE_and_STARS__N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-3.008000__t1_fin_-2.992000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv",
              "../data/p_star_model/t1_scanning/PS__CONNECTANCE_and_STARS__N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-3.024000__t1_fin_-3.008000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv",
              "../data/p_star_model/t1_scanning/PS__CONNECTANCE_and_STARS__N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-3.040000__t1_fin_-3.024000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv"
]

n_nodes = 15
N_links_max = n_nodes*(n_nodes-1)/2

t1_PS_15 = []
L_avrg_PS_15 = []
L_RMS_PS_15 = []
L_err_PS_15_squared = []

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
    t1_PS_15.extend(t1)
    L_avrg_PS_15.extend(L_avrg)
    L_RMS_PS_15.extend(L_RMS)
    L_err_PS_15_squared.extend((np.array(L_var_0)*np.array(N_0) + np.array(L_var_1)*np.array(N_1))/(np.array(N_0)+np.array(N_1)))

##### TRANSITION FORMULA #####
L_ensemble_TF_15 = 1/2*(1 + np.tanh(N_links_max*(t1_ + t2_MF + t3_MF)))

file_names = ["../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_5__N_ITERS_PER_LINK_100__t1_start_-3.500000__t1_fin_-2.500000__t1_step_0.010000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_10__N_AVRGNG_50000__RANDOM_INITIALISATION_0.csv"]

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



file_names = ["../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-3.200000__t1_fin_-3.120000__t1_step_0.004000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_45__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv",
              "../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-3.120000__t1_fin_-3.040000__t1_step_0.004000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_45__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv",
              "../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-3.040000__t1_fin_-2.960000__t1_step_0.004000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_45__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv",
              "../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-2.960000__t1_fin_-2.880000__t1_step_0.004000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_45__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv",
              "../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_10__N_ITERS_PER_LINK_10000__t1_start_-2.880000__t1_fin_-2.796000__t1_step_0.004000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_45__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv"
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



file_names = ["../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-3.040000__t1_fin_-3.024000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv",
              "../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-3.024000__t1_fin_-3.008000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv",
              "../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-3.008000__t1_fin_-2.992000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv",
              "../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-2.992000__t1_fin_-2.976000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv",
              "../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_15__N_ITERS_PER_LINK_10000__t1_start_-2.976000__t1_fin_-2.959000__t1_step_0.002000__t2_MF_2__t3_MF_1__N_DYNAMIC_PAIRS_105__N_AVRGNG_1000__RANDOM_INITIALISATION_0.csv"]

n_nodes = 15
N_links_max = n_nodes*(n_nodes-1)/2

t1_MF_15 = []
L_avrg_MF_15 = []
L_RMS_MF_15 = []
L_err_MF_15_squared = []

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
    t1_MF_15.extend(t1)
    L_avrg_MF_15.extend(L_avrg)
    L_RMS_MF_15.extend(L_RMS)
    L_err_MF_15_squared.extend((np.array(L_var_0)*np.array(N_0) + np.array(L_var_1)*np.array(N_1))/(np.array(N_0)+np.array(N_1)))


# Plotting self-consistency equation
delta_t1 = 0.05
delta_L = 0.000005
t1_range = np.arange(x_min, x_max+delta_t1, delta_t1)
L_range = np.arange(0, 1, delta_L)
L, T1 = np.meshgrid(L_range, t1_range)
SC = L - 1/2*(1 + np.tanh(T1 + 2*t2_MF*L + 3*t3_MF*L**2)) #Self-consistency
CS = plt.contour(T1, L, SC, [0], colors='cyan', zorder=50, linewidths=1.5)
CS.collections[0].set_label('Self-consistency equation')

# plt.title(file_names[0])
plt.axvline(0, color='black')
plt.axhline(0, color='black')

plt.ylim([-.01, 1.01])
plt.xlim([x_min, x_max])

plt.title(title, fontsize=22)
plt.xlabel(r'$t_1$', fontsize=22)
plt.ylabel(r'$\bar L$', fontsize=22)


plt.grid()

plt.plot(x, u5, color="red", label="PDEs numerical solution")
plt.plot(x, u10, color="red")
plt.plot(x, u15, color="red")

# plt.plot(t1_PS_5, L_avrg_PS_5, 'h',markersize=4, label="3-star simulation", zorder=100, color='blue', alpha=.5)
# plt.errorbar(t1_PS_5, L_avrg_PS_5, xerr=None, yerr=[np.minimum(np.sqrt(L_RMS_PS_5[0]**2+(np.array(L_avrg_PS_5)*(1-np.array(L_avrg_PS_5))/100)**2), np.array(L_avrg_PS_5)), np.minimum(np.sqrt(L_RMS_PS_5[0]**2+(np.array(L_avrg_PS_5)*(1-np.array(L_avrg_PS_5))/100)**2), 1-np.array(L_avrg_PS_5))], color="blue", marker='o', markersize=4, ls='none', label="3-star simulation", alpha=.5, elinewidth=1, capsize=2, zorder=50)
plt.errorbar(t1_PS_5, L_avrg_PS_5, xerr=None, yerr=[np.minimum(np.sqrt(L_err_PS_5_squared+(np.array(L_avrg_PS_5)*(1-np.array(L_avrg_PS_5)))**2/50000), np.array(L_avrg_PS_5)), np.minimum(np.sqrt(L_err_PS_5_squared+(np.array(L_avrg_PS_5)*(1-np.array(L_avrg_PS_5)))**2/50000), 1-np.array(L_avrg_PS_5))], color="blue", marker='o', markersize=4, ls='none', label="3-star simulation, n=5", alpha=.5, elinewidth=1, capsize=2, zorder=50)
plt.plot(t1_, L_ensemble_TF_5, color='g', label='Transition formula')

plt.errorbar(t1_PS_10, L_avrg_PS_10, xerr=None, yerr=[np.minimum(np.sqrt(L_err_PS_10_squared+(np.array(L_avrg_PS_10)*(1-np.array(L_avrg_PS_10)))**2/1000), np.array(L_avrg_PS_10)), np.minimum(np.sqrt(L_err_PS_10_squared+(np.array(L_avrg_PS_10)*(1-np.array(L_avrg_PS_10)))**2/1000), 1-np.array(L_avrg_PS_10))], color="royalblue", marker='o', markersize=4, ls='none', label="3-star simulation, n=10", alpha=.5, elinewidth=1, capsize=2, zorder=50)#plot(t1_PS_10, L_avrg_PS_10, 'h',markersize=4, zorder=101, color='blue', alpha=.5)
plt.plot(t1_, L_ensemble_TF_10, color='g')

#plt.plot(t1_PS_15, L_avrg_PS_15, 'h',markersize=4, zorder=101, color='blue', alpha=.5)
plt.errorbar(t1_PS_15, L_avrg_PS_15, xerr=None, yerr=[np.minimum(np.sqrt(L_err_PS_15_squared+(np.array(L_avrg_PS_15)*(1-np.array(L_avrg_PS_15)))**2/1000), np.array(L_avrg_PS_15)), np.minimum(np.sqrt(L_err_PS_15_squared+(np.array(L_avrg_PS_15)*(1-np.array(L_avrg_PS_15)))**2/1000), 1-np.array(L_avrg_PS_15))], color="navy", marker='o', markersize=4, ls='none', label="3-star simulation, n=15", alpha=.5, elinewidth=1, capsize=2, zorder=50)
plt.plot(t1_, L_ensemble_TF_15, color='g')



plt.errorbar(t1_MF_5, L_avrg_MF_5, xerr=None, yerr=[np.minimum(np.sqrt(L_err_MF_5_squared+(np.array(L_avrg_MF_5)*(1-np.array(L_avrg_MF_5)))**2/50000), np.array(L_avrg_MF_5)), np.minimum(np.sqrt(L_err_MF_5_squared+(np.array(L_avrg_MF_5)*(1-np.array(L_avrg_MF_5)))**2/50000), 1-np.array(L_avrg_MF_5))], color="orange", marker='s', markersize=4, ls='none', label="MF3 simulation, n=5", alpha=.5, elinewidth=1, capsize=2, zorder=50)

plt.errorbar(t1_MF_10, L_avrg_MF_10, xerr=None, yerr=[np.minimum(np.sqrt(L_err_MF_10_squared+(np.array(L_avrg_MF_10)*(1-np.array(L_avrg_MF_10)))**2/1000), np.array(L_avrg_MF_10)), np.minimum(np.sqrt(L_err_MF_10_squared+(np.array(L_avrg_MF_10)*(1-np.array(L_avrg_MF_10)))**2/1000), 1-np.array(L_avrg_MF_10))], color="lightcoral", marker='s', markersize=4, ls='none', label="MF3 simulation, n=10", alpha=.5, elinewidth=1, capsize=2, zorder=50)

plt.errorbar(t1_MF_15, L_avrg_MF_15, xerr=None, yerr=[np.minimum(np.sqrt(L_err_MF_15_squared+(np.array(L_avrg_MF_15)*(1-np.array(L_avrg_MF_15)))**2/1000), np.array(L_avrg_MF_15)), np.minimum(np.sqrt(L_err_MF_15_squared+(np.array(L_avrg_MF_15)*(1-np.array(L_avrg_MF_15)))**2/1000), 1-np.array(L_avrg_MF_15))], color="brown", marker='s', markersize=4, ls='none', label="MF3 simulation, n=15", alpha=.7, elinewidth=1, capsize=2, zorder=50)

plt.legend(loc='lower right', prop={'size': 15}).set_zorder(200)

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
axins = inset_axes(ax, width="50%", height="50%", bbox_to_anchor=(0.01, .12, .7, 1.55), bbox_transform=ax.transAxes, loc="lower left") #bbox_to_anchor=(0.64, .05, .7, 1.55)

CS = axins.contour(T1, L, SC, [0], colors='cyan', zorder=100, linewidths=2.5, alpha=.7)

#axins.errorbar(t1_PS_5, L_avrg_PS_5, xerr=None, yerr=[np.minimum(np.sqrt(L_err_PS_5_squared+(np.array(L_avrg_PS_5)*(1-np.array(L_avrg_PS_5)))**2/50000), np.array(L_avrg_PS_5)), np.minimum(np.sqrt(L_err_PS_5_squared+(np.array(L_avrg_PS_5)*(1-np.array(L_avrg_PS_5)))**2/50000), 1-np.array(L_avrg_PS_5))], color="blue", marker='o', markersize=4, ls='none', label="3-star simulation", alpha=.5, elinewidth=1, capsize=2, zorder=50)
axins.plot(t1_PS_5, L_avrg_PS_5, 'h',markersize=4, label="3-star simulation", zorder=100, color='blue', alpha=.5)
#axins.errorbar(t1_PS_5, L_avrg_PS_5, xerr=None, yerr=[np.minimum(np.sqrt(L_RMS_PS_5[0]**2+(np.array(L_avrg_PS_5)*(1-np.array(L_avrg_PS_5))/100)**2), np.array(L_avrg_PS_5)), np.minimum(np.sqrt(L_RMS_PS_5[0]**2+(np.array(L_avrg_PS_5)*(1-np.array(L_avrg_PS_5))/100)**2), 1-np.array(L_avrg_PS_5))], color="blue", marker='o', markersize=4, ls='none', label="3-star simulation", alpha=.5, elinewidth=1, capsize=2, zorder=50)
axins.plot(t1_, L_ensemble_TF_5, color='g')
#axins.errorbar(t1_PS_10, L_avrg_PS_10, xerr=None, yerr=[np.minimum(np.sqrt(L_err_PS_10_squared+(np.array(L_avrg_PS_10)*(1-np.array(L_avrg_PS_10)))**2/1000), np.array(L_avrg_PS_10)), np.minimum(np.sqrt(L_err_PS_10_squared+(np.array(L_avrg_PS_10)*(1-np.array(L_avrg_PS_10)))**2/1000), 1-np.array(L_avrg_PS_10))], color="royalblue", marker='o', markersize=4, ls='none', label="3-star simulation, n=15", alpha=.5, elinewidth=1, capsize=2, zorder=50)
axins.plot(t1_PS_10, L_avrg_PS_10, 'h',markersize=4, zorder=101, color='royalblue', alpha=.5)
axins.plot(t1_, L_ensemble_TF_10, color='g')
#axins.errorbar(t1_PS_15, L_avrg_PS_15, xerr=None, yerr=[np.minimum(np.sqrt(L_err_PS_15_squared+(np.array(L_avrg_PS_15)*(1-np.array(L_avrg_PS_15)))**2/1000), np.array(L_avrg_PS_15)), np.minimum(np.sqrt(L_err_PS_15_squared+(np.array(L_avrg_PS_15)*(1-np.array(L_avrg_PS_15)))**2/1000), 1-np.array(L_avrg_PS_15))], color="navy", marker='o', markersize=4, ls='none', label="3-star simulation, n=15", alpha=.5, elinewidth=1, capsize=2, zorder=50)
axins.plot(t1_PS_15, L_avrg_PS_15, 'h',markersize=4, zorder=101, color='navy', alpha=.5)
# axins.errorbar(t1_PS_15, L_avrg_PS_15, xerr=None, yerr=[np.minimum(np.sqrt(L_RMS_PS_15[0]**2+(np.array(L_avrg_PS_15)*(1-np.array(L_avrg_PS_15))/100)**2), np.array(L_avrg_PS_15)), np.minimum(np.sqrt(L_RMS_PS_15[0]**2+(np.array(L_avrg_PS_15)*(1-np.array(L_avrg_PS_15))/100)**2), 1-np.array(L_avrg_PS_15))], color="blue", marker='o', markersize=4, ls='none', label="3-star simulation", alpha=.5, elinewidth=1, capsize=2, zorder=50)
axins.plot(t1_, L_ensemble_TF_15, color='g')

#axins.errorbar(t1_MF_5, L_avrg_MF_5, xerr=None, yerr=[np.minimum(np.sqrt(L_err_MF_5_squared+(np.array(L_avrg_MF_5)*(1-np.array(L_avrg_MF_5)))**2/50000), np.array(L_avrg_MF_5)), np.minimum(np.sqrt(L_err_MF_5_squared+(np.array(L_avrg_MF_5)*(1-np.array(L_avrg_MF_5)))**2/50000), 1-np.array(L_avrg_MF_5))], color="orange", marker='s', markersize=4, ls='none', label="3-star simulation", alpha=.5, elinewidth=1, capsize=2, zorder=50)
axins.plot(t1_MF_5, L_avrg_MF_5, 's',markersize=4, label="MF3 simulation n=5", zorder=100, color='orange', alpha=.8)
# axins.errorbar(t1_MF_5, L_avrg_MF_5, xerr=None, yerr=[np.minimum(np.sqrt(L_RMS_MF_5[0]**2+(np.array(L_avrg_MF_5)*(1-np.array(L_avrg_MF_5))/100)**2), np.array(L_avrg_MF_5)), np.minimum(np.sqrt(L_RMS_MF_5[0]**2+(np.array(L_avrg_MF_5)*(1-np.array(L_avrg_MF_5))/100)**2), 1-np.array(L_avrg_MF_5))], color="blue", marker='o', markersize=4, ls='none', label="3-star simulation", alpha=.5, elinewidth=1, capsize=2, zorder=50)
#axins.errorbar(t1_MF_10, L_avrg_MF_10, xerr=None, yerr=[np.minimum(np.sqrt(L_err_MF_10_squared+(np.array(L_avrg_MF_10)*(1-np.array(L_avrg_MF_10)))**2/1000), np.array(L_avrg_MF_10)), np.minimum(np.sqrt(L_err_MF_10_squared+(np.array(L_avrg_MF_10)*(1-np.array(L_avrg_MF_10)))**2/1000), 1-np.array(L_avrg_MF_10))], color="lightcoral", marker='s', markersize=4, ls='none', label="3-star simulation", alpha=.5, elinewidth=1, capsize=2, zorder=50)
axins.plot(t1_MF_10, L_avrg_MF_10, 's',markersize=4, zorder=101, color='lightcoral', alpha=.5)
#axins.errorbar(t1_MF_15, L_avrg_MF_15, xerr=None, yerr=[np.minimum(np.sqrt(L_err_MF_15_squared+(np.array(L_avrg_MF_15)*(1-np.array(L_avrg_MF_15)))**2/1000), np.array(L_avrg_MF_15)), np.minimum(np.sqrt(L_err_MF_15_squared+(np.array(L_avrg_MF_15)*(1-np.array(L_avrg_MF_15)))**2/1000), 1-np.array(L_avrg_MF_15))], color="brown", marker='s', markersize=4, ls='none', label="3-star simulation", alpha=.5, elinewidth=1, capsize=2, zorder=50)
axins.plot(t1_MF_15, L_avrg_MF_15, 's',markersize=4, label="MF3 simulation n=15", zorder=100, color='brown', alpha=.5)

axins.plot(x, u5, color="red",label="PDEs numerical solution")
axins.plot(x, u10, color="red")
axins.plot(x, u15, color="red")

axins.axhline(0, color="black")

x1, x2, y1, y2 = -3.495, -2.97, -.0002, .015 # specify the limits
axins.set_xlim(x1, x2) # apply the x-limits
axins.set_ylim(y1, y2) # apply the y-limits
plt.yticks(visible=False)
plt.xticks(visible=False)

from mpl_toolkits.axes_grid1.inset_locator import mark_inset
mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="black", linewidth=1.5, zorder=500, alpha=.6)

plt.tight_layout()

plt.savefig('../figures/ensemble_avrg_with_PDE_solution.png', dpi=300)

plt.show()

fig.clear()
