################################################################
# This script plots phase diagram of MF3 model on t2_t3-plane. #
################################################################
import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(10., 6.))
plt.axvline(0, color="black", alpha=.5)
plt.axhline(0, color="black", alpha=.5)

t1=-1.5
delta_t2 = 0.05
delta_L = 0.00001
t2_range = np.arange(-5, 2.5+delta_t2, delta_t2)
L_range = np.arange(delta_L, 1, delta_L)
L, T2 = np.meshgrid(L_range, t2_range)
dt1_dL = 4*(t1 + T2*L - 1/2.*np.log(L/(1-L)))*(1-L) + 1 # dt1/dL
CS = plt.contour(T2, (-t1-2*T2*L + 1/2.*np.log(L/(1-L)))/(3*L**2), dt1_dL, [0], colors='green', zorder=2, linewidths=2)
CS.collections[0].set_label(r'$t_1=-1.5$'+' singularities')
plt.text(1.62, 3.169, r'$t_1=-1.5$', color='green', zorder=2, size=15)

t1= -.5
delta_t2 = 0.05
delta_L = 0.00001
t2_range = np.arange(-5, 2+delta_t2, delta_t2)
L_range = np.arange(delta_L, 1, delta_L)
L, T2 = np.meshgrid(L_range, t2_range)
dt1_dL = 4*(t1 + T2*L - 1/2.*np.log(L/(1-L)))*(1-L) + 1 # dt1/dL
CS = plt.contour(T2, (-t1-2*T2*L + 1/2.*np.log(L/(1-L)))/(3*L**2), dt1_dL, [0], colors='blue', zorder=2, linewidths=2)
CS.collections[0].set_label(r'$t_1=-0.5$'+' singularities')
plt.text(-1.32, 3.43, r'$t_1=-0.5$', color='blue', zorder=2, size=15)

t1= 0.5
delta_t2 = 0.05
delta_L = 0.00001
t2_range = np.arange(-5, 2+delta_t2, delta_t2)
L_range = np.arange(delta_L, 1, delta_L)
L, T2 = np.meshgrid(L_range, t2_range)
dt1_dL = 4*(t1 + T2*L - 1/2.*np.log(L/(1-L)))*(1-L) + 1 # dt1/dL
CS = plt.contour(T2, (-t1-2*T2*L + 1/2.*np.log(L/(1-L)))/(3*L**2), dt1_dL, [0], colors='navy', zorder=2, linewidths=2)
CS.collections[0].set_label(r'$t_1=0.5$'+' singularities')
plt.text(-3.23, 4.08, r'$t_1=0.5$', color='navy', zorder=2, size=15)

t1= 1.5
delta_t2 = 0.05
delta_L = 0.00001
t2_range = np.arange(-5, 2+delta_t2, delta_t2)
L_range = np.arange(delta_L, 1, delta_L)
L, T2 = np.meshgrid(L_range, t2_range)
dt1_dL = 4*(t1 + T2*L - 1/2.*np.log(L/(1-L)))*(1-L) + 1 # dt1/dL
CS = plt.contour(T2, (-t1-2*T2*L + 1/2.*np.log(L/(1-L)))/(3*L**2), dt1_dL, [0], colors='brown', zorder=2, linewidths=2)
CS.collections[0].set_label(r'$t_1=1.5$'+' singularities')
plt.text(-4.85, 4.634, r'$t_1=1.5$', color='brown', zorder=2, size=15)


#### Critical curve projection ####
L = np.arange(0.001,0.99, 1e-5)
t2 = (2-3*L)/(4*L*(1-L)**2)
t3 = (2*L-1)/(12*L**2*(1-L)**2)
plt.plot(t2, t3, color='red', label='Critical line', zorder=3)

L_ticks = np.arange(0.1,1, .1)
t2_ticks = (2-3*L_ticks)/(4*L_ticks*(1-L_ticks)**2)
t3_ticks = (2*L_ticks-1)/(12*L_ticks**2*(1-L_ticks)**2)
plt.plot(t2_ticks, t3_ticks, marker="|", color='red', linestyle='none')
for i in range(1,8):
    plt.text(t2_ticks[i]-.5, t3_ticks[i]-.25, r"$\langle L\rangle=$" + str(L_ticks[i])[:3], color='r', zorder=3)

plt.xlim([-5, 3.5])
plt.ylim([-2.3, 5])

plt.xlabel(r"$t_2$", fontsize=20)
plt.ylabel(r"$t_3$", fontsize=20)

plt.grid()

plt.legend(loc="lower left", prop={'size': 15})
plt.tight_layout()
plt.savefig("../figures/t2_t3_singularities.png", dpi=300)

plt.show()

