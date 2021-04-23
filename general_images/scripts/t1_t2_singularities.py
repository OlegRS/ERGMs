####################################################################
# This script plots a phase diagram of MF3 model on t1_t2-plane. #
####################################################################
import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(10., 6.))
plt.axvline(0, color="black", alpha=.5)
plt.axhline(0, color="black", alpha=.5)


t3=-1
delta_t2 = 0.05
delta_L = 0.00001
t2_range = np.arange(0, 5+delta_t2, delta_t2)
L_range = np.arange(delta_L, 1, delta_L)
L, T2 = np.meshgrid(L_range, t2_range)
dt1_dL = 12*t3*L**3 - (12*t3-4*T2)*L**2 - 4*T2*L + 1 # dt1/dL
CS = plt.contour(T2, -2*T2*L-3*t3*L**2 + 1/2.*np.log(L/(1-L)), dt1_dL, [0], colors='green', zorder=2, linewidths=2)
CS.collections[0].set_label(r'$t_3=-1$'+' singularities')
plt.text(3.5, -2.9, r'$t_3=-1$', color='green', zorder=2, size=15)

t3=0
delta_t2 = 0.05
delta_L = 0.00001
t2_range = np.arange(0, 5+delta_t2, delta_t2)
L_range = np.arange(delta_L, 1, delta_L)
L, T2 = np.meshgrid(L_range, t2_range)
dt1_dL = 12*t3*L**3 - (12*t3-4*T2)*L**2 - 4*T2*L + 1 # dt1/dL
CS = plt.contour(T2, -2*T2*L-3*t3*L**2 + 1/2.*np.log(L/(1-L)), dt1_dL, [0], colors='blue', zorder=2, linewidths=2)
CS.collections[0].set_label(r'$t_3=0$'+' singularities')
plt.text(2.29, -2.86, r'$t_3=0$', color='blue', zorder=2, size=15)

t3=1
delta_t2 = 0.05
delta_L = 0.00001
t2_range = np.arange(-1, 5+delta_t2, delta_t2)
L_range = np.arange(delta_L, 1, delta_L)
L, T2 = np.meshgrid(L_range, t2_range)
dt1_dL = 12*t3*L**3 - (12*t3-4*T2)*L**2 - 4*T2*L + 1 # dt1/dL
CS = plt.contour(T2, -2*T2*L-3*t3*L**2 + 1/2.*np.log(L/(1-L)), dt1_dL, [0], colors='navy', zorder=2, linewidths=2)
CS.collections[0].set_label(r'$t_3=1$'+' singularities')
plt.text(.6, -2.34, r'$t_3=1$', color='navy', zorder=2, size=15)

t3=2
delta_t2 = 0.05
delta_L = 0.00001
t2_range = np.arange(-4, 5+delta_t2, delta_t2)
L_range = np.arange(delta_L, 1, delta_L)
L, T2 = np.meshgrid(L_range, t2_range)
dt1_dL = 12*t3*L**3 - (12*t3-4*T2)*L**2 - 4*T2*L + 1 # dt1/dL
CS = plt.contour(T2, -2*T2*L-3*t3*L**2 + 1/2.*np.log(L/(1-L)), dt1_dL, [0], colors='brown', zorder=2, linewidths=2)
CS.collections[0].set_label(r'$t_3=2$'+' singularities')
plt.text(-1.29, -1.34, r'$t_3=2$', color='brown', zorder=2, size=15)


#### Critical curve projection ####
L = np.arange(0.001,0.99, 1e-5)
t1 = (4*L-3)/(4*(1-L)**2) + 1/2.*np.log(L/(1-L))
t2 = (2-3*L)/(4*L*(1-L)**2)
plt.plot(t2, t1, color='red', label='Critical line', zorder=3)

L_ticks = np.arange(0.1,1, .1)
t1_ticks = (4*L_ticks-3)/(4*(1-L_ticks)**2) + 1/2.*np.log(L_ticks/(1-L_ticks))
t2_ticks = (2-3*L_ticks)/(4*L_ticks*(1-L_ticks)**2)
plt.plot(t2_ticks, t1_ticks, marker="|", color='red', linestyle='none')
for i in range(1,8):
    plt.text(t2_ticks[i]+.03, t1_ticks[i]+.1, r"$\langle L\rangle=$" + str(L_ticks[i])[:3], color='r', zorder=3)

plt.xlim([-3.5, 5])
plt.ylim([-4, 2.5])

plt.xlabel(r"$t_2$", fontsize=20)
plt.ylabel(r"$t_1$", fontsize=20)

plt.grid()

plt.legend(loc="best", prop={'size': 15})
plt.tight_layout()
plt.savefig("../figures/t1_t2_singularities.png", dpi=300)

plt.show()
