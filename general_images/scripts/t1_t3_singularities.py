################################################################
# This script plots phase diagram of MF3 model on t1_t3-plane. #
################################################################
import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(10., 6.))
plt.axvline(0, color="black", alpha=.5)
plt.axhline(0, color="black", alpha=.5)

t2=-4
delta_t3 = 0.05
delta_L = 0.00001
t3_range = np.arange(0, 5+delta_t3, delta_t3)
L_range = np.arange(delta_L, 1, delta_L)
L, T3 = np.meshgrid(L_range, t3_range)
dt1_dL = 12*T3*L**3 - (12*T3-4*t2)*L**2 - 4*t2*L + 1 # dt1/dL
CS = plt.contour(T3, -2*t2*L-3*T3*L**2 + 1/2.*np.log(L/(1-L)), dt1_dL, [0], colors='orange', zorder=2, linewidths=2)
CS.collections[0].set_label(r'$t_3=-4$'+' singularities')
plt.text(3.1, .8, r'$t_3=-4$', color='orange', zorder=2, size=15)

t2=-2
delta_t3 = 0.05
delta_L = 0.00001
t3_range = np.arange(0, 5+delta_t3, delta_t3)
L_range = np.arange(delta_L, 1, delta_L)
L, T3 = np.meshgrid(L_range, t3_range)
dt1_dL = 12*T3*L**3 - (12*T3-4*t2)*L**2 - 4*t2*L + 1 # dt1/dL
CS = plt.contour(T3, -2*t2*L-3*T3*L**2 + 1/2.*np.log(L/(1-L)), dt1_dL, [0], colors='green', zorder=2, linewidths=2)
CS.collections[0].set_label(r'$t_3=-2$'+' singularities')
plt.text(2.15, -.4, r'$t_3=-2$', color='green', zorder=2, size=15)

t2=0
delta_t3 = 0.05
delta_L = 0.00001
t3_range = np.arange(0, 5+delta_t3, delta_t3)
L_range = np.arange(delta_L, 1, delta_L)
L, T3 = np.meshgrid(L_range, t3_range)
dt1_dL = 12*T3*L**3 - (12*T3-4*t2)*L**2 - 4*t2*L + 1 # dt1/dL
CS = plt.contour(T3, -2*t2*L-3*T3*L**2 + 1/2.*np.log(L/(1-L)), dt1_dL, [0], colors='blue', zorder=2, linewidths=2)
CS.collections[0].set_label(r'$t_3=0$'+' singularities')
plt.text(1.04, -1.26, r'$t_3=0$', color='blue', zorder=2, size=15)

t2=2
delta_t3 = 0.05
delta_L = 0.00001
t3_range = np.arange(-3, 5+delta_t3, delta_t3)
L_range = np.arange(delta_L, 1, delta_L)
L, T3 = np.meshgrid(L_range, t3_range)
dt1_dL = 12*T3*L**3 - (12*T3-4*t2)*L**2 - 4*t2*L + 1 # dt1/dL
CS = plt.contour(T3, -2*t2*L-3*T3*L**2 + 1/2.*np.log(L/(1-L)), dt1_dL, [0], colors='navy', zorder=2, linewidths=2)
CS.collections[0].set_label(r'$t_3=2$'+' singularities')
plt.text(0.096, -2.64, r'$t_3=2$', color='navy', zorder=2, size=15)

t2=4
delta_t3 = 0.05
delta_L = 0.00001
t3_range = np.arange(-5, 5+delta_t3, delta_t3)
L_range = np.arange(delta_L, 1, delta_L)
L, T3 = np.meshgrid(L_range, t3_range)
dt1_dL = 12*T3*L**3 - (12*T3-4*t2)*L**2 - 4*t2*L + 1 # dt1/dL
CS = plt.contour(T3, -2*t2*L-3*T3*L**2 + 1/2.*np.log(L/(1-L)), dt1_dL, [0], colors='brown', zorder=2, linewidths=2)
CS.collections[0].set_label(r'$t_3=4$'+' singularities')
plt.text(-1.5, -2.9, r'$t_3=4$', color='brown', zorder=2, size=15)

#### Critical curve projection ####
L = np.arange(0.001,0.99, 1e-5)
t1 = (4*L-3)/(4*(1-L)**2) + 1/2.*np.log(L/(1-L))
t3 = (2*L-1)/(12*L**2*(1-L)**2)
plt.plot(t3, t1, color='red', label='Critical line', zorder=3)

L_ticks = np.arange(0.1,1, .1)
t1_ticks = (4*L_ticks-3)/(4*(1-L_ticks)**2) + 1/2.*np.log(L_ticks/(1-L_ticks))
t3_ticks = (2*L_ticks-1)/(12*L_ticks**2*(1-L_ticks)**2)
plt.plot(t3_ticks, t1_ticks, marker="|", color='red', linestyle='none')
for i in range(1,8):
    plt.text(t3_ticks[i]-.63, t1_ticks[i]+.048, r"$\langle L\rangle=$" + str(L_ticks[i])[:3], color='r', zorder=3)

plt.xlim([-4.8, 4.8])
plt.ylim([-4, 3])

plt.xlabel(r"$t_3$", fontsize=20)
plt.ylabel(r"$t_1$", fontsize=20)

plt.grid()

plt.legend(loc="upper left", prop={'size': 15})
plt.tight_layout()
plt.savefig("../figures/t1_t3_singularities.png", dpi=300)

plt.show()
