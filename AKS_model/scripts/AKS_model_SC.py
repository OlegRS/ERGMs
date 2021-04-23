##################################################
# This script plots families of self-consistency #
# equations of the AKS model.                    #
##################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

fig, axs = plt.subplots(1, 3, figsize=(15., 6./1.3))

########################################################
x_min, x_max = -6, 7.5
delta_rho = .05
delta_L = .001
rho_range = np.arange(x_min, x_max+.5, delta_rho)
L_range = np.arange(0, 1+2*delta_L, delta_L)
L, RHO = np.meshgrid(L_range, rho_range)

axs[0].grid()
axs[0].set_xlim([x_min,x_max])
axs[0].set_ylim([-.01, 1.01])
axs[0].axvline(0, color="black")
axs[0].axhline(0, color="black")
axs[0].set_title(r"$\tilde\mu\in\{-50, -20, -10, 0, 10, 20, 50\}$" + r"$;\ \tilde\lambda=0.1$", fontsize=15)
axs[0].set_xlabel(r'$\rho$', fontsize=18)
axs[0].set_ylabel(r'$\bar L$', fontsize=18)

mu = -50
lam = .1
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axs[0].contour(RHO, L, SC, [0], colors='red', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = -20
lam = .1
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axs[0].contour(RHO, L, SC, [0], colors='orange', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = -10
lam = .1
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axs[0].contour(RHO, L, SC, [0], colors='lime', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = 0
lam = .1
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axs[0].contour(RHO, L, SC, [0], colors='green', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = 10
lam = .1
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axs[0].contour(RHO, L, SC, [0], colors='cyan', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = 20
lam = .1
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axs[0].contour(RHO, L, SC, [0], colors='blue', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = 50
lam = .1
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axs[0].contour(RHO, L, SC, [0], colors='violet', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

########################################################
x_min, x_max = -30, 50
delta_rho = .05
delta_L = .001
rho_range = np.arange(x_min, x_max+.5, delta_rho)
L_range = np.arange(0, 1+2*delta_L, delta_L)
L, RHO = np.meshgrid(L_range, rho_range)

axs[1].grid()
axs[1].set_xlim([x_min,x_max])
axs[1].set_ylim([-.01, 1.01])
axs[1].axvline(0, color="black")
axs[1].axhline(0, color="black")
axs[1].set_title(r"$\tilde\mu\in\{-50, -20, -10, 0, 10, 20, 50\}$" + r"$;\ \tilde\lambda=1$", fontsize=15)
axs[1].set_xlabel(r'$\rho$', fontsize=18)

mu = -50
lam = 1
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axs[1].contour(RHO, L, SC, [0], colors='red', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = -20
lam = 1
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axs[1].contour(RHO, L, SC, [0], colors='orange', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = -10
lam = 1
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axs[1].contour(RHO, L, SC, [0], colors='lime', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = 1
lam = 1
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axs[1].contour(RHO, L, SC, [0], colors='green', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = 10
lam = 1
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axs[1].contour(RHO, L, SC, [0], colors='cyan', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = 20
lam = 1
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axs[1].contour(RHO, L, SC, [0], colors='blue', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = 50
lam = 1
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axs[1].contour(RHO, L, SC, [0], colors='violet', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

##### ZOOM ######
axins = zoomed_inset_axes(axs[1], 6.9, bbox_to_anchor=(.72, .07, .5, .7), bbox_transform=axs[1].transAxes, loc="lower left")

x1, x2, y1, y2 = -4, -1, -.005, .06 # specify the limits
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.grid()
# plt.yticks(visible=False)
# plt.xticks(visible=False)

mark_inset(axs[1], axins, loc1=2, loc2=3, fc="none", ec="0.5")

mu = -50
lam = 1
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axins.contour(RHO, L, SC, [0], colors='red', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = -20
lam = 1
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axins.contour(RHO, L, SC, [0], colors='orange', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = -10
lam = 1
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axins.contour(RHO, L, SC, [0], colors='lime', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = 1
lam = 1
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axins.contour(RHO, L, SC, [0], colors='green', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = 10
lam = 1
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axins.contour(RHO, L, SC, [0], colors='cyan', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = 20
lam = 1
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axins.contour(RHO, L, SC, [0], colors='blue', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = 50
lam = 1
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axins.contour(RHO, L, SC, [0], colors='violet', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

axins.axhline(0, color='black')

########################################################
x_min, x_max = -40, 50
delta_rho = .05
delta_L = .001
rho_range = np.arange(x_min, x_max+.5, delta_rho)
L_range = np.arange(0, 1+2*delta_L, delta_L)
L, RHO = np.meshgrid(L_range, rho_range)

axs[2].grid()
axs[2].set_xlim([x_min,x_max])
axs[2].set_ylim([-.01, 1.01])
axs[2].axvline(0, color="black")
axs[2].axhline(0, color="black")
axs[2].set_title(r"$\tilde\mu\in\{-50, -20, -10, 0, 10, 20, 50\}$" + r"$;\ \tilde\lambda=2$", fontsize=15)
axs[2].set_xlabel(r'$\rho$', fontsize=18)


mu = -50
lam = 2
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axs[2].contour(RHO, L, SC, [0], colors='red', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = -20
lam = 2
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axs[2].contour(RHO, L, SC, [0], colors='orange', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = -10
lam = 2
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axs[2].contour(RHO, L, SC, [0], colors='lime', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = 1
lam = 2
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axs[2].contour(RHO, L, SC, [0], colors='green', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = 10
lam = 2
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axs[2].contour(RHO, L, SC, [0], colors='cyan', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = 20
lam = 2
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axs[2].contour(RHO, L, SC, [0], colors='blue', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = 50
lam = 2
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axs[2].contour(RHO, L, SC, [0], colors='violet', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

##### ZOOM ######
axins = zoomed_inset_axes(axs[2], 6.9, bbox_to_anchor=(.749, .07, .5, .7), bbox_transform=axs[2].transAxes, loc="lower left")

x1, x2, y1, y2 = -4, -1, -.005, .06 # specify the limits
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.grid()
# plt.yticks(visible=False)
# plt.xticks(visible=False)

mark_inset(axs[2], axins, loc1=2, loc2=3, fc="none", ec="0.5")

mu = -50
lam = 2
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axins.contour(RHO, L, SC, [0], colors='red', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = -20
lam = 2
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axins.contour(RHO, L, SC, [0], colors='orange', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = -10
lam = 2
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axins.contour(RHO, L, SC, [0], colors='lime', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = 1
lam = 2
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axins.contour(RHO, L, SC, [0], colors='green', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = 10
lam = 2
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axins.contour(RHO, L, SC, [0], colors='cyan', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = 20
lam = 2
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axins.contour(RHO, L, SC, [0], colors='blue', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

mu = 50
lam = 2
SC = L - 1/2*(1 + np.tanh(RHO + mu*lam - mu*lam*np.exp(-L/lam))) #Self-consistency
CS = axins.contour(RHO, L, SC, [0], colors='violet', zorder=100, linewidths=2.5, alpha=.7)
CS.collections[0].set_label('MF_theory')

axins.axhline(0, color='black')


##############################################################################

plt.tight_layout()
plt.savefig('../figures/AKS_model_SC.png', dpi=300)
plt.show()
