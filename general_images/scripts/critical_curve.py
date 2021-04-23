##############################################################
# This script plots critical curve on (t2,t3) plane together #
# with positions of other figures from the paper and the     #
# region of narrow connectance bands at low temperatures.    #
##############################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(10., 6.))
L = np.arange(0.001,0.99, 1e-5)
L_ = np.arange(0.1,1, .1)

t2 = (3*L**2-5*L+2)/(4*L*(1-L)**3)
t3 = (2*L-1)/(12*L**2*(1-L)**2)
t2_ = (3*L_**2-5*L_+2)/(4*L_*(1-L_)**3)
t3_ = (2*L_-1)/(12*L_**2*(1-L_)**2)


plt.axhline(0, color="black", alpha=.7)
plt.axvline(0, color="black", alpha=.7)


plt.plot(t2, t3, label="Critical curve", color='red', linewidth=2)
plt.plot(t2_, t3_, marker="|", color='red', linestyle='none')
for i in [0, 1, 7, 8]:
    plt.text(t2_[i], t3_[i], str(L_[i])[:3], color='r', zorder=100)


plt.fill_between(t2, np.maximum(t3, -t2/3,-t2/3), 1000, where = t2>0, color='springgreen', alpha=.3, label=r"Narrow bands bistability")


plt.plot(2, 1, marker='*', markersize=5, color='blue')
# plt.text(2+.2, 1+.2, "1", color='blue')

plt.plot(-75, 45, marker='*', markersize=5, color='blue', linestyle='none', label="Positions of other figures")
plt.text(-75+.5, 45+.5, "4", color='blue')

plt.plot(-7.5, 4.5, marker='*', markersize=5, color='blue')
plt.text(-7.5+.5, 4.5+.5, "3", color='blue')

plt.plot(-.75, .45, marker='*', markersize=5, color='blue')
# plt.text(-.75-3, .45-3, "2", color='blue', zorder=100)

plt.plot(1.8707482993197284, -0.7558578987150416, marker='*', markersize=5, color='blue')

plt.plot(.3, .1, marker='*', markersize=5, color='blue')
# plt.text(.3+.5, .1+.5, "5", color='blue')

plt.plot(-75, 20, marker='*', markersize=5, color='blue')
plt.text(-75+.6, 20-2, "7", color='blue')

plt.plot(-75, 25, marker='*', markersize=5, color='blue')
plt.text(-75+.6, 25-2, "8", color='blue')



plt.grid()

handles, labels = ax.get_legend_handles_labels()
plt.legend([handles[0],handles[2],handles[1]], [labels[0],labels[2],labels[1]], loc="upper right", prop={'size': 18})

plt.xlim([-76, 76])
plt.ylim([-76, 76])
plt.xlabel(r'$t_2$', fontsize=22)
plt.ylabel(r'$t_3$', fontsize=22)


from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

axins = zoomed_inset_axes(ax, 3.6*4, bbox_to_anchor=(0.03, 0.04, .6, .5), bbox_transform=ax.transAxes, loc="lower left")
axins.plot(t2, t3, color='red', linewidth=2)
for i in np.arange(2,7):
    axins.plot(t2_[i], t3_[i], marker="|", color='red', linestyle='none')
    axins.text(t2_[i], t3_[i], str(L_[i])[:3], color='r', zorder=100)

axins.fill_between(t2, np.maximum(t3, -t2/3,-t2/3), 1000, where = t2>0, color='springgreen', alpha=.3)

axins.plot(2, 1, marker='*', markersize=5, color='blue')
axins.text(2+.05, 1+.05, "1", color='blue')

plt.plot(-.75, .45, marker='*', markersize=5, color='blue')
plt.text(-.75+.05, .45+.05, "2", color='blue')

axins.plot(.3, .1, marker='*', markersize=5, color='blue')
axins.text(.3-.11, .1+.06, "5", color='blue')

plt.plot(1.8707482993197284, -0.7558578987150416, marker='*', markersize=5, color='blue')
axins.text(1.8707482993197284-.2*.7, -0.7558578987150416-.4*.7, "6", color='blue')


axins.axhline(0, color='black', alpha=.7)
axins.axvline(0, color='black', alpha=.7)

x1, x2, y1, y2 = -2.1, 2.5, -2.2, 2.2 # specify the limits
axins.set_xlim(x1, x2) # apply the x-limits
axins.set_ylim(y1, y2) # apply the y-limits
plt.grid()

mark_inset(ax, axins, loc1=2, loc2=1, fc="none", ec="0.5")

plt.tight_layout()

plt.savefig('../figures/critical_curve.png', dpi=300)

plt.show()
