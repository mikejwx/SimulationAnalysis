import numpy as np
import matplotlib.pyplot as plt

# Attempting to make an aesthetically appealing schematic for our CT warm plume
# Make some coordinate system
x, y = np.meshgrid(np.arange(1184)*0.1, np.arange(320)*0.1)
x_c, y_c = [108.0, 16.0]
island_area = 50.0
island_radius = np.sqrt(island_area/np.pi)
R = np.sqrt((x - x_c)**2 + (y - y_c)**2)
vp = 2.0/1000. # across-flow anomaly, km/s
u = 10.0/1000. # along-flow wind, km/s
dt = island_radius/vp # time to meet in the middle, s
L = u*dt # warm plume length, km

fig = plt.figure()
axa = fig.add_subplot(2, 1, 1, adjustable = 'box', aspect = 1)
axa.set_ylabel('y (km)')
axa.set_xticklabels([''])
axa.contour(x, y, R, levels = [island_radius], colors = ['k'], linewidths = [2])
# Blue arrow with label: the mean flow
axa.arrow(x = x_c - 1.5*island_radius, y = y_c, dx = - L, dy = 0, head_width = 2, color = 'b', length_includes_head = True)
axa.text(x = x_c - 1.5*island_radius - L, y = y_c, s = '$u \\approx$ 10 m s$^{-1}$', horizontalalignment = 'right', verticalalignment = 'center', color = 'b')
# Purple arrow with across-flow wind anomaly pointing up
axa.arrow(x = x_c - 3*island_radius, y = y_c - vp*dt - 2.5, dy = vp*dt, dx = 0, head_width = 1, color = 'purple', length_includes_head = True)
# Purple arrow with label: across-flow wind anomaly pointing down
axa.arrow(x = x_c - 3*island_radius, y = y_c + vp*dt + 2.5, dy = - vp*dt, dx = 0, head_width = 1, color = 'purple', length_includes_head = True)
axa.text(x = x_c - 3*island_radius - 2.5, y = y_c + 0.5*vp*dt + 2.5, s = '$v^{\prime} \\approx$ 2 m s$^{-1}$', horizontalalignment = 'right', verticalalignment = 'center', color = 'purple')
# Red arrow for island radius
axa.arrow(x = x_c, y = y_c, dx = 0, dy = island_radius, head_width = 1, color = 'r', length_includes_head = True)
axa.text(x = x_c, y = y_c + 1.25*island_radius, s = 'r $\\approx$ 4 km', color = 'r', horizontalalignment = 'center', verticalalignment = 'bottom')
# Figure label
axa.text(x = 62.5, y = 28.5, s = 'a)', verticalalignment = 'center', horizontalalignment = 'left')
axa.set_xlim([60, 118])

axb = fig.add_subplot(2, 1, 2, adjustable = 'box', aspect = 1)
axb.set_ylabel('y (km)')
axb.set_xlabel('x (km)')
axb.contour(x, y, R, levels = [island_radius], colors = ['k'], linewidths = [2])
# Red triangle to highlight the warm plume
axb.fill([x_c, x_c, x_c - L, x_c], [y_c - island_radius, y_c + island_radius, y_c, y_c - island_radius], color = 'r', alpha = 2./3.)
# Blue dashed double-headed arrow 
axb.arrow(x = x_c - 0.5*L, y = y_c - island_radius - 1, dx = - 0.5*L, dy = 0, head_width = 1, color = 'b', ls = ':', length_includes_head = True)
axb.arrow(x = x_c - 0.5*L, y = y_c - island_radius - 1, dx = 0.5*L, dy = 0, head_width = 1, color = 'b', ls = ':', length_includes_head = True)
axb.text(x = x_c - 0.5*L, y = y_c - island_radius - 2.5, s = '$L \\approx u \\times \Delta t$\nLength of Warm Plume', verticalalignment = 'top', horizontalalignment = 'center', color = 'b')
# Purple dashed double-headed arrow
axb.arrow(x = x_c - L, y = y_c, dx = 0 , dy = - island_radius, head_width = 1, color = 'purple', ls = ':', length_includes_head = True)
axb.arrow(x = x_c - L, y = y_c, dx = 0 , dy = island_radius, head_width = 1, color = 'purple', ls = ':', length_includes_head = True)
axb.text(x = x_c - L - 2.5, y = y_c, s = '$\Delta t \\approx \\frac{r}{v^{\prime}}$', color = 'purple', verticalalignment = 'center', horizontalalignment = 'right', fontsize = 16)
# Figure label
axb.text(x = 62.5, y = 28.5, s = 'b)', verticalalignment = 'center', horizontalalignment = 'left')
axb.set_xlim([60, 118])
plt.savefig('../Ch5_Figure05.png', dpi = 250, bbox_inches = 'tight')
plt.show()

