import numpy as np
import matplotlib.pyplot as plt

### A comparison of the balanced initial conditions derived from the spin-up simulations ###
with open('../InitialFields_Moisture_Control.txt') as moisture10:
    moisture10 = moisture10.read()
with open('../InitialFields_Temperature_Control.txt') as temperature10:
    temperature10 = temperature10.read()
with open('../InitialFields_Wind_Control.txt') as wind10:
    wind10 = wind10.read()

mv_z10 = np.array([float(x) for x in moisture10.split('\n')[2].split(' ')[1].split(',')])/1000.
# I didn't put a new line for the moisture data, future versions will include the new line here
mv_data10 = np.array([float(x) for x in moisture10.split('\n')[4].split(' ')[1].split(',')])*100.

theta_z10 = np.array([float(x) for x in temperature10.split('\n')[2].split(' ')[1].split(',')])/1000.
theta_data10 = np.array([float(x) for x in temperature10.split('\n')[4].split(' ')[1].split(',')])

u_z10 = np.array([float(x) for x in wind10.split('\n')[2].split(' ')[1].split(',')])/1000.
u_data10 = np.array([float(x) for x in wind10.split('\n')[3].split(' ')[1].split(',')])
v_data10 = np.array([float(x) for x in wind10.split('\n')[4].split(' ')[1].split(',')])

with open('../InitialFields_Moisture_U05.txt') as moisture05:
    moisture05 = moisture05.read()
with open('../InitialFields_Temperature_U05.txt') as temperature05:
    temperature05 = temperature05.read()
with open('../InitialFields_Wind_U05.txt') as wind05:
    wind05 = wind05.read()

mv_z05 = np.array([float(x) for x in moisture05.split('\n')[2].split(' ')[1].split(',')])/1000.
mv_data05 = np.array([float(x) for x in moisture05.split('\n')[4].split(' ')[1].split(',')])*100.

theta_z05 = np.array([float(x) for x in temperature05.split('\n')[2].split(' ')[1].split(',')])/1000.
theta_data05 = np.array([float(x) for x in temperature05.split('\n')[4].split(' ')[1].split(',')])

u_z05 = np.array([float(x) for x in wind05.split('\n')[2].split(' ')[1].split(',')])/1000.
u_data05 = np.array([float(x) for x in wind05.split('\n')[3].split(' ')[1].split(',')])
v_data05 = np.array([float(x) for x in wind05.split('\n')[4].split(' ')[1].split(',')])

### Make the plot Comparison ###
fig = plt.figure(figsize=(18,6))

ax1 = fig.add_subplot(1, 3, 1)
ax1.plot(u_data10, u_z10, 'b', lw = 2)
ax1.plot(v_data10, u_z10, 'b--', lw = 2)
ax1.plot(u_data05, u_z05, 'g', lw = 2)
ax1.plot(v_data05, u_z05, 'g--', lw = 2)
ax1.set_xlim([-12, 2])
ax1.set_ylim([0, 10])
ax1.set_xlabel('u (m/s) and v (m/s)')
ax1.set_ylabel('height (km)')

ax2 = fig.add_subplot(1, 3, 2, sharey = ax1)
ax2.plot(mv_data10, mv_z10, 'b', lw = 2)
ax2.plot(mv_data05, mv_z05, 'g', lw = 2)
ax2.set_xlabel('RH (%)')
ax2.set_ylim([0, 10])

ax3 = fig.add_subplot(1, 3, 3, sharey = ax1)
ax3.plot(theta_data10, theta_z10, 'b', lw = 2, label = 'Control: U$_g$ = -10')
ax3.plot(theta_data05, theta_z05, 'g', lw = 2, label = 'U$_g$ = -5')
ax3.set_xlabel('$\\theta$ (K)')
ax3.set_ylim([0, 10])
ax3.set_xlim([298, 330])
plt.legend(loc = 2)
plt.savefig('../U10_U05_comparison.png', dpi = 150)
plt.show()

