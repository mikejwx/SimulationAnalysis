import numpy as np
import matplotlib.pyplot as plt

### A comparison of the balanced initial conditions derived from the spin-up simulations ###
<<<<<<< HEAD
### Read the 100 m data ###
with open('../InitialFields_Moisture_Control.txt') as moisture100m:
    moisture100m = moisture100m.read()
with open('../InitialFields_Temperature_Control.txt') as temperature100m:
    temperature100m = temperature100m.read()
with open('../InitialFields_Wind_Control.txt') as wind100m:
    wind100m = wind100m.read()

mv_z100m = np.array([float(x) for x in moisture100m.split('\n')[2].split(' ')[1].split(',')])/1000.
# I didn't put a new line for the moisture data, future versions will include the new line here
mv_data100m = np.array([float(x) for x in moisture100m.split('\n')[4].split(' ')[1].split(',')])*100.

theta_z100m = np.array([float(x) for x in temperature100m.split('\n')[2].split(' ')[1].split(',')])/1000.
theta_data100m = np.array([float(x) for x in temperature100m.split('\n')[4].split(' ')[1].split(',')])

u_z100m = np.array([float(x) for x in wind100m.split('\n')[2].split(' ')[1].split(',')])/1000.
u_data100m = np.array([float(x) for x in wind100m.split('\n')[3].split(' ')[1].split(',')])
v_data100m = np.array([float(x) for x in wind100m.split('\n')[4].split(' ')[1].split(',')])

### Read the 800 m data ###
with open('../InitialFields_Moisture_spinup_800m.txt') as moisture800m:
    moisture800m = moisture800m.read()
with open('../InitialFields_Temperature_spinup_800m.txt') as temperature800m:
    temperature800m = temperature800m.read()
with open('../InitialFields_Wind_spinup_800m.txt') as wind800m:
    wind800m = wind800m.read()

mv_z800m = np.array([float(x) for x in moisture800m.split('\n')[2].split(' ')[1].split(',')])/1000.
mv_data800m = np.array([float(x) for x in moisture800m.split('\n')[4].split(' ')[1].split(',')])*100.

theta_z800m = np.array([float(x) for x in temperature800m.split('\n')[2].split(' ')[1].split(',')])/1000.
theta_data800m = np.array([float(x) for x in temperature800m.split('\n')[4].split(' ')[1].split(',')])

u_z800m = np.array([float(x) for x in wind800m.split('\n')[2].split(' ')[1].split(',')])/1000.
u_data800m = np.array([float(x) for x in wind800m.split('\n')[3].split(' ')[1].split(',')])
v_data800m = np.array([float(x) for x in wind800m.split('\n')[4].split(' ')[1].split(',')])

### Read the 1.6 km data ###
with open('../InitialFields_Moisture_spinup_1600m.txt') as moisture1600m:
    moisture1600m = moisture1600m.read()
with open('../InitialFields_Temperature_spinup_1600m.txt') as temperature1600m:
    temperature1600m = temperature1600m.read()
with open('../InitialFields_Wind_spinup_1600m.txt') as wind1600m:
    wind1600m = wind1600m.read()

mv_z1600m = np.array([float(x) for x in moisture1600m.split('\n')[2].split(' ')[1].split(',')])/1000.
mv_data1600m = np.array([float(x) for x in moisture1600m.split('\n')[4].split(' ')[1].split(',')])*100.

theta_z1600m = np.array([float(x) for x in temperature1600m.split('\n')[2].split(' ')[1].split(',')])/1000.
theta_data1600m = np.array([float(x) for x in temperature1600m.split('\n')[4].split(' ')[1].split(',')])

u_z1600m = np.array([float(x) for x in wind1600m.split('\n')[2].split(' ')[1].split(',')])/1000.
u_data1600m = np.array([float(x) for x in wind1600m.split('\n')[3].split(' ')[1].split(',')])
v_data1600m = np.array([float(x) for x in wind1600m.split('\n')[4].split(' ')[1].split(',')])
=======
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

with open('../InitialFields_Moisture_spinup_1km.txt') as moisture05:
    moisture05 = moisture05.read()
with open('../InitialFields_Temperature_spinup_1km.txt') as temperature05:
    temperature05 = temperature05.read()
with open('../InitialFields_Wind_spinup_1km.txt') as wind05:
    wind05 = wind05.read()

mv_z05 = np.array([float(x) for x in moisture05.split('\n')[2].split(' ')[1].split(',')])/1000.
mv_data05 = np.array([float(x) for x in moisture05.split('\n')[4].split(' ')[1].split(',')])*100.

theta_z05 = np.array([float(x) for x in temperature05.split('\n')[2].split(' ')[1].split(',')])/1000.
theta_data05 = np.array([float(x) for x in temperature05.split('\n')[4].split(' ')[1].split(',')])

u_z05 = np.array([float(x) for x in wind05.split('\n')[2].split(' ')[1].split(',')])/1000.
u_data05 = np.array([float(x) for x in wind05.split('\n')[3].split(' ')[1].split(',')])
v_data05 = np.array([float(x) for x in wind05.split('\n')[4].split(' ')[1].split(',')])
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696

### Make the plot Comparison ###
fig = plt.figure(figsize=(18,6))

ax1 = fig.add_subplot(1, 3, 1)
<<<<<<< HEAD
ax1.plot(u_data100m, u_z100m, 'b', lw = 2)
ax1.plot(v_data100m, u_z100m, 'b--', lw = 2)
ax1.plot(u_data800m, u_z800m, 'g', lw = 2)
ax1.plot(v_data800m, u_z800m, 'g--', lw = 2)
ax1.plot(u_data1600m, u_z1600m, 'y', lw = 2)
ax1.plot(v_data1600m, u_z1600m, 'y--', lw = 2)
=======
ax1.plot(u_data10, u_z10, 'b', lw = 2)
ax1.plot(v_data10, u_z10, 'b--', lw = 2)
ax1.plot(u_data05, u_z05, 'g', lw = 2)
ax1.plot(v_data05, u_z05, 'g--', lw = 2)
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
ax1.set_xlim([-12, 2])
ax1.set_ylim([0, 10])
ax1.set_xlabel('u (m/s) and v (m/s)')
ax1.set_ylabel('height (km)')

ax2 = fig.add_subplot(1, 3, 2, sharey = ax1)
<<<<<<< HEAD
ax2.plot(mv_data100m, mv_z100m, 'b', lw = 2)
ax2.plot(mv_data800m, mv_z800m, 'g', lw = 2)
ax2.plot(mv_data1600m, mv_z1600m, 'y', lw = 2)
=======
ax2.plot(mv_data10, mv_z10, 'b', lw = 2)
ax2.plot(mv_data05, mv_z05, 'g', lw = 2)
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
ax2.set_xlabel('RH (%)')
ax2.set_ylim([0, 10])

ax3 = fig.add_subplot(1, 3, 3, sharey = ax1)
<<<<<<< HEAD
ax3.plot(theta_data100m, theta_z100m, 'b', lw = 2, label = 'Control: (100 m)')
ax3.plot(theta_data800m, theta_z800m, 'g', lw = 2, label = 'Control: (800 m)')
ax3.plot(theta_data1600m, theta_z1600m, 'y', lw = 2, label = 'Control: (1600 m)')
=======
ax3.plot(theta_data10, theta_z10, 'b', lw = 2, label = 'Control: (100 m)')
ax3.plot(theta_data05, theta_z05, 'g', lw = 2, label = 'Coarse: (1 km)')
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
ax3.set_xlabel('$\\theta$ (K)')
ax3.set_ylim([0, 10])
ax3.set_xlim([298, 330])
plt.legend(loc = 2)
<<<<<<< HEAD
plt.savefig('../spinup_plots/100m_800m_1600m_comparison.png', dpi = 150)
=======
plt.savefig('../U10_U05_comparison.png', dpi = 150)
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
plt.show()

