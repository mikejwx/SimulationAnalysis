import numpy as np
import matplotlib.pyplot as plt
from analysis_tools import *
from SkewT_archer import *
from netCDF4 import Dataset
from STASH_keys import *

# Read the data for control and U05 simulations
main_path = '/nerc/n02/n02/xb899100/CloudTrail/'
U10_nc = Dataset(main_path + 'Control/zi_09.nc', 'r')
U05_nc = Dataset(main_path + 'U05/zi_04.nc', 'r')

lsm = Dataset('/work/n02/n02/xb899100/island_masks/lsm50.nc', 'r')
lsm_data = lsm.variables['lsm'][0,0,:,:]*1.
lsm.close()


X, Y = np.meshgrid(np.arange(0., 116., 0.1), np.arange(0., 31.9, 0.1))
my_levels = np.arange(600., 3000.1, 300.)

fig = plt.figure(tight_layout = True)
ax0 = fig.add_subplot(2, 1, 1, adjustable = 'box', aspect = 1)
U10 = ax0.contourf(X, Y, U10_nc.variables[ctz_key][-1,:,:], cmap = 'Greys', vmin = 0, levels = my_levels, extend = 'max')
fig.colorbar(U10, ax = ax0)
ax0.contour(X, Y, lsm_data, colors = ['darkred'], linewidths = [2])
ax0.set_xlabel('x (km)')
ax0.set_ylabel('y (km)')
ax0.set_title('Control')

ax1 = fig.add_subplot(2, 1, 2, adjustable = 'box', aspect = 1)
U05 = ax1.contourf(X, Y, U05_nc.variables[ctz_key][-1,:,:], cmap = 'Greys', vmin = 0, levels = my_levels, extend = 'max')
fig.colorbar(U05, ax = ax1)
ax1.contour(X, Y, lsm_data, colors = ['darkred'], linewidths = [2])
ax1.set_xlabel('x (km)')
ax1.set_ylabel('y (km)')
ax1.set_title('U05')

plt.savefig('../CloudTopHeight.png', dpi = 150)
plt.close('all')

U10_nc.close()
U05_nc.close()


lwp_U10_nc = Dataset(main_path + 'Control/lwp_00.nc', 'r')
lwp_U05_nc = Dataset(main_path + 'U05/lwp_00.nc', 'r')
U10_time_key = [tkey for tkey in lwp_U10_nc.variables.keys() if 'min' in tkey][0]
U05_time_key = [tkey for tkey in lwp_U05_nc.variables.keys() if 'min' in tkey][0]

U10_times = lwp_U10_nc.variables[U10_time_key][:]
U05_times = lwp_U05_nc.variables[U05_time_key][:] + 240.

preday_idx10 = np.where(U10_times == 360.)[0][0]
midday_idx10 = np.where(U10_times == 1080.)[0][0] + 1
preday_idx05 = np.where(U05_times == 360.)[0][0]
midday_idx05 = np.where(U05_times == 1080.)[0][0] + 1

dt10 = U10_times[1] - U10_times[0]
dt05 = U05_times[1] - U05_times[0]

if dt10 < dt05:
    step_10 = int(dt05/dt10)
    step_05 = 1
elif dt10 > dt05:
    step_10 = 1
    step_05 = int(dt10/dt05)
elif dt10 == dt05:
    step_10, step_05 = [1,1]

idx10 = range(preday_idx10, midday_idx10, step_10)
idx05 = range(preday_idx05, midday_idx05, step_05)

lwp_data_10 = np.nanmean(np.where(lwp_U10_nc.variables[lwp_key][idx10,:,:] > 0, 1., 0.), axis = 0)
lwp_data_05 = np.zeros_like(lwp_U10_nc.variables[lwp_key][idx10,:,:])
for it in range(len(idx05)):
    lwp_data_05[it,:,:] = lwp_U05_nc.variables[lwp_key][idx05[it],:,:]

lwp_data_05 = np.nanmean(np.where(lwp_data_05 > 0, 1., 0.), axis = 0)

# Plot of cloud frequency
fig = plt.figure(tight_layout = True)
# Make the plot
ax0 = fig.add_subplot(2, 1, 1, adjustable = 'box', aspect = 1)
lwp10_plt = ax0.contourf(X, Y, lwp_data_10, cmap = 'Greys_r', levels = np.linspace(0.0, 1.0, 11))
fig.colorbar(lwp10_plt, ax = ax0, label = 'Cloud Frequency')
# island mask
ax0.contour(X, Y, lsm_data, colors = ['darkred'], levels = [1e-16], linewidths = [2])
ax0.set_title('Control')

ax1 = fig.add_subplot(2, 1, 2, adjustable = 'box', aspect = 1)
lwp05_plt = ax1.contourf(X, Y, lwp_data_05, cmap = 'Greys_r', levels = np.linspace(0.0, 1.0, 11))
fig.colorbar(lwp05_plt, ax = ax1, label = 'Cloud Frequency')
ax1.contour(X, Y, lsm_data, colors = ['darkred'], levels = [1e-16], linewidths = [2])
ax1.set_title('U05')

plt.savefig('../CloudFrequency.png', dpi = 150.)
plt.show()


