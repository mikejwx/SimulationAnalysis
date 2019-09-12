import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import ctz_key, temp_key, theta_key
from SkewT_archer import skew

# read in the cloud top height data
path_0100m = '/nerc/n02/n02/xb899100/CloudTrail/Control/'
path_0800m = '/nerc/n02/n02/xb899100/CloudTrail/Control_0800m/'

ctz_0100_nc = Dataset(path_0100m + 'zi_09.nc', 'r')
ctz_0100_data = ctz_0100_nc.variables[ctz_key][:]*1.
ctz_0100_data = ctz_0100_data[~np.isnan(ctz_0100_data)].flatten()

ctz_0800_nc = Dataset(path_0800m + 'zi_04.nc', 'r')
ctz_0800_data = ctz_0800_nc.variables[ctz_key][:]*1.
ctz_0800_data = ctz_0800_data[~np.isnan(ctz_0800_data)].flatten()

# grab temperature data
temp_data = {}
temp_data['0100m'] = np.nanmean(Dataset(path_0100m + 'bouy_09.nc', 'r').variables[theta_key][:], axis = (0, 2, 3))
temp_data['0800m'] = np.nanmean(Dataset(path_0800m + 'bouy_04.nc', 'r').variables[theta_key][:], axis = (0, 2, 3))
z_the = Dataset(path_0100m + 'bouy_04.nc', 'r').variables['thlev_zsea_theta'][:]*1.

# plotting parameters
my_bins = np.arange(1, 51)*200.
factor1 = 300
factor2 = 50

# get the histogram data
hist100 = plt.hist(ctz_0100_data, bins = my_bins, normed = True)
hist800 = plt.hist(ctz_0800_data, bins = my_bins, normed = True)
plt.close('all')

y8 = 0.5*(hist800[1][:-1] + hist800[1][1:])

# make the plot

fig = plt.figure(tight_layout = True)
axa = fig.add_subplot(1, 2, 1)
axb = fig.add_subplot(1, 2, 2)

axa.plot(temp_data['0100m'], z_the, color = 'r')
axa.fill_betweenx(y = y8, x1 = np.zeros_like(hist100[0]) + factor1, x2 = hist100[0]*factor2/hist800[0].max()+factor1, color = 'b', alpha = 0.5, edgecolor = '')
axa.set_ylabel('height (m)')
axa.set_xlabel('Temperature (K)')
axa.set_title(u'100 m mean $\\theta$ profile\n histogram of z$_{top}$ \n max(z$_{top}$) = ' + str(int(ctz_0100_data.max())) + ' m')
axa.set_xlim([300, 400])

axb.plot(temp_data['0800m'], z_the, color = 'r')
axb.fill_betweenx(y = y8, x1 = np.zeros_like(hist800[0]) + factor1, x2 = hist800[0]*factor2/hist800[0].max()+factor1, color = 'b', alpha = 0.5, edgecolor = '')
axb.set_xlabel('Temperature (K)')
axb.set_title(u'800 m mean $\\theta$ profile\n histogram of z$_{top}$\n max(z$_{top}$) = ' + str(int(ctz_0800_data.max())) + ' m')
axb.set_xlim([300, 400])

plt.show()

