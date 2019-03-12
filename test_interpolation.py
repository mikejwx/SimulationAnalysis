import numpy as np
from multiprocessing import Pool
from datetime import datetime as dt
from netCDF4 import Dataset
from analysis_tools import bilinear_interpolation, get_cs_coords
import matplotlib.pyplot as plt
"""
Benchmarking and speed improvement for the bilinear_interpolation routine
"""

theta_key = u'STASH_m01s00i004'
bouy_nc = Dataset('/nerc/n02/n02/xb899100/CloudTrail/Control/bouy_09.nc', 'r')
z = bouy_nc.variables['thlev_zsea_theta'][:]
X = np.arange(0., 116000., 100.)
Y = np.arange(0., 31900., 100.)
X, Y = np.meshgrid(X, Y)

x_cs, y_cs = get_cs_coords(108000., 15950., 80., X, Y, max_r = 20000.)

start = dt.now()
theta_interpolated0 = bilinear_interpolation(X, Y, bouy_nc.variables[theta_key][0,:,:,:], x_cs, y_cs, kind = 2)
elapsed_time = [(dt.now() - start).total_seconds()]

print elapsed_time
fig = plt.figure(tight_layout=True)
ax = fig.add_subplot(1, 2, 1)
ax.contourf(theta_interpolated0, levels = np.arange(300., 350.1, 5.))

theta_data = bouy_nc.variables[theta_key][0,:,:,:]
def pool_interpolate(new_coord):
    #print '[' + dt.now().strftime('%H:%M:%S') +']'
    new_x, new_y = new_coord
    return bilinear_interpolation(X, Y, theta_data*1., [new_x], [new_y], kind = 2)

for n_proc in xrange(2, 49):
    print '[' + dt.now().strftime('%H:%M:%S') +'] n_proc = ' + str(n_proc)
    start = dt.now()
    p = Pool(processes = n_proc)
    theta_interpolated = p.map(pool_interpolate, [(x_cs[I], y_cs[I]) for I in xrange(len(x_cs))])#, chunksize = len(x_cs)/n_proc)
    p.close()
    p.join()
    elapsed_time.append((dt.now() - start).total_seconds())

theta_int = np.array([[theta_interpolated[j][i] for i in xrange(len(theta_interpolated[j]))] for j in xrange(len(theta_interpolated))])
ax1 = fig.add_subplot(1, 2, 2)
meow = ax1.contourf(np.transpose(theta_int[:,:,0]), levels = np.arange(300., 350.1, 5.))
fig.colorbar(meow, ax = ax1)
plt.show()

plt.plot([1] + range(2, 49), elapsed_time)
plt.ylim([0, 90])
plt.yticks(range(0, 91, 15))
plt.xlabel('Number of Processors (#)')
plt.ylabel('Time taken to interpolate (seconds)')
plt.show()

