"""
An investigation into the presence of cold pools in our simulated cloud trails
"""
import numpy as np
import matplotlib.pyplot as plt
from SkewT_archer import getThetaE
from netCDF4 import Dataset
from multiprocessing import Pool
from analysis_tools import send_email
import os

# create an animation of surface equivalent potential temperature
# Define keys for the variables needed
pressure_key = u'STASH_m01s00i408'
theta_key = u'STASH_m01s00i004'
temperature_key = u'STASH_m01s16i004'

#Define cartesian coordinate system
x = np.arange(0., 116000., 100.)
y = np.arange(0., 31900., 100.)
X, Y = np.meshgrid(x, y)

lsm = Dataset('/work/n02/n02/xb899100/island_masks/lsm50.nc', 'r')
lsm_data = lsm.variables['lsm'][:]
lsm.close()

for hour in ["{0:02d}".format(h) for h in xrange(0, 24, 3)]:
    # open the netCDFs
    bouy_nc = Dataset('../bouy_' + hour + '.nc', 'r')
    fluxes_nc = Dataset('../fluxes_' + hour + '.nc', 'r')
    
    if hour == '00':
        theta_e = getThetaE(bouy_nc.variables[theta_key][:,0,:,:], bouy_nc.variables[temperature_key][:,0,:,:], fluxes_nc.variables[pressure_key][:,0,:,:], p_units = 'Pa')
        times = bouy_nc.variables['min10_0'][:]
        z = bouy_nc.variables['thlev_zsea_theta'][:]
    else:
        theta_e = np.concatenate((theta_e, getThetaE(bouy_nc.variables[theta_key][:,0,:,:], bouy_nc.variables[temperature_key][:,0,:,:], fluxes_nc.variables[pressure_key][:,0,:,:], p_units = 'Pa')), axis = 0)
        times = np.concatenate((times, bouy_nc.variables['min10_0'][:]), axis = 0)
    bouy_nc.close()
    fluxes_nc.close()

if not os.path.isdir('../theta_e/'):
    os.mkdir('../theta_e/')

# Plot the domain mean equivalent potential temperature
theta_e_mean = np.nanmean(theta_e[1:,:,:], axis = (1, 2))
theta_e_land = np.nanmean(np.where((lsm_data[0,0,:,:] == 1.), theta_e[1:,:,:], np.nan), axis = (1,2))
theta_e_sea = np.nanmean(np.where((lsm_data[0,0,:,:] == 0.), theta_e[1:,:,:], np.nan), axis = (1,2))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(times[1:]/60., theta_e_mean, 'k', lw = 2, label = 'Whole')
ax.plot(times[1:]/60., theta_e_land, 'r', label = 'Land')
ax.plot(times[1:]/60., theta_e_sea, 'b', label = 'Sea')
ax.set_xticks(xrange(0, 24, 3))
ax.set_xlim([0, 24])
ax.set_xlabel('Time (hrs)')
ax.set_ylabel(r'$\theta_e$ (K)')
plt.legend(loc = 2)
plt.savefig('../theta_e/ha_theta_e_timeseries.png', dpi = 100)
plt.show()

X /= 1000.
Y /= 1000.
def my_plot(it):
    # Define some plotting arguments
    my_cmap = 'cool'
    my_island_contour = ['k']
    my_levels = np.linspace(367.5, 371., 15)

    fig = plt.figure(tight_layout = True, figsize = (15, 5))
    ax = fig.add_subplot(1, 1, 1, adjustable = 'box', aspect = 'equal')
    #top left
    ax.contourf(X-np.max(X), Y+np.max(Y), theta_e[it,:,:], cmap = my_cmap, levels = my_levels, extend = 'both')
    ax.contour(X-np.max(X), Y+np.max(Y), lsm_data[0,0,:,:], levels = [1e-08], colors = my_island_contour)
    #top centre
    ax.contourf(X, Y+np.max(Y), theta_e[it,:,:], cmap = my_cmap, levels = my_levels, extend = 'both')
    ax.contour(X, Y+np.max(Y), lsm_data[0,0,:,:], levels = [1e-08], colors = my_island_contour)
    #top right
    ax.contourf(X+np.max(X), Y+np.max(Y), theta_e[it,:,:], cmap = my_cmap, levels = my_levels, extend = 'both')
    ax.contour(X+np.max(X), Y+np.max(Y), lsm_data[0,0,:,:], levels = [1e-08], colors = my_island_contour)
    #centre left
    ax.contourf(X-np.max(X), Y, theta_e[it,:,:], cmap = my_cmap, levels = my_levels, extend = 'both')
    ax.contour(X-np.max(X), Y, lsm_data[0,0,:,:], levels = [1e-08], colors = my_island_contour)
    #centre centre
    ax.contourf(X, Y, theta_e[it,:,:], cmap = my_cmap, levels = my_levels, extend = 'both')
    ax.contour(X, Y, lsm_data[0,0,:,:], levels = [1e-08], colors = my_island_contour)
    #centre right
    the = ax.contourf(X+np.max(X), Y, theta_e[it,:,:], cmap = my_cmap, levels = my_levels, extend = 'both')
    ax.contour(X+np.max(X), Y, lsm_data[0,0,:,:], levels = [1e-08], colors = my_island_contour)
    plt.colorbar(the, ax = ax, label = r'$\theta_e$ (K)')
    #bottom left
    ax.contourf(X-np.max(X), Y-np.max(Y), theta_e[it,:,:], cmap = my_cmap, levels = my_levels, extend = 'both')
    ax.contour(X-np.max(X), Y-np.max(Y), lsm_data[0,0,:,:], levels = [1e-08], colors = my_island_contour)
    #bottom centre
    ax.contourf(X, Y-np.max(Y), theta_e[it,:,:], cmap = my_cmap, levels = my_levels, extend = 'both')
    ax.contour(X, Y-np.max(Y), lsm_data[0,0,:,:], levels = [1e-08], colors = my_island_contour)
    #bottom right
    ax.contourf(X+np.max(X), Y-np.max(Y), theta_e[it,:,:], cmap = my_cmap, levels = my_levels, extend = 'both')
    ax.contour(X+np.max(X), Y-np.max(Y), lsm_data[0,0,:,:], levels = [1e-08], colors = my_island_contour)

    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    # Plot the original domain
    ax.plot([0, 0, np.max(X), np.max(X), 0], [0, np.max(Y), np.max(Y), 0, 0], 'k--')
    ax.set_title(r'$\theta_e$ at ' + str(round(z[0],2)) + ' m and T+' + "{0:04d}".format(int(times[it]))+' mins')
    plt.savefig('../theta_e/theta_e_T'+ "{0:04d}".format(int(times[it])) + '.png', dpi = 100)
    plt.close('all')

p = Pool()
p.map(my_plot, range(1,len(times)))
plt.show()
p.close()

my_command = 'convert -delay 30 -loop 0 ../theta_e/*.png ../theta_e/theta_e_anim.gif'
os.system(my_command)

send_email(message = 'Finished creating gif.', subject = 'cold_pools_thetaE.py', attachments = ['../theta_e/theta_e_anim.gif'], isAttach = True)

