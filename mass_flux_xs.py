import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from analysis_tools import fromComponents, bilinear_interpolation, get_cs_coords
from scipy import integrate, interpolate
from datetime import datetime as dt
#plt.switch_backend('agg')

"""
Compute the mass flux across the cloud trail for our experiments, how does it
change for the changing environment at different heights and different 
distances downwind of the island.
"""

print '[' + dt.now().strftime('%H:%M:%S') + '] Defining the paths to data...'
# Define the paths in which the data are stored

paths = {'Control'       : '/nerc/n02/n02/xb899100/CloudTrail/Control/',
         'B1/3'          : '/nerc/n02/n02/xb899100/CloudTrail/H125/',
         'B3.0'          : '/nerc/n02/n02/xb899100/CloudTrail/H375/',
         'BL_RHm25'      : '/nerc/n02/n02/xb899100/CloudTrail/RH_BLm25/',
         'FA_RHm25'      : '/nerc/n02/n02/xb899100/CloudTrail/RH_FAm25/',
         'U05'           : '/nerc/n02/n02/xb899100/CloudTrail/U05/'}

# 'Control_short' : '/nerc/n02/n02/xb899100/CloudTrail/Control_short/',

my_colors = {'Control'       : 'black',
             'B1/3'          : 'cyan',
             'B3.0'          : 'orange',
             'Control_short' : 'grey',
             'BL_RHm25'      : 'olive',
             'FA_RHm25'      : 'magenta',
             'U05'           : 'green'}

print '[' + dt.now().strftime('%H:%M:%S') + '] Determining short sims...'
short_sims = ['Control_short', 'BL_RHm25', 'FA_RHm25', 'U05']

print '[' + dt.now().strftime('%H:%M:%S') + '] Defining the required variable keys...'
# Define the keys for rho and winds
from STASH_keys import rho_key, u_key, v_key, w_key

print '[' + dt.now().strftime('%H:%M:%S') + '] Starting the routine for each experiment...'
mass_flux_xs = {}
for key in paths.keys():
    print '[' + dt.now().strftime('%H:%M:%S') + '] Starting experiment ' + key + '...'
    if key in short_sims:
        print '[' + dt.now().strftime('%H:%M:%S') + '] ' + key + ' is a short experiment...'
        hour = '04'
        start_t = 300.
        end_t = 480.
    else:
        print '[' + dt.now().strftime('%H:%M:%S') + '] ' + key + ' is a long experiment...'
        hour = '09'
        start_t = 540.
        end_t = 720.
    
    print '[' + dt.now().strftime('%H:%M:%S') + '] Opening the netCDF...'
    # Read in swath netCDF
    rho_swath_nc  = Dataset(paths[key] + 'rho_swath_' + hour + '.nc', 'r')
    w_swath_nc    = Dataset(paths[key] + 'w_swath_' + hour + '.nc', 'r')
    
    # Find the start and end time idxs
    rho_time_key = [tkey for tkey in rho_swath_nc.variables.keys() if 'min' in tkey][0]
    
    rho_times = rho_swath_nc.variables[rho_time_key][:]
    
    start_idx = np.where(np.abs(rho_times - start_t) == np.min(np.abs(rho_times - start_t)))[0][0]
    end_idx   = np.where(np.abs(rho_times - end_t)   == np.min(np.abs(rho_times - end_t)))[0][0]
    
    t_idx = range(start_idx, end_idx + 1)
    print '[' + dt.now().strftime('%H:%M:%S') + '] Reading the required dimensions...'
    # Want the cross section across the trail averaged in time for 3-hours prior
    # to the peak island heating
    
    # Define our height dimension
    z_theta = w_swath_nc.variables['thlev_zsea_theta'][:]*1.
    
    # We know from our domain definition where the centre of the island should be
    R_i = 1000.0*(50.0/np.pi)**0.5 # island radius
    x_c = 100000.0 + R_i
    y_c = 4*R_i
    
    # Define the distance away from the island centre
    x_prime = rho_swath_nc.variables['x_prime'][:]
    y_prime = rho_swath_nc.variables['y_prime'][:]
    x_R = x_prime - x_c
    y_R = y_prime - y_c
    R = -np.sign(x_R)*np.sqrt(x_R**2 + y_R**2)
    
    # Define the coordinate system in terms of along the wind and across the wind
    i_along, i_across = np.where(R == np.min(np.abs(R)))
    R_along = R[i_along[0],:]
    R_across = R[:,i_across[0]]
    
    # Select a chunk from the 'along the wind' direction
    r_target0 = R_i + 00000.
    r_target1 = R_i + 72000.
    
    idx0 = np.where(np.abs(R_along - r_target0) == np.min(np.abs(R_along - r_target0)))[0][0]
    idx1 = np.where(np.abs(R_along - r_target1) == np.min(np.abs(R_along - r_target1)))[0][0] + 1
    
    print '[' + dt.now().strftime('%H:%M:%S') + '] Computing mass flux...'
    # Average in time
    mf = np.nanmean(rho_swath_nc.variables[rho_key][t_idx,:,:,:]*w_swath_nc.variables[w_key][t_idx,:,:,:], axis = 0)
    # Average along the CT in the x' direction
    mass_flux_xs[key] = np.nanmean(mf[:,:,idx0:idx1], axis = 2)
    
    # Grab the landmask before closing
    lsm = rho_swath_nc.variables['lsm'][:]
    
    rho_swath_nc.close()
    w_swath_nc.close()

### Plot the mass flux filled contour, and the mass flux at chosen heights ###
fig = plt.figure(figsize = (12, 18))
counter = 1
# Choose a height to plot
z_targetBL = 350.
print '[' + dt.now().strftime('%H:%M:%S') + '] Plotting the mass flux for ' + key + '...'
iz_targetBL = np.where(np.abs(z_theta - z_targetBL) == np.min(np.abs(z_theta - z_targetBL)))[0][0]

# Choose a height to plot
z_targetCL = 1500.
iz_targetCL = np.where(np.abs(z_theta - z_targetCL) == np.min(np.abs(z_theta - z_targetCL)))[0][0]

# Choose a height for plot top
z_targetPT = 4000.
iz_targetPT = np.where(np.abs(z_theta - z_targetPT) == np.min(np.abs(z_theta - z_targetPT)))[0][0] + 1

for key in paths.keys():
    # Plot the mass flux at that height
    ax0 = fig.add_subplot(len(paths.keys()), 2, counter)
    ax0.plot(R_across/1000., mass_flux_xs[key][iz_targetBL,:], 'k', label = str(int(z_theta[iz_targetBL])) + ' m', lw = 2)
    ax0.plot(R_across/1000., mass_flux_xs[key][iz_targetCL,:], 'k--', label = str(int(z_theta[iz_targetCL])) + ' m', lw = 2)
    ax0.set_xlim([-5, 5])
    ax0.set_ylim([-0.3, 0.6])
    ax0.plot([-5,5],[0,0], color = 'grey', linestyle = ':')
    ax0.set_ylabel('Mass Flux ($kg m^{-2} s^{-1}$)')
    ax0.set_xlabel('Distance from along wind centrepoint (i.e. y$^{\prime}$, km)')
    ax0.legend(loc = 'upper left', frameon = False)
    ax0.set_title(key)
    
    counter += 1
    ax1 = fig.add_subplot(len(paths.keys()), 2, counter, adjustable = 'box', aspect = 1)
    xs = ax1.contourf(R_across/1000., z_theta[:iz_targetPT]/1000., mass_flux_xs[key][:iz_targetPT,:], cmap = 'bwr', levels = [-0.3, -0.2, -0.1, 0.1, 0.2, 0.3], extend = 'both')
    ax1.plot(R_across/1000., np.zeros_like(R_across) + z_theta[iz_targetBL]/1000., 'k')
    ax1.plot(R_across/1000., np.zeros_like(R_across) + z_theta[iz_targetCL]/1000., 'k--')
    fig.colorbar(xs, ax = ax1)
    ax1.set_xlabel('Distance from along wind centrepoint (i.e. y$^{\prime}$, km)')
    ax1.set_ylabel('Height (km)')
    ax1.set_title(key)
    
    counter += 1

fig.suptitle('Mass flux (' + str(int((r_target0 - R_i)/1000.)) + ' - ' + str(int((r_target1 - R_i)/1000.)) + ' km) downwind of lee-side, (9am - 12pm), \nfor each experiment', fontsize = 20)
fig.tight_layout(rect=[0.0, 0.05, 1.0, 0.95])
plt.savefig('../mass_flux_' + str(int((r_target0 - R_i)/1000.)) + '-' + str(int((r_target1 - R_i)/1000.)) + 'km_downwind_all_exp.png', dpi = 150)
plt.show()



# Plot the LWP at midday for each case
from STASH_keys import lwp_key

X, Y = np.meshgrid(np.arange(0., 116000., 100.), np.arange(0., 31900., 100.))
x_prime_p = x_prime%(X.max() + 100.)
y_prime_p = y_prime%(Y.max() + 100.)
swath_mask = np.zeros_like(x_prime)
swath_mask[:,idx0:idx1] += 1.
fig = plt.figure()
from math import ceil
ncol = 2
for key in paths.keys():
    lwp_nc = Dataset(paths[key] + 'lwp_00.nc', 'r')
    lwp_time_key = [tkey for tkey in lwp_nc.variables.keys() if 'min' in tkey][0]
    
    # Check if we're dealing with a short simulation
    if key in short_sims:    
        l_short = True
    else:
        l_short = False
    
    # If it is a short simulation, make sure to look for the correct times
    if l_short:
        lwp_times = lwp_nc.variables[lwp_time_key][:] + 240.
    else:
        lwp_times = lwp_nc.variables[lwp_time_key][:]
    
    preday_idx = np.where(lwp_times == 360.)[0][0]
    midday_idx = np.where(lwp_times == 1080.)[0][0] + 1
    
    lwp_data = np.nanmean(lwp_nc.variables[lwp_key][preday_idx:midday_idx,:,:], axis = 0)*1000.
    
    n_plot = paths.keys().index(key) + 1
    ax = fig.add_subplot(ceil(len(paths.keys())/float(ncol)), ncol, n_plot, adjustable = 'box', aspect = 1)
    lwp_plt = ax.contourf(X/1000., Y/1000., lwp_data, cmap = 'Greys_r', levels = [0, 25, 50, 75, 100, 125, 150, 175, 200], extend = 'max')
    fig.colorbar(lwp_plt, ax = ax, label = 'LWP (g/m$^{2}$)')
    # island mask
    ax.contour(x_prime/1000., y_prime/1000., lsm, colors = ['r'], levels = [1e-16])
    # main bit
    ax.contourf(x_prime/1000., y_prime/1000., swath_mask, colors = ['r'], levels = [1e-16, 1e2], alpha = 0.25)
    # overlap left
    ax.contourf((x_prime + (X.max() + 100.))/1000., y_prime/1000., swath_mask, colors = ['r'], levels = [1e-16, 1e2], alpha = 0.25)
    # overlap right
    ax.contourf((x_prime - (X.max() + 100.))/1000., y_prime/1000., swath_mask, colors = ['r'], levels = [1e-16, 1e2], alpha = 0.25)
    # overlap top
    ax.contourf(x_prime/1000., (y_prime + (Y.max() + 100.))/1000., swath_mask, colors = ['r'], levels = [1e-16, 1e2], alpha = 0.25)
    # overlap bottom
    ax.contourf(x_prime/1000., (y_prime - (Y.max() + 100.))/1000., swath_mask, colors = ['r'], levels = [1e-16, 1e2], alpha = 0.25)
    ax.set_xlim([0, X.max()/1000.])
    ax.set_ylim([0, Y.max()/1000.])
    ax.set_title('Experiment: ' + key)

plt.savefig('../corresponding_lwp_' + str(int((r_target0 - R_i)/1000.)) + '-' + str(int((r_target1 - R_i)/1000.)) + 'km.png', dpi = 150.)
plt.show()



