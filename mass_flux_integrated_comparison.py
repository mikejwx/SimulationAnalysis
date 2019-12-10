import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from analysis_tools import fromComponents, bilinear_interpolation, get_cs_coords, lcl as get_lcl, downwind_rectangle
from scipy import integrate, interpolate
from datetime import datetime as dt
#plt.switch_backend('agg')
from STASH_keys import rho_key, u_key, v_key, w_key, pthe_key, q_key, theta_key
from SkewT_archer import PTtoTemp

"""
Compute the updraught mass flux over a range of distances downwind of the island.
Compare between the heat flux experiments.
"""
l_testing = True

# Define a function to do the integration
def get_Mc(y_prime, M):
    """
    Takes the rotated coordinate system (x_prime, y_prime) and the mass flux (M)
    and averages it in the x_prime direction, integrates in the y_prime direction.
    
    To do this, it must first do some housekeeping on the arrays to make sure
    that we are integrating along the right cartesian axis, rather than the 
    array axis.
    
    y_prime -> the rotated y coordinate, increasing to the right of the line
               directly downwind of the island (i.e. positive to thr right, and
               negative to the left).
    M       -> the time-mean mass flux array at a given height
    
    All above arrays are 2-D.
    """
    
    # First squeeze in the x_prime direction
    M_x         = np.array([np.nanmean(np.where((lower <= y_prime)*(y_prime <= lower + 100.0), M, np.nan)) for lower in np.arange(-5050.0, 5000.0, 100.0)])
    new_y_prime = np.arange(-5000.0, 5000.1, 100.0)
    
    # Then integrate it in y_prime
    M_int = integrate.trapz(x = new_y_prime, y = M_x)
    
    return M_int

# read in the data
print '[' + dt.now().strftime('%H:%M:%S') + '] Defining the paths to data...'
# Define the paths in which the data are stored

paths = {'H250E250' : '/nerc/n02/n02/xb899100/CloudTrail/Control/',
         'H125E375' : '/nerc/n02/n02/xb899100/CloudTrail/H125/',
         'H375E125' : '/nerc/n02/n02/xb899100/CloudTrail/H375/',
         'H400E250' : '/nerc/n02/n02/xb899100/CloudTrail/H400E250/',
         'H350E250' : '/nerc/n02/n02/xb899100/CloudTrail/H350E250/',
         'H300E250' : '/nerc/n02/n02/xb899100/CloudTrail/H300E250/',
         'H200E250' : '/nerc/n02/n02/xb899100/CloudTrail/H200E250/',
         'H150E250' : '/nerc/n02/n02/xb899100/CloudTrail/H150E250/',
         'H100E250' : '/nerc/n02/n02/xb899100/CloudTrail/H100E250/',
         'H050E250' : '/nerc/n02/n02/xb899100/CloudTrail/H050E250/'}

print '[' + dt.now().strftime('%H:%M:%S') + '] Starting the routine for each experiment...'
mass_flux_xs = {}
lcl = {}
# For each experiment
for key in paths.keys():
    print '[' + dt.now().strftime('%H:%M:%S') + '] Starting experiment ' + key + '...'
    hour = '09'    # Which time chunk of output are we analysing?
    start_t = 540. # What time are we starting the analysis?
    end_t = 720.   # What time are we finishing the analysis?
    
    print '[' + dt.now().strftime('%H:%M:%S') + '] Opening the netCDF...'
    # Read in swath netCDF
    rho_nc    = Dataset(paths[key] + 'fluxes_' + hour + '.nc', 'r')
    u_nc      = Dataset(paths[key] + 'u_' + hour + '.nc', 'r')
    v_nc      = Dataset(paths[key] + 'v_' + hour + '.nc', 'r')
    thermo_nc = Dataset(paths[key] + 'bouy_' + hour + '.nc', 'r')
    
    # Find the start and end time idxs
    rho_time_key = [tkey for tkey in rho_nc.variables.keys() if 'min' in tkey][0]
    
    rho_times = rho_nc.variables[rho_time_key][:]
    start_idx = np.where(np.abs(rho_times - start_t) == np.min(np.abs(rho_times - start_t)))[0][0]
    end_idx   = np.where(np.abs(rho_times - end_t)   == np.min(np.abs(rho_times - end_t)))[0][0]
    
    t_idx = range(start_idx, end_idx + 1)
    
    print '[' + dt.now().strftime('%H:%M:%S') + '] Reading the required dimensions...'
    # Want the cross section across the trail averaged in time between start_t 
    # and end_t
    
    ### Set-up the coordinate system ###
    # Define our height dimension
    z_theta = thermo_nc.variables['thlev_zsea_theta'][:]*1.
    
    # We know from our domain definition where the centre of the island should be
    R_i = 1000.0*(50.0/np.pi)**0.5 # island radius, metres
    x_c = 100000.0 + [R_i if key in ['H250E250', 'H125E375', 'H375E125'] else 2*R_i][0]
    y_c = 4*R_i
    
    mass_flux_xs[key] = {}
    
    # use our new function to define a rectangle downwind
    # Define the existing cartesian coordinate system
    X, Y = np.meshgrid(np.arange(thermo_nc.variables[theta_key][0,0,0,:].shape[0])*100.0, np.arange(thermo_nc.variables[theta_key][0,0,:,0].shape[0])*100.0)
    
    # Compute and store the horizontal mean LCL
    lcl[key] = [np.nanmean(np.array([get_lcl(temp = PTtoTemp(theta = thermo_nc.variables[theta_key][idx,0,:,:], PIN = thermo_nc.variables[pthe_key][idx,0,:,:], t_units = 'K', p_units = 'Pa'), q = thermo_nc.variables[q_key][idx,0,:,:], pres = thermo_nc.variables[pthe_key][idx,:,:,:], z = z_theta) for idx in t_idx])) if pthe_key in thermo_nc.variables.keys() else 650.0][0]
    
    # Get the mean wind speed and direction below LCL
    iz_LCL = np.where(np.abs(z_theta - lcl[key]) == np.min(np.abs(z_theta - lcl[key])))[0][0]
    speed, wind_dir = fromComponents(u = np.nanmean(u_nc.variables[u_key][t_idx,:iz_LCL,:,:]), v = np.nanmean(v_nc.variables[v_key][t_idx,:iz_LCL,:,:]))
    
    # Get the distance away from the island-edge of the surface warm plume.
<<<<<<< HEAD
    R_wp = np.nanmax(np.where(np.nanmean(thermo_nc.variables[theta_key][start_idx:end_idx,0,:,:], axis = 0) - np.nanmean(thermo_nc.variables[theta_key][start_idx:end_idx,0,:,:]) > 0.1, np.sqrt((X-x_c)**2 + (Y-y_c)**2) - R_i, np.nan))
=======
    R_wp = 0#np.nanmax(np.where(np.nanmean(thermo_nc.variables[theta_key][start_idx:end_idx,0,:,:], axis = 0) - np.nanmean(thermo_nc.variables[theta_key][start_idx:end_idx,0,:,:]) > 0.1, np.sqrt((X-x_c)**2 + (Y-y_c)**2) - R_i, np.nan))
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
    # Compute the rectangle and the new cartesian coordinate system
    mask, mass_flux_xs[key]['y_prime'], mass_flux_xs[key]['x_prime'] = downwind_rectangle(wind_dir = wind_dir, x_c = x_c, y_c = y_c, X = X, Y = Y, island_radius = R_i, dist_0 = 0+R_wp, dist_1 = (end_t - start_t)*60.0*speed+R_wp, half_width = 5000.0)
    
    print '[' + dt.now().strftime('%H:%M:%S') + '] Computing mass flux...'
    # Convert rho from rho levels to theta levels
    rho_theta_levels = np.array([interpolate.interp1d(x = rho_nc.variables['rholev_zsea_rho'][:], y = rho_nc.variables[rho_key][idx,:,:,:], axis = 0, fill_value = 'extrapolate')(z_theta[iz_LCL]) for idx in t_idx])
    # Average in time
    mf = np.nanmean(rho_theta_levels*thermo_nc.variables[w_key][t_idx,iz_LCL,:,:], axis = 0)
    mass_flux_xs[key]['mass_flux'] = mf*mask
    
    rho_nc.close()
    u_nc.close()
    v_nc.close()
    thermo_nc.close()

# do an integration...
circulation_MF = {}
for key in paths.keys():
    circulation_MF[key] = get_Mc(mass_flux_xs[key]['y_prime'], mass_flux_xs[key]['mass_flux']*1.0)

import matplotlib
my_cmap = matplotlib.cm.get_cmap('Reds')
my_colors = {}
for key in circulation_MF.keys():
    my_colors[key] = my_cmap(float(key[1:4])/500.0)

# do plot
offset = 4
fig = plt.figure()
<<<<<<< HEAD
axa = fig.add_subplot(1, 1, 1, adjustable = 'box')
=======
axa = fig.add_subplot(1, 1, 1, adjustable = 'box', aspect = 1)
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
[axa.plot(int(key[1:4]), circulation_MF[key], color = my_colors[key], 
     marker = ['o' if int(key[5:8]) == 250 else '^' if int(key[5:8]) == 375 else 'v'][0], ls = 'None',
     markersize = 10, label = [u'E = 250 W m$^{-2}$' if (int(key[5:8]) == int(key[1:4])) else u'E = 375 W m$^{-2}$' if int(key[5:8]) == 375 else u'E = 125 W m$^{-2}$' if int(key[5:8]) == 125 else None][0]) for key in paths.keys()]
axa.set_ylabel(u'Circulation Mass Flux (kg m$^{-1}$ s$^{1}$)')
<<<<<<< HEAD
#axa.set_ylabel(u'Circulation Updraught (m s$^{-1}$)')
=======
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
axa.set_xlabel(u'Max Island Surface Sensible Heat Flux (W m$^{-2}$)')
# Plot a line of best fit
my_x = np.array([int(key[1:4]) for key in paths.keys()])
my_y = np.array([circulation_MF[key] for key in paths.keys()])
my_fit = np.polyfit(my_x, my_y, 1)
axa.plot([0,500], [np.sum([my_fit[i]*(x**(len(my_fit) - i - 1)) for i in range(len(my_fit))]) for x in [0,500]], 'k--')
axa.text(250-offset, circulation_MF['H250E250']+offset, u'y = ' + str(round(my_fit[0], 4)) + 'x ' +['+ ' if my_fit[1] >= 0 else ''][0] + str(round(my_fit[1], 4)), rotation = 180*np.arctan(my_fit[0])/np.pi, color = 'k', ha = 'left', va = 'bottom')
<<<<<<< HEAD
#my_max = (max(100*my_y) + 2)/100.
# Plot a one to one line
my_max = 50*(max(my_y)/50 + 1)
#axa.plot([0,500], [0,500], 'r--')
#axa.text(100-offset, 100+offset, 'one-to-one', rotation = 45, color = 'red', ha = 'left', va = 'bottom')
=======
# Plot a one to one line
my_max = 50*(max(my_y)/50 + 1)
axa.plot([0,500], [0,500], 'r--')
axa.text(100-offset, 100+offset, 'one-to-one', rotation = 45, color = 'red', ha = 'left', va = 'bottom')
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
# Domain shape parameters
axa.set_xlim([0, 500])
axa.set_ylim([0, my_max])
axa.text(250+offset, circulation_MF['H250E250']-offset, 'Control', va = 'top')
axa.legend(loc = 0, numpoints = 1)
<<<<<<< HEAD
plt.savefig('../circulationM_Heat_cloudbandonly.png', dpi = 150, bbox_inches = 'tight')
=======
plt.savefig('../circulationM_Heat.png', dpi = 150, bbox_inches = 'tight')
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
plt.show()

# Against peak total energy flux
fig = plt.figure()
<<<<<<< HEAD
axa = fig.add_subplot(1, 1, 1, adjustable = 'box')
[axa.plot(int(key[1:4]) + int(key[5:8]), circulation_MF[key], color = my_colors[key], 
     marker = ['o' if int(key[5:8]) == 250 else 'v' if int(key[5:8]) == 375 else '^'][0], ls = 'None',
     markersize = 10, label = [u'$\\beta = 1$' if (int(key[5:8]) == int(key[1:4])) else u'$\\beta = 1/3$' if int(key[5:8]) == 375 else u'$\\beta = 3$' if int(key[5:8]) == 125 else None][0]) for key in paths.keys()]
axa.set_ylabel(u'Circulation Mass Flux (kg m$^{-1}$ s$^{-1}$)')
=======
axa = fig.add_subplot(1, 1, 1, adjustable = 'box', aspect = 1)
[axa.plot(int(key[1:4]) + int(key[5:8]), circulation_MF[key], color = my_colors[key], 
     marker = ['o' if int(key[5:8]) == 250 else 'v' if int(key[5:8]) == 375 else '^'][0], ls = 'None',
     markersize = 10, label = [u'$\\beta = 1$' if (int(key[5:8]) == int(key[1:4])) else u'$\\beta = 1/3$' if int(key[5:8]) == 375 else u'$\\beta = 3$' if int(key[5:8]) == 125 else None][0]) for key in paths.keys()]
axa.set_ylabel(u'Circulation Mass Flux (kg m$^{-1}$ s$^{1}$)')
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
axa.set_xlabel(u'Total Energy Flux (W m$^{-2}$)')
# Plot a line of best fit
my_x = np.array([int(key[1:4]) + int(key[5:8]) for key in paths.keys()])
my_y = np.array([circulation_MF[key] for key in paths.keys()])
my_fit = np.polyfit(my_x, my_y, 1)
axa.plot([0,700], [np.sum([my_fit[i]*(x**(len(my_fit) - i - 1)) for i in range(len(my_fit))]) for x in [0,700]], 'k--')
axa.text(500-offset, circulation_MF['H250E250']+offset, u'y = ' + str(round(my_fit[0], 4)) + 'x ' +['+ ' if my_fit[1] >= 0 else ''][0] + str(round(my_fit[1], 4)), rotation = 180*np.arctan(my_fit[0])/np.pi, color = 'k', ha = 'left', va = 'bottom')
# Domain shape parameters
axa.set_xlim([250, 700])
axa.set_ylim([0, my_max])
<<<<<<< HEAD
axa.text(500+offset, circulation_MF['H250E250']-offset, 'Control', va = 'top', ha = 'right')
axa.legend(loc = 0, numpoints = 1)
plt.savefig('../circulationM_Heat_Rnet_cloudbandonly.png', dpi = 150, bbox_inches = 'tight')
=======
axa.text(500+offset, circulation_MF['H250E250']-offset, 'Control', va = 'top')
axa.legend(loc = 0, numpoints = 1)
plt.savefig('../circulationM_Heat_Rnet.png', dpi = 150, bbox_inches = 'tight')
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
plt.show()

if l_testing:
    # Run the following to demonstrate the steps of the calculation
    # Step one calculate mass flux
    R = np.sqrt((X - x_c)**2 + (Y - y_c)**2)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, adjustable = 'box', aspect = 1)
    my_plt = ax.contourf(X/1000, Y/1000, mass_flux_xs[key]['mass_flux'], cmap = 'bwr', levels = [level for level in np.arange(-1.0, 1.1, 0.1) if level != 0], extend = 'both')
    plt.colorbar(my_plt, ax = ax, orientation = 'horizontal', label = 'Mass Flux')
    ax.contour(X/1000, Y/1000, R, levels = [R_i], colors = ['k'])
    plt.show()
    
    # Step two crop
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, adjustable = 'box', aspect = 1)
    my_plt = ax.contourf(X/1000, Y/1000, mass_flux_xs[key]['mass_flux'], cmap = 'bwr', levels = [level for level in np.arange(-1.0, 1.1, 0.1) if level != 0], extend = 'both')
    plt.colorbar(my_plt, ax = ax, orientation = 'horizontal', label = 'Mass Flux')
    ax.contour(X/1000, Y/1000, R, levels = [R_i], colors = ['k'])
    ax.contour(X/1000, Y/1000, np.where(mask== mask, 1.0, 0.0), levels = [0.5], colors = ['purple'])
    plt.show()
    
    # step three average in x_prime
    plt.plot(np.arange(-5000.0, 5000.1, 100.0)/1000, np.array([np.nanmean(np.where((lower <= mass_flux_xs[key]['y_prime'])*(mass_flux_xs[key]['y_prime'] <= lower + 100.0), mass_flux_xs[key]['mass_flux'], np.nan)) for lower in np.arange(-5050.0, 5000.0, 100.0)]), 'k')
    plt.xlabel('y_prime (km)')
    plt.ylabel('Mass Flux')
    plt.show()

