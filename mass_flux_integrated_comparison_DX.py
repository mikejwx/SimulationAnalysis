import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from analysis_tools import fromComponents, bilinear_interpolation, get_cs_coords, lcl as get_lcl, downwind_rectangle
from scipy import integrate, interpolate
from datetime import datetime as dt
from STASH_keys import rho_key, u_key, v_key, w_key, pthe_key, q_key, theta_key
from SkewT_archer import PTtoTemp

"""
Compute the updraught mass flux over a range of distances downwind of the island.
Compare between the grid spacing flux experiments.
"""

l_testing = False

# Define a function to do the integration
def get_Mc(y_prime, M, dx = 500.):
    """
    Takes the rotated coordinate system (x_prime, y_prime) and the mass flux (M)
    and averages it in the x_prime direction, integrates in the y_prime direction.
    
    To do this, it must first do some housekeeping on the arrays to make sure
    that we are integrating along the right cartesian axis, rather than the 
    array axis.
    
    y_prime -> the rotated y coordinate, increasing to the right of the line
               directly downwind of the island (i.e. positive to the right, and
               negative to the left).
    M       -> the time-mean mass flux array at a given height
    
    All above arrays are 2-D.
    """
    
    # First squeeze in the x_prime direction
    M_x         = np.array([np.nanmean(np.where((lower <= y_prime)*(y_prime <= lower + dx), M, np.nan)) for lower in np.arange(-5050.0, 5000.0, dx)])
    new_y_prime = np.arange(-5000.0, 5000.1, dx)
    
    # Then integrate it in y_prime
    M_int = integrate.trapz(x = new_y_prime, y = M_x)
    
    return M_int

# read in the data
print '[' + dt.now().strftime('%H:%M:%S') + '] Defining the paths to data...'
# Define the paths in which the data are stored

paths = {'DX0100'  : '/nerc/n02/n02/xb899100/CloudTrail/Control2/',
         'DX0200'  : '/nerc/n02/n02/xb899100/CloudTrail/Control_0200m_HRIC_INV/',
         'DX0400'  : '/nerc/n02/n02/xb899100/CloudTrail/Control_0400m_HRIC_INV/',
         'DX0800'  : '/nerc/n02/n02/xb899100/CloudTrail/Control_0800m_HRIC_INV/',
         'DX1600'  : '/nerc/n02/n02/xb899100/CloudTrail/Control_1600m_HRIC_INV/'}

print '[' + dt.now().strftime('%H:%M:%S') + '] Starting the routine for each experiment...'
mass_flux_xs = {}
lcl = {}
# For each experiment
for key in paths.keys():
    # Initialise a dictionary entry for this experiment
    mass_flux_xs[key] = {}
    print '[' + dt.now().strftime('%H:%M:%S') + '] Starting experiment ' + key + '...'
    hours = ['06', '09'] # Which time chunks of output are we analysing?
    for hour in hours:
        print '[' + dt.now().strftime('%H:%M:%S') + '] Opening the ' + hour + ' netCDF...'
        # Read in swath netCDF
        with Dataset(paths[key] + 'fluxes_' + hour + '.nc', 'r') as rho_nc:
            if hour == hours[0]:
                # Find the time key
                rho_time_key = [tkey for tkey in rho_nc.variables.keys() if 'min' in tkey][0]
                # Grab the times
                mass_flux_xs[key]['rho_times'] = rho_nc.variables[rho_time_key][:]*1.
                # Grab the rho heights
                z_rho = rho_nc.variables['rholev_zsea_rho'][:]*1.
                # Grab the data
                mass_flux_xs[key][rho_key] = rho_nc.variables[rho_key][:]*1.
            else:
                mass_flux_xs[key]['rho_times'] = np.concatenate((mass_flux_xs[key]['rho_times'], rho_nc.variables[rho_time_key][:]), axis = 0)
                mass_flux_xs[key][rho_key] = np.concatenate((mass_flux_xs[key][rho_key], rho_nc.variables[rho_key][:]*1.), axis = 0)
        with Dataset(paths[key] + 'u_' + hour + '.nc', 'r') as u_nc:
            if hour == hours[0]:
                # Grab the u-wind component
                mass_flux_xs[key][u_key] = u_nc.variables[u_key][:]*1.
            else:
                mass_flux_xs[key][u_key] = np.concatenate((mass_flux_xs[key][u_key], u_nc.variables[u_key][:]), axis = 0)
        
        with Dataset(paths[key] + 'v_' + hour + '.nc', 'r') as v_nc:
            if hour == hours[0]:
                # Grab the v-wind component
                mass_flux_xs[key][v_key] = v_nc.variables[v_key][:]*1.
            else:
                mass_flux_xs[key][v_key] = np.concatenate((mass_flux_xs[key][v_key], v_nc.variables[v_key][:]), axis = 0)
        
        with Dataset(paths[key] + 'bouy_' + hour + '.nc', 'r') as thermo_nc:
            if hour == hours[0]:
                # Grab the height dimension
                z_theta = thermo_nc.variables['thlev_zsea_theta'][:]*1.
                # Grab the w-wind component
                mass_flux_xs[key][w_key] = thermo_nc.variables[w_key][:]*1.
                # Grab the potential temperature
                mass_flux_xs[key][theta_key] = thermo_nc.variables[theta_key][:]*1.
                # Grab the specific humidity
                mass_flux_xs[key][q_key] = thermo_nc.variables[q_key][:]*1.
                # Grab the pressure on theta levels if it exists
                if pthe_key in thermo_nc.variables.keys():
                    mass_flux_xs[key][pthe_key] = thermo_nc.variables[pthe_key][:]*1.
            else:
                mass_flux_xs[key][w_key] = np.concatenate((mass_flux_xs[key][w_key], thermo_nc.variables[w_key][:]*1.), axis = 0)
                mass_flux_xs[key][theta_key] = np.concatenate((mass_flux_xs[key][theta_key], thermo_nc.variables[theta_key][:]*1.), axis = 0)
                mass_flux_xs[key][q_key] = np.concatenate((mass_flux_xs[key][q_key], thermo_nc.variables[q_key][:]*1.), axis = 0)
                if pthe_key in thermo_nc.variables.keys():
                    mass_flux_xs[key][pthe_key] = np.concatenate((mass_flux_xs[key][pthe_key], thermo_nc.variables[pthe_key][:]*1.), axis = 0)
        
    start_t = 540. # What time are we starting the analysis?
    end_t   = 720. # What time are we finishing the analysis?
    start_idx = np.where(np.abs(mass_flux_xs[key]['rho_times'] - start_t) == np.min(np.abs(mass_flux_xs[key]['rho_times'] - start_t)))[0][0]
    end_idx   = np.where(np.abs(mass_flux_xs[key]['rho_times'] - end_t)   == np.min(np.abs(mass_flux_xs[key]['rho_times'] - end_t)))[0][0]
    t_idx = range(start_idx, end_idx + 1)
    
    print '[' + dt.now().strftime('%H:%M:%S') + '] Reading the required dimensions...'
    # Want the cross section across the trail averaged in time between start_t 
    # and end_t
    
    ### Set-up the coordinate system ###
    # We know from our domain definition where the centre of the island should be
    R_i = 1000.0*(50.0/np.pi)**0.5 # island radius, metres
    x_c = [108000.0 if key != 'DX0100' else 100000.0 + R_i][0]
    y_c = [16000.0 if key != 'DX0100' else 4*R_i][0]
    
    # use our new function to define a rectangle downwind
    # Define the existing cartesian coordinate system
    grid_spacing = float(key[2:])
    X, Y = np.meshgrid(np.arange(mass_flux_xs[key][theta_key][0,0,0,:].shape[0])*grid_spacing, np.arange(mass_flux_xs[key][theta_key][0,0,:,0].shape[0])*grid_spacing)
    
    # Compute and store the horizontal mean LCL if the required variables are present
    lcl[key] = [np.nanmean(np.array([get_lcl(temp = PTtoTemp(theta = mass_flux_xs[key][theta_key][idx,0,:,:], PIN = mass_flux_xs[key][pthe_key][idx,0,:,:], t_units = 'K', p_units = 'Pa'), q = mass_flux_xs[key][q_key][idx,0,:,:], pres = mass_flux_xs[key][pthe_key][idx,:,:,:], z = z_theta) for idx in t_idx])) if pthe_key in mass_flux_xs[key].keys() else 650.0][0]
    
    # Get the mean wind speed and direction below LCL
    iz_LCL = np.where(np.abs(z_theta - lcl[key]) == np.min(np.abs(z_theta - lcl[key])))[0][0]
    speed, wind_dir = fromComponents(u = np.nanmean(mass_flux_xs[key][u_key][t_idx,:iz_LCL,:,:]), v = np.nanmean(mass_flux_xs[key][v_key][t_idx,:iz_LCL,:,:]))
    # Get the distance away from the island-edge of the surface warm plume.
    R_wp = np.nanmax(np.where(np.nanmean(mass_flux_xs[key][theta_key][start_idx:end_idx,0,:,:], axis = 0) - np.nanmean(mass_flux_xs[key][theta_key][start_idx:end_idx,0,:,:]) > 0.1, np.sqrt((X-x_c)**2 + (Y-y_c)**2) - R_i, np.nan))
    # Compute the rectangle and the new cartesian coordinate system
    mask, mass_flux_xs[key]['y_prime'], mass_flux_xs[key]['x_prime'] = downwind_rectangle(wind_dir = wind_dir, x_c = x_c, y_c = y_c, X = X, Y = Y, island_radius = R_i, dist_0 = 0 + R_wp, dist_1 = (end_t - start_t)*60.0*speed + R_wp, half_width = 5000.0)
    print '[' + dt.now().strftime('%H:%M:%S') + '] Computing mass flux...'
    # Convert rho from rho levels to theta levels
    rho_theta_levels = np.array([interpolate.interp1d(x = z_rho, y = mass_flux_xs[key][rho_key][idx,:,:,:], axis = 0, fill_value = 'extrapolate')(z_theta[iz_LCL]) for idx in t_idx])
    # Average in time
    mf = np.nanmean(rho_theta_levels*mass_flux_xs[key][w_key][t_idx,iz_LCL,:,:], axis = 0)
    mass_flux_xs[key]['mass_flux_raw'] = mf
    mass_flux_xs[key]['mass_flux']     = mf*mask

# do an integration...
circulation_MF = {}
for key in paths.keys():
    circulation_MF[key] = get_Mc(mass_flux_xs[key]['y_prime'], mass_flux_xs[key]['mass_flux']*1.0, dx = float(key[2:]))

import matplotlib
my_cmap = matplotlib.cm.get_cmap('Reds')
my_colors = {}
for key in circulation_MF.keys():
    my_colors[key] = my_cmap(float(key[2:])/1600.0)

# do plot
offset = 0.04
keys = paths.keys()
keys.sort()
fig = plt.figure()
axa = fig.add_subplot(1, 1, 1, adjustable = 'box')
# plot a point with a different colour for each dx key
[axa.plot(float(key[2:])/1000., circulation_MF[key], color = my_colors[key], marker = ['*' if key == 'DX0100' else 'o'][0], ls = 'None', markersize = 10, label = [key if key != 'DX0100' else 'Ch05'][0]) for key in keys]
axa.set_ylabel(u'Circulation Mass Flux (kg m$^{-1}$ s$^{1}$)')
axa.set_xlabel(u'Horizontal Grid Length (km)')
# Plot a one to one line
my_y = np.array([circulation_MF[key] for key in circulation_MF.keys()])
my_max = 50*(max(my_y)/50 + 1)
# Domain shape parameters
axa.set_xlim([0, 2.])
axa.set_ylim([0, my_max])
axa.legend(loc = 'lower right', numpoints = 1, frameon = 0)
plt.savefig('../circulationM_Heat_cloudbandonly_DX.png', dpi = 150, bbox_inches = 'tight')
plt.show()
crash
################################################################################
# Introduce the hourly-mean time evolution of Mc for each experiment           #
################################################################################
circulation_MF_hourly = {}
for key in paths.keys():
    print key
    ### Set-up the coordinate system ###
    # We know from our domain definition where the centre of the island should be
    R_i = 1000.0*(50.0/np.pi)**0.5 # island radius, metres
    x_c = [108000.0 if key != 'DX0100' else 100000.0 + R_i][0]
    y_c = [16000.0 if key != 'DX0100' else 4*R_i][0]
    
    # Define the existing cartesian coordinate system
    grid_spacing = float(key[2:])
    X, Y = np.meshgrid(np.arange(mass_flux_xs[key][theta_key][0,0,0,:].shape[0])*grid_spacing, np.arange(mass_flux_xs[key][theta_key][0,0,:,0].shape[0])*grid_spacing)
    
    for hour in range(6, 12):
        start_t = hour*60.       # What time are we starting the analysis?
        end_t   = (hour + 1)*60. # What time are we finishing the analysis?
        start_idx = np.where(np.abs(mass_flux_xs[key]['rho_times'] - start_t) == np.min(np.abs(mass_flux_xs[key]['rho_times'] - start_t)))[0][0]
        end_idx   = np.where(np.abs(mass_flux_xs[key]['rho_times'] - end_t)   == np.min(np.abs(mass_flux_xs[key]['rho_times'] - end_t)))[0][0]
        t_idx = range(start_idx, end_idx + 1)
        # Want the cross section across the trail averaged in time between start_t 
        # and end_t
        # Compute and store the horizontal mean LCL if the required variables are present
        lcl[key] = [np.nanmean(np.array([get_lcl(temp = PTtoTemp(theta = mass_flux_xs[key][theta_key][idx,0,:,:], PIN = mass_flux_xs[key][pthe_key][idx,0,:,:], t_units = 'K', p_units = 'Pa'), q = mass_flux_xs[key][q_key][idx,0,:,:], pres = mass_flux_xs[key][pthe_key][idx,:,:,:], z = z_theta) for idx in t_idx])) if pthe_key in mass_flux_xs[key].keys() else 650.0][0]
        
        # Get the mean wind speed and direction below LCL
        iz_LCL = np.where(np.abs(z_theta - lcl[key]) == np.min(np.abs(z_theta - lcl[key])))[0][0]
        speed, wind_dir = fromComponents(u = np.nanmean(mass_flux_xs[key][u_key][t_idx,:iz_LCL,:,:]), v = np.nanmean(mass_flux_xs[key][v_key][t_idx,:iz_LCL,:,:]))
        # Get the distance away from the island-edge of the surface warm plume.
        R_wp = np.nanmax(np.where(np.nanmean(mass_flux_xs[key][theta_key][start_idx:end_idx,0,:,:], axis = 0) - np.nanmean(mass_flux_xs[key][theta_key][start_idx:end_idx,0,:,:]) > 0.1, np.sqrt((X-x_c)**2 + (Y-y_c)**2) - R_i, np.nan))
        # Compute the rectangle and the new cartesian coordinate system
        mask, mass_flux_xs[key]['y_prime'], mass_flux_xs[key]['x_prime'] = downwind_rectangle(wind_dir = wind_dir, x_c = x_c, y_c = y_c, X = X, Y = Y, island_radius = R_i, dist_0 = 0 + R_wp, dist_1 = (end_t - start_t)*60.0*speed + R_wp, half_width = 5000.0)
        
        print '[' + dt.now().strftime('%H:%M:%S') + '] Computing mass flux...'
        # Convert rho from rho levels to theta levels
        rho_theta_levels = np.array([interpolate.interp1d(x = z_rho, y = mass_flux_xs[key][rho_key][idx,:,:,:], axis = 0, fill_value = 'extrapolate')(z_theta[iz_LCL]) for idx in t_idx])
        # Average in time
        mf = np.nanmean(rho_theta_levels*mass_flux_xs[key][w_key][t_idx,iz_LCL,:,:], axis = 0)
        if hour == 6:
            mass_flux_xs[key]['mass_flux_raw_hourly'] = np.array([mf])
            mass_flux_xs[key]['mass_flux_hourly']     = np.array([mf*mask])
        else:
            mass_flux_xs[key]['mass_flux_raw_hourly'] = np.concatenate((mass_flux_xs[key]['mass_flux_raw_hourly'], np.array([mf])), axis = 0)
            mass_flux_xs[key]['mass_flux_hourly']     = np.concatenate((mass_flux_xs[key]['mass_flux_hourly'], np.array([mf*mask])), axis = 0)
        # do the integration to get Mc
        if hour == hours[0]:
            circulation_MF_hourly[key] = np.array([get_Mc(mass_flux_xs[key]['y_prime'], mass_flux_xs[key]['mass_flux_hourly'][hours.index(hour),:,:]*1.0, dx = float(key[2:]))])
        else:
            circulation_MF_hourly[key] = np.concatenate((circulation_MF_hourly[key], np.array([get_Mc(mass_flux_xs[key]['y_prime'], mass_flux_xs[key]['mass_flux_hourly'][hours.index(hour),:,:]*1.0, dx = float(key[2:]))])), axis = 0)

my_keys = ['DX0100', 'DX0200', 'DX0400', 'DX0800', 'DX1600']
my_colors = {'DX0100' : 'yellow', 'DX0200' : 'gold', 'DX0400' : 'orange', 'DX0800' : 'red', 'DX1600' : 'darkred'}
# Make a plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
[ax.plot(hours, circulation_MF_hourly[key], color = my_colors[key], lw = 2, label = key) for key in my_keys]
plt.legend(loc = 0)
plt.savefig('../circulationM_Heat_cloudbandonly_DX_timeseries.png', dpi = 150, bbox_inches = 'tight')
plt.show()


################################################################################
### Testing Section ###
################################################################################
if l_testing:
    key = 'DX0400'
    X, Y = np.meshgrid(np.arange(mass_flux_xs[key]['mass_flux'].shape[1])*400.0, np.arange(mass_flux_xs[key]['mass_flux'].shape[0])*400.0)
    # Run the following to demonstrate the steps of the calculation
    # Step one calculate mass flux
    R = np.sqrt((X - x_c)**2 + (Y - y_c)**2)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, adjustable = 'box', aspect = 1)
    my_plt = ax.contourf(X/1000, Y/1000, mass_flux_xs[key]['mass_flux_raw'], cmap = 'bwr', levels = [level for level in np.arange(-1.0, 1.1, 0.1) if level != 0], extend = 'both')
    plt.colorbar(my_plt, ax = ax, orientation = 'horizontal', label = 'Mass Flux')
    ax.contour(X/1000, Y/1000, R, levels = [R_i], colors = ['k'])
    plt.show()
    
    # Step two crop
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, adjustable = 'box', aspect = 1)
    my_plt = ax.contourf(X/1000, Y/1000, mass_flux_xs[key]['mass_flux'], cmap = 'bwr', levels = [level for level in np.arange(-1.0, 1.1, 0.1) if level != 0], extend = 'both')
    plt.colorbar(my_plt, ax = ax, orientation = 'horizontal', label = 'Mass Flux')
    ax.contour(X/1000, Y/1000, R, levels = [R_i], colors = ['k'])
    ax.contour(X/1000, Y/1000, np.where(mass_flux_xs[key]['mass_flux'] == mass_flux_xs[key]['mass_flux'], 1.0, 0.0), levels = [0.5], colors = ['purple'])
    plt.show()

    # step three average in x_prime
    plt.plot(np.arange(-5000.0, 5000.1, 100.0)/1000, np.array([np.nanmean(np.where((lower <= mass_flux_xs[key]['y_prime'])*(mass_flux_xs[key]['y_prime'] <= lower + 100.0), mass_flux_xs[key]['mass_flux'], np.nan)) for lower in np.arange(-5050.0, 5000.0, 100.0)]), 'k')
    plt.xlabel('y_prime (km)')
    plt.ylabel('Mass Flux')
    plt.show()

