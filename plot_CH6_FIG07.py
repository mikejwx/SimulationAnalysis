import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from analysis_tools import fromComponents, bilinear_interpolation, get_cs_coords, lcl as get_lcl, downwind_rectangle
from scipy import integrate, interpolate
from datetime import datetime as dt
from STASH_keys import rho_key, u_key, v_key, w_key, pthe_key, q_key, theta_key, lcl_key, zi_new_key
from SkewT_archer import PTtoTemp
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

"""
Compute the updraught mass flux over a range of distances downwind of the island.
Compare between the heat flux experiments.
"""
l_testing = False

# Define a function to do the integration
def get_Mc(y_prime, M, dx):
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
    # M_up = positive M
    M = np.where(M > 0, M, 0.0)
    # Average along trail
    M_x         = np.array([np.nanmean(np.where((lower <= y_prime)*(y_prime <= lower + dx), M, np.nan)) for lower in np.arange(-3050.0, 3000.0, dx)])
    new_y_prime = np.arange(-3000.0, 3000.1, dx)
    # Integrate across trail
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
zi = {}
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
    zi_nc     = Dataset(paths[key] + 'zi_' + hour + '.nc', 'r')
    
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
    x_c = 100000.0 + [R_i if key in ['DX0100'] else 8000.0][0]
    y_c = [4*R_i if key in ['DX0100'] else 16000.0][0]
    
    mass_flux_xs[key] = {}
    
    # use our new function to define a rectangle downwind
    # Define the existing cartesian coordinate system
    grid_spacing = float(key[2:])
    X, Y = np.meshgrid(np.arange(thermo_nc.variables[theta_key][0,0,0,:].shape[0])*grid_spacing, np.arange(thermo_nc.variables[theta_key][0,0,:,0].shape[0])*grid_spacing)
    
    # Compute and store the horizontal mean LCL
    lcl[key] = zi_nc.variables[lcl_key][t_idx,:,:].mean()
    zi[key]  = zi_nc.variables[zi_new_key][t_idx,:,:].mean()
    # Get the mean wind speed and direction below LCL
    iz_LCL = np.where(np.abs(z_theta - lcl[key]) == np.min(np.abs(z_theta - lcl[key])))[0][0]
    iz_zi  = np.where(np.abs(z_theta - zi[key]) == np.min(np.abs(z_theta - zi[key])))[0][0]
    u_mean = u_nc.variables[u_key][t_idx,:iz_zi,:,:].mean(axis = (0, 2, 3))
    v_mean = v_nc.variables[v_key][t_idx,:iz_zi,:,:].mean(axis = (0, 2, 3))
    U_zi = integrate.trapz(x = z_theta[:iz_zi], y = u_mean)/zi[key]
    V_zi = integrate.trapz(x = z_theta[:iz_zi], y = v_mean)/zi[key]
    speed, wind_dir = fromComponents(u = U_zi, v = V_zi)
    
    # Get the distance away from the island-edge of the surface warm plume.
    mass_flux_xs[key]['theta_anom'] = np.nanmean(np.array([thermo_nc.variables[theta_key][idx,1,:,:] - thermo_nc.variables[theta_key][idx,1,:,:].mean() for idx in t_idx]), axis = 0)
    R_wp = 4*R_i#np.nanmax(np.where(mass_flux_xs[key]['theta_anom'] > 0.1, np.sqrt((X-x_c)**2 + (Y-y_c)**2) - R_i, np.nan))
    
    # Compute the rectangle and the new cartesian coordinate system
    # distance travelled during time averaging period
    mass_flux_xs[key]['mask'], mass_flux_xs[key]['y_prime'], mass_flux_xs[key]['x_prime'] = downwind_rectangle(wind_dir = wind_dir, x_c = x_c, y_c = y_c, X = X, Y = Y, island_radius = R_i, dist_0 = R_wp, dist_1 = R_wp + ((end_t-start_t)*60.0*speed - R_wp), half_width = 3000.0)
    
    print '[' + dt.now().strftime('%H:%M:%S') + '] Computing mass flux...'
    # Convert rho from rho levels to theta levels
    rho_theta_levels = np.array([interpolate.interp1d(x = rho_nc.variables['rholev_zsea_rho'][:], y = rho_nc.variables[rho_key][idx,:,:,:], axis = 0, fill_value = 'extrapolate')(z_theta[iz_LCL]) for idx in t_idx])
    # Average in time
    mf = np.nanmean(rho_theta_levels*thermo_nc.variables[w_key][t_idx,iz_LCL,:,:], axis = 0)
    mass_flux_xs[key]['mass_flux'] = mf
    mass_flux_xs[key]['mass_flux_masked'] = mf*mass_flux_xs[key]['mask']
    
    rho_nc.close()
    u_nc.close()
    v_nc.close()
    thermo_nc.close()
    zi_nc.close()

# do an integration...
circulation_MF = {}
for key in paths.keys():
    circulation_MF[key] = get_Mc(mass_flux_xs[key]['y_prime'], mass_flux_xs[key]['mass_flux_masked']*1.0, float(key[2:]))

import matplotlib
my_cmap = matplotlib.cm.get_cmap('Reds')
my_colors = {}
for key in circulation_MF.keys():
    my_colors[key] = my_cmap((float(key[2:]) + 100.)/2000.0)

DX = {}
for key in paths.keys():
    DX[key] = float(key[2:])/1000.

keys = ['DX0100', 'DX0200', 'DX0400', 'DX0800', 'DX1600']
### Make the plot ###
fig = plt.figure()
# Against peak total energy flux
axa = fig.add_subplot(1, 1, 1)
[axa.plot(DX[key], circulation_MF[key], color = my_colors[key], 
     marker = ['s' if key == 'DX0100' else 'o'][0], ls = 'None',
     markersize = 10, markeredgecolor = 'none', label = key + [' (Ch05)' if key == 'DX0100' else ''][0]) for key in keys]
axa.set_ylabel(u'$M_{c}$ (kg m$^{-1}$ s$^{-1}$)')
axa.set_xlabel(u'Horizontal Grid Spacing (km)')
# Domain shape parameters
axa.set_xlim([0, 1.75])
axa.set_ylim([150, 400])
axa.set_yticks(np.arange(150, 401, 50))
axa.legend(loc = 'lower right', frameon = False, numpoints = 1, fontsize = 12)

plt.subplots_adjust(hspace = 0.4, bottom = 0.2)
plt.savefig('../Ch6_Figure07.png', dpi = 250, bbox_inches = 'tight')
plt.show()

