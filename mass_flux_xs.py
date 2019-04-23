import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from analysis_tools import fromComponents, bilinear_interpolation, get_cs_coords
from scipy import integrate, interpolate
from datetime import datetime as dt
plt.switch_backend('agg')

"""
Compute the mass flux across the cloud trail for our experiments, how does it
change for the changing environment at different heights and different 
distances downwind of the island.
"""

my_h = 100.
print '[' + dt.now().strftime('%H:%M') + '] Defining the paths to data...'
# Define the paths in which the data are stored
paths = {'Control'       : '/nerc/n02/n02/xb899100/CloudTrail/Control/',
         'B_1_3'         : '/nerc/n02/n02/xb899100/CloudTrail/H125/',
         'B_3.0'         : '/nerc/n02/n02/xb899100/CloudTrail/H375/',
         'Control_short' : '/nerc/n02/n02/xb899100/CloudTrail/Control_short/',
         'BL_RHm25'      : '/nerc/n02/n02/xb899100/CloudTrail/RH_BLm25/',
         'FA_RHm25'      : '/nerc/n02/n02/xb899100/CloudTrail/RH_FAm25/',
         'U05'           : '/nerc/n02/n02/xb899100/CloudTrail/U05/'}

print '[' + dt.now().strftime('%H:%M') + '] Determining short sims...'
short_sims = ['Control_short', 'BL_RHm25', 'FA_RHm25', 'U05']

print '[' + dt.now().strftime('%H:%M') + '] Defining the required variable keys...'
# Define the keys for rho and winds
rho_key = u'STASH_m01s00i389'
u_key   = u'STASH_m01s00i002'
v_key   = u'STASH_m01s00i003'
w_key   = u'STASH_m01s00i150'

print '[' + dt.now().strftime('%H:%M') + '] Starting the routine for each experiment...'
mass_flux_xs = {}
for key in paths.keys():
    print '[' + dt.now().strftime('%H:%M') + '] Starting experiment ' + key + '...'
    if key in short_sims:
        print '[' + dt.now().strftime('%H:%M') + '] ' + key + ' is a short experiment...'
        hour = '04'
        start_t = 300.
        end_t = 480.
    else:
        print '[' + dt.now().strftime('%H:%M') + '] ' + key + ' is a long experiment...'
        hour = '09'
        start_t = 540.
        end_t = 720.
    
    print '[' + dt.now().strftime('%H:%M') + '] Opening the netCDF...'
    # Read in rho and winds
    rho_nc  = Dataset(paths[key] + 'fluxes_' + hour + '.nc', 'r')
    wind_nc = Dataset(paths[key] + 'wind_' + hour + '.nc', 'r')
    
    print '[' + dt.now().strftime('%H:%M') + '] Reading the required dimensions...'
    # Want the cross section across the trail averaged in time for 3-hours prior
    # to the start of the cloud trail, find the appropriate times
    # Find the time keys for the respective netCDF
    rho_times_key = [tkey for tkey in rho_nc.variables.keys() if 'min' in tkey][0]
    wind_times_key = [tkey for tkey in wind_nc.variables.keys() if 'min' in tkey][0]
    # Read the times
    rho_times = rho_nc.variables[rho_times_key][:]
    wind_times = wind_nc.variables[wind_times_key][:]
    # Find the indexes for all the times in the preceeding period of interest
    if len(wind_times) > len(rho_times):
        # do the rho times first, and only do the matching wind times
        its_rho = [it for it in xrange(rho_times.shape[0]) if start_t <= rho_times[it] <= end_t]
        its_wind = [np.where(wind_times == rho_times[it])[0][0] for it in its_rho]
    else:
        # do the wind times first and only do the matching rho times
        its_wind = [it for it in xrange(wind_times.shape[0]) if start_t <= wind_times[it] <= end_t]
        its_rho = [np.where(rho_times == wind_times[it])[0][0] for it in its_wind]
    
    print '[' + dt.now().strftime('%H:%M') + '] Determining the mean wind direction...'
    # Find the average wind for that time period in u- and v- components
    u_0 = np.nanmean(wind_nc.variables[u_key][its_wind,:,:,:], axis = (0, 2, 3))
    v_0 = np.nanmean(wind_nc.variables[v_key][its_wind,:,:,:], axis = (0, 2, 3))
    # Integrate in height (1st get the height dimension)
    z_theta = wind_nc.variables['thlev_zsea_theta'][:]
    z_rho = rho_nc.variables['rholev_zsea_rho'][:]
    # Define upper height limit for the integration
    z_top = 500.
    iz_top = np.where(np.min(np.abs(z_theta - z_top)) == np.abs(z_theta - z_top))[0][0]
    # Calculate the wind direction
    u_mean = integrate.trapz(x = z_theta[:iz_top], y = u_0[:iz_top])/z_top
    v_mean = integrate.trapz(x = z_theta[:iz_top], y = v_0[:iz_top])/z_top
    wind_direction = fromComponents(u_mean, v_mean)[1]
    
    print '[' + dt.now().strftime('%H:%M') + '] The mean wind direction is ' + str(int(wind_direction)) + '...'
    # Find the cross sectional slice through the domain along the wind
    # This requires us to define where we want to originate the slice from
    # e.g. centre of the island, requires land mask
    print '[' + dt.now().strftime('%H:%M') + '] Defining the along-wind cross section...'
    landseamask = Dataset('/work/n02/n02/xb899100/island_masks/lsm50.nc', 'r')
    
    lsm = landseamask.variables['lsm'][0,0,:,:]*1.
    y = landseamask.variables['latitude'][:]*1.
    x = landseamask.variables['longitude'][:]*1.
    
    # create 2D coordinate mesh
    X, Y = np.meshgrid(x, y)
    
    landseamask.close()
    
    # We know from our domain definition where the centre of the island should be
    R_i = 1000.0*(50.0/np.pi)**0.5 # island radius
    x_c = 100000.0 + R_i
    y_c = 4*R_i
    
    # From this center point, draw line in either direction parallel to wind_dir_0
    # get the coordinates of this line at 100 m resolution.
    x_cs, y_cs = get_cs_coords(x_c, y_c, wind_direction, X, Y, h = my_h)
    
    R = -np.sign(x_cs - x_c)*np.sqrt((x_cs - x_c)**2 + (y_cs - y_c)**2)
    
    print '[' + dt.now().strftime('%H:%M') + '] Defining the across-wind cross section, 10km downwind of the leeward island edge...'
    # Define a target distance downwind of the island mid-point
    r_target = 0. + R_i # in metres, the plus R_i makes this the distance from island edge
    
    # Find where the distance to the island is nearest to the target distance
    R_pos = np.where((R > 0), R, np.nan)
    iR = np.where(np.abs(R_pos - r_target) == np.nanmin(np.abs(R_pos - r_target)))[0][0]
    
    x_cs2, y_cs2 = get_cs_coords(x_cs[iR], y_cs[iR], wind_direction+90., X, Y, h = my_h, max_r = 5000.)
    # get the distance from away from the along-flow cross section
    R_across = np.sign(y_cs2 - y_cs[iR])*np.sqrt((x_cs2 - x_cs[iR])**2 + (y_cs2 - y_cs[iR])**2)
    
    # now that we have these coordinates for our across trail cross section, we need to interpolate to them
    # Take the time mean
    print '[' + dt.now().strftime('%H:%M') + '] Reading the data data...'
    #rho0 = rho_nc.variables[rho_key][its_rho,:,:,:]
    w = wind_nc.variables[w_key][its_wind,:,:,:]
    print '[' + dt.now().strftime('%H:%M') + '] Interpolating rho to theta levels...'
    rho = np.zeros_like(w)
    for it in its_rho:
        rho[its_rho.index(it),:,:,:] = interpolate.interp1d(x = z_rho, y = rho_nc.variables[rho_key][it,:,:,:], axis = 0, fill_value = 'extrapolate')(z_theta)
    
    print '[' + dt.now().strftime('%H:%M') + '] Closing netCDF...'
    rho_nc.close()
    wind_nc.close()
    
    print '[' + dt.now().strftime('%H:%M') + '] Computing mass flux...'
    mf = np.nanmean(rho*w, axis = 0)
    
    print '[' + dt.now().strftime('%H:%M') + '] Interpolating to the across-wind cross section coordinates n = 1...'
    # interpolate rho to w coordinates
    mass_flux_xs[key] = bilinear_interpolation(X, Y, mf, x_cs2, y_cs2, kind = 2)
    count = 1
    # take the mean of island diameter sized chunk, means ~80 points
    n_max = 5*int((R_i*2)/my_h)
    for ir in xrange(1, n_max + 1):
        count += 1
        print '[' + dt.now().strftime('%H:%M') + '] Interpolating to the across-wind cross section coordinates n = ' + str(count) + '...'
        x_cs2, y_cs2 = get_cs_coords(x_cs[iR+ir], y_cs[iR+ir], wind_direction+90., X, Y, h = my_h, max_r = 5000.)
        mass_flux_xs[key] += bilinear_interpolation(X, Y, mf, x_cs2, y_cs2, kind = 2)
    
    mass_flux_xs[key] /= float(count)
    
    z_target = 350.
    print '[' + dt.now().strftime('%H:%M') + '] Plotting the mass flux for ' + key + ' near ' + str(int(z_target)) + ' m...'
    iz_target = np.where(np.abs(z_theta - z_target) == np.min(np.abs(z_theta - z_target)))[0][0]
    plt.plot(R_across, mass_flux_xs[key][iz_target,:], label = key)

plt.xlim([-5000, 5000])
plt.plot([-5000, 5000], [0,0], 'grey', ls = '--')
plt.title('Mass flux ' + str(int((r_target - R_i)/1000.)) + ' km downwind of leeward side averaged over the 3hrs prior\n to peak heating (9-12pm) for each experiment')
plt.ylabel('Mass Flux ($kg m^{-2} s^{-1}$)')
plt.xlabel('Distance from along wind centrepoint (i.e. y$^{\prime}$, m)')
plt.legend(loc = 0)
plt.savefig('../mass_flux_' + str(int((r_target - R_i)/1000.)) + 'km_downwind_all_exp.png', dpi = 150)
#plt.show()

