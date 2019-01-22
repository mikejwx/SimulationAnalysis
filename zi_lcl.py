import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import interpolate
from datetime import datetime as dt
from analysis_tools import zi, lcl

"""
Code to calculate the mixed layer depth from a simple level of neutral buoyancy
method, and calculate the surface based lifting condensation level.

Then compute the domain mean of these quantities and plot their time series.
"""

# Read in the data
hours = ["{0:02d}".format(hour) for hour in xrange(0,24,3)]
nhours = len(hours)

# for our 1-day simulation there are 1440 minutes
# with 10 minute output, there are 144 time steps
times  = np.arange(0., 1440.1, 10.)
ntimes = times.shape[0]
zi_ts  = np.zeros(1)
lcl_ts = np.zeros(1)

Rd = 287.05
cp = 1005.
p0 = 1e5
create_netCDF = True
pres_key  = 'STASH_m01s00i408'
theta_key = 'STASH_m01s00i004'
q_key     = 'STASH_m01s00i010' # specific humidity returns a masked array where values are below the reset (e.g. at the first time step)
temp_key  = 'STASH_m01s16i004'
z_i_key = 'boundary layer depth'
lcl_key = 'lifting condensation level'
for hour in hours:
    if create_netCDF:
        print '[' + dt.now().strftime("%H:%M:%S") + '] Open the netCDF on hour ' + hour
        print '[' + dt.now().strftime("%H:%M:%S") + '] Open pres_nc'
        pres_nc   = Dataset('fluxes_'+hour+'.nc', 'r')

        print '[' + dt.now().strftime("%H:%M:%S") + '] Open theta_nc'
        theta_nc  = Dataset('bouy_'+hour+'.nc', 'r')

        if hour == '00':
            print '[' + dt.now().strftime("%H:%M:%S") + '] Grabbing the height coordinates'
            z_theta = theta_nc.variables['thlev_zsea_theta'][:]
            it_start = 1 # because we have moisture resetting, the first time step is filled with missing data for q
        else:
            it_start = 0

        print '[' + dt.now().strftime("%H:%M:%S") + '] Creating zi netCDF'
        # Create a new netcdf for the boundary layer height each day
        zi_nc = Dataset('zi_'+hour+'.nc', 'w')
        print '[' + dt.now().strftime("%H:%M:%S") + '] Creating dimensions for the netCDF'
        # create dimensions for that netcdf
        time_dim = zi_nc.createDimension('time', theta_nc.variables['min10_0'][it_start:].shape[0])
        lat_dim  = zi_nc.createDimension('lat', 319)
        lon_dim  = zi_nc.createDimension('lon', 1160)
        print '[' + dt.now().strftime("%H:%M:%S") + '] Creating variables for the netCDF'
        # create variables for each dimension
        times_var = zi_nc.createVariable('time', np.float64, ('time',))
        lats_var  = zi_nc.createVariable('latitude', np.float64, ('lat','lon'))
        lons_var  = zi_nc.createVariable('longitude', np.float64, ('lat','lon'))
        # create the variable to store the boundary layer heights
        zi_var    = zi_nc.createVariable(z_i_key, np.float64, ('time', 'lat', 'lon'))
        lcl_var   = zi_nc.createVariable(lcl_key, np.float64, ('time', 'lat', 'lon'))
        print '[' + dt.now().strftime("%H:%M:%S") + '] Populating the dimension variables'
        # populate the dimension variables
        times_var[:] = theta_nc.variables['min10_0'][it_start:]*1.
        lats_var[:]  = theta_nc.variables['latitude_t'][:]*1.
        lons_var[:]  = theta_nc.variables['longitude_t'][:]*1.
        
        for it in xrange(it_start, len(theta_nc.variables['min10_0'][:])):
            print '[' + dt.now().strftime("%H:%M:%S") + '] Working on time slice number ' + str(it) + ', for netCDF hour ' + str(hour)
            theta_v = theta_nc.variables[theta_key][it,:,:,:]*(1. + 0.608*theta_nc.variables[q_key][it,:,:,:])
            
            # populate the boundary layer height variable
            print '[' + dt.now().strftime("%H:%M:%S") + '] Calculating zi...'
            zi_var[it-it_start,:,:] = zi(theta_v, z_theta)[:]
            print '[' + dt.now().strftime("%H:%M:%S") + '] Calculating lcl...'
            lcl_var[it-it_start,:,:] = lcl(theta_nc.variables[temp_key][it,0,:,:]*1.,theta_nc.variables[q_key][it,0,:,:]*1., pres_nc.variables[pres_key][it,:,:,:]*1., z_theta)
            
            print '[' + dt.now().strftime("%H:%M:%S") + '] Concatenating the z_i for hour ' + hour
            zi_ts = np.concatenate((zi_ts, np.array([np.nanmean(zi_var[it-it_start,:,:])])), axis = 0)
            
            print '[' + dt.now().strftime("%H:%M:%S") + '] Concatenating the lcl for hour ' + hour
            lcl_ts = np.concatenate((lcl_ts, np.array([np.nanmean(lcl_var[it-it_start,:,:])])), axis = 0)
        
        # close the new netcdf
        zi_nc.close()
        print '[' + dt.now().strftime("%H:%M:%S") + '] Closing the netCDF'
        theta_nc.close()
        pres_nc.close()
    else:
        print '[' + dt.now().strftime("%H:%M:%S") + '] Open the netCDF for hour ' + hour
        print '[' + dt.now().strftime("%H:%M:%S") + '] Opening zi_'+hour+'.nc'
        zi_lcl = Dataset('zi_' + hour + '.nc', 'r')

        z_i_data = zi_lcl.variables[z_i_key][:]*1.
        lcl_data = zi_lcl.variables[lcl_key][:]*1.
        zi_lcl.close()

        if hour == '00':
            print '[' + dt.now().strftime("%H:%M:%S") + '] This is the first time chunk, initialise the time series arrays'
            zi_ts = np.zeros(1)
            lcl_ts = np.zeros(1)

        print '[' + dt.now().strftime("%H:%M:%S") + '] Concatenating zi and lcl'
        zi_ts = np.concatenate((zi_ts, np.mean(z_i_data, axis = (1, 2))))
        lcl_ts = np.concatenate((lcl_ts, np.mean(lcl_data, axis = (1, 2))))

# cloud cover and liquid water path stuff
print '[' + dt.now().strftime("%H:%M:%S") + '] Opening lwp netCDF...'
lwp_nc = Dataset('lwp_00.nc', 'r')
lwp_key   = 'STASH_m01s30i405'

print '[' + dt.now().strftime("%H:%M:%S") + '] Creating the mean lwp and cloud cover time series...'
lwp_ts = np.zeros_like(times)
cc_ts = np.zeros_like(times)
for i in xrange(times.shape[0]):
    lwp_ts[i] = np.mean(lwp_nc.variables[lwp_key][2*i,:,:])
    cloud_mask = np.where((lwp_nc.variables[lwp_key][2*i,:,:] > 1e-08), 1., 0.0)
    cc_ts[i] = np.mean(cloud_mask)

# reset the missing data in the first index to be 0
lwp_ts[0] = 0.0
cc_ts[0] = 0.0

lwp_nc.close()

print '[' + dt.now().strftime("%H:%M:%S") + '] Making the plot...'
plt.subplot(221)
item1, = plt.plot(times/60., zi_ts, 'k', lw = 2)
plt.ylim([500, 900])
plt.ylabel('Height (m)')
plt.xlim([0, 24])
plt.xticks(range(0, 25, 6))
plt.xlabel('Time (hours)')

plt.subplot(222)
item2, = plt.plot(times/60., lcl_ts, 'k--', lw = 2)
plt.ylabel('Height (m)')
plt.ylim([500, 900])
plt.xlim([0, 24])
plt.xticks(range(0, 25, 6))
plt.xlabel('Time (hours)')

plt.subplot(223)
item3, = plt.plot(times/60., lwp_ts*1000., 'k:', lw = 2)
plt.ylabel('Liquid Water Path (g/m$^2$)')
plt.ylim([0, 30])
plt.xlim([0, 24])
plt.xticks(range(0, 25, 6))
plt.xlabel('Time (hours)')

plt.subplot(224)
item4, = plt.plot(times/60., cc_ts, 'k')
plt.ylabel('Cloud Cover')
plt.ylim([0, 0.4])
plt.xlim([0, 24])
plt.xticks(range(0, 25, 6))
plt.xlabel('Time (hours)')

lgd = plt.legend([item1, item2, item3, item4], ['$z_i$ (m)','$lcl$ (m)','liquid water path (g/m$^2$)','total cloud cover'], loc = 'center left', bbox_to_anchor = (1., 0.5), ncol = 1)
plt.tight_layout()
plt.savefig('boundary_layer_changes_sfc_new.png', bbox_extra_artists=(lgd,), dpi = 150, bbox_inches = 'tight')
plt.show()


