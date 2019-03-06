import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from netCDF4 import Dataset
from scipy import interpolate
from datetime import datetime as dt
from analysis_tools import zi, lcl, find_h, get_CTZ

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

# Define some constants
Rd = 287.05
cp = 1005.
p0 = 1e5

# Do we want to create zi netCDF or read from zi netCDF?
create_netCDF = True

# Define keys to read variables from UM output
pres_key  = u'STASH_m01s00i408'
theta_key = u'STASH_m01s00i004'
q_key     = u'STASH_m01s00i010' # specific humidity returns a masked array where values are below the reset at the first time step
temp_key  = u'STASH_m01s16i004'
mcl_key   = u'STASH_m01s00i392'
u_key     = u'STASH_m01s00i002'
v_key     = u'STASH_m01s00i003'

# Define keys to create variables in our zi netCDFs
zi0_key = u'boundary layer depth'
zi1_key = u'new boundary layer depth'
lcl_key = u'lifting condensation level'
ctz_key = u'cloud top height'
paths = ['/work/n02/n02/xb899100/cylc-run/u-bg023/share/data/history/', '/work/n02/n02/xb899100/cylc-run/u-bg113/share/data/history/']
for path in paths:
    for hour in hours:
        if create_netCDF:
            print '[' + dt.now().strftime("%H:%M:%S") + '] Open the netCDF on hour ' + hour
            print '[' + dt.now().strftime("%H:%M:%S") + '] Open pres_nc'
            pres_nc   = Dataset(path + 'fluxes_'+hour+'.nc', 'r')

            print '[' + dt.now().strftime("%H:%M:%S") + '] Open theta_nc'
            theta_nc  = Dataset(path + 'bouy_'+hour+'.nc', 'r')
            
            print '[' + dt.now().strftime("%H:%M:%S") + '] Open mr_nc'
            mr_nc  = Dataset(path + 'mr_'+hour+'.nc', 'r')
            
            print '[' + dt.now().strftime("%H:%M:%S") + '] Open wind_nc'
            wind_nc  = Dataset(path + 'wind_'+hour+'.nc', 'r')
            
            if hour == '00':
                print '[' + dt.now().strftime("%H:%M:%S") + '] Grabbing the height coordinates'
                time_key = [i for i in wind_nc.variables.keys() if 'min' in i][0]
                z_theta = theta_nc.variables['thlev_zsea_theta'][:]
                it_start = 1 # because we have moisture resetting, the first time step is filled with missing data for q
            else:
                it_start = 0
            
            print '[' + dt.now().strftime("%H:%M:%S") + '] Creating zi netCDF'
            # Create a new netcdf for the boundary layer height each day
            zi_nc = Dataset(path + 'zi_'+hour+'.nc', 'w')
            
            print '[' + dt.now().strftime("%H:%M:%S") + '] Creating dimensions for the netCDF'
            # create dimensions for that netcdf
            time_dim = zi_nc.createDimension('time', wind_nc.variables[time_key][it_start:].shape[0])
            lat_dim  = zi_nc.createDimension('lat', 319)
            lon_dim  = zi_nc.createDimension('lon', 1160)
            
            print '[' + dt.now().strftime("%H:%M:%S") + '] Creating variables for the netCDF'
            # create variables for each dimension
            times_var = zi_nc.createVariable('time', np.float64, ('time',))
            lats_var  = zi_nc.createVariable('latitude', np.float64, ('lat','lon'))
            lons_var  = zi_nc.createVariable('longitude', np.float64, ('lat','lon'))
            # create the variable to store the boundary layer heights
            zi0_var   = zi_nc.createVariable(zi0_key, np.float64, ('time', 'lat', 'lon'))
            zi1_var   = zi_nc.createVariable(zi1_key, np.float64, ('time', 'lat', 'lon'))
            lcl_var   = zi_nc.createVariable(lcl_key, np.float64, ('time', 'lat', 'lon'))
            ctz_var   = zi_nc.createVariable(ctz_key, np.float64, ('time', 'lat', 'lon'))
            
            print '[' + dt.now().strftime("%H:%M:%S") + '] Populating the dimension variables'
            # populate the dimension variables
            times_var[:] = theta_nc.variables[time_key][it_start:]*1.
            lats_var[:]  = theta_nc.variables['latitude_t'][:]*1.
            lons_var[:]  = theta_nc.variables['longitude_t'][:]*1.
            
            for it in xrange(it_start, len(theta_nc.variables[time_key][:])):
                print '[' + dt.now().strftime("%H:%M:%S") + '] Working on time slice number ' + str(it) + ', for netCDF hour ' + str(hour)
                theta_v = theta_nc.variables[theta_key][it,:,:,:]*(1. + 0.608*theta_nc.variables[q_key][it,:,:,:])
                
                # populate the boundary layer height variable
                print '[' + dt.now().strftime("%H:%M:%S") + '] Calculating zi0...'
                zi0_var[it-it_start,:,:] = zi(theta_v, z_theta)[:]
                print '[' + dt.now().strftime("%H:%M:%S") + '] Calculating zi1...'
                zi1_var[it-it_start,:,:] = find_h(theta_v, wind_nc.variables[u_key][it,:,:,:], wind_nc.variables[v_key][it,:,:,:], z_theta)
                print '[' + dt.now().strftime("%H:%M:%S") + '] Calculating lcl...'
                lcl_var[it-it_start,:,:] = lcl(theta_nc.variables[temp_key][it,0,:,:]*1.,theta_nc.variables[q_key][it,0,:,:]*1., pres_nc.variables[pres_key][it,:,:,:]*1., z_theta)
                print '[' + dt.now().strftime("%H:%M:%S") + '] Calculating ctz...'
                ctz_var[it-it_start,:,:] = get_CTZ(mr_nc.variables[mcl_key][it,:,:,:], z_theta)
            # close the new netcdf
            zi_nc.close()
            print '[' + dt.now().strftime("%H:%M:%S") + '] Closing the netCDF'
            theta_nc.close()
            pres_nc.close()
            mr_nc.close()
            wind_nc.close()
        
        print '[' + dt.now().strftime("%H:%M:%S") + '] Opening zi_'+hour+'.nc'
        zi_lcl = Dataset(path + 'zi_' + hour + '.nc', 'r')

        zi0_data = zi_lcl.variables[zi0_key][:]
        zi1_data = zi_lcl.variables[zi1_key][:]
        lcl_data = zi_lcl.variables[lcl_key][:]
        zi_lcl.close()

        if hour == '00':
            print '[' + dt.now().strftime("%H:%M:%S") + '] This is the first time chunk, initialise the time series arrays'
            zi0_ts = np.zeros(1)
            zi1_ts = np.zeros(1)
            lcl_ts = np.zeros(1)

        print '[' + dt.now().strftime("%H:%M:%S") + '] Concatenating zi and lcl'
        zi0_ts = np.concatenate((zi0_ts, np.mean(zi0_data, axis = (1, 2))))
        zi1_ts = np.concatenate((zi1_ts, np.mean(zi1_data, axis = (1, 2))))
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

    fig = plt.figure()

    ax = fig.add_subplot(1, 3, 1)
    ax.plot(times/60., zi0_ts, 'r', lw = 2)
    ax.plot(times/60., zi1_ts, 'b', lw = 1)
    ax.plot(times/60., lcl_ts, 'k--')
    ax.set_ylim([400, 1000])
    ax.set_ylabel('Height (m)')
    ax.set_xlim([0, 24])
    ax.set_xticks(range(0, 25, 6))
    ax.set_xlabel('Time (hours)')
    plt.legend(loc = 2)

    ax1 = fig.add_subplot(1, 3, 2)
    ax1.plot(times/60., lwp_ts*1000., 'k:', lw = 2)
    ax1.set_ylabel('Liquid Water Path (g m$^{-2}$)')
    ax1.set_ylim([0, 30])
    ax1.set_xlim([0, 24])
    ax1.set_xticks(range(0, 25, 6))
    ax1.set_xlabel('Time (hours)')

    ax2 = fig.add_subplot(1, 3, 3)
    ax2.plot(times/60., cc_ts, 'k')
    ax2.set_ylabel('Cloud Cover')
    ax2.set_ylim([0, 0.4])
    ax2.set_xlim([0, 24])
    ax2.set_xticks(range(0, 25, 6))
    ax2.set_xlabel('Time (hours)')

    fig.tight_layout()
    plt.savefig('../zi_lcl_cc_lwp_exp0'+ str(paths.index(path)) + '.png', dpi = 150)
    plt.close('all')


