import numpy as np
import matplotlib.pyplot as plt
#plt.switch_backend('agg')
from netCDF4 import Dataset
from scipy import interpolate
from datetime import datetime as dt
from analysis_tools import zi, lcl, find_h, get_CTZ
from SkewT_archer import PTtoTemp

"""
Code to calculate the mixed layer depth from a simple level of neutral buoyancy
method, and calculate the surface based lifting condensation level.

Then compute the domain mean of these quantities and plot their time series.
"""
l_short = 0
l_spinup = 1
# Read in the data
if l_short:
    hours = ["{0:02d}".format(hour) for hour in xrange(0,13,4)]
elif l_spinup:
    hours = ["{0:02d}".format(d) for d in xrange(1, 11)] # days for the spinup simulation
else:
    hours = ["{0:02d}".format(hour) for hour in xrange(0, 24, 3)]

nhours = len(hours)

# for our 1-day simulation there are 1440 minutes
# with 10 minute output, there are 144 time steps
zi_ts  = np.zeros(1)
lcl_ts = np.zeros(1)

# Define some constants
Rd = 287.05
cp = 1005.
p0 = 1e5

# Do we want to create zi netCDF or read from zi netCDF?
create_netCDF = 0

# Define keys to read variables from UM output
pres_key  = u'STASH_m01s00i408' # Pressure on theta levels
pres_key1 = u'STASH_m01s00i407' # Pressure on rho levels (to be used for spin-up sims
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
exp01_path        = '/work/n02/n02/xb899100/cylc-run/u-bg023/share/data/history/'
exp02_path        = '/work/n02/n02/xb899100/cylc-run/u-bg113/share/data/history/'
control_path      = '/work/n02/n02/xb899100/cylc-run/u-bd527/share/data/history/'
U10_spinup_path   = '/nerc/n02/n02/xb899100/CloudTrail/Control_Spinup/'
U05_spinup_path   = '/nerc/n02/n02/xb899100/CloudTrail/U05_Spinup/'
RHm25_path        = '/work/n02/n02/xb899100/cylc-run/u-bg665/share/data/history/'
U05v2_spinup_path = '/work/n02/n02/xb899100/cylc-run/u-bg952/share/data/history/'

#paths = [U10_spinup_path, U05v2_spinup_path, U05_spinup_path]
#ID = ['Control_spinup', 'U05v2_spinup', 'U05_spinup']

paths = [U05_spinup_path]
ID = ['U05_spinup']

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
            
            if hour == hours[0]:
                print '[' + dt.now().strftime("%H:%M:%S") + '] Grabbing the height coordinates'
                time_key = [i for i in wind_nc.variables.keys() if 'min' in i][0]
                z_theta = theta_nc.variables['thlev_zsea_theta'][:]
                z_rho = pres_nc.variables['rholev_zsea_rho'][:]
                it_start = 1 # because we have moisture resetting, the first time step is filled with missing data for q
                dim_y, dim_x = theta_nc.variables[theta_key][0,0,:,:].shape
            else:
                it_start = 0
            
            print '[' + dt.now().strftime("%H:%M:%S") + '] Creating zi netCDF'
            # Create a new netcdf for the boundary layer height each day
            zi_nc = Dataset(path + 'zi_' + hour + '.nc', 'w')
            
            print '[' + dt.now().strftime("%H:%M:%S") + '] Creating dimensions for the netCDF'
            # create dimensions for that netcdf
            time_dim = zi_nc.createDimension('time', wind_nc.variables[time_key][it_start:].shape[0])
            lat_dim  = zi_nc.createDimension('lat', dim_y)
            lon_dim  = zi_nc.createDimension('lon', dim_x)
            
            print '[' + dt.now().strftime("%H:%M:%S") + '] Creating variables for the netCDF'
            # create variables for each dimension
            times_var = zi_nc.createVariable('time', np.float64, ('time',))
            lats_var  = zi_nc.createVariable('latitude', np.float64, ('lat','lon'))
            lons_var  = zi_nc.createVariable('longitude', np.float64, ('lat','lon'))
            # create the variable to store the boundary layer heights
            zi0_var   = zi_nc.createVariable(zi0_key, np.float64, ('time', 'lat', 'lon'))
            zi1_var   = zi_nc.createVariable(zi1_key, np.float64, ('time', 'lat', 'lon'))
            lcl_var   = zi_nc.createVariable(lcl_key, np.float64, ('time', 'lat', 'lon'))
            if mcl_key in mr_nc.variables.keys():
                ctz_var   = zi_nc.createVariable(ctz_key, np.float64, ('time', 'lat', 'lon'))
            
            print '[' + dt.now().strftime("%H:%M:%S") + '] Populating the dimension variables'
            # populate the dimension variables
            times_var[:] = theta_nc.variables[time_key][it_start:]*1.
            lats_var[:]  = theta_nc.variables['latitude_t'][:]*1.
            lons_var[:]  = theta_nc.variables['longitude_t'][:]*1.
            
            for it in xrange(it_start, len(theta_nc.variables[time_key][:])):
                print '[' + dt.now().strftime("%H:%M:%S") + '] Working on time slice number ' + str(it) + ', for netCDF hour ' + str(hour)
                if not l_spinup:
                    if type(theta_nc.variables[q_key][it,:,:,:]) == np.ma.core.MaskedArray:
                        q_data = theta_nc.variables[q_key][it,:,:,:].data
                    else:
                        q_data = theta_nc.variables[q_key][it,:,:,:]
                    theta_v = theta_nc.variables[theta_key][it,:,:,:]*(1. + 0.608*q_data)
                else:
                    theta_v = theta_nc.variables[theta_key][it,:,:,:]*(1. + 0.608*mr_nc.variables[q_key][it,:,:,:])
                # populate the boundary layer height variable
                print '[' + dt.now().strftime("%H:%M:%S") + '] Calculating zi0...'
                zi0_var[it-it_start,:,:] = zi(theta_v, z_theta)[:]
                print '[' + dt.now().strftime("%H:%M:%S") + '] Calculating zi1...'
                zi1_var[it-it_start,:,:] = find_h(theta_v, wind_nc.variables[u_key][it,:,:,:], wind_nc.variables[v_key][it,:,:,:], z_theta)
                print '[' + dt.now().strftime("%H:%M:%S") + '] Calculating lcl...'
                if not l_spinup:
                    lcl_var[it-it_start,:,:] = lcl(theta_nc.variables[temp_key][it,0,:,:]*1.,theta_nc.variables[q_key][it,0,:,:]*1., pres_nc.variables[pres_key][it,:,:,:]*1., z_theta)
                else:
                    pressure = interpolate.interp1d(x = z_rho, y = pres_nc.variables[pres_key1][it,:,:,:], fill_value = 'extrapolate', axis = 0)(z_theta)
                    temperature = PTtoTemp(theta_nc.variables[theta_key][it,0,:,:]*1., pressure[0,:,:], t_units = 'K', p_units = 'Pa')
                    lcl_var[it-it_start,:,:] = lcl(temperature,mr_nc.variables[q_key][it,0,:,:]*1., pressure, z_theta)
                if mcl_key in mr_nc.variables.keys():
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
        if hour == hours[0]:
            zi0_data = zi_lcl.variables[zi0_key][:]
            zi1_data = zi_lcl.variables[zi1_key][:]
            lcl_data = zi_lcl.variables[lcl_key][:]
            times    = zi_lcl.variables['time'][:]
        else:
            zi0_data = np.concatenate((zi0_data, zi_lcl.variables[zi0_key][:]), axis = 0)
            zi1_data = np.concatenate((zi1_data, zi_lcl.variables[zi1_key][:]), axis = 0)
            lcl_data = np.concatenate((lcl_data, zi_lcl.variables[lcl_key][:]), axis = 0)
            times    = np.concatenate((times, zi_lcl.variables['time'][:]), axis = 0)
        zi_lcl.close()
        
    print '[' + dt.now().strftime("%H:%M:%S") + '] Creating time series'
    zi0_ts = np.mean(zi0_data, axis = (1, 2))
    zi1_ts = np.mean(zi1_data, axis = (1, 2))
    lcl_ts = np.mean(lcl_data, axis = (1, 2))
    
    # cloud cover and liquid water path stuff
    print '[' + dt.now().strftime("%H:%M:%S") + '] Opening lwp netCDF...'
    lwp_key   = 'STASH_m01s30i405'
    if l_spinup:
        for hour in hours:
            lwp_nc = Dataset(path + 'lwp_' + hour + '.nc', 'r')
            if hour == hours[0]:
                times_lwp = lwp_nc.variables[[key for key in lwp_nc.variables.keys() if 'min' in key][0]][:]
                factor = int((times[1] - times[0])/(times_lwp[1] - times_lwp[0]))
                print '[' + dt.now().strftime("%H:%M:%S") + '] Creating the mean lwp and cloud cover time series...'
                lwp_ts = np.zeros_like(times)
                cc_ts = np.zeros_like(times)
                lwp_data = lwp_nc.variables[lwp_key][:]*1.
            else:
                lwp_data = np.concatenate((lwp_data, lwp_nc.variables[lwp_key][:]*1.), axis = 0)
        
        for i in xrange(times.shape[0]):
            lwp_ts[i] = np.mean(lwp_data[factor*i,:,:])
            cloud_mask = np.where((lwp_data[factor*i,:,:] > 1e-08), 1., 0.0)
            cc_ts[i] = np.mean(cloud_mask)
    else:
        lwp_nc = Dataset(path + 'lwp_00.nc', 'r')
        times_lwp = lwp_nc.variables[[key for key in lwp_nc.variables.keys() if 'min' in key][0]][:]
        factor = int((times[1] - times[0])/(times_lwp[1] - times_lwp[0]))
        print '[' + dt.now().strftime("%H:%M:%S") + '] Creating the mean lwp and cloud cover time series...'
        lwp_ts = np.zeros_like(times)
        cc_ts = np.zeros_like(times)
        
        for i in xrange(times.shape[0]):
            lwp_ts[i] = np.mean(lwp_nc.variables[lwp_key][factor*i,:,:])
            cloud_mask = np.where((lwp_nc.variables[lwp_key][factor*i,:,:] > 1e-08), 1., 0.0)
            cc_ts[i] = np.mean(cloud_mask)
    
    # reset the missing data in the first index to be 0
    lwp_ts[0] = 0.0
    cc_ts[0] = 0.0
    
    lwp_nc.close()
    
    print '[' + dt.now().strftime("%H:%M:%S") + '] Making the plot...'
    if l_short:
        times = times/60. + 4.
    if l_spinup:
        times = times/1440.
    else:
        times /= 60.
    
    if l_spinup:
        fig = plt.figure(tight_layout = True, figsize = (10, 12))
        ax = fig.add_subplot(3, 1, 1)
        ax1 = fig.add_subplot(3, 1, 2)
        ax2 = fig.add_subplot(3, 1, 3)
        ax.set_xlim([0, 10])
        ax.set_xlabel('Time (days)')
        ax1.set_xlim([0, 10])
        ax1.set_xlabel('Time (days)')
        ax2.set_xlim([0, 10])
        ax2.set_xlabel('Time (days)')
        ax.set_title(ID[paths.index(path)])
    else:
        fig = plt.figure(tight_layout=True,figsize=(15,5))
        ax = fig.add_subplot(1, 3, 1)
        ax1 = fig.add_subplot(1, 3, 2)
        ax2 = fig.add_subplot(1, 3, 3)
        ax.set_xlim([0, 24])
        ax.set_xticks(range(0, 25, 6))
        ax.set_xlabel('Time (hours)')
        ax1.set_xlim([0, 24])
        ax1.set_xticks(range(0, 25, 6))
        ax1.set_xlabel('Time (hours)')
        ax2.set_xlim([0, 24])
        ax2.set_xticks(range(0, 25, 6))
        ax2.set_xlabel('Time (hours)')
    
    ax.plot(times, zi0_ts, 'r', lw = 2, label = 'z$_i$ Parcel Method')
    ax.plot(times, zi1_ts, 'b', lw = 1, label = 'z$_i$ Rib Method')
    ax.plot(times, lcl_ts, 'k--', label = 'LCL')
    ax.set_ylabel('Height (m)')
    ax.legend(loc = 0)
    
    ax1.plot(times, lwp_ts*1000., 'k:', lw = 2)
    ax1.set_ylabel('Liquid Water Path (g m$^{-2}$)')
    
    ax2.plot(times, cc_ts, 'k')
    ax2.set_ylabel('Cloud Cover')
    
    plt.savefig('../zi_lcl_cc_lwp_' + ID[paths.index(path)] + '.png', dpi = 150)
    plt.show()
    plt.close('all')


