"""
Regrid the wind netCDF from staggered horizontal grid points and rho 
vertical levels onto theta grid points and theta levels.

Store the regridded winds to one netCDF for each wind hour containing both
wind components.

Creates a netCDF for three-hourly time chunks containing u, v, and w wind
components in the theta coordinate system.
"""

def main(path, l_spinup, l_short):
    import numpy as np
    import matplotlib.pyplot as plt
    from netCDF4 import Dataset
    from analysis_tools import regrid, transform_winds
    from datetime import datetime as dt
    import os
    
    # Use the naming convention for the input files
    if l_short:
        hours = ["{0:02d}".format(h) for h in xrange(0, 13, 4)]
    elif l_spinup:
        hours = ["{0:02d}".format(d) for d in xrange(1, 11)] # days for the spinup simulations
    else:
        hours = ["{0:02}".format(h) for h in xrange(0, 24, 3)]
    
    from STASH_keys import u_key, v_key, w_key, s_key, n_key
    
    for hour in hours:
        ########################################################################
        #                                                                      #
        # Section One: Read in the wind netCDF                                 #
        #                                                                      #
        ########################################################################
        
        print '[' + dt.now().strftime("%H:%M:%S") + '] Reading winds for ' + hour
        if 'u_' + hour + '.nc' in os.listdir(path):
            u_nc = Dataset(path + 'u_' + hour + '.nc', 'r')
            v_nc = Dataset(path + 'v_' + hour + '.nc', 'r')
            # Define a target coordinate system
            w_nc = Dataset(path + 'bouy_' + hour + '.nc', 'r')
            
            time_key = [i for i in w_nc.variables.keys() if 'min' in i][0]
            
            if l_diagnostics: print '[' + dt.now().strftime("%H:%M:%S") + '] Creating new winds netCDF'
            # Create a new netcdf for the regridded wind components
            winds_nc = Dataset(path + 'wind_' + hour + '.nc', 'w')
            
            if l_diagnostics: print '[' + dt.now().strftime("%H:%M:%S") + '] Creating dimensions for the netCDF'
            # Create dimensions for that netcdf
            time_dim = winds_nc.createDimension(time_key, w_nc.variables[w_key][:].shape[0])
            z_dim    = winds_nc.createDimension('thlev_zsea_theta', w_nc.variables[w_key][:].shape[1])
            lat_dim  = winds_nc.createDimension('latitude_t', w_nc.variables[w_key][:].shape[2])
            lon_dim  = winds_nc.createDimension('longitude_t', w_nc.variables[w_key][:].shape[3])
            
            if l_diagnostics: print '[' + dt.now().strftime("%H:%M:%S") + '] Creating variables for the netCDF'
            # Create variables for each dimension
            times_var = winds_nc.createVariable(time_key, np.float64, (time_key,))
            z_var     = winds_nc.createVariable('thlev_zsea_theta', np.float64, ('thlev_zsea_theta',))
            lats_var  = winds_nc.createVariable('latitude_t', np.float64, ('latitude_t','longitude_t'))
            lons_var  = winds_nc.createVariable('longitude_t', np.float64, ('latitude_t','longitude_t'))
            # Create the variables to store the wind components
            u_var   = winds_nc.createVariable(u_key, np.float64, (time_key, 'thlev_zsea_theta', 'latitude_t', 'longitude_t'))
            v_var   = winds_nc.createVariable(v_key, np.float64, (time_key, 'thlev_zsea_theta', 'latitude_t', 'longitude_t'))
            w_var   = winds_nc.createVariable(w_key, np.float64, (time_key, 'thlev_zsea_theta', 'latitude_t', 'longitude_t'))
            s_var   = winds_nc.createVariable(s_key, np.float64, (time_key, 'thlev_zsea_theta', 'latitude_t', 'longitude_t'))
            n_var   = winds_nc.createVariable(n_key, np.float64, (time_key, 'thlev_zsea_theta', 'latitude_t', 'longitude_t'))
            
            if l_diagnostics: print '[' + dt.now().strftime("%H:%M:%S") + '] Populating the dimension variables'
            # populate the dimension variables
            times_var[:] = w_nc.variables[time_key][:]
            z_var[:]     = w_nc.variables['thlev_zsea_theta'][:]
            lats_var[:]  = w_nc.variables['latitude_t'][:]
            lons_var[:]  = w_nc.variables['longitude_t'][:]
            
            if l_diagnostics: print '[' + dt.now().strftime("%H:%M:%S") + '] Populating the wind variables'
            if l_diagnostics: print '[' + dt.now().strftime("%H:%M:%S") + '] Regridding u-winds'
            u_regrid = regrid(w_nc, u_nc, u_key)[:]
            u_var[:] = u_regrid[:]*1.
            
            if l_diagnostics: print '[' + dt.now().strftime("%H:%M:%S") + '] Regridding v-winds'
            v_regrid = regrid(w_nc, v_nc, v_key)[:]
            v_var[:] = v_regrid[:]*1.
            
            if l_diagnostics: print '[' + dt.now().strftime("%H:%M:%S") + '] Copying over w-winds'
            w_var[:] = w_nc.variables[w_key][:]
            
            if l_diagnostics: print '[' + dt.now().strftime("%H:%M:%S") + '] Closing the old winds ncs'
            u_nc.close()
            v_nc.close()
            w_nc.close()
            
            if l_diagnostics: print '[' + dt.now().strftime("%H:%M:%S") + '] Transforming the winds'
            s_var[:], n_var[:] = transform_winds(u_regrid[:]*1., v_regrid[:]*1.)
            
            if l_diagnostics: print '[' + dt.now().strftime("%H:%M:%S") + '] Closing the new winds_nc'
            winds_nc.close()
            
            print '[' + dt.now().strftime("%H:%M:%S") + '] Regridding hour ' + hour + ' complete.\n'

