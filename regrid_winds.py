import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from analysis_tools import regrid
from datetime import datetime as dt

def main():
    """
    Regrid the wind netCDF from staggered horizontal grid points and rho 
    vertical levels onto theta grid points and theta levels.

    Store the regridded winds to one netCDF for each wind hour containing both
    wind components.

    Creates a netCDF for three-hourly time chunks containing u, v, and w wind
    components in the theta coordinate system.
    """

    # Use the naming convention for the input files
    hours = ["{0:02d}".format(h) for h in xrange(0, 24, 3)]
    u_key = u'STASH_m01s00i002'
    v_key = u'STASH_m01s00i003'
    w_key = u'STASH_m01s00i150'

    for hour in hours:
        ########################################################################
        #                                                                      #
        # Section One: Read in the wind netCDF                                 #
        #                                                                      #
        ########################################################################
        
        print '[' + dt.now().strftime("%H:%M:%S") + '] Reading winds_' + hour + '.nc'
        u_nc = Dataset('../u_' + hour + '.nc', 'r')
        v_nc = Dataset('../v_' + hour + '.nc', 'r')
        # Define a target coordinate system
        w_nc = Dataset('../bouy_' + hour + '.nc', 'r')
        
        print '[' + dt.now().strftime("%H:%M:%S") + '] Creating new winds netCDF'
        # Create a new netcdf for the regridded wind components
        winds_nc = Dataset('../wind_' + hour + '.nc', 'w')
        
        print '[' + dt.now().strftime("%H:%M:%S") + '] Creating dimensions for the netCDF'
        # Create dimensions for that netcdf
        time_dim = winds_nc.createDimension('min10_0', w_nc.variables[w_key][:].shape[0])
        z_dim    = winds_nc.createDimension('thlev_zsea_theta', w_nc.variables[w_key][:].shape[1])
        lat_dim  = winds_nc.createDimension('latitude_t', w_nc.variables[w_key][:].shape[2])
        lon_dim  = winds_nc.createDimension('longitude_t', w_nc.variables[w_key][:].shape[3])
        
        print '[' + dt.now().strftime("%H:%M:%S") + '] Creating variables for the netCDF'
        # Create variables for each dimension
        times_var = winds_nc.createVariable('min10_0', np.float64, ('min10_0',))
        z_var     = winds_nc.createVariable('thlev_zsea_theta', np.float64, ('thlev_zsea_theta',))
        lats_var  = winds_nc.createVariable('latitude_t', np.float64, ('latitude_t','longitude_t'))
        lons_var  = winds_nc.createVariable('longitude_t', np.float64, ('latitude_t','longitude_t'))
        # Create the variables to store the wind components
        u_var   = winds_nc.createVariable(u_key, np.float64, ('min10_0', 'thlev_zsea_theta', 'latitude_t', 'longitude_t'))
        v_var   = winds_nc.createVariable(v_key, np.float64, ('min10_0', 'thlev_zsea_theta', 'latitude_t', 'longitude_t'))
        w_var   = winds_nc.createVariable(w_key, np.float64, ('min10_0', 'thlev_zsea_theta', 'latitude_t', 'longitude_t'))
        
        print '[' + dt.now().strftime("%H:%M:%S") + '] Populating the dimension variables'
        # populate the dimension variables
        times_var[:] = w_nc.variables['min10_0'][:]
        z_var[:]     = w_nc.variables['thlev_zsea_theta'][:]
        lats_var[:]  = w_nc.variables['latitude_t'][:]
        lons_var[:]  = w_nc.variables['longitude_t'][:]
        
        print '[' + dt.now().strftime("%H:%M:%S") + '] Populating the wind variables'
        print '[' + dt.now().strftime("%H:%M:%S") + '] Regridding u-winds'
        u_var[:] = regrid(w_nc, u_nc, u_key)[:]
        
        print '[' + dt.now().strftime("%H:%M:%S") + '] Regridding v-winds'
        v_var[:] = regrid(w_nc, v_nc, v_key)[:]
        
        print '[' + dt.now().strftime("%H:%M:%S") + '] Copying over w-winds'
        w_var[:] = w_nc.variables[w_key][:]
        
        print '[' + dt.now().strftime("%H:%M:%S") + '] Closing the new winds_nc'
        winds_nc.close()
        
        print '[' + dt.now().strftime("%H:%M:%S") + '] Closing the old winds ncs'
        u_nc.close()
        v_nc.close()
        w_nc.close()
        
        print '[' + dt.now().strftime("%H:%M:%S") + '] Regridding hour ' + hour + ' complete.'

main()
