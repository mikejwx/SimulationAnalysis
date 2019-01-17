 # script to plot skew-T from the simulations
from netCDF4 import Dataset
from scipy import interpolate
execfile('SkewT.py')
def horizontalMean(data):
    """
    Function to compute and return the horizontal mean of a given array.
    Will check if the array is 3- or 4- dimensional. If 3D, then assume that 
    array = f(z, y, x), and if 4D, assume array = f(t, z, y, x).
    """
    
    if len(data.shape) == 3:
        h_mean = np.zeros(data.shape[0])
        for k in xrange(data.shape[0]):
            h_mean[k] = np.mean(data[k,:,:])
    elif len(data.shape) == 4:
        h_mean = np.zeros((data.shape[0], data.shape[1]))
        for t in xrange(data.shape[0]):
            for k in xrange(data.shape[1]):
                h_mean[t, k] = np.mean(data[t, k, :, :])
    else:
        raise Exception("Unrecognised array shape in horizontalMean")
    
    return h_mean

def RDP(x, y, e):
    """
    Algorithm to choose the minimum number of points that retains the shape of 
    y with maximum error, "error". Returns the values of y, and the coordinates 
    x for those values.
    """
    # find the end of the data
    end = len(x)
    index0 = 0
    index1 = end -1
    # The first point of our data is the same as the first point in our simplification
    x_simplified = [x[index0]]
    y_simplified = [y[index0]]
    #plt.plot(x, y, 'k-', marker = 'o', lw = 2)
    #plt.show()
    while (index0 < (end - 1)):
        #first estimate is a straight line between the first and last point of the data
        y_est = (y[index1] - y[index0])*(x - x[index0])/(x[index1] - x[index0]) + y[index0]
        d = np.abs(y - y_est)[index0:index1]
        dmax = np.max(d)
        while (dmax > e)*(index0 < index1):
            #plt.plot(x, y, 'k-', marker = 'o', lw = 2)
            #plt.plot(x[index0:(index1+1)], y_est[index0:(index1+1)], 'b--')
            index1 = index0 + np.where(d == dmax)[0][0] # find the index of distance_max > error
            #plt.plot(x[index1], y[index1], 'ro')
            #plt.show()
            y_est = (y[index1] - y[index0])*(x - x[index0])/(x[index1] - x[index0]) + y[index0] # draw straight line to that point
            d = np.abs(y[index0:index1] - y_est[index0:index1]) # recalculate the differences
            dmax = np.max(d) # find the maximum differences
        x_simplified.append(x[index1])
        y_simplified.append(y[index1])
        index0 = index1
        index1 = end - 1
    return np.array(x_simplified), np.array(y_simplified)

def regrid(target_grid, current_grid, current_data_key):
    """
    Function to attempt to regrid an already read-in netCDF B (current_grid) to 
    the same grid as already read-in netCDF A (target_grid).
    Should only work if there is some overlap between the range of the functions
    otherwise the interpolation won't make sense.
    
    Regrids in all four dimensions, although going from e.g. 1 minute data to
    10 or 15 minute data shouldn't really require regridding. But going from
    10 minute to 1 minute data definitely requires regridding and the linear 
    interpolation approach used here is probably not appropriate.
    """
    # Read the netCDF that contains the target grid
    target_x_key     = [key for key in target_grid.variables.keys() if ('longitude' in key) and ('_l' not in key)][0]
    target_x = target_grid.variables[target_x_key][:]
    target_y_key     = [key for key in target_grid.variables.keys() if ('latitude' in key) and ('_l' not in key)][0]
    target_y = target_grid.variables[target_y_key][:]
    target_z_key     = [key for key in target_grid.variables.keys() if ('zsea' in key) and ('bounds' not in key)][0]
    target_z = target_grid.variables[target_z_key][:]
    target_t_key     = [key for key in target_grid.variables.keys() if 'min' in key][0]
    target_t = target_grid.variables[target_t_key][:]
    
    # Read the netCDF that needs to be transformed
    current_x_key     = [key for key in current_grid.variables.keys() if ('longitude' in key) and ('_l' not in key)][0]
    current_x = current_grid.variables[current_x_key][:]
    current_y_key     = [key for key in current_grid.variables.keys() if ('latitude' in key) and ('_l' not in key)][0]
    current_y = current_grid.variables[current_y_key][:]
    current_z_key     = [key for key in current_grid.variables.keys() if ('zsea' in key) and ('bounds' not in key)][0]
    current_z = current_grid.variables[current_z_key][:]
    current_t_key     = [key for key in current_grid.variables.keys() if 'min' in key][0]
    current_t = current_grid.variables[current_t_key][:]
    current_data = current_grid.variables[current_data_key][:]
     
    ### Procedure ###
    # Will usually only need to interpolate in z, sometimes in x or y, probably
    # never required to interpolate in more than two dimensions
    # 1. interpolate in time
    if target_t_key != current_t_key:
        print 'Regridding in time...'
        # Interpolate in time
        current_data = interpolate.interp1d(current_t, current_data, axis = 0, fill_value = 'extrapolate')(target_t)
        print 'Complete.'
    # 2. interpolate in z
    if target_z_key != current_z_key:
        print 'Regridding in height...'
        # Interpolate in z
        current_data = interpolate.interp1d(current_z, current_data, axis = 1, fill_value = 'extrapolate')(target_z)
        print 'Complete.'
    # 3. interpolate in y
    if target_y_key != current_y_key:
        print 'Regridding in latitude...'
        # Interpolate in y
        for i in xrange(current_data.shape[3]):
            current_data[:,:,:,i] = interpolate.interp1d(current_y[:,i], current_data[:,:,:,i], axis = 2, fill_value = 'extrapolate')(target_y[:,i])
        print 'Complete.'
    # 4. interpolate in x
    if target_x_key != current_x_key:
        print 'Regridding in longitude...'
        # Interpolate in x
        for j in xrange(current_data.shape[2]):
            current_data[:,:,j,:] = interpolate.interp1d(current_x[j,:], current_data[:,:,j,:], axis = 2, fill_value = 'extrapolate')(target_x[j,:])
        print 'Complete.'
    
    return current_data

#read the data
days = ["{0:02d}".format(d) for d in xrange(1, 11)]
Rd = 287.05
cp = 1005.
g = 9.81

theta_key = u'STASH_m01s00i004'
pressure_key = u'STASH_m01s00i407'
q_key = u'STASH_m01s00i010'
mv_key = u'STASH_m01s00i391'
u_key = u'STASH_m01s00i002'
v_key = u'STASH_m01s00i003'
for day in days:
    print 'starting day ' + day
    
    bouy = Dataset('bouy_'+day+'.nc', 'r')
    fluxes = Dataset('fluxes_'+day+'.nc', 'r')
    mr = Dataset('mr_'+day+'.nc', 'r')
    u = Dataset('u_'+day+'.nc', 'r')
    v = Dataset('v_'+day+'.nc', 'r')
    
    theta = bouy.variables[theta_key][:]
    #regrid the pressures and specific humidity
    pressure_regrid = regrid(bouy, fluxes, pressure_key)
    q_regrid = regrid(bouy, mr, q_key)
    mv_regrid = regrid(bouy, mr, mv_key)
    u_regrid = regrid(bouy, u, u_key)
    v_regrid = regrid(bouy, v, v_key)
    
    if day == '01':
        temperature = np.zeros((theta.shape[0]*len(days), theta.shape[1]))
        pressure_rg = np.zeros((theta.shape[0]*len(days), theta.shape[1]))
        dewpoint_rg = np.zeros((theta.shape[0]*len(days), theta.shape[1]))
        u_rg        = np.zeros((theta.shape[0]*len(days), theta.shape[1]))
        v_rg        = np.zeros((theta.shape[0]*len(days), theta.shape[1]))
        q_rg        = np.zeros((theta.shape[0]*len(days), theta.shape[1]))
        mv_rg       = np.zeros((theta.shape[0]*len(days), theta.shape[1]))
        rh_rg       = np.zeros((theta.shape[0]*len(days), theta.shape[1]))
        theta_rg    = np.zeros((theta.shape[0]*len(days), theta.shape[1]))
        z           = bouy.variables['thlev_zsea_theta'][:]
    
    bouy.close()
    fluxes.close()
    mr.close()
    u.close()
    v.close()
    
    #convert from potential temperature to temperature
    temp = theta/((100000./pressure_regrid)**(Rd/cp))
    #horizontally average the temperatures
    temperature[theta.shape[0]*days.index(day):theta.shape[0]*(days.index(day)+1),:] = horizontalMean(temp)[:]
    pressure_rg[theta.shape[0]*days.index(day):theta.shape[0]*(days.index(day)+1),:] = horizontalMean(pressure_regrid)[:]
    dewpoint_rg[theta.shape[0]*days.index(day):theta.shape[0]*(days.index(day)+1),:] = horizontalMean(getDew(q_regrid, pressure_regrid/100.))[:]
    u_rg[theta.shape[0]*days.index(day):theta.shape[0]*(days.index(day)+1),:] = horizontalMean(u_regrid)[:]
    v_rg[theta.shape[0]*days.index(day):theta.shape[0]*(days.index(day)+1),:] = horizontalMean(v_regrid)[:]
    q_rg[theta.shape[0]*days.index(day):theta.shape[0]*(days.index(day)+1),:] = horizontalMean(q_regrid)[:]
    mv_rg[theta.shape[0]*days.index(day):theta.shape[0]*(days.index(day)+1),:] = horizontalMean(mv_regrid)[:]
    rh_rg[theta.shape[0]*days.index(day):theta.shape[0]*(days.index(day)+1),:] = horizontalMean(q_regrid/getQ(temp, 100., pressure_regrid))
    theta_rg[theta.shape[0]*days.index(day):theta.shape[0]*(days.index(day)+1),:] = horizontalMean(theta)[:]

# Time output every ten minutes
times = np.arange(1., 14400.*len(days), 10.)/60.

### Make skew-T log-p every three hours
#for it in xrange(0, temperature.shape[0], 18):
#    plotSkewT(temperature[it,:],dewpoint_rg[it,:], pressure_rg[it,:]/100., u = u_rg[it,:]/0.5144, v = v_rg[it,:]/0.5144, my_title = 'T+'+"{0:02d}".format(int(times[it]))+'hours')
#    plt.savefig('skewT_'+"{0:03d}".format(int(times[it]))+'.png', dpi = 150)
#    plt.close('all')

dt_i = theta.shape[0]
with open('balanced_settings.txt', 'a') as my_file:
    my_file.write('The following are the thermodynamic profiles averaged over \n the last four days of simulation, when it is roughly in equilibrium.\n')
    my_file.write('Potential Temperature (Theta)\n')
    # Find minimum number of required levels to reproduce theta profile
    theta_levels, theta_init = RDP(z, np.mean(theta_rg[-4*dt_i:,:], axis = 0), 0.1)
    n_thlev = len(theta_levels)
    my_file.write('Number of Theta levels:\n')
    my_file.write(str(n_thlev) + '\n')
    my_file.write('Theta Levels (m):\n')
    for k in xrange(n_thlev-1):
        my_file.write(str(theta_levels[k]) + ',')
    my_file.write(str(theta_levels[-1]) + '\n')
    my_file.write('Theta (K):\n')
    for k in xrange(n_thlev-1):
        my_file.write(str(theta_init[k]) + ',')
    my_file.write(str(theta_init[-1]) + '\n')
    
    # Find the minimum number of required levels to reproduce the mv profile
    mv_levels, mv_init = RDP(z, np.mean(rh_rg[-4*dt_i:, :], axis = 0), 0.01)
    n_mvlev = len(mv_levels)
    my_file.write('Relative Humidity (RH)\n')
    my_file.write('Number of mv levels:\n')
    my_file.write(str(n_mvlev) + '\n')
    my_file.write('mv Levels (m):\n')
    for k in xrange(n_mvlev-1):
        my_file.write(str(mv_levels[k]) + ',')
    my_file.write(str(mv_levels[-1]) + '\n')
    my_file.write('mv (RH decimal):\n')
    for k in xrange(n_mvlev-1):
        my_file.write(str(mv_init[k]) + ',')
    my_file.write(str(mv_init[-1]) + '\n')
plt.plot(np.array(mv_init)*100., np.array(mv_levels)/1000.)
plt.xlabel('RH (%)')
plt.ylabel('Height (km)')
plt.show()

