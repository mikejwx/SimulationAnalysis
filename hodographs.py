import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import interpolate

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

def toComponents(speed, direction):
    "takes the wind speed and wind direction and breaks the wind speed into"
    "u- and v- components to be returned."
    "speed = wind speed"
    "direction = wind direction"
    u = - speed*np.sin(direction*np.pi/180)
    v = - speed*np.cos(direction*np.pi/180)
    
    return u, v

def fromComponents(u, v, isList = False):
    "Takes the two horizontal wind components and combines them into a wind"
    "speed and direction to be returned."
    "u = the zonal wind component"
    "v = the meridional wind component"
    speed = np.sqrt(np.array(u)**2 + np.array(v)**2)
    if isList:
        direction = np.arcsin(v/speed)*180/np.pi
        for iV in range(len(u)):
            if speed[iV] > 0:
                if (u[iV] >= 0) and (v[iV] > 0):
                    #if u > 0 and v > 0 SW quadrant
                    direction[iV] = 270.0 - direction[iV]
                elif (u[iV] > 0) and (v[iV] <= 0):
                    #if u > 0 and v < 0 NW quadrant
                    direction[iV] = 270.0 - direction[iV]
                elif (u[iV] <= 0) and (v[iV] > 0):
                    #if u < 0 and v > 0 SE quadrant
                    direction[iV] = 90.0 + direction[iV]
                elif (u[iV] < 0) and (v[iV] <= 0):
                    #if u < 0 and v < 0 NE quadrant
                    direction[iV] = 90.0 + direction[iV]
            else:
                direction[iV] = 0.
    else:
        direction = np.arcsin(v/speed)*180/np.pi
        if speed > 0:
            if (u >= 0) and (v > 0):
                #if u > 0 and v > 0 SW quadrant
                direction = 270.0 - direction
            elif (u > 0) and (v <= 0):
                #if u > 0 and v < 0 NW quadrant
                direction = 270.0 - direction
            elif (u <= 0) and (v > 0):
                #if u < 0 and v > 0 SE quadrant
                direction = 90.0 + direction
            elif (u < 0) and (v <= 0):
                #if u < 0 and v < 0 NE quadrant
                direction = 90.0 + direction
        else:
            direction = 0.
    return speed, direction

def main():
    """
    Read in all of the wind data and compute horizontally averaged wind
    profiles. The produce several plots.
    """

    days = ["{0:02d}".format(x) for x in xrange(0,1,3)]#24,3)]
    ndays = len(days)
    # on the first day
    print 'Working on day ' + days[0]
    u = Dataset('u_' + days[0] + '.nc', 'r') # read in netcdf for u on day 01
    v = Dataset('v_' + days[0] + '.nc', 'r') # read in netcdf for v on day 01
    z = u.variables['rholev_zsea_rho'][:]    # read in the heights from u netcdf
    target_z = 200                           # level at which we want to analyse the winds
    iz = np.where(abs(z - target_z) == np.min(abs(z - target_z)))[0][0] # index for that level
    dt_i = u.variables['min15'].shape[0]     # length of the time dimension
    zi = np.nanmean(Dataset('zi_' + days[0] + '.nc', 'r').variables['boundary layer depth'][:], axis = (1,2))
    u_mean = np.zeros((u.variables['STASH_m01s00i002'].shape[0]*ndays, u.variables['STASH_m01s00i002'].shape[1])) # create an array for the uwinds = u_mean[t,z]
    v_mean = np.zeros_like(u_mean) # create counterpart array for the vwinds
    times = np.zeros(u.variables['min15'].shape[0]*ndays) # create a time dimension

    # populate our uwinds, vwinds, and time dimension with data from day 01
    u_mean[:dt_i,:] = horizontalMean(u.variables['STASH_m01s00i002'][:])
    v_mean[:dt_i,:] = horizontalMean(v.variables['STASH_m01s00i003'][:])
    times[:dt_i] = u.variables['min15'][:]

    # close the netcdfs
    u.close()
    v.close()

    for it in xrange(1, ndays):
        # day[0] = day 01 is already done, start from day[1] = day 02
        print 'Working on day ' + days[it]
        u = Dataset('u_' + days[it] + '.nc', 'r')
        v = Dataset('v_' + days[it] + '.nc', 'r')
        zi = np.concatenate((zi, np.nanmean(Dataset('zi_' + days[it] + '.nc', 'r').variables['boundary layer depth'][:], axis = (1,2))), axis = 0)
        u_mean[(dt_i*it):(dt_i*(it+1)),:] = horizontalMean(u.variables['STASH_m01s00i002'][:])
        v_mean[(dt_i*it):(dt_i*(it+1)),:] = horizontalMean(v.variables['STASH_m01s00i003'][:])
        times[(dt_i*it):(dt_i*(it+1))] = u.variables['min15'][:]
        
        u.close()
        v.close()
        
    print 'making the plot'
    fig = plt.figure()
    fig.set_size_inches(12, 10)
    plt.subplot(221)
    #plt.scatter(u1, v1, c = np.arange(len(u1)), cmap = 'hot')
    #plt.colorbar()
    plt.plot(u_mean[:,iz], v_mean[:,iz], 'b')
    plt.plot(u_mean[0,iz], v_mean[0,iz], 'b', marker = '^', markersize = 10)
    plt.xlabel('u-wind (m/s)')
    plt.ylabel('v-wind (m/s)')
    plt.title('Hodograph at z='+str(round(z[iz], 1))+' m')

    plt.subplot(222)
    times = times/1440.
    plt.plot(times, u_mean[:,iz], label = 'u-wind')
    plt.plot(times, v_mean[:,iz], label = 'v-wind')
    plt.legend()
    plt.xlabel('Time (days)')
    plt.ylabel('Wind Speed (m/s)')
    plt.title('Timeseries at z = '+ str(round(z[iz], 1))+' m\nMean over last 4 days (u = ' + str(round(np.mean(u_mean[(-4*dt_i):,iz]),2)) + ', v = ' + str(round(np.mean(v_mean[(-4*dt_i):,iz]),2)) + ')')

    ax3= plt.subplot(223)
    # dt_i is the length of a day
    ax3.plot(np.mean(u_mean[-4*dt_i:,:], axis = 0), z/1000., 'k', lw = 2, label = 'u-wind')
    ax3.plot(np.mean(v_mean[-4*dt_i:,:], axis = 0), z/1000., 'k--', lw = 2, label = 'v-wind')
    u_g = np.array([-10.0, -10.0, 0.0, 0.0])
    u_g_z = np.array([0, 9000., 15000., 40000.])/1000.
    v_g = np.zeros_like(u_g)
    ax3.plot(u_g, u_g_z, 'r', label = 'u$_g$-wind = u$_0$')
    ax3.plot(v_g, u_g_z, 'r', ls = '--', label = 'v$_g$-wind = v$_0$')
    ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.ylim([0, 10])
    plt.ylabel('Height (km)')
    plt.xlabel('wind speed (m/s)')
    plt.title('Mean wind profile in lowest 10 km over last 4 days')
    
    plt.savefig('Inertial_'+str(int(round(z[iz],0)))+'.png', dpi = 100, bbox_inches = 'tight')
    plt.show()
    
    # Making the wind speed and wind direction time series
    wind_spd, wind_dir = fromComponents(u_mean[:,iz], v_mean[:,iz], isList = True)
    plt.subplot(211)
    plt.plot(times, wind_spd, 'k', lw = 2, label = 'Wind Speed')
    plt.xlabel('Time (days)')
    plt.ylabel('Wind Speed (m s$^{-1}$)')
    plt.legend(loc = 2)
    
    plt.subplot(212)
    plt.plot(times, wind_dir, 'k', lw = 2, label = 'Wind Direction')
    plt.xlabel('Time (days)')
    plt.ylabel('Wind Direction ($^{\circ}$)')
    plt.legend(loc = 2)
    plt.suptitle('Timeseries at z = '+ str(round(z[iz], 1))+' m')
    
    plt.savefig('wind_character.png', dpi = 150, bbox_inches = 'tight')
    plt.show()
    
    # write the mean wind speed over the last 4 days into our balanced 
    # conditions txt file
    ### Get the mean wind speed profiles for the last 4 days
    # Initialise the balanced wind arrays
    u_balanced = np.zeros_like(z)
    v_balanced = np.zeros_like(z)
    target_zi  = np.mean(zi[-4*dt_i:])
    count = 0
    for timestep in xrange(-4*dt_i, u_mean.shape[0]):
        u_balanced += interpolate.interp1d(z/zi[timestep], u_mean[timestep,:], fill_value = 'extrapolate')(z/target_zi)
        v_balanced += interpolate.interp1d(z/zi[timestep], v_mean[timestep,:], fill_value = 'extrapolate')(z/target_zi)
        count += 1.
    u_balanced /= count
    v_balanced /= count
    
    # Plot the mean wind speed for the last four days
    plt.plot(u_balanced, z/1000., 'k', lw = 2, label = 'uwind')
    plt.plot(v_balanced, z/1000., 'k--', lw = 2, label = 'vwind')
    plt.plot(u_g, u_g_z, 'r', label = 'u$_g$-wind = u$_0$')
    plt.plot(v_g, u_g_z, 'r', ls = '--', label = 'v$_g$-wind = v$_0$')
    plt.ylim([0, 3])
    plt.title('Mean winds over the last four days')
    plt.legend()
    plt.savefig('originalProf.png', dpi = 150, bbox_inches = 'tight')
    plt.show()
    
    u_g_interpolated = interpolate.interp1d(u_g_z*1000., u_g, fill_value = 'extrapolate')(z)
    v_g_interpolated = np.zeros_like(z)
    z_balanced = z.copy()
    
    u_balanced_new = np.zeros_like(z)
    v_balanced_new = np.zeros_like(z)
    
    # only change below the boundary layer
    for k in xrange(len(z)):
        factor = z[k]/target_zi
        if factor < 1.0:
            u_balanced_new[k] = u_balanced[k]
            v_balanced_new[k] = v_balanced[k]
        elif factor < 2.0:
            factor = factor - 1.
            u_balanced_new[k] = factor*u_g_interpolated[k] + (1. - factor)*u_balanced[k]
            v_balanced_new[k] = factor*v_g_interpolated[k] + (1. - factor)*v_balanced[k] 
        else:
            u_balanced_new[k] = u_g_interpolated[k]
            v_balanced_new[k] = v_g_interpolated[k]
    
    # Minimize required points to reproduce the profiles
    z_u, u_u = RDP(z_balanced, u_balanced_new, 0.001)
    v_u = interpolate.interp1d(z_balanced, v_balanced_new, fill_value = 'extrapolate')(z_u)
    
    plt.plot(u_u, z_u/1000., 'k', lw = 2, label = 'uwind')
    plt.plot(v_u, z_u/1000., 'k--', lw = 2, label = 'vwind')
    plt.plot(u_g, u_g_z, 'r', label = 'u$_g$-wind = u$_0$')
    plt.plot(v_g, u_g_z, 'r', ls = '--', label = 'v$_g$-wind = v$_0$')
    plt.ylim([0, 3])
    plt.legend()
    plt.title('Profile now tends to the geostrophic wind above boundary layer')
    plt.savefig('adjustedProf.png', dpi = 150, bbox_inches = 'tight')
    plt.show()
    
    n_uvlev = len(z_u)
    if n_uvlev < 100:
        print 'n = ' + str(n_uvlev) + '\nYAY!!!!!'
        with open('balanced_settings.txt', 'a') as my_new_file:
            my_new_file.write('The following are the wind profiles averaged over \n the last day of simulation, when it is roughly in equilibrium.\n')
            my_new_file.write('Number of uv levels:\n')
            my_new_file.write(str(n_uvlev) + '\n')
            my_new_file.write('uv levels (m):\n')
            for k in xrange(n_uvlev-1):
                my_new_file.write(str(z_u[k])+',')
            my_new_file.write(str(z_u[k])+'\n')
            my_new_file.write('\nu init wind (m/s):\n')
            for k in xrange(n_uvlev-1):
                my_new_file.write(str(round(u_u[k],2))+',')
            my_new_file.write(str(round(u_u[-1],2))+'\n')
            my_new_file.write('v init wind (m/s):\n')
            for k in xrange(n_uvlev-1):
                my_new_file.write(str(round(v_u[k],2))+',')
            my_new_file.write(str(round(v_u[-1],2)))
    else:
        print 'BOOO!!!!'
   
main()


