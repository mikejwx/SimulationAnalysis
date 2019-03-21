import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate,integrate
from netCDF4 import Dataset
from SkewT_archer import Lv

hour = '09'
print 'Defining keys for variables and dimensions required.'
# Need to use the boundary layer heat flux (W/m2) variable where the z = 0 data should == the surface flux but at higher temporal resolution
shf_key  = u'STASH_m01s03i216' # need for the surface sensible heat flux
mcl_key  = u'STASH_m01s00i392' # need for the total cloud liquid water
mr_key   = u'STASH_m01s00i394' # need for the total suspended rain liquid
rain_key = u'STASH_m01s04i201' # need this for the surface rainfall accumulation
rho_key  = u'STASH_m01s00i389' # need for the vertical integral
zth_key  = 'thlev_zsea_theta' # theta levels
zrh_key  = 'rholev_zsea_rho' # rho levels. ugh.

u_key = u'STASH_m01s00i002'
v_key = u'STASH_m01s00i003'

print 'Defining a coordinate system and taking advantage of periodic boundary conditions.'
# Define a coordinate system
x = np.arange(0., 116000., 100.)
y = np.arange(0., 31900., 100.)
X, Y = np.meshgrid(x, y)

# define a periodic coordinate system
X_per = np.concatenate((X-116000., X, X+116000.), axis = 1)
X_per = np.concatenate((X_per, X_per, X_per), axis = 0)
Y_per = np.concatenate((Y, Y, Y), axis = 1)
Y_per = np.concatenate((Y_per - 31900., Y_per, Y_per +  31900.), axis = 0)

# Pick a starting point over the island and advect it forward
# e.g. pick the island centrepoint
# Define the coordinates of the island centre point
print 'Defining the starting point.'
x_c, y_c = [108000., 15950.]

IDs = ['Control', 'H125', 'H375']
#IDs = ['H125', 'H375']
for ID in IDs:
    base_path = '/nerc/n02/n02/xb899100/CloudTrail/' + ID + '/'
    print 'base_path = ' + base_path

    """
    We want to compare the total boundary layer moisture fluxes with the ls cloud moisture fluxes
    -> in qinc there is a variable that contains the bl and ls cld
        This is the tendency due to the boundary layer scheme and the large scale cloud scheme kg/kg/timestep
    -> in fluxes there is a variable that contains bl
        This is the boundary layer moisture fluxes kg/m2/s
    To get the ls cloud, simply subtract the bl from bl+ls cld, we want to consider vertically integrated quantities, and follow the flow (Lagrangian perspective)
    """

    ### Start the lagrangian routine ###
    # Open the netCDF
    print 'Opening the relevant netCDF'
    mr_nc     = Dataset(base_path + 'mr_' + hour + '.nc', 'r')
    fluxes_nc = Dataset(base_path + 'fluxes_' + hour + '.nc', 'r')
    wind_nc   = Dataset(base_path + 'wind_' + hour + '.nc', 'r')
    lwp_nc    = Dataset(base_path + 'lwp_00.nc', 'r')
    
    print 'Reading the relevant dimensions.'
    # Read the height dimensions
    z_the   = mr_nc.variables[zth_key][:]*1.
    z_rho   = fluxes_nc.variables[zrh_key][:]*1.
    # Read in the times
    times = mr_nc.variables['min10_0'][:]
    train_key = 'accum15'
    if train_key in lwp_nc.variables.keys():
        times_rain = lwp_nc.variables[train_key][:]*1.
    
    # Initialise some arrays for the terms of our partial budget
    print 'Initialising arrays to store our lagrangian heating.'
    SHF = np.zeros_like(times)
    CON = np.zeros_like(times)
    PRE = np.zeros(1) # need to start an array with one zero - the first time will always be zero and then concatenate on subsequent times

    print 'Defining the wind with which to advect the parcel.'
    # Advect our parcel with the boundary layer winds lowest ~ 500 m
    iz500 = np.where(np.abs(z_the - 500.) == np.min(np.abs(z_the - 500.)))[0][0]+1

    print 'Defining the parcel area.'
    # Compute the distance from the parcel
    R = np.sqrt((X_per-x_c)**2. + (Y_per-y_c)**2.)
    # Only consider grid cells that are less than 2 km away from the centrepoint, using the periodic coordinate system
    mask = np.where((R < 2000.), 1., 0.)
    mask = np.max([mask[(0+j*319):(319+j*319),(0+i*1160):(1160+i*1160)] for j in xrange(3) for i in xrange(3)], axis = 0)
    
    # Define the land-sea mask
    R = np.sqrt((X-x_c)**2. + (Y-y_c)**2.)
    lsm  = np.where((R < 1000*np.sqrt(50./np.pi)), 1., 0)
    
    print 'Starting the advection routine.'
    # Initialise the parcel coordinates at the current time, i.e. upwind of the island centrepoint
    x_p = (x_c - np.nanmean(np.where(mask, integrate.trapz(y = wind_nc.variables[u_key][0,1:(iz500+1),:,:], x = z_the[:iz500], axis = 0)/np.max(z_the[:iz500]), np.nan))*(times[1] - times[0])*60.)%116000. 
    # N.B. wind nc starts at z = 0, while other nc on theta levels start at z = 2 m
    y_p = (y_c - np.nanmean(np.where(mask, integrate.trapz(y = wind_nc.variables[v_key][0,1:(iz500+1),:,:], x = z_the[:iz500], axis = 0)/np.max(z_the[:iz500]), np.nan))*(times[1] - times[0])*60.)%31900.

    # Start lists of each parcel centrepoint coordinate to track the parcel, for the condensate track
    x_pcon = [x_p]
    y_pcon = [y_p]
    # And for the rain swaths
    x_prain = [x_p]
    y_prain = [y_p]
    
    ############################################################################
    # Do the lagrangian accumulation for the condensate
    ############################################################################
    print 'Defining the parcel area.'
    # Compute the distance from the parcel
    R = np.sqrt((X_per-x_p)**2. + (Y_per-y_p)**2.)
    # Only consider grid cells that are less than 2 km away from the centrepoint, using the periodic coordinate system
    mask = np.where((R < 2000.), 1., 0.)
    mask = np.max([mask[(0+j*319):(319+j*319),(0+i*1160):(1160+i*1160)] for j in xrange(3) for i in xrange(3)], axis = 0)
    
    for it in xrange(len(times)):
        # Define our point at this time, starting from the island centre point
        if it != 0:
            print 'Calculating the new location of our parcel.'
            # if it == 0, then this is the starting point and we don't need to advect it along or redefine which grid points are in the parcel
            u_mean = np.nanmean(np.where(mask, integrate.trapz(y = wind_nc.variables[u_key][it,1:(iz500+1),:,:], x = z_the[:iz500], axis = 0)/np.max(z_the[:iz500]), np.nan))
            v_mean = np.nanmean(np.where(mask, integrate.trapz(y = wind_nc.variables[v_key][it,1:(iz500+1),:,:], x = z_the[:iz500], axis = 0)/np.max(z_the[:iz500]), np.nan))
            # N.B. wind nc starts at z = 0, while other nc on theta levels start at z = 2 m
            x_p = (x_pcon[-1] + u_mean*(times[it] - times[it-1])*60.)%116000. 
            y_p = (y_pcon[-1] + v_mean*(times[it] - times[it-1])*60.)%31900.
            x_pcon.append(x_p)
            y_pcon.append(y_p)
        print 'Our parcel is currently centred at:'
        print 'x_p = ' + str(int(x_p)) + ', y_p = ' + str(int(y_p))
        # Redefine our mask
        R = np.sqrt((X_per - x_p)**2. + (Y_per - y_p)**2.)
        mask = np.where((R < 2000.), 1., 0.)
        mask = np.max([mask[(0+j*319):(319+j*319),(0+i*1160):(1160+i*1160)] for j in xrange(3) for i in xrange(3)], axis = 0)
        
        print 'Computing the sensible heat flux and the latent heating in the cloud at time ' + str(int(times[it]))
        # Compute the surface sensible heat fluxes in that ring.
        SHF[it]  = np.nanmean(np.where(mask, fluxes_nc.variables[shf_key][it,0,:,:], np.nan))
        
        # For the condensate latent heating part, this requires a density weighting
        # Density needs to be interpolated onto theta levels
        rho_data = interpolate.interp1d(x = z_rho, y = fluxes_nc.variables[rho_key][it,:,:,:]*1., axis = 0, fill_value = 'extrapolate')(z_the)
        q_con    = mr_nc.variables[mcl_key][it,:,:,:]
        # If rain mixing ratio is in there, include it
        if mr_key in mr_nc.variables.keys():
            q_con += mr_nc.variables[mr_key][it,:,:,:]
            mr_flag = False
        else:
            mr_flag = True
        CON[it]  = np.nanmean(np.where(mask, integrate.trapz(y = Lv*rho_data*q_con, x = z_the, axis = 0), np.nan))
    
    ############################################################################
    # Do the lagrangian accumulation for the surface rainfall
    ############################################################################
    # If surface rainfall is in there, include it    
    if (rain_key in lwp_nc.variables.keys()):
        l_rain = True
        itr0 = np.where(np.abs(times_rain - times[0]) == np.min(np.abs(times_rain - times[0])))[0][0]
        times_r = np.array([540.])
        # If there is rainfall accumulation available
        for itr in xrange(itr0, times_rain.shape[0]):
            if (540. < times_rain[itr] <= 720.):
                times_r = np.concatenate((times_r, np.array([times_rain[itr]])), axis = 0)
                # If the time in accumulation is within the time range of the other netCDF
                # find the indexes in wind_nc that are nearest to the time in accum15
                it0 = np.where(np.abs(times - times_rain[itr]) == np.min(np.abs(times - times_rain[itr])))[0][0]
                it1 = it0 + 2
                print 'Calculating the new location of our parcel.'
                # if it == 0, then this is the starting point and we don't need to advect it along or redefine which grid points are in the parcel
                u_mean = np.nanmean(np.where(mask, integrate.trapz(y = np.nanmean(wind_nc.variables[u_key][it0:it1,1:(iz500+1),:,:], axis = 0), x = z_the[:iz500], axis = 0)/np.max(z_the[:iz500]), np.nan))
                v_mean = np.nanmean(np.where(mask, integrate.trapz(y = np.nanmean(wind_nc.variables[v_key][it0:it1,1:(iz500+1),:,:], axis = 0), x = z_the[:iz500], axis = 0)/np.max(z_the[:iz500]), np.nan))
                # N.B. wind nc starts at z = 0, while other nc on theta levels start at z = 2 m
                x_p = (x_prain[-1] + u_mean*(times_rain[itr] - times_rain[itr-1])*60.)%116000. 
                y_p = (y_prain[-1] + v_mean*(times_rain[itr] - times_rain[itr-1])*60.)%31900.
                x_prain.append(x_p)
                y_prain.append(y_p)
                # Define the swath between the previous two time steps [-1] and [-2]
                swath = np.zeros_like(X)
                for x_coord in np.arange(x_prain[-2], x_prain[-1], -250.):
                    x_pp = x_coord
                    if np.sign(y_prain[-1] - y_prain[-2]) == np.sign(v_mean):
                        y_pp = (y_prain[-1] - y_prain[-2])*(x_coord - x_prain[-2])/(x_prain[-1] - x_prain[-2]) + y_prain[-2]
                    else:
                        y_pp = ((y_prain[-1]-31900.) - y_prain[-2])*(x_coord - x_prain[-2])/(x_prain[-1] - x_prain[-2]) + y_prain[-2]
                    R = np.sqrt((X_per - x_pp)**2. + (Y_per - y_pp)**2.)
                    mask = np.where((R < 2000.), 1., 0.)
                    mask = np.max([mask[(0+j*319):(319+j*319),(0+i*1160):(1160+i*1160)] for j in xrange(3) for i in xrange(3)], axis = 0)
                    swath += mask
                swath = np.where((swath > 0), 1., 0.)
                PRE = np.concatenate((PRE, np.array([np.nanmean(Lv*np.where((swath == 1.), lwp_nc.variables[rain_key][itr,:,:] - lwp_nc.variables[rain_key][itr-1,:,:], np.nan))])), axis = 0) # do something with the rain
    else:
        l_rain = False
    
    # Close our netCDFs
    mr_nc.close()
    fluxes_nc.close()
    wind_nc.close()
    lwp_nc.close()
    
    CON_change = np.concatenate((np.array([CON[0]]), CON[1:] - CON[:-1]), axis = 0) # don't need to multiply by Dtime because we haven't divided by it here
    if l_rain:
        PRE_i = interpolate.interp1d(x = times_r, y = PRE, fill_value = 'extrapolate')(times)
        PRE_change = np.concatenate((np.array([PRE_i[0]]), PRE_i[1:] - PRE_i[:-1]), axis = 0)
        CON_change += PRE_change
    # instead of plotting against time, plot against distance from origin
    dist = np.array([np.sqrt((x_pcon[i] - x_pcon[0])**2. + (y_pcon[i] - y_pcon[0])**2.) for i in xrange(len(x_pcon))])/1000.
    dist_rain = np.array([np.sqrt((x_prain[i] - x_prain[0])**2. + (y_prain[i] - y_prain[0])**2.) for i in xrange(len(x_prain))])/1000.
    
    ### Make the plot ###
    print 'Making the plot.'
    # Top Panel : the time series following the flow of the heating 
    fig = plt.figure(tight_layout = True, figsize = (12,6))
    ax = fig.add_subplot(1, 2, 1)
    ax.plot(dist, np.cumsum(SHF)*(times[1] - times[0])*60./1000., '--b*', label = '$\int_{}^{} \\rho c_{pd} (\overline{w^{\prime} \\theta^{\prime}})_{sfc} Dt$')
    ax.plot(dist, np.cumsum(CON_change)/1000., '-b*', label = '$\int_{}^{} L_{v} \int_{}^{} \\rho \Delta q_{con} dz Dt$')
    #ax.plot(dist_rain, np.cumsum(PRE_change)/1000., '-r*', label = 'RAIN')
    #ax.plot(times/60., np.cumsum(SHF), 'b--', label = '$\\rho c_{pd} (\overline{w^{\prime} \\theta^{\prime}})_{sfc}$')
    #ax.plot(times/60., np.cumsum(CON), 'b', label = '$L_{v} \int \\rho q_{con} dz$')
    plt.legend(loc = 0)
    #ax.set_xlabel('Time (hrs)')
    ax.set_xlabel('Distance from start point (km)')
    ax.set_ylabel('Flux (kJ/m$^{2}$)')
    my_title = 'Lagrangian time integral'
    if mr_flag:
        my_title += '\n(not including rain mixing ratio)'
    if not l_rain:
        my_title += '\n (not including surface rain)'
    ax.set_title(my_title)

    # Bottom Panel: Where the parcel has gone as it follows the flow
    ax = fig.add_subplot(1, 2, 2, adjustable = 'box', aspect = 1)
    R = np.sqrt((X-x_c)**2. + (Y-y_c)**2.)
    # plot the land-sea mask
    ax.contourf(X/1000., Y/1000., lsm, levels = [0., 1e-16, 1.], colors = ['w', 'k'])
    # plot the footprints of the parcel
    for i in xrange(len(x_pcon)):
        R = np.sqrt((X_per - x_pcon[i])**2. + (Y_per - y_pcon[i])**2.)
        mask = np.where((R < 2000.), 1., 0.)
        mask = np.max([mask[(0+j*319):(319+j*319),(0+i*1160):(1160+i*1160)] for j in xrange(3) for i in xrange(3)], axis = 0)
        ax.contourf(X/1000., Y/1000., mask, levels = [1e-16, 1.], colors = ['r'], alpha = 0.5)

    # plot the track of the parcel
    ax.plot(np.array(x_pcon)/1000., np.array(y_pcon)/1000., '--r*')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_title('This is for the ' + ID + ' experiment')
    plt.savefig('../lagrangian_latent_heat_fluxes_' + ID + '.png', dpi = 150)
    plt.pause(10)
    plt.close('all')
    #plt.show()

