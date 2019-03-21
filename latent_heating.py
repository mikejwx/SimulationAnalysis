import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate,integrate
from netCDF4 import Dataset
from SkewT_archer import Lv

hour = '09'
print 'Defining keys for variables and dimensions required.'
# Need to use the boundary layer heat flux (W/m2) variable where the z = 0 data should == the surface flux but at higher temporal resolution
shf_key = u'STASH_m01s03i216' # need for the surface sensible heat flux
mcl_key = u'STASH_m01s00i392' # need for the total cloud liquid water
mr_key = u'STASH_m01s00i394' # need for the total suspended rain liquid
rho_key = u'STASH_m01s00i389' # need for the vertical integral
zth_key = 'thlev_zsea_theta' # theta levels
zrh_key = 'rholev_zsea_rho' # rho levels. ugh.

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

    print 'Reading the relevant dimensions.'
    # Read the height dimensions
    z_the   = mr_nc.variables[zth_key][:]*1.
    z_rho   = fluxes_nc.variables[zrh_key][:]*1.
    # Read in the times
    times = mr_nc.variables['min10_0'][:]

    # Initialise some arrays for the terms of our partial budget
    print 'Initialising arrays to store our lagrangian heating.'
    SHF = np.zeros_like(times)
    CON = np.zeros_like(times)
    #PRE = np.zeros_like(times) # precipitation is currently slightly different so we will come back to this later.

    print 'Defining the wind with which to advect the parcel.'
    # Advect our parcel with the boundary layer winds lowest ~ 500 m
    iz500 = np.where(np.abs(z_the - 500.) == np.min(np.abs(z_the - 500.)))[0][0]+1

    print 'Defining the parcel area.'
    # Compute the distance from the parcel
    R = np.sqrt((X_per-x_c)**2. + (Y_per-y_c)**2.)
    # Only consider grid cells that are less than 2 km away from the centrepoint, using the periodic coordinate system
    mask = np.where((R < 2000.), 1., 0.)
    mask = np.max([mask[(0+j*319):(319+j*319),(0+i*1160):(1160+i*1160)] for j in xrange(3) for i in xrange(3)], axis = 0)

    print 'Starting the advection routine.'
    # Initialise the parcel coordinates at the current time, i.e. upwind of the island centrepoint
    it = 1
    x_p = (x_c - np.nanmean(np.where(mask, integrate.trapz(y = wind_nc.variables[u_key][it,1:(iz500+1),:,:], x = z_the[:iz500], axis = 0)/np.max(z_the[:iz500]), np.nan))*(times[it] - times[it-1])*60.)%116000. 
    # N.B. wind nc starts at z = 0, while other nc on theta levels start at z = 2 m
    y_p = (y_c - np.nanmean(np.where(mask, integrate.trapz(y = wind_nc.variables[v_key][it,1:(iz500+1),:,:], x = z_the[:iz500], axis = 0)/np.max(z_the[:iz500]), np.nan))*(times[it] - times[it-1])*60.)%31900.

    # Start lists of each parcel centrepoint coordinate to track the parcel
    x_ps = [x_p]
    y_ps = [y_p]

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
            x_p = (x_p + np.nanmean(np.where(mask, integrate.trapz(y = wind_nc.variables[u_key][it,1:(iz500+1),:,:], x = z_the[:iz500], axis = 0)/np.max(z_the[:iz500]), np.nan))*(times[it] - times[it-1])*60.)%116000. 
            # N.B. wind nc starts at z = 0, while other nc on theta levels start at z = 2 m
            y_p = (y_p + np.nanmean(np.where(mask, integrate.trapz(y = wind_nc.variables[v_key][it,1:(iz500+1),:,:], x = z_the[:iz500], axis = 0)/np.max(z_the[:iz500]), np.nan))*(times[it] - times[it-1])*60.)%31900.
            x_ps.append(x_p)
            y_ps.append(y_p)
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
        if mr_key in mr_nc.variables.keys():
            q_con += mr_nc.variables[mr_key][it,:,:,:]
            mr_flag = False
        else:
            mr_flag = True
        CON[it]  = np.nanmean(np.where(mask, integrate.trapz(y = Lv*rho_data*q_con, x = z_the, axis = 0), np.nan))

    # Close our netCDFs
    mr_nc.close()
    fluxes_nc.close()
    wind_nc.close()

    CON_change = np.concatenate((np.array([CON[0]]), CON[1:] - CON[:-1]), axis = 0) # don't neet to multiply by Dtime because we haven't divided by it here
    # instead of plotting against time, plot against distance from origin
    dist = np.array([np.sqrt((x_ps[i] - x_ps[0])**2. + (y_ps[i] - y_ps[0])**2.) for i in xrange(len(x_ps))])/1000.
    ### Make the plot ###
    print 'Making the plot.'
    # Top Panel : the time series following the flow of the heating 
    fig = plt.figure(tight_layout = True, figsize = (12,6))
    ax = fig.add_subplot(1, 2, 1)
    ax.plot(dist, np.cumsum(SHF)*(times[1] - times[0])*60., '--b*', label = '$\int_{}^{} \\rho c_{pd} (\overline{w^{\prime} \\theta^{\prime}})_{sfc} Dt$')
    ax.plot(dist, np.cumsum(CON_change), '-b*', label = '$\int_{}^{} L_{v} \int_{}^{} \\rho \Delta q_{con} dz Dt$')
    #ax.plot(times/60., np.cumsum(SHF), 'b--', label = '$\\rho c_{pd} (\overline{w^{\prime} \\theta^{\prime}})_{sfc}$')
    #ax.plot(times/60., np.cumsum(CON), 'b', label = '$L_{v} \int \\rho q_{con} dz$')
    plt.legend(loc = 0)
    #ax.set_xlabel('Time (hrs)')
    ax.set_xlabel('Distance from start point (km)')
    ax.set_ylabel('Flux (J/m$^{2}$)')
    if mr_flag:
        ax.set_title('Lagrangian time integral (not including rain)')
    else:
        ax.set_title('Lagrangian time integral')

    # Bottom Panel: Where the parcel has gone as it follows the flow
    ax = fig.add_subplot(1, 2, 2, adjustable = 'box', aspect = 1)
    R = np.sqrt((X-x_c)**2. + (Y-y_c)**2.)
    # plot the island location
    ax.contourf(X/1000., Y/1000., np.where((R < 1000*np.sqrt(50./np.pi)), 1., 0), levels = [0., 1e-16, 1.], colors = ['w', 'k'])
    # plot the footprints of the parcel
    for i in xrange(len(x_ps)):
        R = np.sqrt((X_per - x_ps[i])**2. + (Y_per - y_ps[i])**2.)
        mask = np.where((R < 2000.), 1., 0.)
        mask = np.max([mask[(0+j*319):(319+j*319),(0+i*1160):(1160+i*1160)] for j in xrange(3) for i in xrange(3)], axis = 0)
        ax.contourf(X/1000., Y/1000., mask, levels = [1e-16, 1.], colors = ['r'], alpha = 0.5)

    # plot the track of the parcel
    ax.plot(np.array(x_ps)/1000., np.array(y_ps)/1000., '--r*')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_title('This is for the ' + ID + ' experiment')
    plt.savefig('../lagrangian_latent_heat_fluxes_' + ID + '.png', dpi = 150)
    plt.close('all')
    #plt.show()


