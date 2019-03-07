"""
Analysis of a 10-day spin up simulation to arrive at balanced winds and surface
fluxes for our cloud trail simulations.

This script replaces 'temperature_profiles.py' and 'surface_fluxes.py' which are
now depreciated.
"""
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import integrate, interpolate
from analysis_tools import RDP, find_h, lcl, get_CC, get_CTZ, regrid

### Copy over from 'temperature_profiles.py' ###

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

### Copy over from 'hodographs.py' ###

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
    plt.savefig('../adjustedProf.png', dpi = 150, bbox_inches = 'tight')
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

### Copy over from 'surface_fluxes.py' ###

days = ["{0:02d}".format(day) for day in xrange(1, 11)]
ndays = len(days)
lhf_key = u'STASH_m01s03i234'
shf_key = u'STASH_m01s03i217'
E = np.zeros(1)
H = np.zeros(1)
for day in days:
    print 'STARTING: day ' + day
    fluxes = Dataset('fluxes_' + day + '.nc', 'r')
    E = np.concatenate((E, np.mean(fluxes.variables[lhf_key][:], axis = (1, 2))))
    H = np.concatenate((H, np.mean(fluxes.variables[shf_key][:], axis = (1, 2))))
    if day == '01':
        dt_i = fluxes.variables[lhf_key].shape[0]
    fluxes.close()

times = np.arange(1., 24*ndays*60., 60.)/1440.
ax1 = plt.subplot(211)
ax1.plot(times, H[1:], lw = 2, color = 'red')
plt.ylabel('Sensible Heat Flux (Wm$^{-2}$)')

ax2 = plt.subplot(212, sharex = ax1)
ax2.plot(times, E[1:], lw = 2, color = 'blue')
plt.ylabel('Latent Heat Flux (Wm$^{-2}$)')
plt.xlabel('Time (days)')

plt.savefig('Flux_equilibrium.png', dpi = 100)
plt.show()

with open('balanced_settings.txt', 'a') as my_file:
    my_file.write('\n')
    my_file.write('Initial Conditions and Forcing for a Balanced Simulation.\n')
    my_file.write('At least the last three days of a 10-day simulation were in approximate balance between the prescribed forcing and the environmental response.\n')
    my_file.write('\n')
    my_file.write('Surface Sensible Heat Flux = ' + str(round(np.mean(H[-4*dt_i:]), 4)) + ' W/m2\n')
    my_file.write('Surface Latent Heat Flux = ' + str(round(np.mean(E[-4*dt_i:]), 4)) + ' W/m2\n')

