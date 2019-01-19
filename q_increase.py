### investigation into the increase in water vapour mixing ratio aloft
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import interpolate

# Read in the data
hours = ["{0:02d}".format(h) for h in xrange(0, 24, 3)]
nhours = len(hours)

# keys for the variables we're investigating
q_lsrain_key = u'STASH_m01s04i182'
q_blpcld_key = u'STASH_m01s09i182'
q_advect_key = u'STASH_m01s12i184'
q_diffus_key = u'STASH_m01s13i182'
q_qtbcld_key = u'STASH_m01s15i182'
q_totali_key = u'STASH_m01s30i182'
q_ideali_key = u'STASH_m01s53i182'

zi_key       = u'boundary layer depth'
theta_key    = u'STASH_m01s00i004'
pressure_key = u'STASH_m01s00i407'

for hour in hours:
    print 'Starting hour: ' + hour
    # import the q_incs
    q_inc   = Dataset('../qinc_' + hour + '.nc', 'r')
    bouy    = Dataset('../bouy_' + hour + '.nc', 'r')
    fluxes  = Dataset('../fluxes_' + hour + '.nc', 'r')
    zi_hour = Dataset('../zi_' + hour + '.nc', 'r')
    
    if hour == '00':
        print 'Initialising arrays'
        # initialize some arrays for our diagnostics
        z     = q_inc.variables['thlev_zsea_theta'][:]*1.
        z_rho = fluxes.variables['rholev_zsea_rho'][:]*1.
        t_len = q_inc.variables[q_lsrain_key].shape[0]*nhours
        times = np.arange(0., 1440.1, 10.)
        rain       = np.zeros((t_len, z.shape[0]))
        bl_cld     = np.zeros((t_len, z.shape[0]))
        advection  = np.zeros((t_len, z.shape[0]))
        diffusion  = np.zeros((t_len, z.shape[0]))
        qt_bal_cld = np.zeros((t_len, z.shape[0]))
        total      = np.zeros((t_len, z.shape[0]))
        idealised  = np.zeros((t_len, z.shape[0]))
        zi = np.zeros(t_len) # Boundary Layer Depth
        zf = np.zeros(t_len) # Freezing Level
        zt = np.zeros(t_len) # Tropopause Level
    
    print 'Store the horizontally averaged profiles'
    # store the horizontally averaged profiles
    rain[hours.index(hour)*18:(hours.index(hour)+1)*18, :] = np.mean(q_inc.variables[q_lsrain_key][:], axis = (2,3))
    
    bl_cld[hours.index(hour)*18:(hours.index(hour)+1)*18, :] = np.mean(q_inc.variables[q_blpcld_key][:], axis = (2,3))
    
    advection[hours.index(hour)*18:(hours.index(hour)+1)*18, :] = np.mean(q_inc.variables[q_advect_key][:], axis = (2,3))
    
    diffusion[hours.index(hour)*18:(hours.index(hour)+1)*18, :] = np.mean(q_inc.variables[q_diffus_key][:], axis = (2,3))
    
    qt_bal_cld[hours.index(hour)*18:(hours.index(hour)+1)*18, :] = np.mean(q_inc.variables[q_qtbcld_key][:], axis = (2,3))
    
    total[hours.index(hour)*18:(hours.index(hour)+1)*18, :] = np.mean(q_inc.variables[q_totali_key][:], axis = (2,3))
    
    idealised[hours.index(hour)*18:(hours.index(hour)+1)*18, :] = np.mean(q_inc.variables[q_ideali_key][:], axis = (2,3))
    
    zi[hours.index(hour)*18:(hours.index(hour)+1)*18] = np.mean(zi_day.variables[zi_key][:], axis = (1,2))
    
    temperature = bouy.variables[theta_key][:, 1:, :, :]/(100000./interpolate.interp1d(z_rho, fluxes.variables[pressure_key][:], axis = 1, fill_value = 'extrapolate')(z))**(287.05/1005.)
    freezing_level = np.zeros_like(temperature[:,0,:,:])
    for it in xrange(temperature.shape[0]):
        for iy in xrange(temperature.shape[2]):
            for ix in xrange(temperature.shape[3]):
                freezing_level[it, iy, ix] = interpolate.interp1d(temperature[it, :, iy, ix], z, fill_value = 'extrapolate')(273.15)
    
    zf[hours.index(hour)*18:(hours.index(hour)+1)*18] = np.mean(freezing_level, axis = (1,2))
    
    tropopause_level = np.zeros_like(temperature[:,0,:,:])
    for it in xrange(temperature.shape[0]):
        for iy in xrange(temperature.shape[2]):
            for ix in xrange(temperature.shape[3]):
                tropopause_level[it, iy, ix] = interpolate.interp1d(temperature[it,:,iy,ix], z, fill_value = 'extrapolate')(np.max(np.min(temperature[it,:,:,:], axis = 0)))

    zt[hours.index(hour)*18:(hours.index(hour)+1)*18] = np.mean(tropopause_level, axis = (1,2))
    
    print 'Closing the netCDF'
    q_inc.close()
    bouy.close()
    fluxes.close()
    zi_day.close()
    

def my_plot_style(variable, my_title, file_name, show = True):
    """
    Plot each thing with my plot style.
    -> Log y-axis with my labels
    -> blue-white-red color bar centred on zero
    -> time x-axis the bottom in days
    """
    fig = plt.figure()
    ax = plt.subplot()
    CL = ax.contourf(times/60., z, np.transpose(variable), levels = np.linspace(np.min(variable), np.max(variable), 21), cmap = 'bwr', vmin = -np.max(np.abs(variable)), vmax = np.max(np.abs(variable)))
    ax.plot(times/60., zi, color = 'k', lw = 2, label = 'Boundary Layer Depth')
    ax.plot(times/60., zf, color = 'b', lw = 2, ls = '--', label = 'Freezing Level')
    ax.plot(times/60., zt, color = 'r', lw = 2, ls = ':', label = 'Tropopause Level')
    ax.legend(loc = 2)
    ax.set_yscale('log', basey = 10, nonposy = 'clip', subsy = [2, 3, 4, 5, 6, 7, 8, 9])
    ax.set_ylim([100, 40000])
    ax.set_ylabel('Height (m)')
    ax.set_xlim([0, 10])
    ax.set_xlabel('Time (hours)')
    ax.set_title(my_title)
    plt.colorbar(mappable = CL)
    plt.savefig(file_name, dpi = 150)
    if show:
        plt.show()


print 'Plotting the tendency due to large scale rain scheme'
my_plot_style(rain, 'Tendency due to large scale rain (kg/kg/tstep)', '../tendency_lsrain.png')

print 'Plotting the tendency due to boundary layer and large scale cloud schemes'
my_plot_style(bl_cld, 'Tendency due to boundary layer and large scale cloud schemes (kg/kg/tstep)', '../tendency_blpcld.png')

print 'Plotting the tendency due to advection'
my_plot_style(advection, 'Tendency due to advection (kg/kg/tstep)', '../tendency_advect.png')

print 'Plotting the tendency due to diffusion'
my_plot_style(diffusion, 'Tendency due to diffusion (kg/kg/tstep)', '../tendency_diffus.png')

print 'Plotting the tendency due to qt_bal_cld'
my_plot_style(qt_bal_cld, 'Tendency due to qt_bal_cld scheme (kg/kg/tstep)',  '../tendency_qtbcld.png')

print 'Plotting the total tendency'
my_plot_style(total, 'Total tendency due to all schemes (kg/kg/tstep)', '../tendency_totali.png')

print 'Plotting the tendency due to idealised scheme'
my_plot_style(idealised, 'Tendency due to idealised scheme (kg/kg/tstep)', '../tendency_ideali.png')


