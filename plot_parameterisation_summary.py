import numpy as np
import matplotlib as mpl
<<<<<<< HEAD
=======
mpl.use('Agg')
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import *
from SkewT_archer import *
import os
from scipy import integrate
<<<<<<< HEAD
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
=======
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696

qlconv_key = u'STASH_m01s05i213'
conv_rain_amt_key =  u'STASH_m01s05i214'
################################################################################
############################ read in the data ##################################
################################################################################
experiments = ['1600mS']
paths = {'1600mS' : '/nerc/n02/n02/xb899100/CloudTrail/Control_1600m_HRIC_INV_S/'}

hours = ["{0:02d}".format(hour) for hour in range(0, 24, 3)]
my_data = {}
for exp in experiments:
    base_path = paths[exp]
    my_data[exp] = {}
    for hour in hours:
        with Dataset(base_path + 'bouy_' + hour + '.nc', 'r') as bouy_nc:
            if hour == hours[0]:
                z_theta = bouy_nc.variables['thlev_zsea_theta'][:]*1.
                # find the time key for the main data
                time_key = [tkey for tkey in bouy_nc.variables.keys() if 'min' in tkey][0]
                my_data[exp][time_key] = bouy_nc.variables[time_key][1:]*1.
                my_data[exp][theta_key] = bouy_nc.variables[theta_key][1:,:,:,:]*1.
                my_data[exp][q_key] = bouy_nc.variables[q_key][1:,:,:,:]*1.
            else:
                my_data[exp][theta_key] = np.concatenate((my_data[exp][theta_key], bouy_nc.variables[theta_key][:]*1.), axis = 0)
                my_data[exp][q_key] = np.concatenate((my_data[exp][q_key], bouy_nc.variables[q_key][:]*1.), axis = 0)
                my_data[exp][time_key] = np.concatenate((my_data[exp][time_key], bouy_nc.variables[time_key][:]*1.), axis = 0)
        
        with Dataset(base_path + 'wind_' + hour + '.nc', 'r') as wind_nc:
            if hour == hours[0]:
                my_data[exp][n_key] = wind_nc.variables[n_key][1:,:,:,:]*1.
            else:
                my_data[exp][n_key] = np.concatenate((my_data[exp][n_key], wind_nc.variables[n_key][:]*1.), axis = 0)
        
        with Dataset(base_path + 'fluxes_' + hour + '.nc', 'r') as fluxes_nc:
            if hour == hours[0]:
                z_rho = fluxes_nc.variables['rholev_zsea_rho'][:]*1.
                my_data[exp][rho_key] = fluxes_nc.variables[rho_key][1:,:,:,:]*1.
                my_data[exp][pthe_key] = fluxes_nc.variables[pthe_key][1:,:,:,:]*1.
            else:
                my_data[exp][rho_key] = np.concatenate((my_data[exp][rho_key], fluxes_nc.variables[rho_key][:]*1.), axis = 0)
                my_data[exp][pthe_key] = np.concatenate((my_data[exp][pthe_key], fluxes_nc.variables[pthe_key][:]*1.), axis = 0)
        
        with Dataset(base_path + 'conv_' + hour + '.nc', 'r') as conv_nc:
            if hour == hours[0]:
                z_conv = conv_nc.variables['thlev_zsea_theta'][:]*1.
                my_data[exp][qlconv_key] = conv_nc.variables[qlconv_key][:]*1.
                my_data[exp][conv_rain_amt_key] = conv_nc.variables[conv_rain_amt_key][:]*1.
                # find the time key for the convective rain amount
                conv_accum_time_key = [tkey for tkey in conv_nc.variables.keys() if 'accum' in tkey][0]
                my_data[exp][conv_accum_time_key] = conv_nc.variables[conv_accum_time_key][:]*1.
            else:
                my_data[exp][qlconv_key] = np.concatenate((my_data[exp][qlconv_key], conv_nc.variables[qlconv_key][:]*1.), axis = 0)
                my_data[exp][conv_rain_amt_key] = np.concatenate((my_data[exp][conv_rain_amt_key][:], conv_nc.variables[conv_rain_amt_key][:]*1.), axis = 0)
                my_data[exp][conv_accum_time_key] = np.concatenate((my_data[exp][conv_accum_time_key][:], conv_nc.variables[conv_accum_time_key][:]*1.), axis = 0)
    
    with Dataset(base_path + 'lwp_00.nc', 'r') as lwp_nc:
        my_data[exp][lwp_key] = lwp_nc.variables[lwp_key][1:,:,:]*1.
        my_data[exp][ls_rain_amt_key] = lwp_nc.variables[ls_rain_amt_key][:]*1.
        # find the key corresponding to the lwp times
        lwp_time_key = [tkey for tkey in lwp_nc.variables.keys() if 'min' in tkey][0]
        my_data[exp][lwp_time_key] = lwp_nc.variables[lwp_time_key][1:]*1.
        # find the key corresponding to the rain accumulation times
        rain_time_key = [tkey for tkey in lwp_nc.variables.keys() if 'accum' in tkey][0]
        my_data[exp][rain_time_key] = lwp_nc.variables[rain_time_key][:]*1.

# manufacture some coordinate system and get distance from island centre
x, y = np.meshgrid(np.arange(my_data[exp][lwp_key].shape[2])*1.6, np.arange(my_data[exp][lwp_key].shape[1])*1.6)
x_c, y_c = 108.0, 16.0
R = np.sqrt((x - x_c)**2 + (y - y_c)**2)
island_area = 50.0
R_i = np.sqrt(island_area/np.pi)

# index of 10 m height
iz10 = np.where(np.abs(z_theta - 10.0) == np.min(np.abs(z_theta - 10.0)))[0][0]

################################################################################
############################ make some plots  ##################################
################################################################################

# get the time indexes for daytime
# for lwp -> cloud frequency
lwp_tidx = [i for i in range(len(my_data[exp][lwp_time_key])) if (360.0 < my_data[exp][lwp_time_key][i])*(my_data[exp][lwp_time_key][i] < 1080.0)]
lwp_tidx2 = [i for i in range(len(my_data[exp][time_key])) if (360.0 < my_data[exp][time_key][i])*(my_data[exp][time_key][i] < 1080.0)]
# for total daytime rainfall
idx0 = np.where(my_data[exp][rain_time_key] == 360.0)[0][0]
idx1 = np.where(my_data[exp][rain_time_key] == 1080.0)[0][0]
# for midday
idx12 = np.where(my_data[exp][time_key] == 720.0)[0][0]

<<<<<<< HEAD
=======
# compute equivalent potential temperature
thetae_key = 'equivalent potential temperature'
thetae_anom_key = 'equivalent potential temperature anomaly'
my_data[exp][thetae_key] = getThetaE(my_data[exp][theta_key][idx12,0,:,:], PTtoTemp(my_data[exp][theta_key][idx12,0,:,:], my_data[exp][pthe_key][idx12,0,:,:], t_units = 'K', p_units = 'Pa'), my_data[exp][pthe_key][idx12,0,:,:], t_units = 'K', p_units = 'Pa')
my_data[exp][thetae_anom_key] = my_data[exp][thetae_key] - my_data[exp][thetae_key].mean()

>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
# compute the convective liquid water path
convective_lwp = 'convective lwp'
my_data[exp][convective_lwp] = np.empty_like(my_data[exp][qlconv_key][:,0,:,:])
for it in range(my_data[exp][qlconv_key].shape[0]):
    my_data[exp][convective_lwp][it,:,:] = integrate.trapz(x = z_conv, y = my_data[exp][qlconv_key][it,:,:,:], axis = 0)

# merge the resolved and convective lwp
total_lwp_key = 'total lwp'
my_data[exp][total_lwp_key] = np.empty_like(my_data[exp][convective_lwp])
for it in range(my_data[exp][time_key].size):
    it_match = np.where(my_data[exp][lwp_time_key] == my_data[exp][time_key][it])[0][0]
    my_data[exp][total_lwp_key][it,:,:] = my_data[exp][lwp_key][it_match,:,:] + my_data[exp][convective_lwp][it,:,:]

# set the levels for plotting
<<<<<<< HEAD
my_temp_levels = np.array([-2.0, -1.0, -0.5, -0.25, -0.1, 0.1, 0.25, 0.5, 1.0, 2.0])
n_temp_levels = float(len(my_temp_levels))
my_temp_cmap = mpl.cm.get_cmap('bwr')
my_temp_colors = my_temp_cmap((np.arange(n_temp_levels)+1.0)/(n_temp_levels))

my_rain_levels = np.array([0.01, 0.5, 1., 2., 4., 8., 16., 32.])
n_rain_levels = float(len(my_rain_levels))
my_rain_cmap = mpl.cm.get_cmap('YlGnBu')
my_rain_colors = my_rain_cmap((np.arange(n_rain_levels)+1.0)/(n_rain_levels))


my_cldfreq_levels = np.arange(0, 0.61, 0.075)
my_wind_levels = [level for level in np.arange(-3.5, 3.51, 0.5) if level != 0.0]

fig = plt.figure(figsize = (12, 7))
=======
my_cldfreq_levels = np.arange(0, 0.61, 0.075)
my_wind_levels = [level for level in np.arange(-3.5, 3.51, 0.5) if level != 0.0]
my_temp_levels = [level for level in np.arange(-5., 5.01, 0.5) if level != 0.0]
my_rain_levels = [0.01, 0.5, 1., 2., 4., 8., 16., 32.]
my_rain_colors = ['blue', 'cornflowerblue', 'olive', 'gold', 'orange', 'red', 'magenta', 'ghostwhite']

fig = plt.figure(figsize = (8.5, 5))
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
axa = fig.add_subplot(2, 2, 1, adjustable = 'box', aspect = 1)
axb = fig.add_subplot(2, 2, 2, adjustable = 'box', aspect = 1)
axc = fig.add_subplot(2, 2, 3, adjustable = 'box', aspect = 1)
axd = fig.add_subplot(2, 2, 4, adjustable = 'box', aspect = 1)

cld_freq_plt = axa.contourf(x, y, np.nanmean(np.where(my_data[exp][total_lwp_key][lwp_tidx2,:,:] > 0, 1.0, 0.0), axis = 0), levels = my_cldfreq_levels, cmap = 'Greys_r', extend = 'max')
cbax = plt.colorbar(cld_freq_plt, ax = axa, orientation = 'horizontal', label = 'Cloud Frequency')
cbax.set_ticks(my_cldfreq_levels[::2])
axa.contour(x, y, R, levels = [R_i], colors = ['r'])
axa.set_xticklabels([''])
axa.set_ylabel('y (km)')
<<<<<<< HEAD
axa.set_title('a) Cloud Frequency')
# Index of the 2 m height level
iz2 = np.where(np.abs(z_theta - 2.0) == np.min(np.abs(z_theta - 2.0)))[0][0]
warm_plume_plt = axb.contourf(x, y, my_data[exp][theta_key][idx12,iz2,:,:] - my_data[exp][theta_key][idx12,iz2,:,:].mean(), colors = my_temp_colors, levels = my_temp_levels, extend = 'both')
warm_plume_plt.cmap.set_under('navy')
warm_plume_plt.cmap.set_over('firebrick')
cbax = plt.colorbar(warm_plume_plt, ax = axb, orientation = 'horizontal', label = u'$\\theta^{\prime}$ (K)')
cbax.set_ticks([-2, -1, -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 1, 2])
cbax.ax.set_xticklabels([-2, -1, -0.5, -0.25, -0.1, '', 0.1, 0.25, 0.5, 1, 2])
axb.contour(x, y, R, levels = [R_i], colors = ['k'])
axb.set_yticklabels([''])
axb.set_xticklabels([''])
axb.set_title('b) Surface Warm Plume')

wind_plt = axc.contourf(x, y, my_data[exp][n_key][idx12,iz10,:,:], cmap = 'bwr', levels = my_wind_levels, extend = 'both')
wind_plt.cmap.set_under('navy')
wind_plt.cmap.set_over('firebrick')
axins = inset_axes(axc, width = "100%", height = "20%", loc = 3, bbox_to_anchor = (0.0, -0.725, 1.0, 1.0), bbox_transform = axc.transAxes, borderpad = 0.0)
cbax = plt.colorbar(wind_plt, cax = axins, orientation = 'horizontal', label = 'Across-flow wind (m s$^{-1}$)')
=======

warm_plume_plt = axb.contourf(x, y, my_data[exp][thetae_anom_key], cmap = 'bwr', levels = my_temp_levels, extend = 'both')
cbax = plt.colorbar(warm_plume_plt, ax = axb, orientation = 'horizontal', label = u'$\\theta_{e,sfc}^{\prime}$ (K)')
cbax.set_ticks([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])
axb.contour(x, y, R, levels = [R_i], colors = ['k'])
axb.set_yticklabels([''])
axb.set_xticklabels([''])

wind_plt = axc.contourf(x, y, my_data[exp][n_key][idx12,iz10,:,:], cmap = 'bwr', levels = my_wind_levels, extend = 'both')
cbax = plt.colorbar(wind_plt, ax = axc, orientation = 'horizontal', label = 'Across-flow wind (m s$^{-1}$)')
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
cbax.set_ticks([-3, -2, -1, 0, 1, 2, 3])
axc.contour(x, y, R, levels = [R_i], colors = ['k'])
axc.set_ylabel('y (km)')
axc.set_xlabel('x (km)')
<<<<<<< HEAD
axc.set_title('c) Across-flow wind')

precip = (my_data[exp][ls_rain_amt_key][idx1,:,:] + my_data[exp][conv_rain_amt_key][idx1,:,:]) - (my_data[exp][ls_rain_amt_key][idx0,:,:] - my_data[exp][conv_rain_amt_key][idx0,:,:])
rain_plt = axd.contourf(x, y, precip, colors = my_rain_colors, levels = my_rain_levels, extend = 'max')
axins = inset_axes(axd, width = "100%", height = "20%", loc = 3, bbox_to_anchor = (0.0, -0.725, 1.0, 1.0), bbox_transform = axd.transAxes, borderpad = 0.0)
cbax = plt.colorbar(rain_plt, cax = axins, orientation = 'horizontal', label = 'Precipitation Total (mm)')
cbax.set_ticklabels([int(level) if level > 1 else level for level in my_rain_levels])
axd.contour(x, y, R, levels = [R_i], colors = ['k'])
axd.set_title('d) P$_{max}$ = ' + str(round(precip.max(), 2)) + ' mm')
axd.set_yticklabels([''])
axd.set_xlabel('x (km)')

plt.subplots_adjust(bottom = 0.3, wspace = 0.2, hspace = 0.1)
plt.savefig('../parametrisedSummary_1600m.png', dpi = 250, bbox_inches = 'tight')
plt.show()
=======

precip = (my_data[exp][ls_rain_amt_key][idx1,:,:] + my_data[exp][conv_rain_amt_key][idx1,:,:]) - (my_data[exp][ls_rain_amt_key][idx0,:,:] - my_data[exp][conv_rain_amt_key][idx0,:,:])
rain_plt = axd.contourf(x, y, precip, colors = my_rain_colors, levels = my_rain_levels, extend = 'max')
cbax = plt.colorbar(rain_plt, ax = axd, orientation = 'horizontal', label = 'Precipitation Total (mm)')
cbax.set_ticklabels([int(level) if level > 1 else level for level in my_rain_levels])
axd.contour(x, y, R, levels = [R_i], colors = ['k'])
axd.set_title('Maximum Precip = ' + str(round(precip.max(), 2)) + ' mm')
axd.set_yticklabels([''])
axd.set_xlabel('x (km)')

plt.subplots_adjust(wspace = 0.2, hspace = 0.05)
plt.savefig('../parametrisedSummary_1600m.png', dpi = 150, bbox_inches = 'tight')
plt.close('all')
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696

