"""
Script to compare the three boundary layer options

For grid spacing of 200 m, the weighting used to blend the 1-D Lock(2000) 
Boundary Layer scheme with the 3-D Smagorinsky scheme is around 50/50. Here we 
compare boundary layer structure for 1-D Lock, blended and 3-D Smag.
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import lcl_key, zi_new_key, ctz_key, u_key, v_key, mcl_key, mci_key, q_key, theta_key, sshf_key, slhf_key, lwp_key

### Read the data ###
exps  = ['1-D Lock', 'Blended', '3-D Smag']
paths = {'1-D Lock' : '/work/n02/n02/xb899100/cylc-run/u-bl874/share/data/history/',
         'Blended'  : '/work/n02/n02/xb899100/cylc-run/u-bl875/share/data/history/',
         '3-D Smag' : '/work/n02/n02/xb899100/cylc-run/u-bl786/share/data/history/'}

my_data = {}
days = ["{0:02d}".format(d) for d in range(7, 11)]
for exp in exps:
    print 'Reading data for experiment: ' + exp
    my_data[exp] = {}
    for day in days:
        print ' Reading data for day: ' + day
        zi_nc   = Dataset(paths[exp] + 'zi_' + day + '.nc', 'r')
        qinc_nc = Dataset(paths[exp] + 'qinc_' + day + '.nc', 'r')
        tinc_nc = Dataset(paths[exp] + 'tempinc_' + day + '.nc', 'r')
        wind_nc = Dataset(paths[exp] + 'wind_' + day + '.nc','r')
        mr_nc   = Dataset(paths[exp] + 'mr_' + day + '.nc', 'r')
        lwp_nc  = Dataset(paths[exp] + 'lwp_' + day + '.nc', 'r')
        bouy_nc = Dataset(paths[exp] + 'bouy_' + day + '.nc', 'r')
        if day == '07':
            # This is the first day and we want to initialise the array
            # Boundary layer stuff
            my_data[exp][lcl_key]    = zi_nc.variables[lcl_key][:]*1.
            my_data[exp][zi_new_key] = zi_nc.variables[zi_new_key][:]*1.
            my_data[exp][ctz_key]    = zi_nc.variables[ctz_key][:]*1.
            my_data[exp][sshf_key]   = tinc_nc.variables[sshf_key][:]*1.
            my_data[exp][slhf_key]   = qinc_nc.variables[slhf_key][:]*1.
            # Winds
            my_data[exp][u_key]      = np.nanmean(wind_nc.variables[u_key][:], axis = (2, 3))
            my_data[exp][v_key]      = np.nanmean(wind_nc.variables[v_key][:], axis = (2, 3))
            # Cloud
            my_data[exp][mcl_key]    = mr_nc.variables[mcl_key][:]*1.
            my_data[exp][mci_key]    = mr_nc.variables[mci_key][:]*1.
            my_data[exp][lwp_key]    = lwp_nc.variables[lwp_key][:]*1.
            # Thermodynamics
            my_data[exp][q_key]      = np.nanmean(mr_nc.variables[q_key][:], axis = (2, 3))
            my_data[exp][theta_key]  = np.nanmean(bouy_nc.variables[theta_key][:], axis = (2, 3))
            # Get the height coordinates
            z_theta = wind_nc.variables['thlev_zsea_theta'][:]*1.
            times   = wind_nc.variables['min10'][:]*1.
            times_lwp = lwp_nc.variables['min1'][:]*1.
        else:
            # This is not the first day and we want to concatenate the array along the time axis
            # Boundary Layer stuff
            my_data[exp][lcl_key]    = np.concatenate((my_data[exp][lcl_key], zi_nc.variables[lcl_key][:]*1.), axis = 0)
            my_data[exp][zi_new_key] = np.concatenate((my_data[exp][zi_new_key], zi_nc.variables[zi_new_key][:]*1.), axis = 0)
            my_data[exp][ctz_key]    = np.concatenate((my_data[exp][ctz_key], zi_nc.variables[ctz_key][:]*1.), axis = 0)
            my_data[exp][sshf_key]   = np.concatenate((my_data[exp][sshf_key], tinc_nc.variables[sshf_key][:]*1.), axis = 0)
            my_data[exp][slhf_key]   = np.concatenate((my_data[exp][slhf_key], qinc_nc.variables[slhf_key][:]*1.), axis = 0)
            # Winds
            my_data[exp][u_key]      = np.concatenate((my_data[exp][u_key], np.nanmean(wind_nc.variables[u_key][:], axis = (2, 3))), axis = 0)
            my_data[exp][v_key]      = np.concatenate((my_data[exp][v_key], np.nanmean(wind_nc.variables[v_key][:], axis = (2, 3))), axis = 0)
            # Cloud
            my_data[exp][mcl_key]    = np.concatenate((my_data[exp][mcl_key], mr_nc.variables[mcl_key][:]*1.), axis = 0)
            my_data[exp][mci_key]    = np.concatenate((my_data[exp][mci_key], mr_nc.variables[mci_key][:]*1.), axis = 0)
            my_data[exp][lwp_key]    = np.concatenate((my_data[exp][lwp_key], lwp_nc.variables[lwp_key][:]*1.), axis = 0)
            # Thermodynamics
            my_data[exp][q_key]      = np.concatenate((my_data[exp][q_key], np.nanmean(mr_nc.variables[q_key][:], axis = (2, 3))), axis = 0)
            my_data[exp][theta_key]  = np.concatenate((my_data[exp][theta_key], np.nanmean(bouy_nc.variables[theta_key][:], axis = (2, 3))), axis = 0)
            times     = np.concatenate((times, wind_nc.variables['min10'][:]*1.), axis = 0)
            times_lwp = np.concatenate((times_lwp, lwp_nc.variables['min1'][:]*1.), axis = 0)
        zi_nc.close()
        qinc_nc.close()
        tinc_nc.close()
        wind_nc.close()
        mr_nc.close()
        lwp_nc.close()
        bouy_nc.close()
<<<<<<< HEAD
"""
=======

>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
### Box/whisker plots of the boundary layer stuff from zi_DD.nc files ###
fig = plt.figure()
axa = fig.add_subplot(2, 2, 1)
lcl_bp = axa.boxplot([np.nanmean(my_data[exp][lcl_key], axis = (1, 2)) for exp in exps], labels = exps, widths = 0.8)
axa.set_title('LCL')
axa.set_ylabel('Height (m)')

axb = fig.add_subplot(2, 2, 2)
zin_bp = axb.boxplot([np.nanmean(my_data[exp][zi_new_key], axis = (1, 2)) for exp in exps], labels = exps, widths = 0.8)
axb.set_title(u'z$_{i}$')
axb.set_ylabel('Height (m)')

axc = fig.add_subplot(2, 2, 3)
ctz_bp = axc.boxplot([np.nanmean(my_data[exp][ctz_key], axis = (1, 2)) for exp in exps], labels = exps, widths = 0.8)
axc.set_title(u'z$_{cld top}$')
axc.set_ylabel('Height (m)')

axd = fig.add_subplot(2, 2, 4)
E_0_bp = axd.boxplot([np.nanmean(my_data[exp][slhf_key], axis = (1, 2)) for exp in exps], labels = exps, widths = 0.8)
axd.set_title(u'E$_{sfc}$')
axd.set_ylabel(u'Latent Heat Flux (W m$^{-2}$)')

# boxplot stylisations
from matplotlib.patches import Polygon
my_boxplots = [lcl_bp, zin_bp, ctz_bp, E_0_bp]
my_axes = [axa, axb, axc, axd]
exp_cols = {'1-D Lock' : 'red',
            'Blended'  : 'purple',
            '3-D Smag' : 'blue'}
for bp in my_boxplots:
    for i in range(len(exps)):
        mylw = 2.
        myalpha = 0.5
        box = bp['boxes'][i]
        boxX = [box.get_xdata()[j] for j in range(5)]
        boxY = [box.get_ydata()[j] for j in range(5)]
        box_coords = np.column_stack([boxX, boxY])
        
        my_axes[my_boxplots.index(bp)].add_patch(Polygon(box_coords, facecolor=exp_cols[exps[i]], alpha = myalpha))
        plt.setp(bp['whiskers'][2*i:(2*i+2)], color = exp_cols[exps[i]], ls = '-', lw = mylw)
        plt.setp(bp['boxes'][i], color = exp_cols[exps[i]], ls = '-', lw = mylw)
        plt.setp(bp['caps'][2*i:(2*i+2)], color = exp_cols[exps[i]], ls = '-', lw = mylw)
        plt.setp(bp['medians'][i], color = exp_cols[exps[i]], lw = 2)
    
    plt.setp(bp['fliers'], color = 'k', marker = '.', markeredgecolor = 'k')

fig.suptitle('Values over the last 4 days of the Spinup simulations, dx = 200 m')
plt.subplots_adjust(top = 0.88, hspace = 0.3, wspace = 0.3)
plt.savefig('../BLSchemeComparison_metrics.png', dpi = 150)
plt.show()

### Profiles ###
fig = plt.figure(figsize = (18, 9))
axa = fig.add_subplot(1, 2, 1)
axa.set_ylim([0, 4000])
axa.set_ylabel('height (m)')

axb = fig.add_subplot(1, 2, 2)
axb.set_ylim([0, 4000])
axb.set_xlabel(u'$\\theta$ (K)')
axb.set_xticks([300, 304, 308, 312])
axb.set_xlim([300, 314])
axb.set_yticklabels([''])

axc = axb.twiny()
axc.set_xlim(axb.get_xlim())
axc.set_ylim([0, 4000])
axc.set_xticks([4, 8, 12, 16, 20])
axc.set_xlim([2, 20])
axc.set_xlabel('q$_{v}$ (g kg$^{-1}$)')
axc.set_yticklabels([''])

for exp in exps:
    # Plot u-wind
    axa.plot(np.nanmean(my_data[exp][u_key], axis = 0), z_theta, color = exp_cols[exp], lw = 1, label = 'u-wind ' + exp)
    # Plot v-wind
    axa.plot(np.nanmean(my_data[exp][v_key], axis = 0), z_theta, color = exp_cols[exp], ls = '--', label = 'v-wind ' + exp)
    # Plot theta
    axb.plot(np.nanmean(my_data[exp][theta_key], axis = 0), z_theta, color = exp_cols[exp], lw = 1, label = u'$\\theta$ ' + exp)
    # Plot q
    axc.plot(np.nanmean(my_data[exp][q_key], axis = 0)*1000.0, z_theta, color = exp_cols[exp], ls = '--', label = u'q$_{v}$ ' + exp)

axa.plot([0,0], [0, 40000], lw = 0.5, ls = ':', color = 'grey')
axa.legend(loc = 'upper center')
axb.legend(loc = 'center left')
axc.legend(loc = 'right')
plt.savefig('../BLSchemeComparison_profiles.png', dpi = 150)
plt.show()

### cloud field and properties? ###
x = np.arange(32)*0.2
y = np.arange(32)*0.2
lwp_levels = np.linspace(0.0001, 5000.0, 11)
mcl_levels = np.linspace(0.0001, 5., 11)

def do_plot(ax1, ax2, ax3, exp_num):
    it_max_lwp = np.where(my_data[exps[exp_num]][lwp_key] == np.nanmax(my_data[exps[exp_num]][lwp_key]))[0][0]
    test = np.abs(times - times_lwp[it_max_lwp])
    it_max = np.where(test == np.nanmin(test))[0][0]
    test2 = np.abs(times[it_max] - times_lwp)
    it_max_lwp = np.where(test2 == np.nanmin(test2))[0][0]
    total_water_xs = my_data[exps[exp_num]][lwp_key][it_max_lwp,:,:]
    iy_max, ix_max = np.where(total_water_xs == np.nanmax(total_water_xs))
    ix_max = ix_max[0]
    iy_max = iy_max[0]
    ax1.contourf(x, y, total_water_xs*1000.0, cmap = 'Greys', levels = lwp_levels, extend = 'max')
    ax1.plot(x, np.zeros_like(y)+y[iy_max], color = 'r')
    ax1.plot(np.zeros_like(x)+x[ix_max], y, color = 'r')
    
    ax2.contourf(y, z_theta/1000.0, my_data[exps[exp_num]][mcl_key][it_max,:,:,ix_max]*1000.0, cmap = 'Greys', levels = mcl_levels)
    ax2.contourf(y, z_theta/1000.0, my_data[exps[exp_num]][mci_key][it_max,:,:,ix_max]*1000.0, cmap = 'Blues', levels = mcl_levels)
    
    ax3.contourf(x, z_theta/1000.0, my_data[exps[exp_num]][mcl_key][it_max,:,iy_max,:]*1000.0, cmap = 'Greys', levels = mcl_levels)
    ax3.contourf(x, z_theta/1000.0, my_data[exps[exp_num]][mci_key][it_max,:,iy_max,:]*1000.0, cmap = 'Blues', levels = mcl_levels)


fig = plt.figure(figsize = (12, 12))
axa = fig.add_subplot(3, 3, 1, adjustable = 'box', aspect = 1)
axa.set_ylabel('y (km)')
axa.set_xticklabels([''])

axb = fig.add_subplot(3, 3, 2, adjustable = 'box', aspect = 1)
axb.set_ylim([0, 6.4])
axb.set_title(exps[0])
axb.set_ylabel('height (km)')
axb.set_xticklabels([''])

axc = fig.add_subplot(3, 3, 3, adjustable = 'box', aspect = 1)
axc.set_ylim([0, 6.4])
axc.set_yticklabels([''])
axc.set_xticklabels([''])

axd = fig.add_subplot(3, 3, 4, adjustable = 'box', aspect = 1)
axd.set_ylabel('y (km)')
axd.set_xticklabels([''])

axe = fig.add_subplot(3, 3, 5, adjustable = 'box', aspect = 1)
axe.set_ylim([0, 6.4])
axe.set_title(exps[1])
axe.set_ylabel('height (km)')
axe.set_xticklabels([''])

axf = fig.add_subplot(3, 3, 6, adjustable = 'box', aspect = 1)
axf.set_ylim([0, 6.4])
axf.set_yticklabels([''])
axf.set_xticklabels([''])

axg = fig.add_subplot(3, 3, 7, adjustable = 'box', aspect = 1)
axg.set_ylabel('y (km)')
axg.set_xlabel('x (km)')

axh = fig.add_subplot(3, 3, 8, adjustable = 'box', aspect = 1)
axh.set_ylim([0, 6.4])
axh.set_title(exps[2])
axh.set_ylabel('height (km)')
axh.set_xlabel('y (km)')

axi = fig.add_subplot(3, 3, 9, adjustable = 'box', aspect = 1)
axi.set_ylim([0, 6.4])
axi.set_yticklabels([''])
axi.set_xlabel('x (km)')
# 1-D Lock
do_plot(axa, axb, axc, 0)
# Blended
do_plot(axd, axe, axf, 1)
# 3-D Smag
do_plot(axg, axh, axi, 2)
plt.savefig('../BLSchemeComparison_Clouds.png', dpi = 150)
plt.show()

### Mean cloud hovmoller (?) ###
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
mean_mcl_levels = np.linspace(0.0001, 0.125, 11)
fig = plt.figure()
axa = fig.add_subplot(3, 1, 1)
axa.set_ylim([0, 6.4])
axa.set_ylabel('height (km)')
axa.set_xticklabels([''])
axa.set_title(exps[0])
axb = fig.add_subplot(3, 1, 2)
axb.set_ylim([0, 6.4])
axb.set_ylabel('height (km)')
axb.set_xticklabels([''])
axb.set_title(exps[1])
axc = fig.add_subplot(3, 1, 3)
axc.set_ylim([0, 6.4])
axc.set_ylabel('height (km)')
axc.set_xlabel('time (days)')
axc.set_title(exps[2])

axa.contourf(times/60./24., z_theta/1000.0, 1000.0*np.transpose(np.nanmean(my_data[exps[0]][mcl_key], axis = (2, 3))), cmap = 'Greys', levels = mean_mcl_levels, extend = 'max')
plt_mcl = axb.contourf(times/60./24., z_theta/1000.0, 1000.0*np.transpose(np.nanmean(my_data[exps[1]][mcl_key], axis = (2, 3))), cmap = 'Greys', levels = mean_mcl_levels, extend = 'max')
axbins = inset_axes(axb, width = "5%", height = "350%", loc = 6, bbox_to_anchor = (1.05, 0., 1, 1), bbox_transform = axb.transAxes, borderpad = 0.)
plt.colorbar(plt_mcl, cax = axbins, label = u'm$_{cl}$ g kg$^{-1}$')
axc.contourf(times/60./24., z_theta/1000.0, 1000.0*np.transpose(np.nanmean(my_data[exps[2]][mcl_key], axis = (2, 3))), cmap = 'Greys', levels = mean_mcl_levels, extend = 'max')
plt.suptitle('Comparison of Cloud Liquid Water Content')
plt.subplots_adjust(top = 0.88, right = 0.8)
plt.savefig('../BLSchemeComparison_CloudHovmoller.png', dpi = 150)
plt.show()

### Mean cloud fraction hovmoller (?) ###
mean_mcl_levels = np.linspace(0.0001, 0.25, 11)
fig = plt.figure()
axa = fig.add_subplot(3, 1, 1)
axa.set_ylim([0, 6.4])
axa.set_ylabel('height (km)')
axa.set_xticklabels([''])
axa.set_title(exps[0])
axb = fig.add_subplot(3, 1, 2)
axb.set_ylim([0, 6.4])
axb.set_ylabel('height (km)')
axb.set_xticklabels([''])
axb.set_title(exps[1])
axc = fig.add_subplot(3, 1, 3)
axc.set_ylim([0, 6.4])
axc.set_ylabel('height (km)')
axc.set_xlabel('time (days)')
axc.set_title(exps[2])

axa.contourf(times/60./24., z_theta/1000.0, np.transpose(np.nanmean(np.where(my_data[exps[0]][mcl_key] > 0, 1., 0.), axis = (2, 3))), cmap = 'Greys', levels = mean_mcl_levels, extend = 'max')
plt_mcf = axb.contourf(times/60./24., z_theta/1000.0, np.transpose(np.nanmean(np.where(my_data[exps[1]][mcl_key] > 0, 1., 0.), axis = (2, 3))), cmap = 'Greys', levels = mean_mcl_levels, extend = 'max')
axbins = inset_axes(axb, width = "5%", height = "350%", loc = 6, bbox_to_anchor = (1.05, 0., 1, 1), bbox_transform = axb.transAxes, borderpad = 0.)
plt.colorbar(plt_mcf, cax = axbins, label = u'Cloud Fraction')
axc.contourf(times/60./24., z_theta/1000.0, np.transpose(np.nanmean(np.where(my_data[exps[2]][mcl_key] > 0, 1., 0.), axis = (2, 3))), cmap = 'Greys', levels = mean_mcl_levels, extend = 'max')
plt.suptitle('Comparison of Cloud Fraction')
plt.subplots_adjust(top = 0.88, right = 0.8)
plt.savefig('../BLSchemeComparison_CloudFracHovmoller.png', dpi = 150)
plt.show()
<<<<<<< HEAD
"""

### Plot box and whisker plots of cloud cover and liquid water path over the last four days ###
for exp in exps:
    my_data[exp]['cloud_cover'] = np.nanmean(np.where(my_data[exp][lwp_key] > 0., 1., 0.), axis = (1, 2))

fig = plt.figure()
axa = fig.add_subplot(1, 2, 1)
ccc_bp = axa.boxplot([my_data[exp]['cloud_cover']*100. for exp in exps], labels = exps, widths = 0.8)
axa.set_title('Cloud Cover')
axa.set_ylabel('Cloud Cover (%)')

axb = fig.add_subplot(1, 2, 2)
lwp_bp = axb.boxplot([1000.0*np.nanmean(my_data[exp][lwp_key], axis = (1, 2)) for exp in exps], labels = exps, widths = 0.8)
axb.set_title('LWP')
axb.set_ylabel(u'LWP g kg$^{-1}$')

# boxplot stylisations
from matplotlib.patches import Polygon
my_boxplots = [ccc_bp, lwp_bp]
my_axes = [axa, axb]
exp_cols = {'1-D Lock' : 'red',
            'Blended'  : 'purple',
            '3-D Smag' : 'blue'}
for bp in my_boxplots:
    for i in range(len(exps)):
        mylw = 2.
        myalpha = 0.5
        box = bp['boxes'][i]
        boxX = [box.get_xdata()[j] for j in range(5)]
        boxY = [box.get_ydata()[j] for j in range(5)]
        box_coords = np.column_stack([boxX, boxY])
        
        my_axes[my_boxplots.index(bp)].add_patch(Polygon(box_coords, facecolor=exp_cols[exps[i]], alpha = myalpha))
        plt.setp(bp['whiskers'][2*i:(2*i+2)], color = exp_cols[exps[i]], ls = '-', lw = mylw)
        plt.setp(bp['boxes'][i], color = exp_cols[exps[i]], ls = '-', lw = mylw)
        plt.setp(bp['caps'][2*i:(2*i+2)], color = exp_cols[exps[i]], ls = '-', lw = mylw)
        plt.setp(bp['medians'][i], color = exp_cols[exps[i]], lw = 2)
    
    plt.setp(bp['fliers'], color = 'k', marker = '.', markeredgecolor = 'k')

fig.suptitle('Values over the last 4 days of the Spinup simulations, dx = 200 m')
plt.subplots_adjust(top = 0.88, hspace = 0.3, wspace = 0.3)
plt.savefig('../BLSchemeComparison_cloud_stuff.png', dpi = 150)
plt.show()
=======

>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696

