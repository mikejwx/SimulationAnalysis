import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from SkewT_archer import getThetaE, PTtoTemp
from netCDF4 import Dataset
from multiprocessing import Pool
from analysis_tools import send_email
import os
from STASH_keys import theta_key
import string

### Heat flux experiments ###
paths = {'H250E250' : '/nerc/n02/n02/xb899100/CloudTrail/Control2/',
         #'H400E250' : '/nerc/n02/n02/xb899100/CloudTrail/H400E250/',
         #'H350E250' : '/nerc/n02/n02/xb899100/CloudTrail/H350E250/',
         #'H300E250' : '/nerc/n02/n02/xb899100/CloudTrail/H300E250/',
         #'H200E250' : '/nerc/n02/n02/xb899100/CloudTrail/H200E250/',
         #'H150E250' : '/nerc/n02/n02/xb899100/CloudTrail/H150E250/',
         #'H100E250' : '/nerc/n02/n02/xb899100/CloudTrail/H100E250/',
         #'H050E250' : '/nerc/n02/n02/xb899100/CloudTrail/H050E250/',
### Wind Speed experiments ###
         #'U05' : '/nerc/n02/n02/xb899100/CloudTrail/U05/',
         #'U06' : '/nerc/n02/n02/xb899100/CloudTrail/U06/',
         #'U07' : '/nerc/n02/n02/xb899100/CloudTrail/U07/',
         #'U08' : '/nerc/n02/n02/xb899100/CloudTrail/U08/',
         #'U09' : '/nerc/n02/n02/xb899100/CloudTrail/U09/',
### Resolution experiments ###
         'DX0200' : '/nerc/n02/n02/xb899100/CloudTrail/Control_0200m_HRIC_INV/',
         'DX0400' : '/nerc/n02/n02/xb899100/CloudTrail/Control_0400m_HRIC_INV/',
         'DX0800' : '/nerc/n02/n02/xb899100/CloudTrail/Control_0800m_HRIC_INV/',
         'DX1600' : '/nerc/n02/n02/xb899100/CloudTrail/Control_1600m_HRIC_INV/',
         'DX1600S': '/nerc/n02/n02/xb899100/CloudTrail/Control_1600m_HRIC_INV_S/'}

# define constants to compute the surface equivalent potential temperature
p_sfc = 101700.0
#keys = ['H050E250', 'H100E250', 'H150E250', 'H200E250', 'H250E250', 'H300E250', 'H350E250', 'H400E250', 'U05', 'U06', 'U07', 'U08', 'U09', 'DX0200', 'DX0400', 'DX0800', 'DX1600', 'DX1600S']
keys = ['H250E250', 'DX0200', 'DX0400', 'DX0800', 'DX1600']
hours = ["{0:02d}".format(hour) for hour in range(0, 24, 3)]
hours_U05 = ["{0:02d}".format(hour) for hour in range(0, 16, 4)]
thetaE = {}
for key in keys:
    print key
    # get cheap and dirty temperature using the initial surface pressure and theta
    if key == 'U05':
        for hour in hours_U05:
            bouy_nc = Dataset(paths[key] + 'bouy_' + hour + '.nc', 'r')
            if hour == hours_U05[0]:
                theta       = bouy_nc.variables[theta_key][:,0,:,:]*1.
                temp        = PTtoTemp(theta, p_sfc, t_units = 'K', p_units = 'Pa')
                thetaE[key] = getThetaE(theta, temp, p_sfc, t_units = 'K', p_units = 'Pa')
            else:
                theta       = bouy_nc.variables[theta_key][:,0,:,:]*1.
                temp        = PTtoTemp(theta, p_sfc, t_units = 'K', p_units = 'Pa')
                thetaE[key] = np.concatenate((thetaE[key], getThetaE(theta, temp, p_sfc, t_units = 'K', p_units = 'Pa')), axis = 0)
    else:
        for hour in hours:
            bouy_nc = Dataset(paths[key] + 'bouy_' + hour + '.nc', 'r')
            if hour == hours[0]:
                theta       = bouy_nc.variables[theta_key][:,0,:,:]*1.
                temp        = PTtoTemp(theta, p_sfc, t_units = 'K', p_units = 'Pa')
                thetaE[key] = getThetaE(theta, temp, p_sfc, t_units = 'K', p_units = 'Pa')
            else:
                theta       = bouy_nc.variables[theta_key][:,0,:,:]*1.
                temp        = PTtoTemp(theta, p_sfc, t_units = 'K', p_units = 'Pa')
                thetaE[key] = np.concatenate((thetaE[key], getThetaE(theta, temp, p_sfc, t_units = 'K', p_units = 'Pa')), axis = 0)
"""
### Plot how the minimum and maximum surface thetae changes with U and H ###
heat_cmap = mpl.cm.get_cmap('Reds')
wind_cmap = mpl.cm.get_cmap('Blues')
my_colors = {}
for key in keys:
    if 'E250' in key:
        my_colors[key] = heat_cmap(float(key[1:4])/400.0)
    elif 'U' in key:
        my_colors[key] = wind_cmap((float(key[1:])-4)/10.0)

fig = plt.figure()
axa = fig.add_subplot(2, 2, 1)
axa.set_xticklabels([''])
axa.set_xlim([4, 11])
axa.set_ylim([0, 80])
axa.set_ylabel(u'$\\theta_{e, max}^{\prime}$ (K)')

axb = fig.add_subplot(2, 2, 2)
axb.set_xticklabels([''])
axb.set_ylim([0, 80])
axb.set_yticklabels([''])
axb.set_xlim([25, 425])

axc = fig.add_subplot(2, 2, 3)
axc.set_xlim([4, 11])
axc.set_ylim([-6, 0])
axc.set_xlabel('Wind Speed (m s$^{-1}$)')
axc.set_ylabel(u'$\\theta_{e, min}^{\prime}$ (K)')

axd = fig.add_subplot(2, 2, 4)
axd.set_xlim([25, 425])
axd.set_yticklabels([''])
axd.set_ylim([-6, 0])
axd.set_xlabel('Island SHF (W m$^{-2}$)')
for key in keys:
    # Do the maximum pertubation
    if ('U' in key) or (key == 'H250E250'):
        if key not in ['U05', 'H250E250']:
            print key + ' axa'
            thetaE_mean_temp = np.nanmean(thetaE[key], axis = (1, 2))[23:97]
            axa.plot(float(key[1:]), np.nanmax([thetaE[key][it,:,:] - thetaE_mean_temp[it] for it in range(len(thetaE_mean_temp))]), color = my_colors[key], marker = 'o', markeredgecolor = 'None')
        elif key == 'H250E250':
            print key + ' axa'
            thetaE_mean_temp = np.nanmean(thetaE[key], axis = (1, 2))
            axa.plot(10.0, np.nanmax([thetaE[key][it,:,:] - thetaE_mean_temp[it] for it in range(len(thetaE_mean_temp))]), color = wind_cmap(1.0), marker = 'o', markeredgecolor = 'None')
        else:
            print key + ' axa'
            thetaE_mean_temp = np.nanmean(thetaE[key], axis = (1, 2))
            axa.plot(float(key[1:]), np.nanmax([thetaE[key][it,:,:] - thetaE_mean_temp[it] for it in range(len(thetaE_mean_temp))]), color = my_colors[key], marker = 'o', markeredgecolor = 'None')
    if 'E250' in key:
        print key + ' axb'
        thetaE_mean_temp = np.nanmean(thetaE[key], axis = (1, 2))
        axb.plot(float(key[1:4]), np.nanmax([thetaE[key][it,:,:] - thetaE_mean_temp[it] for it in range(len(thetaE_mean_temp))]), color = my_colors[key], marker = 'o', markeredgecolor = 'None')
    # Do the minimum perturbation
    if ('U' in key) or (key == 'H250E250'):
        if key not in ['U05', 'H250E250']:
            print key + ' axc'
            thetaE_mean_temp = np.nanmean(thetaE[key], axis = (1, 2))[23:97]
            axc.plot(float(key[1:]), np.nanmin([thetaE[key][it,:,:] - thetaE_mean_temp[it] for it in range(len(thetaE_mean_temp))]), color = my_colors[key], marker = 'o', markeredgecolor = 'None')
        elif key == 'H250E250':
            print key + ' axc'
            thetaE_mean_temp = np.nanmean(thetaE[key], axis = (1, 2))
            axc.plot(10.0, np.nanmin([thetaE[key][it,:,:] - thetaE_mean_temp[it] for it in range(len(thetaE_mean_temp))]), color = wind_cmap(1.0), marker = 'o', markeredgecolor = 'None')
        else:
            print key + ' axc'
            thetaE_mean_temp = np.nanmean(thetaE[key], axis = (1, 2))
            axc.plot(float(key[1:]), np.nanmin([thetaE[key][it,:,:] - thetaE_mean_temp[it] for it in range(len(thetaE_mean_temp))]), color = my_colors[key], marker = 'o', markeredgecolor = 'None')
    if 'E250' in key:
        print key + ' axd'
        thetaE_mean_temp = np.nanmean(thetaE[key], axis = (1, 2))
        axd.plot(float(key[1:4]), np.nanmin([thetaE[key][it,:,:] - thetaE_mean_temp[it] for it in range(len(thetaE_mean_temp))]), color = my_colors[key], marker = 'o', markeredgecolor = 'None')

plt.savefig('../warmplume_coldpool_intensity.png', dpi = 150, bbox_inches = 'tight')
plt.show()

### Plot the surface thetae for the heat flux experiments ###
my_levels = np.array([-2.0, -1.0, -0.5, -0.25, -0.1, 0.1, 0.25, 0.5, 1.0, 2.0])
n_levels = float(len(my_levels))
my_cmap = mpl.cm.get_cmap('bwr')
my_colors = my_cmap((np.arange(n_levels)+1.0)/(n_levels))

heat_exp =  ['H050E250', 'H100E250', 'H150E250', 'H200E250', 'H250E250', 'H300E250', 'H350E250', 'H400E250']
fig, axes = plt.subplots(nrows = 4, ncols = 2)
for key in heat_exp:
    ax = axes.flat[heat_exp.index(key)]
    ax.set_adjustable('box')
    ax.set_aspect(1)
    idx = thetaE[key].shape[0]/2
    x, y = np.meshgrid(np.arange(thetaE[key][idx,:,:].shape[1])*0.1, np.arange(thetaE[key][idx,:,:].shape[0])*0.1)
    
    im = ax.contourf(x, y, thetaE[key][idx,:,:] - np.nanmean(thetaE[key][idx,:,:]), levels = my_levels, colors = my_colors, extend = 'both')
    ax.set_title(key + ' at T+720 mins', fontsize = 12)
    
    if (heat_exp.index(key)+1)%2 == 0:
        ax.set_yticklabels([''])
    else:
        ax.set_ylabel('y (km)')
    if (heat_exp.index(key)+1) in [1, 2, 3, 4, 5, 6]:
        ax.set_xticklabels([''])
    else:
        ax.set_xlabel('x (km)')

cbar_ax = fig.add_axes([0.1, 0.1, 0.85, 0.02])
cbar = fig.colorbar(im, cax = cbar_ax, orientation = 'horizontal', label = '$\\theta_{e}^{\prime}$ (K)')
cbar.ax.set_xticklabels(['-2', '-1', '-0.5', '-0.25', '-0.1', '0.1', '0.25', '0.5', '1', '2'])
plt.subplots_adjust(left = 0.1, bottom = 0.2, right = 0.95, top = 0.9, wspace = 0.1, hspace = 0.1)
plt.savefig('../heat_thetaE_sfc.png', dpi = 150, bbox_inches = 'tight')
plt.show()

### plot the surface thetae for the wind experiments ###
wind_exp =  ['U05', 'U06', 'U07', 'U08', 'U09', 'H250E250']
fig, axes = plt.subplots(nrows = 3, ncols = 2)
for key in wind_exp:
    ax = axes.flat[wind_exp.index(key)]
    ax.set_adjustable('box')
    ax.set_aspect(1)
    idx = thetaE[key].shape[0]/2
    x, y = np.meshgrid(np.arange(thetaE[key][idx,:,:].shape[1])*0.1, np.arange(thetaE[key][idx,:,:].shape[0])*0.1)
    
    im = ax.contourf(x, y, thetaE[key][idx,:,:] - np.nanmean(thetaE[key][idx,:,:]), levels = my_levels, colors = my_colors, extend = 'both')
    ax.set_title([key if 'H' not in key else 'U10'][0] + ' at T+720 mins', fontsize = 12)
    
    if (wind_exp.index(key)+1)%2 == 0:
        ax.set_yticklabels([''])
    else:
        ax.set_ylabel('y (km)')
    if (wind_exp.index(key)+1) in [1, 2, 3, 4]:
        ax.set_xticklabels([''])
    else:
        ax.set_xlabel('x (km)')

cbar_ax = fig.add_axes([0.1, 0.1, 0.85, 0.02])
cbar = fig.colorbar(im, cax = cbar_ax, orientation = 'horizontal', label = '$\\theta_{e}^{\prime}$ (K)')
cbar.ax.set_xticklabels(['-2', '-1', '-0.5', '-0.25', '-0.1', '0.1', '0.25', '0.5', '1', '2'])
plt.subplots_adjust(left = 0.1, bottom = 0.2, right = 0.95, top = 0.9, wspace = 0.1, hspace = 0.1)
plt.savefig('../wind_thetaE_sfc.png', dpi = 150, bbox_inches = 'tight')
plt.show()
"""

### plot the surface thetae for the grid spacing experiments ###
my_levels = np.array([-2.0, -1.0, -0.5, -0.25, -0.1, 0.1, 0.25, 0.5, 1.0, 2.0])
n_levels = float(len(my_levels))
my_cmap = mpl.cm.get_cmap('bwr')
my_colors = my_cmap((np.arange(n_levels)+1.0)/(n_levels))
dx_exp =  ['H250E250', 'DX0200', 'DX0400', 'DX0800', 'DX1600']
A = string.maketrans('','')
nodigs = A.translate(A, string.digits)
fig = plt.figure(figsize = (6, 10))
for key in dx_exp:
    ax = fig.add_subplot(5, 1, dx_exp.index(key) + 1, adjustable = 'box', aspect = 1)
    idx = thetaE[key].shape[0]/2
    grid_spacing = [0.1 if key == 'H250E250' else float(key.translate(A, nodigs))/1000.][0]
    x, y = np.meshgrid(np.arange(thetaE[key][idx,:,:].shape[1])*grid_spacing, np.arange(thetaE[key][idx,:,:].shape[0])*grid_spacing)
    
    im = ax.contourf(x, y, thetaE[key][idx,:,:] - np.nanmean(thetaE[key][idx,:,:]), levels = my_levels, colors = my_colors, extend = 'both')
    ax.set_title([key if 'H' not in key and 'U' not in key else 'DX0100'][0] + ' at T+720 mins', fontsize = 12)
    
    ax.set_ylabel('y (km)')
    if (dx_exp.index(key)+1) in [1, 2, 3, 4]:
        ax.set_xticklabels([''])
    else:
        ax.set_xlabel('x (km)')

cbar_ax = fig.add_axes([0.12, 0.1, 0.78, 0.02])
cbar = fig.colorbar(im, cax = cbar_ax, orientation = 'horizontal', label = '$\\theta_{e}^{\prime}$ (K)')
cbar.ax.set_xticklabels(['-2', '-1', '-0.5', '-0.25', '-0.1', '0.1', '0.25', '0.5', '1', '2'])
plt.subplots_adjust(bottom = 0.17, wspace = 0.17, hspace = 0.30)
plt.savefig('../dx_thetaE_sfc.png', dpi = 150, bbox_inches = 'tight')
plt.show()


