import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import theta_key, rho_key, q_key, pthe_key, prho_key, mv_key, mr_key, zi_new_key, temp_key
from SkewT_archer import getCAPE
from scipy import interpolate, integrate

### read in all of this data into dictionaries
#   for the control cloud trail simulation, and for the 800 m same-domain-size, 
#   cloud trail-free simulation.

path_100 = '/nerc/n02/n02/xb899100/CloudTrail/Control/'
path_800 = '/work/n02/n02/xb899100/cylc-run/u-bn261/share/data/history/'

hours = ["{0:02d}".format(hour) for hour in range(0, 24, 3)]

data_100m = {}
data_800m = {}
l_CAPE = True
for hour in hours:
    print 'Reading the 100 m data for hour ' + hour
    bouy_nc   = Dataset(path_100 + 'bouy_' + hour + '.nc', 'r')
    fluxes_nc = Dataset(path_100 + 'fluxes_' + hour + '.nc', 'r')
    mr_nc     = Dataset(path_100 + 'mr_' + hour + '.nc', 'r')
    zi_nc     = Dataset(path_100 + 'zi_' + hour + '.nc', 'r')
    
    if hour == hours[0]:
        z_theta = bouy_nc.variables['thlev_zsea_theta'][:]*1.
        z_rho   = fluxes_nc.variables['rholev_zsea_rho'][:]*1.
        data_100m['time']      = zi_nc.variables['time'][:]
        data_100m['mean_temp'] = np.nanmean(bouy_nc.variables[temp_key][1:,:,:,:], axis = (1, 2, 3))
        data_100m['pwat']      = np.nanmean(np.array([integrate.trapz(x = z_theta, y = interpolate.interp1d(x = z_rho, y = fluxes_nc.variables[rho_key][it,:,:,:], fill_value = 'extrapolate', axis = 0)(z_theta)*mr_nc.variables[mv_key][it,:,:,:], axis = 0) for it in range(1, mr_nc.variables[mv_key].shape[0])]), axis = (1, 2))
        if l_CAPE:
            data_100m['CAPE']  = np.zeros_like(bouy_nc.variables[temp_key][1:,0,0,0])
            data_100m['CIN']   = np.zeros_like(bouy_nc.variables[temp_key][1:,0,0,0])
            for it in range(1, bouy_nc.variables[temp_key].shape[0]):
                data_100m['CAPE'][it-1], data_100m['CIN'][it-1], T_parcel, P_parcel, LCL_pressure, LFC_pressure = getCAPE(np.nanmean(bouy_nc.variables[temp_key][it,:-1,:,:], axis = (1, 2)), np.nanmean(bouy_nc.variables[q_key][it,:-1,:,:], axis = (1, 2)), np.nanmean(fluxes_nc.variables[pthe_key][it,:-1,:,:], axis = (1, 2)), parcel_type = 1, t_units = 'K', q_units = 'kg/kg', p_units = 'Pa')
        data_100m['zi']        = np.nanmean(zi_nc.variables[zi_new_key][:], axis = (1, 2))
    else:
        data_100m['time']      = np.concatenate((data_100m['time'], zi_nc.variables['time'][:]), axis = 0)
        data_100m['mean_temp'] = np.concatenate((data_100m['mean_temp'], np.nanmean(bouy_nc.variables[temp_key][:], axis = (1, 2, 3))), axis = 0)
        data_100m['pwat']      = np.concatenate((data_100m['pwat'], np.nanmean(np.array([integrate.trapz(x = z_theta, y = interpolate.interp1d(x = z_rho, y = fluxes_nc.variables[rho_key][it,:,:,:], fill_value = 'extrapolate', axis = 0)(z_theta)*mr_nc.variables[mv_key][it,:,:,:], axis = 0) for it in range(mr_nc.variables[mv_key].shape[0])]), axis = (1, 2))), axis = 0)
        if l_CAPE:
            CAPE_temp = np.zeros_like(bouy_nc.variables[temp_key][:,0,0,0])
            CIN_temp  = np.zeros_like(bouy_nc.variables[temp_key][:,0,0,0])
            for it in range(bouy_nc.variables[temp_key].shape[0]):
                CAPE_temp[it], CIN_temp[it], T_parcel, P_parcel, LCL_pressure, LFC_pressure = getCAPE(np.nanmean(bouy_nc.variables[temp_key][it,:-1,:,:], axis = (1, 2)), np.nanmean(bouy_nc.variables[q_key][it,:-1,:,:], axis = (1, 2)), np.nanmean(fluxes_nc.variables[pthe_key][it,:-1,:,:], axis = (1, 2)), parcel_type = 1, t_units = 'K', q_units = 'kg/kg', p_units = 'Pa')
            data_100m['CAPE'] = np.concatenate((data_100m['CAPE'], CAPE_temp), axis = 0)
            data_100m['CIN']  = np.concatenate((data_100m['CIN'], CIN_temp), axis = 0)
        data_100m['zi']        = np.concatenate((data_100m['zi'], np.nanmean(zi_nc.variables[zi_new_key][:], axis = (1, 2))), axis = 0)
    
    bouy_nc.close()
    fluxes_nc.close()
    mr_nc.close()
    zi_nc.close()
    print 'Reading the 800 m data for hour ' + hour
    bouy_nc   = Dataset(path_800 + 'bouy_' + hour + '.nc', 'r')
    fluxes_nc = Dataset(path_800 + 'fluxes_' + hour + '.nc', 'r')
    mr_nc     = Dataset(path_800 + 'mr_' + hour + '.nc', 'r')
    zi_nc     = Dataset(path_800 + 'zi_' + hour + '.nc', 'r')
    
    if hour == hours[0]:
        z_theta = bouy_nc.variables['thlev_zsea_theta'][:]*1.
        z_rho   = fluxes_nc.variables['rholev_zsea_rho'][:]*1.
        data_800m['time']      = zi_nc.variables['time'][:]
        data_800m['mean_temp'] = np.nanmean(bouy_nc.variables[temp_key][1:,:,:,:], axis = (1, 2, 3))
        data_800m['pwat']      = np.nanmean(np.array([integrate.trapz(x = z_theta, y = interpolate.interp1d(x = z_rho, y = fluxes_nc.variables[rho_key][it,:,:,:], fill_value = 'extrapolate', axis = 0)(z_theta)*mr_nc.variables[mv_key][it,:,:,:], axis = 0) for it in range(1, mr_nc.variables[mv_key].shape[0])]), axis = (1, 2))
        if l_CAPE:
            data_800m['CAPE']  = np.zeros_like(bouy_nc.variables[temp_key][1:,0,0,0])
            data_800m['CIN']   = np.zeros_like(bouy_nc.variables[temp_key][1:,0,0,0])
            for it in range(1, bouy_nc.variables[temp_key].shape[0]):
                data_800m['CAPE'][it-1], data_800m['CIN'][it-1], T_parcel, P_parcel, LCL_pressure, LFC_pressure = getCAPE(np.nanmean(bouy_nc.variables[temp_key][it,:-1,:,:], axis = (1, 2)), np.nanmean(bouy_nc.variables[q_key][it,:-1,:,:], axis = (1, 2)), np.nanmean(fluxes_nc.variables[pthe_key][it,:-1,:,:], axis = (1, 2)), parcel_type = 1, t_units = 'K', q_units = 'kg/kg', p_units = 'Pa')
        data_800m['zi']        = np.nanmean(zi_nc.variables[zi_new_key][:], axis = (1, 2))
    else:
        data_800m['time']      = np.concatenate((data_800m['time'], zi_nc.variables['time'][:]), axis = 0)
        data_800m['mean_temp'] = np.concatenate((data_800m['mean_temp'], np.nanmean(bouy_nc.variables[temp_key][:], axis = (1, 2, 3))), axis = 0)
        data_800m['pwat']      = np.concatenate((data_800m['pwat'], np.nanmean(np.array([integrate.trapz(x = z_theta, y = interpolate.interp1d(x = z_rho, y = fluxes_nc.variables[rho_key][it,:,:,:], fill_value = 'extrapolate', axis = 0)(z_theta)*mr_nc.variables[mv_key][it,:,:,:], axis = 0) for it in range(mr_nc.variables[mv_key].shape[0])]), axis = (1, 2))), axis = 0)
        if l_CAPE:
            CAPE_temp = np.zeros_like(bouy_nc.variables[temp_key][:,0,0,0])
            CIN_temp  = np.zeros_like(bouy_nc.variables[temp_key][:,0,0,0])
            for it in range(bouy_nc.variables[temp_key].shape[0]):
                CAPE_temp[it], CIN_temp[it], T_parcel, P_parcel, LCL_pressure, LFC_pressure = getCAPE(np.nanmean(bouy_nc.variables[temp_key][it,:-1,:,:], axis = (1, 2)), np.nanmean(bouy_nc.variables[q_key][it,:-1,:,:], axis = (1, 2)), np.nanmean(fluxes_nc.variables[pthe_key][it,:-1,:,:], axis = (1, 2)), parcel_type = 1, t_units = 'K', q_units = 'kg/kg', p_units = 'Pa')
            data_800m['CAPE'] = np.concatenate((data_800m['CAPE'], CAPE_temp), axis = 0)
            data_800m['CIN']  = np.concatenate((data_800m['CIN'], CIN_temp), axis = 0)
        data_800m['zi']        = np.concatenate((data_800m['zi'], np.nanmean(zi_nc.variables[zi_new_key][:], axis = (1, 2))), axis = 0)
    bouy_nc.close()
    fluxes_nc.close()
    mr_nc.close()
    zi_nc.close()

# Creating time series of four different metrics to demonstrate drift in our 
# simulations
# Convert time from minutes to hours
data_100m['time'] /= 60.
data_800m['time'] /= 60.

print 'Creating the plot'
fig = plt.figure(figsize = (10, 10), tight_layout = True)
axa = fig.add_subplot(2, 2, 1)
axa.set_xlim([0, 24])
axa.set_xticks(range(0, 24, 3))
axa.set_xticklabels([''])

axb = fig.add_subplot(2, 2, 2)
axb.set_xlim([0, 24])
axb.set_xticks(range(0, 24, 3))
axb.set_xticklabels([''])

axc = fig.add_subplot(2, 2, 3)
axc.set_xlim([0, 24])
axc.set_xticks(range(0, 24, 3))

axd = fig.add_subplot(2, 2, 4)
axd.set_xlim([0, 24])
axd.set_xticks(range(0, 24, 3))

axa.plot(data_100m['time'], data_100m['mean_temp'], label = '100 m w/ island')
axa.set_ylim([276, 277])
axa.set_title('Mean Temperature (K)')
axa.text(2, 276.9, 'a)')

axb.plot(data_100m['time'], data_100m['pwat'])
axb.text(2, 43.25, 'b)')
axb.set_title('Column Integrated Water Vapour (g m$^{-2}$)')
if l_CAPE:
    axc.plot(data_100m['time'], data_100m['CAPE'])
    axc.set_title('CAPE (J kg$^{-1}$)')
    axc.text(2, 3680, 'c)')

axd.plot(data_100m['time'], data_100m['zi'])
axd.set_title('Boundary Layer Height (m)')
axd.text(2, 810, 'd)')

axa.plot(data_800m['time'], data_800m['mean_temp'], label = '800 m w/o island')
axa.legend(loc = 4, frameon = False)
axb.plot(data_800m['time'], data_800m['pwat'])
if l_CAPE:
    axc.plot(data_800m['time'], data_800m['CAPE'])

axd.plot(data_800m['time'], data_800m['zi'])
plt.savefig('../drift_timeseries_100vs800.png', dpi = 150, bbox_inches = 'tight')
plt.show()

