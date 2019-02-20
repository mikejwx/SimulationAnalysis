"""
Plot the diurnal cycle of temperature, dew point, relatice humidity, wind speed,
and wind direction over the island.
"""
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from SkewT_archer import getDew
# get a land mask
landseamask = Dataset('/work/n02/n02/xb899100/island_masks/lsm50.nc', 'r')

lsm = landseamask.variables['lsm'][0,0,:,:]*1.
y = landseamask.variables['latitude'][:]*1.
x = landseamask.variables['longitude'][:]*1.

# create 2D coordinate mesh
X, Y = np.meshgrid(x, y)

landseamask.close()

# read the data and store the mean 2 m temperature and 10 m winds
hours = ["{0:02d}".format(h) for h in xrange(0, 24, 3)]

temp_key = u'STASH_m01s16i004'
p_key = u'STASH_m01s00i408'
q_key = u'STASH_m01s00i010'
u_key = u'STASH_m01s00i002'
v_key = u'STASH_m01s00i003'

for hour in hours:
    # Read the data for the 
    bouy_nc = Dataset('../bouy_'+hour+'.nc', 'r')
    wind_nc = Dataset('../wind_'+hour+'.nc', 'r')
    fluxes_nc = Dataset('../fluxes_'+hour+'.nc', 'r')
    
    if hour == '00':
        z = bouy_nc.variables['thlev_zsea_theta'][:]
        iz2m   = np.where(np.abs(z - 2.) == np.min(np.abs(z - 2.)))[0][0]
        iz10m  = np.where(np.abs(z - 10.) == np.min(np.abs(z - 10.)))[0][0]
        times  = bouy_nc.variables['min10_0'][:]
        # Initialise our arrays to store the diurnal cycle
        T2m = np.nanmean(np.where((lsm == 1.), bouy_nc.variables[temp_key][1:,iz2m,:,:], np.nan), axis = (1, 2))
        Td2m = np.nanmean(np.where((lsm == 1.), getDew(bouy_nc.variables[q_key][1:,iz2m,:,:], fluxes_nc.variables[p_key][1:,iz2m,:,:]/100.), np.nan), axis = (1, 2))
        spd10m = np.nanmean(np.where((lsm == 1.), np.sqrt(wind_nc.variables[u_key][1:,iz10m,:,:]**2. + wind_nc.variables[v_key][1:,iz10m,:,:]**2.), np.nan), axis = (1, 2))
        times = bouy_nc.variables['min10_0'][1:]
    else:
        T2m = np.concatenate((T2m, np.nanmean(np.where((lsm == 1.), bouy_nc.variables[temp_key][1:,iz2m,:,:], np.nan), axis = (1, 2))), axis = 0)
        Td2m = np.concatenate((Td2m, np.nanmean(np.where((lsm == 1.), getDew(bouy_nc.variables[q_key][1:,iz2m,:,:], fluxes_nc.variables[p_key][1:,iz2m,:,:]/100.), np.nan), axis = (1, 2))), axis = 0)
        spd10m = np.concatenate((spd10m, np.nanmean(np.where((lsm == 1.), np.sqrt(wind_nc.variables[u_key][1:,iz10m,:,:]**2. + wind_nc.variables[v_key][1:,iz10m,:,:]**2.), np.nan), axis = (1, 2))), axis = 0)
        times = np.concatenate((times, bouy_nc.variables['min10_0'][1:]), axis = 0)
    bouy_nc.close()
    wind_nc.close()
    fluxes_nc.close()

# Convert times from minutes to hours
times /= 60.

fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(1, 2, 1)
ax.plot(times, T2m - 273.15, 'r', lw = 2)
ax.plot(times, Td2m - 273.15, 'g', lw = 2)
ax.set_xlabel('Time (hrs)')
ax.set_xticks(range(0, 25, 3))
ax.set_xlim([0, 24])
ax.set_ylabel('Temperature ($^{\circ}$C)')
ax.set_title('2m Temperature/Dew Point averaged over island')

ax = fig.add_subplot(1, 2, 2)
ax.plot(times, spd10m, 'k', lw = 2)
ax.set_xlabel('Time (hrs)')
ax.set_xticks(range(0, 25, 3))
ax.set_xlim([0, 24])
ax.set_ylabel('Wind Speed (m s$^{-1}$)')
ax.set_title('10m wind speed averaged over island')
plt.show()

lwp_key = u'STASH_m01s30i405'
lwp_nc = Dataset('../lwp_00.nc', 'r')
lwp_island = np.nanmean(np.where((lsm == 1.), lwp_nc.variables[lwp_key][1:,:,:], np.nan), axis = (1, 2))*1000.
lwp_times = lwp_nc.variables['min5_0'][1:]/60.
lwp_nc.close()

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(lwp_times, lwp_island, 'k', lw = 2)
ax.set_xticks(range(0, 24, 3))
ax.set_xlim([0, 24])
ax.set_xlabel('Time (hrs)')
ax.set_ylabel(r'LWP (g kg$^{-1}$)')
ax.set_title('LWP averaged over island')

plt.show()
