### Find the cloud weighted wind direction ###
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import interpolate, integrate
from analysis_tools import toComponents, fromComponents, regrid

days = ["{0:02d}".format(x) for x in xrange(1, 11)]
ndays = len(days)
u_key = u'STASH_m01s00i002'
v_key = u'STASH_m01s00i003'
q_key = u'STASH_m01s00i392'
for day in days:
    print 'Starting day ' + day
    u = Dataset('u_'+day+'.nc', 'r')
    v = Dataset('v_'+day+'.nc', 'r')
    mr = Dataset('mr_'+day+'.nc', 'r')
    
    if day == '01':
        dt_i = u.variables[u_key].shape[0]
        u_cld = np.zeros(dt_i*ndays)
        v_cld = np.zeros_like(u_cld)
        z = mr.variables['thlev_zsea_theta'][:]*1.
    # get the cloud liquid water mixing ratio, u- and v- wind components
    cld = mr.variables[q_key][:]*1.
    u_rg = regrid(mr, u, u_key)
    v_rg = regrid(mr, v, v_key)
    ## weight wind components by the cloud liquid water mixing ratio
    # 1. multiply by the cld liquid water mixing ratio
    u_weight = u_rg*cld
    v_weight = v_rg*cld
    # 2. vertically integrate the weighted profiles
    u_weight = integrate.trapz(u_weight, z, axis = 1)
    v_weight = integrate.trapz(v_weight, z, axis = 1)
    # 3. divide by the vertically integrated cloud liquid water mixing ratio
    total_cld = integrate.trapz(cld, z, axis = 1)
    u_weight /= total_cld
    v_weight /= total_cld
    # 4. take the horizontal mean
    u_weight = np.nanmean(np.where((total_cld > 0), u_weight, np.nan), axis = (1, 2))
    v_weight = np.nanmean(np.where((total_cld > 0), v_weight, np.nan), axis = (1, 2))
    # 5. store to u_cld and v_cld
    u_cld[dt_i*days.index(day):dt_i*(days.index(day)+1)] = u_weight
    v_cld[dt_i*days.index(day):dt_i*(days.index(day)+1)] = v_weight
    u.close()
    v.close()
    mr.close()
# convert from components to wind direction in degrees
wind_direction = fromComponents(u_cld, v_cld, 1)[1]

plt.plot(wind_direction)
plt.show()
print 'The mean wind direction over the last 4 days is: ' + str(round(np.mean(wind_direction[-4*dt_i:]), 1))

