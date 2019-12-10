import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import q_key, temp_key, pthe_key
from SkewT_archer import getQ, Lv
from scipy import integrate

# create the variable storage infrastructure, using dictionaries
my_RH_data = {}

paths = ['/nerc/n02/n02/xb899100/CloudTrail/RH_BLm25/', '/nerc/n02/n02/xb899100/CloudTrail/RH_FAm25/', '/nerc/n02/n02/xb899100/CloudTrail/Control_short/']
for path in paths:
    path_key = path.split('/')[-2]
    my_RH_data[path_key] = {}
    with Dataset(path + 'bouy_00.nc', 'r') as bouy_nc:
        print 'reading bouy_nc for ' + path_key
        my_RH_data[path_key][q_key]    = bouy_nc.variables[q_key][1,:-1,:,:].mean(axis = (1,2))
        my_RH_data[path_key][temp_key] = bouy_nc.variables[temp_key][1,:-1,:,:].mean(axis = (1,2))
        z = bouy_nc.variables['thlev_zsea_theta'][:-1]*1.
    
    with Dataset(path + 'fluxes_00.nc', 'r') as fluxes_nc:
        my_RH_data[path_key][pthe_key] = fluxes_nc.variables[pthe_key][1,:-1,:,:].mean(axis = (1,2))
    
    my_RH_data[path_key]['q_sat'] = getQ(my_RH_data[path_key][temp_key], [100.], my_RH_data[path_key][pthe_key], t_units = 'K', p_units = 'Pa')
    my_RH_data[path_key]['RH']    = 100.*my_RH_data[path_key][q_key]/my_RH_data[path_key]['q_sat']
    
# compute delta column water vapours
CIWV_RHmBL25 = integrate.trapz(x = z, y = my_RH_data['RH_BLm25'][q_key])
CIWV_RHmFA25 = integrate.trapz(x = z, y = my_RH_data['RH_FAm25'][q_key])
CIWV_Control_Short = integrate.trapz(x = z, y = my_RH_data['Control_short'][q_key])

lhf = 167.2327 #W/m2
rho_sfc = 1.1771 #kg/m3
lhf = lhf/Lv/rho_sfc

dt = ((CIWV_Control_Short - CIWV_RHmBL25)/lhf)/86400.
print 'BLm25: ' + str(round(dt, 1)) + ' days'
dt = ((CIWV_Control_Short - CIWV_RHmFA25)/lhf)/86400.
print 'FAm25: ' + str(round(dt, 1)) + ' days'
