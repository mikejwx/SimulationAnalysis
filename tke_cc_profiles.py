import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import interpolate, integrate
import os
"""
Code to read in the u, v, w, and mcl data and plot vertical profiles of resolved
scale turbulent kinetic energy and cloud cover. These plots can be made at
several times during the duration of the simulation, or for mean over given
time intervals.
"""

def horizontalMean(data, height = True):
    """
    Function to compute and return the horizontal mean of a given array.
    Will check if the array is 3- or 4- dimensional. If 3D, then assume that 
    array = f(z, y, x), and if 4D, assume array = f(t, z, y, x).
    
    If the logical 'height' is true, then input arrays are assumed to be array 
    = f(t, y, x), and time series of horizontal mean for single level are 
    returned.
    """
    
    if (len(data.shape) == 3) and height:
        h_mean = np.zeros(data.shape[0])
        for k in xrange(data.shape[0]):
            h_mean[k] = np.mean(data[k,:,:])
    elif (len(data.shape) == 4) and height:
        h_mean = np.zeros((data.shape[0], data.shape[1]))
        for t in xrange(data.shape[0]):
            for k in xrange(data.shape[1]):
                h_mean[t, k] = np.mean(data[t, k, :, :])
    elif (len(data.shape) == 3) and not height:
        h_mean = np.zeros(data.shape[0])
        for t in xrange(data.shape[0]):
            h_mean[t] = np.mean(data[t,:,:])
    else:
        raise Exception("Unrecognised array shape in horizontalMean")
    
    return h_mean

def get_TKE(U, V, W, rho = 1., start = 0, end = None):
    """
    Function to compute the vertical profile of turbulent kinetic energy from 
    given 4D U, V, and W fields.
    
    First, calculates the perturbations from horizontal mean profiles of U, V, 
    and W. Then, calculates the turbulent kinetic energy, and then the 
    horizontal mean so that there is a time series of profiles.
    
    U, V, and W must all be on the same grid.
    """
    U_mean = horizontalMean(U)
    V_mean = horizontalMean(V)
    W_mean = horizontalMean(W)
    U_p = np.zeros_like(U)
    V_p = np.zeros_like(V)
    W_p = np.zeros_like(W)
    for t in xrange(U.shape[0]):
        for k in xrange(U.shape[1]):
            U_p[t,k,:,:] = U[t,k,:,:] - U_mean[t,k]
            V_p[t,k,:,:] = V[t,k,:,:] - V_mean[t,k]
            W_p[t,k,:,:] = W[t,k,:,:] - W_mean[t,k]
    if type(rho) != float:
        TKE_profile = 0.5*horizontalMean(rho)*(horizontalMean(U_p**2)+ horizontalMean(V_p**2) + horizontalMean(W_p**2))
    else:
        TKE_profile = 0.5*rho*(horizontalMean(U_p**2)+ horizontalMean(V_p**2) + horizontalMean(W_p**2))
    
    return TKE_profile

def get_CC(m_cl, start = 0, end = None, threshold = 0.):
    """
    Function to compute the vertical profile of cloud cover from given 4D m_cl 
    field.
    
    First, finds the number of cloudy pixels, i.e. those that contain non-zero 
    m_cl. This is a binary, 0 = no cloud, 1 = cloud. Then calculates the cloud 
    cover on each vertical level as the mean of the binary field at each level. 
    Finally, does a time mean over a given period. Returns a vertical profile of
    this cloud cover estimate.
    """
    cloud_mask = m_cl[:]
    cloud_mask[m_cl > threshold] = 1.
    cloud_mask[m_cl <= threshold] = 0.
    cc_profile = np.mean(horizontalMean(cloud_mask)[start:end,:], axis = 0)
    
    return cc_profile

def regrid(target_grid, current_grid, current_data_key):
    """
    Function to attempt to regrid an already read-in netCDF B (current_grid) to 
    the same grid as already read-in netCDF A (target_grid).
    Should only work if there is some overlap between the range of the functions
    otherwise the interpolation won't make sense.
    
    Regrids in all four dimensions, although going from e.g. 1 minute data to
    10 or 15 minute data shouldn't really require regridding. But going from
    10 minute to 1 minute data definitely requires regridding and the linear 
    interpolation approach used here is probably not appropriate.
    """
    # Read the data that contains the target grid
    target_x_key     = [key for key in target_grid.variables.keys() if ('longitude' in key) and ('_l' not in key)][0]
    target_x = target_grid.variables[target_x_key][:]
    target_y_key     = [key for key in target_grid.variables.keys() if ('latitude' in key) and ('_l' not in key)][0]
    target_y = target_grid.variables[target_y_key][:]
    target_z_key     = [key for key in target_grid.variables.keys() if ('zsea' in key) and ('bounds' not in key)][0]
    target_z = target_grid.variables[target_z_key][:]
    target_t_key     = [key for key in target_grid.variables.keys() if 'min' in key][0]
    target_t = target_grid.variables[target_t_key][:]
    
    # Read the netCDF that needs to be transformed
    current_x_key     = [key for key in current_grid.variables.keys() if ('longitude' in key) and ('_l' not in key)][0]
    current_x = current_grid.variables[current_x_key][:]
    current_y_key     = [key for key in current_grid.variables.keys() if ('latitude' in key) and ('_l' not in key)][0]
    current_y = current_grid.variables[current_y_key][:]
    current_z_key     = [key for key in current_grid.variables.keys() if ('zsea' in key) and ('bounds' not in key)][0]
    current_z = current_grid.variables[current_z_key][:]
    current_t_key     = [key for key in current_grid.variables.keys() if 'min' in key][0]
    current_t = current_grid.variables[current_t_key][:]
    current_data = current_grid.variables[current_data_key][:]
     
    ### Procedure ###
    # Will usually only need to interpolate in z, sometimes in x or y, probably
    # never required to interpolate in more than two dimensions
    # 1. interpolate in time
    if target_t_key != current_t_key:
        print 'Regridding in time...'
        # Interpolate in time
        current_data = interpolate.interp1d(current_t, current_data, axis = 0, fill_value = 'extrapolate')(target_t)
        print 'Complete.'
    # 2. interpolate in z
    if target_z_key != current_z_key:
        print 'Regridding in height...'
        # Interpolate in z
        current_data = interpolate.interp1d(current_z, current_data, axis = 1, fill_value = 'extrapolate')(target_z)
        print 'Complete.'
    # 3. interpolate in y
    if target_y_key != current_y_key:
        print 'Regridding in latitude...'
        # Interpolate in y
        for i in xrange(current_data.shape[3]):
            current_data[:,:,:,i] = interpolate.interp1d(current_y[:,i], current_data[:,:,:,i], axis = 2, fill_value = 'extrapolate')(target_y[:,i])
        print 'Complete.'
    # 4. interpolate in x
    if target_x_key != current_x_key:
        print 'Regridding in longitude...'
        # Interpolate in x
        for j in xrange(current_data.shape[2]):
            current_data[:,:,j,:] = interpolate.interp1d(current_x[j,:], current_data[:,:,j,:], axis = 2, fill_value = 'extrapolate')(target_x[j,:])
        print 'Complete.'
    
    return current_data

# Read in the data
days = ["{0:02d}".format(day) for day in xrange(1,11)]
ndays = len(days)

# initialise a tke time series variable
tke_ts = np.zeros(144*ndays)
mcl = np.zeros((144*ndays, 141))
mcf = np.zeros((144*ndays, 141))

for day in days:
    print 'Open the netCDF on day ' + day
    u_nc = Dataset('u_'+day+'.nc', 'r')
    u_key   = 'STASH_m01s00i002'
    v_nc = Dataset('v_'+day+'.nc', 'r')
    v_key   = 'STASH_m01s00i003'
    w_nc = Dataset('bouy_'+day+'.nc', 'r')
    w_key   = 'STASH_m01s00i150'
    mr_nc = Dataset('mr_'+day+'.nc', 'r')
    mcl_key = 'STASH_m01s00i392'
    mcf_key = 'STASH_m01s00i393'
    rho_nc = Dataset('fluxes_'+day+'.nc', 'r')
    rho_key = 'STASH_m01s00i389'
    if day == '01':
        print 'Grabbing the height coordinates'
        z_rho = u_nc.variables['rholev_zsea_rho'][:]
        z_theta = w_nc.variables['thlev_zsea_theta'][:]
        iz_3000 = np.where(abs(z_theta - 3000) == np.min(abs(z_theta - 3000)))[0][0]+1


    # u and v need to be regridded onto theta levels
    print 'Regridding the winds...'
    u_regrid = regrid(w_nc, u_nc, u_key)
    u_nc.close()
    v_regrid = regrid(w_nc, v_nc, v_key)
    v_nc.close()
    rho_regrid = regrid(w_nc, rho_nc, rho_key)
    rho_nc.close()
    
    TKE = get_TKE(u_regrid, v_regrid, w_nc.variables[w_key][:], rho_regrid)
    for it in xrange(TKE.shape[0]):
        tke_ts[(int(day)-1)*144 + it] = integrate.trapz(y = TKE[it,:iz_3000], x = z_theta[:iz_3000])
    
    tke_profile = np.mean(get_TKE(u_regrid, v_regrid, w_nc.variables[w_key][:]), axis = 0)
    
    w_nc.close()
    
    cc_profile = get_CC(mr_nc.variables[mcl_key][:]+mr_nc.variables[mcf_key][:])
    mcl[days.index(day)*144:(days.index(day)+1)*144,:] = np.max(mr_nc.variables[mcl_key][:], axis = (2,3))
    mcf[days.index(day)*144:(days.index(day)+1)*144,:] = np.max(mr_nc.variables[mcf_key][:], axis = (2,3))
    mr_nc.close()
    
    print 'Making the plot...'
    fig = plt.figure()
    plt.subplot(121)
    plt.plot(tke_profile, z_theta, 'k', lw = 2, label = 'TKE Profile')
    plt.ylabel('Height (m)')
    plt.xlabel('TKE ($m^2 s^{-2}$)')
    plt.ylim([0, 15000])
    plt.xlim([0, 0.5])
    
    plt.subplot(122)
    plt.plot(cc_profile, z_theta, 'k--', lw = 2, label = 'Cloud Cover Profile')
    plt.xlabel('Cloud Cover')
    plt.xlim([0, 0.1])
    plt.ylim([0, 15000])
    plt.suptitle('Profiles averages over day ' + day, fontsize = 14)
    lgd = plt.legend(loc = 'upper center', bbox_to_anchor = (-0.25, -0.1))
    plt.savefig('tke_cc_profiles_day_' + day + '.png', dpi = 150, bbox_extra_artists=(lgd,), bbox_inches = 'tight')
    plt.close('all')
    #plt.show()

times = np.arange(1., 1440.*ndays, 10.)
plt.plot(times/1440., tke_ts, 'k', lw = 2)
plt.xlabel('Time (days)')
plt.ylabel('Vertically Integrated TKE (kg s$^{-2}$)')
plt.title('Integrated below ' + str(z_theta[iz_3000-1]) + ' m')
plt.savefig('VITKE_10-day.png', dpi = 150)
plt.show()

fig = plt.figure()
fig.set_size_inches(15, 8)
ax = plt.subplot()
CL = ax.contourf(times/1440., z_theta, np.transpose(mcl), levels = np.linspace(1e-08, 5e-03, 21), cmap = 'Greys')
ax.set_yscale('log', basey = 10, nonposy = 'clip', subsy = [2, 3, 4, 5, 6, 7, 8, 9])
ax.set_ylim([100, 40000])
ax.set_xlim([0, 10])
plt.colorbar(mappable = CL, label = 'Cloud Liquid')
CI = ax.contourf(times/1440., z_theta, np.transpose(mcf), levels = np.linspace(1e-08, 5e-05, 21), cmap = 'Blues')
ax.set_yscale('log', basey = 10, nonposy = 'clip', subsy = [2, 3, 4, 5, 6, 7, 8, 9])
ax.set_ylim([100, 40000])
ax.set_xlim([0, 10])
plt.colorbar(CI, label = 'Cloud Ice')
ax.set_ylabel('Height (m)')
ax.set_yticks([100, 500, 1000, 2500, 5000, 10000, 15000, 40000])
ax.set_yticklabels([100, 500, 1000, 2500, 5000, 10000, 15000, 40000])
ax.set_xlabel('Time (days)')
plt.savefig('cloud_timeseries.png', dpi = 150)
plt.show()



