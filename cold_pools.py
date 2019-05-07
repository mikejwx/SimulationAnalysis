"""
Cold pool phase speeds for our U05 and Control_short simulation

Comparison to the background flow
"""
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import interpolate, integrate
from STASH_keys import theta_key, mv_key, mcl_key, mr_key, u_key, v_key, w_key
from SkewT_archer import g
#plt.switch_backend('agg')
from analysis_tools import send_email

# Function to compute the cold pool phase speed
def C(theta, theta_bar, qv, qv_bar, qcl, qr, z):
    """
    Computes the vertical integral of buoyancy where surface is negatively
    buoyant. The square root of this is the cold pool phase speed which is 
    returned by this function.
    ----------------------------------------------------------------------------
    INPUT:
    theta     = 4D array of potential temperature (K) with theta[t, z, y, x]
    theta_bar = 1D array of the potential temperature base state (K) with theta[z]
    qv        = 4D array of specific humidity (kg/kg) with qv[t,z,y,x]
    qv_bar    = 1D array of the specific humidity base state (kg/kg) with qv[z]
    qcl       = 4D array of the cloud liquid specific humidity (kg/kg) with qcl[t,z,y,x]
    qr        = 4D array of the rain water specific humidity (kg/kg) with qr[t,z,y,x]
    z         = 1D height coordinates
    
    OUTPUT:
    C         = 3D array of the cold pool phase speed C[t,y,x]
                N.B. where there is no negative surface buoyancy, np.nan is used
                as the fill_value
    """
    B = np.zeros_like(theta)
    C2 = np.zeros_like(theta[:,0,:,:])
    for it in xrange(theta.shape[0]):
        print 'C -> Working on time step ' + str(it)
        # Compute the vertical buoyancy profile for each time step
        B[it,:,:,:] = g*(np.transpose((np.transpose(theta[it,:,:,:]) - theta_bar)/theta_bar) + 0.61*np.transpose(np.transpose(qv[it,:,:,:]) - qv_bar) - qcl[it,:,:,:] - qr[it,:,:,:])
        for j in xrange(theta.shape[2]):
            for i in xrange(theta.shape[3]):
                # For each horizontal coordinate...
                if B[it,0,j,i] < 0:
                    # If there is a negative surface buoyancy anomaly, assume it is a cold pool
                    # Find the index for the depth of the cold pool
                    H_idx = 0
                    while B[it,H_idx,j,i] <= 0:
                        H_idx += 1
                    
                    # Integrate to the height where B increases to 0
                    C2[it,j,i] = 2*integrate.trapz(y = -B[it,:H_idx,j,i], x = z[:H_idx])
                else:
                    C2[it,j,i] = np.nan
    C = np.sqrt(C2)
    return C

# Define paths to data
REF_path = '/nerc/n02/n02/xb899100/CloudTrail/Control_short/'
U05_path = '/nerc/n02/n02/xb899100/CloudTrail/U05/'
l_testing = True

if l_testing:
    hour = '08'
    print 'Working on hour ' + hour
    # Open data for the Control_short experiment which did not have significant cold pools
    bouy_REF_nc = Dataset(REF_path + 'bouy_' + hour + '.nc', 'r')
    mr_REF_nc   = Dataset(REF_path + 'mr_' + hour + '.nc', 'r')
    wind_REF_nc = Dataset(REF_path + 'wind_' + hour + '.nc', 'r')
    
    # Open data for the U05 experiment which did have significant cold pools
    bouy_U05_nc = Dataset(U05_path + 'bouy_' + hour + '.nc', 'r')
    mr_U05_nc   = Dataset(U05_path + 'mr_' + hour + '.nc', 'r')
    wind_U05_nc = Dataset(U05_path + 'wind_' + hour + '.nc', 'r')
    
    print 'Reading the data'
    # Get the height coordinate
    z = bouy_REF_nc.variables['thlev_zsea_theta'][:]
    
    # Read data for REF
    theta_REF = bouy_REF_nc.variables[theta_key][:]
    mv_REF    = mr_REF_nc.variables[mv_key][:]
    mcl_REF   = mr_REF_nc.variables[mv_key][:]*0
    mr_REF    = mr_REF_nc.variables[mv_key][:]*0
    u_REF     = wind_REF_nc.variables[u_key][:]
    v_REF     = wind_REF_nc.variables[v_key][:]
    w_REF     = wind_REF_nc.variables[w_key][:]
    
    # Read data for U05
    theta_U05 = bouy_U05_nc.variables[theta_key][:]
    mv_U05    = mr_U05_nc.variables[mv_key][:]
    mcl_U05   = mr_U05_nc.variables[mcl_key][:]
    mr_U05    = mr_U05_nc.variables[mr_key][:]
    u_U05     = wind_U05_nc.variables[u_key][:]
    v_U05     = wind_U05_nc.variables[v_key][:]
    w_U05     = wind_U05_nc.variables[w_key][:]
    
    # Get the background state for REF
    theta_bar_REF = np.nanmean(theta_REF, axis = (2, 3))
    mv_bar_REF    = np.nanmean(mv_REF, axis = (2, 3))
    
    # Get the background state for U05
    theta_bar_U05 = np.nanmean(theta_U05, axis = (2, 3))
    mv_bar_U05    = np.nanmean(mv_U05, axis = (2, 3))
    
    # Define the buoyancy variables for each experiment
    theta_p_REF = np.zeros_like(theta_REF)
    theta_p_U05 = np.zeros_like(theta_U05)
    B_REF = np.zeros_like(theta_REF)
    B_U05 = np.zeros_like(theta_U05)
    for j in xrange(theta_REF.shape[2]):
        for i in xrange(theta_REF.shape[3]):
            theta_p_REF[:,:,j,i] = (theta_REF[:,:,j,i] - theta_bar_REF)
            theta_p_U05[:,:,j,i] = (theta_U05[:,:,j,i] - theta_bar_U05)
            B_REF[:,:,j,i] = (theta_REF[:,:,j,i] - theta_bar_REF)/theta_bar_REF + 0.61*(mv_REF[:,:,j,i] - mv_bar_REF) - mcl_REF[:,:,j,i] - mr_REF[:,:,j,i]
            B_U05[:,:,j,i] = (theta_U05[:,:,j,i] - theta_bar_U05)/theta_bar_U05 + 0.61*(mv_U05[:,:,j,i] - mv_bar_U05) - mcl_U05[:,:,j,i] - mr_U05[:,:,j,i]
    
# Plot some vertical profiles through some particularly cold surface spots
it, iy, ix = np.where(B_REF[:,0,:,:] == np.min(B_REF[:,0,:,:]))
it, iy, ix = [it[0], iy[0], ix[0]]
fig = plt.figure()
ax0 = fig.add_subplot(2, 2, 1)
ax0.plot(theta_p_REF[it,:,iy,ix], z)
ax0.set_title('$\\theta^{\prime}_{REF}$')

# Get the density current propagation speed
iz = 0
while B_REF[it,iz,iy,ix] < 0:
    iz += 1

C2 = - 2*integrate.trapz(B_REF[it,:iz,iy,ix], z[:iz])
C_REF = np.sqrt(C2)

ax1 = fig.add_subplot(2, 2, 2)
ax1.plot(B_REF[it,:,iy,ix], z)
ax1.set_title('B$_{REF}$, C$_{REF}$ = ' + str(round(C_REF, 1)) + ' m/s')

it, iy, ix = np.where(B_U05[:,0,:,:] == np.min(B_U05[:,0,:,:]))
it, iy, ix = [it[0], iy[0], ix[0]]
ax2 = fig.add_subplot(2, 2, 3)
ax2.plot(theta_p_U05[it,:,iy,ix], z)
ax2.set_title('$\\theta^{\prime}_{U05}$')

# Get the density current propagation speed
iz = 0
while B_U05[it,iz,iy,ix] < 0:
    iz += 1

C2 = - 2*integrate.trapz(B_U05[it,:iz,iy,ix], z[:iz])
C_U05 = np.sqrt(C2)

ax3 = fig.add_subplot(2, 2, 4)
ax3.plot(B_U05[it,:,iy,ix], z)
ax3.set_title('B$_{U05}$, C$_{U05}$ = ' + str(round(C_U05, 1)) + ' m/s')
plt.show()
    
else:
    hours = ["{0:02d}".format(h) for h in xrange(0, 16, 4)]
    for hour in hours:
        print 'Working on hour ' + hour
        # Open data for the Control_short experiment which did not have significant cold pools
        bouy_REF_nc = Dataset(REF_path + 'bouy_' + hour + '.nc', 'r')
        mr_REF_nc   = Dataset(REF_path + 'mr_' + hour + '.nc', 'r')
        wind_REF_nc = Dataset(REF_path + 'wind_' + hour + '.nc', 'r')
        
        # Open data for the U05 experiment which did have significant cold pools
        bouy_U05_nc = Dataset(U05_path + 'bouy_' + hour + '.nc', 'r')
        mr_U05_nc   = Dataset(U05_path + 'mr_' + hour + '.nc', 'r')
        wind_U05_nc = Dataset(U05_path + 'wind_' + hour + '.nc', 'r')
        
        if hour == hours[0]:
            print 'Reading the data'
            z = bouy_REF_nc.variables['thlev_zsea_theta'][:]
            
            # Read data for REF
            times_key_REF = [key for key in bouy_REF_nc.variables.keys() if 'min' in key][0]
            theta_REF = bouy_REF_nc.variables[theta_key][1:,:,:,:]
            mv_REF    = mr_REF_nc.variables[mv_key][1:,:,:,:]
            mcl_REF   = mr_REF_nc.variables[mv_key][1:,:,:,:]*0
            mr_REF    = mr_REF_nc.variables[mv_key][1:,:,:,:]*0
            u_REF     = wind_REF_nc.variables[u_key][1:,:,:,:]
            v_REF     = wind_REF_nc.variables[v_key][1:,:,:,:]
            w_REF     = wind_REF_nc.variables[w_key][1:,:,:,:]
            times_REF = bouy_REF_nc.variables[times_key_REF][:]

            # Read data for U05
            times_key_U05 = [key for key in bouy_U05_nc.variables.keys() if 'min' in key][0]
            theta_U05 = bouy_U05_nc.variables[theta_key][1:,:,:,:]
            mv_U05    = mr_U05_nc.variables[mv_key][1:,:,:,:]
            mcl_U05   = mr_U05_nc.variables[mcl_key][1:,:,:,:]
            mr_U05    = mr_U05_nc.variables[mr_key][1:,:,:,:]
            u_U05     = wind_U05_nc.variables[u_key][1:,:,:,:]
            v_U05     = wind_U05_nc.variables[v_key][1:,:,:,:]
            w_U05     = wind_U05_nc.variables[w_key][1:,:,:,:]
            times_U05 = bouy_U05_nc.variables[times_key_U05][:]

            # Get the background state
            theta_bar_REF = np.nanmean(theta_REF, axis = (2, 3))
            mv_bar_REF    = np.nanmean(mv_REF, axis = (2, 3))

            # Get the background state
            theta_bar_U05 = np.nanmean(theta_U05, axis = (2, 3))
            mv_bar_U05    = np.nanmean(mv_U05, axis = (2, 3))
            
            ### Do analysis ###
            print 'Computing the speed'
            C_REF = C(theta_REF, theta_bar_REF, mv_REF, mv_bar_REF, mcl_REF, mr_REF, z)
            C_U05 = C(theta_U05, theta_bar_U05, mv_U05, mv_bar_U05, mcl_U05, mr_U05, z)
        else:
            print 'Reading the data'
            # Read data for REF
            theta_REF = bouy_REF_nc.variables[theta_key][:]
            mv_REF    = mr_REF_nc.variables[mv_key][:]
            mcl_REF   = mr_REF_nc.variables[mv_key][:]*0
            mr_REF    = mr_REF_nc.variables[mv_key][:]*0
            u_REF     = wind_REF_nc.variables[u_key][:]
            v_REF     = wind_REF_nc.variables[v_key][:]
            w_REF     = wind_REF_nc.variables[w_key][:]
            times_REF = np.concatenate((times_REF, bouy_REF_nc.variables[times_key_REF][:]), axis = 0)
            
            # Read data for U05
            theta_U05 = bouy_U05_nc.variables[theta_key][:]
            mv_U05    = mr_U05_nc.variables[mv_key][:]
            mcl_U05   = mr_U05_nc.variables[mcl_key][:]
            mr_U05    = mr_U05_nc.variables[mr_key][:]
            u_U05     = wind_U05_nc.variables[u_key][:]
            v_U05     = wind_U05_nc.variables[v_key][:]
            w_U05     = wind_U05_nc.variables[w_key][:]
            times_U05 = np.concatenate((times_U05, bouy_U05_nc.variables[times_key_U05][:]), axis = 0)
            
            # Get the background state
            theta_bar_REF = np.nanmean(theta_REF, axis = (2, 3))
            mv_bar_REF    = np.nanmean(mv_REF, axis = (2, 3))
            
            # Get the background state
            theta_bar_U05 = np.nanmean(theta_U05, axis = (2, 3))
            mv_bar_U05    = np.nanmean(mv_U05, axis = (2, 3))
            
            ### Do analysis ###
            print 'Computing the speed'
            C_REF = np.concatenate((C_REF, C(theta_REF, theta_bar_REF, mv_REF, mv_bar_REF, mcl_REF, mr_REF, z)), axis = 0)
            C_U05 = np.concatenate((C_U05, C(theta_U05, theta_bar_U05, mv_U05, mv_bar_U05, mcl_U05, mr_U05, z)), axis = 0)
        
        # Close netCDF
        bouy_REF_nc.close()
        mr_REF_nc.close()
        wind_REF_nc.close()
        bouy_U05_nc.close()
        mr_U05_nc.close()
        wind_U05_nc.close()
        
        C_REF_ts = np.nanmean(C_REF, axis = (1, 2))
        C_U05_ts = np.nanmean(C_U05, axis = (1, 2))
    
    ### Do plotting ###
    print 'Plotting'
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    if (len(times_REF) == len(C_REF_ts)) and (len(times_U05) == len(C_U05_ts)):
        ax.plot(times_REF+240., C_REF_ts, color = 'k', marker = 'o', markerfacecolor = 'white', label = 'REF')
        ax.plot([240., 1200.], [10., 10.], color = 'grey', ls = '--', label = 'REF:spd')
        ax.plot([240., 1200.], [5., 5.], color = 'blue', ls = '--', label = 'U05:spd')
        ax.plot(times_U05+240., C_U05_ts, color = 'k', marker = 'o', markerfacecolor = 'k', label = 'U05')
    else:
        ax.plot(np.linspace(240., 1200., len(C_REF_ts)), C_REF_ts, color = 'k', marker = 'o', markerfacecolor = 'white', label = 'REF')
        ax.plot(np.linspace(240., 1200., len(C_U05_ts)), C_U05_ts, color = 'k', marker = 'o', markerfacecolor = 'k', label = 'U05')
        ax.plot([240., 1200.], [10., 10.], color = 'grey', ls = '--', label = 'REF:spd')
        ax.plot([240., 1200.], [5., 5.], color = 'blue', ls = '--', label = 'U05:spd')
    plt.legend(loc = 0)
    plt.savefig('../cold_pool_stuff.png')
    
    send_email(message = 'Cold pools is complete.', subject = 'cold_pools.py -> complete', attachments = ['../cold_pool_stuff.png'], isAttach = True)



