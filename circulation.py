import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

"""
Trying to compare the equations from Kirshbaum and Fairman for the 'seabreeze'
strength. These are empirically derived equations for the boundary layer vertical
velocity and cross flow perturbation in cloud trails.
"""

def V_sb_prediction(H):
    """
    Compute the expected sea breeze flow from a given H for our island dimensions
    which are hard coded here.
    
    Equation is originally from Kirshbaum and Fairman 2015.
    """
    island_area = 50. # km2
    l_x = np.sqrt(island_area/np.pi)*1000. # island radius (m)
    T0 = 300.0 # (K)
    rho = 1.1771 # kg/m3 - taken from the surface values
    cpd = 1005. # J/kg/K
    Rd = 287.05 # J/kg/K
    U_sc = 9.5 # m/s
    g = 9.81 # m/s2
    
    V_sb = 0.78*(g*H*l_x/(rho*cpd*T0*U_sc))**0.5
    
    return V_sb

def W_sb_prediction(H):
    """
    Compute the expected sea breeze vertical velocity from a given H for our
    island dimensions which are hard coded here.
    
    Equation is originally from Kirshbaum and Fairman (2015).
    """
    island_area = 50. # km2
    l_x = np.sqrt(island_area/np.pi)*1000. # island radius (m)
    l_y = l_x*1. # (m)
    z_cb = 700. # (m)
    T0 = 300.0 # (K)
    rho = 1.1771 # kg/m3 - taken from the surface values
    cpd = 1005. # J/kg/K
    Rd = 287.05 # J/kg/K
    U_sc = 9.5 # m/s
    g = 9.81 # m/s2
    
    W_sb = 1.25*(z_cb*H*l_x/(l_y*rho*cpd*T0*U_sc))**0.5
    
    return W_sb

#### First let's explore the time evolution of these parameters in our experiments. ####
# Use the original netcdf
path_pt1 = '/work/n02/n02/xb899100/cylc-run/'
path_pt2 = '/share/data/history/'
experiments = ['control', 'exp01', 'exp02']
suites = {'control' : 'u-bd527', 
          'exp01' : 'u-bg023', 
          'exp02' : 'u-bg113'}
V_sb_sim = {}
W_sb_sim = {}

# Define keys for the wind components
w_key = u'STASH_m01s00i150'
n_key = u'Flow-Perpendicular'
hours = ["{0:02d}".format(h) for h in xrange(0, 24, 3)]

for experiment in experiments:
    for hour in hours:
        # Read the data
        wind_nc = Dataset(path_pt1 + suites[experiment] + path_pt2 + 'wind_' + hour + '.nc', 'r')
        if hour == '00':
            time_key = [i for i in wind_nc.variables.keys() if 'min' in i][0]
            times = wind_nc.variables[time_key][:]
            if experiment == 'control':
                z = wind_nc.variables['thlev_zsea_theta'][:]
                iz_n = np.where(np.abs(z - 10.) == np.min(np.abs(z - 10.)))[0][0]
                iz_w = np.where(np.abs(z - 350.) == np.min(np.abs(z - 350.)))[0][0]
            # Assume the maximum w and n are the values corresponding to the CT
            V_sb_sim[experiment] = [i for i in np.max(wind_nc.variables[n_key][:,iz_n,:,:], axis = (1, 2))]
            W_sb_sim[experiment] = [i for i in np.max(wind_nc.variables[w_key][:,iz_w,:,:], axis = (1, 2))]
        else:
            [V_sb_sim[experiment].append(i) for i in np.max(wind_nc.variables[n_key][:,iz_n,:,:], axis = (1, 2))]
            [W_sb_sim[experiment].append(i) for i in np.max(wind_nc.variables[w_key][:,iz_w,:,:], axis = (1, 2))]
            times = np.concatenate((times, wind_nc.variables[time_key][:]*1.), axis = 0)
        # Want to add in a land sea mask, then plot dots for where the maximum in vertical velocity is occurring (need to collect this from above, likely using np.where), and color them for when the heating is on vs. off.
    #plt.show()    
    """
    range_max = 10.
    fig = plt.figure(tight_layout = True, figsize=(10,10))
    ax = fig.add_subplot(2, 2, 1, adjustable='box', aspect=1)
    # Plot the simulated horizontal ct component timeseries
    ax.plot(times/60., V_sb_sim[experiment], 'b', lw = 2)
    ax.set_xlabel('Time (hrs)')
    ax.set_ylabel('V$_{sb,sim}$ (m s$^{-1}$)')
    ax.plot([6., 6.], [0., range_max], 'k--')
    ax.plot([18., 18.], [0., range_max], 'k--')
    ax.set_xlim([0, 24])
    ax.set_ylim([0, range_max])
    ax.set_xticks(xrange(0, 24, 3))
    ax.set_title('Horizontal CT component for exp=' + experiment)
    # Plot the simulated vertical ct component timeseries
    ax = fig.add_subplot(2, 2, 2, adjustable='box', aspect=1)
    ax.plot(times/60., W_sb_sim[experiment], 'k', lw = 2)
    ax.set_xlabel('Time (hrs)')
    ax.set_ylabel('W$_{sb,sim}$ (m s$^{-1}$)')
    ax.plot([6., 6.], [0., range_max], 'k--')
    ax.plot([18., 18.], [0., range_max], 'k--')
    ax.set_xlim([0, 24])
    ax.set_ylim([0, range_max])
    ax.set_xticks(xrange(0, 24, 3))
    ax.set_title('Vertical CT component for exp=' + experiment)
    # Plot the relationship between the two, we expect W to be directly proportional to V
    ax = fig.add_subplot(2, 2, 3, adjustable='box', aspect=1)
    ax.plot([V_sb_sim[experiment][i] for i in xrange(len(times)) if not (6 <= times[i]/60. <= 18)], [W_sb_sim[experiment][i] for i in xrange(len(times)) if not (6 <= times[i]/60. <= 18)], 'ko', markersize = 10)
    ax.plot([V_sb_sim[experiment][i] for i in xrange(len(times)) if 6 <= times[i]/60. <= 18], [W_sb_sim[experiment][i] for i in xrange(len(times)) if 6 <= times[i]/60. <= 18], 'ro', markersize = 10)
    ax.set_xlabel('V$_{sb,sim}$')
    ax.set_ylabel('W$_{sb,sim}$')
    ax.set_title('Horizontal vs. Vertical CT component for exp=' + experiment)
    ax.set_xlim([0, range_max])
    ax.set_ylim([0, range_max])
    plt.show()
    """
#### Choose the maximum V_sb_sim and W_sb_sim for each experiment ####

H_control = 250.
H_exp01 = 125.
H_exp02 = 375.
H_sim = [H_control, H_exp01, H_exp02]
V_sb_sim_list = [np.max(V_sb_sim[experiment]) for experiment in experiments]
W_sb_sim_list = [np.max(W_sb_sim[experiment]) for experiment in experiments]

H_prediction = np.arange(50., 500.1, 5.)

#### Make our plots ####
fig = plt.figure(figsize = (15, 6))
ax0 = fig.add_subplot(1, 2, 1)
ax0.plot(H_prediction, V_sb_prediction(H_prediction), 'k--', label = 'V$_{sb}$ (KF15)')
ax0.plot(H_sim, V_sb_sim_list, 'b*', markersize = 15, label = 'V$_{sb,sim}$')
ax0.set_title('Maximum Cross-Flow Component')
ax0.set_xlabel('Sensible Heat Flux (W m$^{-2}$)')
ax0.set_ylabel('Cross-Flow Wind Speed (m s$^{-1}$)')
plt.legend(loc = 2)

ax1 = fig.add_subplot(1, 2, 2)
ax1.plot(H_prediction, W_sb_prediction(H_prediction), 'k--', label = 'W$_{sb}$ (KF15)')
ax1.plot(H_sim, W_sb_sim_list, 'b*', markersize = 15, label = 'W$_{sb,sim}$')
ax1.set_title('Maximum Vertical Velocity')
ax1.set_xlabel('Sensible Heat Flux (W m$^{-2}$)')
ax1.set_ylabel('Vertical Velocity (m s$^{-1}$)')
plt.legend(loc = 2)

plt.savefig('../KF15_comparison.png', dpi = 100)
plt.show()


