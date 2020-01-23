import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from analysis_tools import fromComponents, bilinear_interpolation, get_cs_coords, lcl as get_lcl, downwind_rectangle
from scipy import integrate, interpolate
from datetime import datetime as dt
from STASH_keys import rho_key, u_key, v_key, w_key, pthe_key, q_key, theta_key, lcl_key, zi_new_key
from SkewT_archer import PTtoTemp
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

"""
Compute the updraught mass flux over a range of distances downwind of the island.
Compare between the heat flux experiments.
"""
l_testing = True

# Define a function to do the integration
def get_Mc(y_prime, M):
    """
    Takes the rotated coordinate system (x_prime, y_prime) and the mass flux (M)
    and averages it in the x_prime direction, integrates in the y_prime direction.
    
    To do this, it must first do some housekeeping on the arrays to make sure
    that we are integrating along the right cartesian axis, rather than the 
    array axis.
    
    y_prime -> the rotated y coordinate, increasing to the right of the line
               directly downwind of the island (i.e. positive to thr right, and
               negative to the left).
    M       -> the time-mean mass flux array at a given height
    
    All above arrays are 2-D.
    """
    # M_up = positive M
    M = np.where(M > 0, M, 0.0)
    # Average along trail
    M_x         = np.array([np.nanmean(np.where((lower <= y_prime)*(y_prime <= lower + 100.0), M, np.nan)) for lower in np.arange(-3050.0, 3000.0, 100.0)])
    new_y_prime = np.arange(-3000.0, 3000.1, 100.0)
    # Integrate across trail
    M_int = integrate.trapz(x = new_y_prime, y = M_x)
    
    return M_int

# read in the data
print '[' + dt.now().strftime('%H:%M:%S') + '] Defining the paths to data...'
# Define the paths in which the data are stored

paths = {'Control_short' : '/nerc/n02/n02/xb899100/CloudTrail/Control_short/',
         'BLm25' : '/nerc/n02/n02/xb899100/CloudTrail/RH_BLm25/',
         'FAm25' : '/nerc/n02/n02/xb899100/CloudTrail/RH_FAm25/',
         'Control' : '/nerc/n02/n02/xb899100/CloudTrail/Control2/',
         'Control_predawn' : '/nerc/n02/n02/xb899100/CloudTrail/Control2/'}

print '[' + dt.now().strftime('%H:%M:%S') + '] Starting the routine for each experiment...'
mass_flux_xs = {}
lcl = {}
zi = {}
# For each experiment
for key in paths.keys():
    print '[' + dt.now().strftime('%H:%M:%S') + '] Starting experiment ' + key + '...'
    hour = ['09' if key == 'Control' else '03' if key == 'Control_predawn' else '04'][0]    # Which time chunk of output are we analysing?
    start_t = [180. if key == 'Control_predawn' else 540.][0] # What time are we starting the analysis?
    end_t = [360. if key == 'Control_predawn' else 720.][0]   # What time are we finishing the analysis?
    
    print '[' + dt.now().strftime('%H:%M:%S') + '] Opening the netCDF...'
    # Read in swath netCDF
    rho_nc    = Dataset(paths[key] + 'fluxes_' + hour + '.nc', 'r')
    u_nc      = Dataset(paths[key] + 'u_' + hour + '.nc', 'r')
    v_nc      = Dataset(paths[key] + 'v_' + hour + '.nc', 'r')
    thermo_nc = Dataset(paths[key] + 'bouy_' + hour + '.nc', 'r')
    zi_nc = Dataset(paths[key] + 'zi_' + hour + '.nc', 'r')
    
    # Find the start and end time idxs
    rho_time_key = [tkey for tkey in rho_nc.variables.keys() if 'min' in tkey][0]
    zi_time_key = 'time'
    rho_times = rho_nc.variables[rho_time_key][:]+[0.0 if key in ['Control', 'Control_predawn'] else 120.0][0]
    zi_times = zi_nc.variables[zi_time_key][:]+[0.0 if key in ['Control', 'Control_predawn'] else 120.0][0]
    start_idx = np.where(np.abs(rho_times - start_t) == np.min(np.abs(rho_times - start_t)))[0][0]
    end_idx   = np.where(np.abs(rho_times - end_t)   == np.min(np.abs(rho_times - end_t)))[0][0]
    print '[' + dt.now().strftime('%H:%M:%S') + '] Finding the indexes'
    t_idx = [it for it in range(rho_times.size) if (start_t <= rho_times[it])*(rho_times[it] <= end_t)*(rho_times[it]%30. == 0)]
    t_idx2 = [it for it in range(zi_times.size) if (start_t <= zi_times[it])*(zi_times[it] <= end_t)*(zi_times[it]%30. == 0)]
    print '[' + dt.now().strftime('%H:%M:%S') + '] Reading the required dimensions...'
    # Want the cross section across the trail averaged in time between start_t 
    # and end_t
    
    ### Set-up the coordinate system ###
    # Define our height dimension
    z_theta = thermo_nc.variables['thlev_zsea_theta'][:]*1.
    
    # We know from our domain definition where the centre of the island should be
    R_i = 1000.0*(50.0/np.pi)**0.5 # island radius, metres
    x_c = 100000.0 + R_i
    y_c = 4*R_i
    
    mass_flux_xs[key] = {}
    
    # use our new function to define a rectangle downwind
    # Define the existing cartesian coordinate system
    X, Y = np.meshgrid(np.arange(thermo_nc.variables[theta_key][0,0,0,:].shape[0])*100.0, np.arange(thermo_nc.variables[theta_key][0,0,:,0].shape[0])*100.0)
    
    # Compute and store the horizontal mean LCL
    lcl[key] = np.nanmean(zi_nc.variables[lcl_key][t_idx2,:,:])
    zi[key] = np.nanmean(zi_nc.variables[zi_new_key][t_idx2,:,:])
    
    # Get the mean wind speed and direction below LCL
    iz_LCL = np.where(np.abs(z_theta - lcl[key]) == np.min(np.abs(z_theta - lcl[key])))[0][0]
    iz_zi  = np.where(np.abs(z_theta - zi[key]) == np.min(np.abs(z_theta - zi[key])))[0][0]
    u_mean = np.array([u_nc.variables[u_key][idx,:iz_zi,:,:] for idx in t_idx]).mean(axis = (0, 2, 3))
    v_mean = np.array([v_nc.variables[v_key][idx,:iz_zi,:,:] for idx in t_idx]).mean(axis = (0, 2, 3))
    U_zi = integrate.trapz(x = z_theta[:iz_zi], y = u_mean)/zi[key]
    V_zi = integrate.trapz(x = z_theta[:iz_zi], y = v_mean)/zi[key]
    speed, wind_dir = fromComponents(u = U_zi, v = V_zi)
    
    # Get the distance away from the island-edge of the surface warm plume.
    mass_flux_xs[key]['theta_anom'] = np.nanmean(np.array([thermo_nc.variables[theta_key][idx,1,:,:] - thermo_nc.variables[theta_key][idx,1,:,:].mean() for idx in t_idx]), axis = 0)
    R_wp = 4*R_i#np.nanmax(np.where(mass_flux_xs[key]['theta_anom'] > 0.1, np.sqrt((X-x_c)**2 + (Y-y_c)**2) - R_i, np.nan))
    
    # Compute the rectangle and the new cartesian coordinate system
    # distance travelled during time averaging period
    mass_flux_xs[key]['mask'], mass_flux_xs[key]['y_prime'], mass_flux_xs[key]['x_prime'] = downwind_rectangle(wind_dir = wind_dir, x_c = x_c, y_c = y_c, X = X, Y = Y, island_radius = R_i, dist_0 = R_wp, dist_1 = R_wp + ((end_t-start_t)*60.0*speed - R_wp), half_width = 3000.0)
    
    print '[' + dt.now().strftime('%H:%M:%S') + '] Computing mass flux...'
    # Convert rho from rho levels to theta levels
    rho_theta_levels = np.array([interpolate.interp1d(x = rho_nc.variables['rholev_zsea_rho'][:], y = rho_nc.variables[rho_key][idx,:,:,:], axis = 0, fill_value = 'extrapolate')(z_theta[iz_zi]) for idx in t_idx])
    # Average in time
    mf = np.nanmean(rho_theta_levels*thermo_nc.variables[w_key][t_idx,iz_zi,:,:], axis = 0)
    mass_flux_xs[key]['mass_flux'] = mf
    mass_flux_xs[key]['mass_flux_masked'] = mf*mass_flux_xs[key]['mask']
    
    rho_nc.close()
    u_nc.close()
    v_nc.close()
    thermo_nc.close()
    zi_nc.close()

# do an integration...
circulation_MF = {}
for key in paths.keys():
    circulation_MF[key] = get_Mc(mass_flux_xs[key]['y_prime'], mass_flux_xs[key]['mass_flux_masked']*1.0)

print zi
print lcl
print circulation_MF
"""
import matplotlib
my_cmap = matplotlib.cm.get_cmap('Reds')
my_colors = {}
for key in circulation_MF.keys():
    my_colors[key] = my_cmap(float(key[1:4])/500.0)

g = 9.81
T0 = 302.3
q0 = 0.017
rho0 = 1.1771
cvirt = 0.61
H = {}
Q = {}
B = {}
t0 = 720.0
dt = 720.0
for key in paths.keys():
    H[key] = float(key[1:4])
    Q[key] = H[key] + float(key[5:])
    wpqvp = (Q[key] - H[key])/(1.17*2.501e6)
    wptp = H[key]/(1.17*1005.)
    B[key] = (g/T0)*((wptp)*(1. + cvirt*q0) + cvirt*T0*(wpqvp))

### Make the plot ###
fig = plt.figure(figsize = (10,7.5))
# Against peak total energy flux
axa = fig.add_subplot(2, 3, 1, adjustable = 'box')
[axa.plot(Q[key], circulation_MF[key], color = my_colors[key], 
     marker = ['o' if int(key[5:8]) == 250 else 'v' if int(key[5:8]) == 375 else '^'][0], ls = 'None',
     markersize = 7.5, markeredgecolor = 'none') for key in paths.keys()]
axa.set_ylabel(u'Circulation Mass Flux (kg m$^{-1}$ s$^{-1}$)')
axa.set_xlabel(u'$Q_{net} = H_{0} + E_{0}$ (W m$^{-2}$)')
# Plot a line of best fit
my_x = np.array([int(key[1:4]) + int(key[5:8]) for key in paths.keys()])
my_y = np.array([circulation_MF[key] for key in paths.keys()])
my_fit = np.polyfit(my_x, my_y, 1)
axa.plot([0,750], [np.sum([my_fit[i]*(x**(len(my_fit) - i - 1)) for i in range(len(my_fit))]) for x in [0,750]], 'k--')
axa.set_title(u'a) y = ' + str(round(my_fit[0], 2)) + ' x ' +['+ ' if my_fit[1] >= 0 else ''][0] + str(round(my_fit[1], 4)), fontsize = 12)
# Domain shape parameters
axa.set_xlim([250, 750])
axa.set_ylim([0, 400])
axa.set_yticks(np.arange(0, 401, 100))
# point to the control point
axa.annotate('Control', xy = (500, circulation_MF['H250E250']), xytext = (600, 0.75*circulation_MF['H250E250']), arrowprops = dict(color = 'k', width = 0.5, headwidth = 7.5, headlength = 7.5))
# point to the bowen ratio experiments
axa.annotate('$\\beta$ = 1/3', xy = (500, circulation_MF['H125E375']), xytext = (500, 0.5*circulation_MF['H125E375']), arrowprops = dict(color = 'k', width = 0.5, headwidth = 7.5, headlength = 7.5), horizontalalignment = 'center')
axa.annotate('$\\beta$ = 3', xy = (500, circulation_MF['H375E125']), xytext = (500, 1.25*circulation_MF['H375E125']), arrowprops = dict(color = 'k', width = 0.5, headwidth = 7.5, headlength = 7.5), horizontalalignment = 'center')

# Against peak suface sensible heat flux
axb = fig.add_subplot(2, 3, 2, adjustable = 'box')
[axb.plot(H[key], circulation_MF[key], color = my_colors[key], 
     marker = ['o' if int(key[5:8]) == 250 else '^' if int(key[5:8]) == 375 else 'v'][0], ls = 'None',
     markersize = 7.5, markeredgecolor = 'none') for key in paths.keys()]
axb.set_xlabel(u'$H_{0}$ (W m$^{-2}$)')
# Plot a line of best fit
my_x = np.array([int(key[1:4]) for key in paths.keys()])
my_y = np.array([circulation_MF[key] for key in paths.keys()])
my_fit = np.polyfit(my_x, my_y, 1)
axb.plot([0,500], [np.sum([my_fit[i]*(x**(len(my_fit) - i - 1)) for i in range(len(my_fit))]) for x in [0,500]], 'k--')
axb.set_title(u'b) y = ' + str(round(my_fit[0], 2)) + ' x ' +['+ ' if my_fit[1] >= 0 else ''][0] + str(round(my_fit[1], 4)), fontsize = 12)
# Domain shape parameters
axb.set_xlim([0, 500])
axb.set_ylim([0, 400])
axb.set_yticks(np.arange(0, 401, 100))
axb.set_yticklabels([''])
axb.annotate('Control', (250, circulation_MF['H250E250']), (250, 0.5*circulation_MF['H250E250']), horizontalalignment = 'center', arrowprops = dict(color = 'k', width = 0.5, headwidth = 7.5, headlength = 7.5))

# Against peak land buoyancy flux (thetav'w')
axc = fig.add_subplot(2, 3, 3, adjustable = 'box')
[axc.plot(B[key], circulation_MF[key], color = my_colors[key], marker = ['o' if int(key[5:8]) == 250 else '^' if int(key[5:8]) == 375 else 'v'][0], ls = 'None', markersize = 7.5, markeredgecolor = 'none') for key in paths.keys()]
axc.set_xlabel(u'$\\frac{g}{\\theta_{v,0}}$ ($\overline{w^{\\prime} \\theta_{v}^{\\prime}}$)$_{0}$ (m$^{2}$ s$^{-3}$)')
# Plot a line of best fit
my_x = np.array([B[key] for key in paths.keys()])
my_y = np.array([circulation_MF[key] for key in paths.keys()])
my_fit = np.polyfit(my_x, my_y, 1)
axc.plot([0,0.015], [np.sum([my_fit[i]*(x**(len(my_fit) - i - 1)) for i in range(len(my_fit))]) for x in [0,0.015]], 'k--')
axc.set_title(u'c) y = ' + str(round(my_fit[0], 3)) + ' x ' +['+ ' if my_fit[1] >= 0 else ''][0] + str(round(my_fit[1], 4)), fontsize = 12)
# Domain shape parameters
axc.set_xlim([0, 0.015])
axc.set_ylim([0, 400])
axc.set_yticks(np.arange(0, 401, 100))
axc.set_xticks(np.arange(0, 0.0151, 0.005))
axc.set_xticklabels(['0', '5e-3', '10e-3', '15e-3'])
axc.set_yticklabels([''])
axc.annotate('Control', (B['H250E250'], circulation_MF['H250E250']), (0.01, 0.75*circulation_MF['H250E250']), arrowprops = dict(color = 'k', width = 0.5, headwidth = 7.5, headlength = 7.5))
"""

fig = plt.figure()
my_keys = ['Control_short', 'BLm25', 'FAm25']
my_letters = ['a)', 'b)', 'c)']
for key in my_keys:
    ax = fig.add_subplot(3, 1, my_keys.index(key) + 1, adjustable = 'box', aspect = 1)
    R = np.sqrt((X - x_c)**2 + (Y - y_c)**2)
    my_plt = ax.contourf(X/1000, Y/1000, mass_flux_xs[key]['mass_flux'], cmap = 'bwr', levels = [level for level in np.arange(-1.0, 1.1, 0.125) if round(level,1) != 0], extend = 'both')
    ax.contour(X/1000., Y/1000., mass_flux_xs[key]['mask'], levels = [0.5], colors = ['k'], linestyles = ['--'])
    my_plt.cmap.set_over('firebrick')
    my_plt.cmap.set_under('navy')
    ax.contour(X/1000, Y/1000, R, levels = [R_i], colors = ['k'], linewidths = [2])

my_cb_ticks = [-1, -2./3., -1./3., 0, 1./3., 2./3., 1]
my_cb_ticks = [round(tick, 2) for tick in my_cb_ticks]
cbar_ax = inset_axes(ax, width = "100%", height = "7.5%", loc = 3, bbox_to_anchor = (0.0, -0.35, 1.0, 1.0), bbox_transform = ax.transAxes, borderpad = 0.0)
cbar = fig.colorbar(my_plt, cax = cbar_ax, orientation = 'horizontal', label = 'Mass Flux (kg m$^{-2}$ s$^{-1}$)')
cbar.set_ticks(my_cb_ticks)
cbar.ax.set_xticklabels([str(tick) for tick in my_cb_ticks])
ax.set_title(my_letters[my_keys.index(key)] + ' ' + key + ' Mass Flux at z = LCL', fontsize = 12) 
ax.set_xlabel('x (km)')
ax.set_ylabel('y (km)')

plt.subplots_adjust(hspace = 0.4, bottom = 0.2)
#plt.savefig('../Ch5_Figure14.png', dpi = 250, bbox_inches = 'tight')
plt.show()

