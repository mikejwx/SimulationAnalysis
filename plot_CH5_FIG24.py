import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import *
from scipy import ndimage
import matplotlib as mpl

### Define functions to compute the analytical solutions ###
def vct(U):
    g = 9.81
    H = 250.0
    lx = island_radius*1.
    ly = island_radius*1.
    rho = 1.17
    cpd = 1005.0
    zh = 600.0
    T0 = 300.0
    vct = (g*H*lx*lx*ly)/(4.*rho*cpd*zh*zh*T0*U*U)
    return vct

def lwp(U):
    g = 9.81
    H = 250.0
    lx = island_radius*1.
    ly = island_radius*1.
    rho = 1.17
    cpd = 1005.0
    zh = 600.0
    T0 = 300.0
    lwp = (4.*rho*cpd*zh*zh*T0*U*U*U)/(g*H*lx*lx)
    return lwp

### Read in theta p and n at T = 720 for wind experiments ###
paths = {'U10' : '/nerc/n02/n02/xb899100/CloudTrail/Control2/',
         'U09' : '/nerc/n02/n02/xb899100/CloudTrail/U09/',
         'U08' : '/nerc/n02/n02/xb899100/CloudTrail/U08/',
         'U07' : '/nerc/n02/n02/xb899100/CloudTrail/U07/',
         'U06' : '/nerc/n02/n02/xb899100/CloudTrail/U06/',
         'U05' : '/nerc/n02/n02/xb899100/CloudTrail/U05/'}

keys = ['U05', 'U06', 'U07', 'U08', 'U09', 'U10']
my_data = {}
for key in keys:
    my_data[key] = {}
    with Dataset(paths[key] + 'bouy_' + ['04' if key in ['U05'] else '09'][0] + '.nc', 'r') as bouy_nc:
        time_key = [tkey for tkey in bouy_nc.variables.keys() if 'min' in tkey][0]
        times = bouy_nc.variables[time_key][:] + [240.0 if key in ['U05'] else 0.0][0]
        t_idx = [it for it in range(times.size) if (660. <= times[it])*(times[it] <= 720.)]
        z = bouy_nc.variables['thlev_zsea_theta'][:]*1.
        iz = np.where(np.abs(z-2) == np.min(np.abs(z-2)))[0][0]
        my_data[key]['thetap'] = np.array([bouy_nc.variables[theta_key][it, iz, :,:] - bouy_nc.variables[theta_key][it, iz, :,:].mean() for it in t_idx]).mean(axis = 0)
    
    with Dataset(paths[key] + 'wind_' + ['04' if key in ['U05'] else '09'][0] +'.nc', 'r') as wind_nc:
        time_key = [tkey for tkey in wind_nc.variables.keys() if 'min' in tkey][0]
        times = wind_nc.variables[time_key][:] + [240.0 if key in ['U05'] else 0.0][0]
        t_idx = [it for it in range(times.size) if (660. <= times[it])*(times[it] <= 720.)]
        z = wind_nc.variables['thlev_zsea_theta'][:]*1.
        iz = np.where(np.abs(z-10) == np.min(np.abs(z-10)))[0][0]
        my_data[key][n_key] = wind_nc.variables[n_key][t_idx, iz, :,:].mean(axis = 0)

### find the maximum vct ###
island_radius = np.sqrt(50.0/np.pi)*1000.0
my_sims = {}
for key in keys:
    my_sims[key] = {}
    my_sims[key]['vct_sim'] = np.abs(my_data[key][n_key]).max()
    my_sims[key]['vct_est'] = vct(float(key[1:]))
    X, Y = np.meshgrid(np.arange(my_data[key]['thetap'].shape[1])*100.0, np.arange(my_data[key]['thetap'].shape[0])*100.0)
    x_c = 100000. + [island_radius if key in ['U05', 'U10'] else + 8000.0][0]
    y_c = [4*island_radius if key in ['U05', 'U10'] else 16000.0][0]
    R = np.sqrt((X - x_c)**2 + (Y - y_c)**2)
    
    warm_patches, n_patches = ndimage.label(np.where(my_data[key]['thetap'] > 0.1, 1.0, 0.0))
    bp_no = 0
    for patch_no in range(n_patches+1):
        # Assuming patch 0 is the null patch
        # Assuming the warm plume is overlapping the island
        if np.max(np.where(R < island_radius, 1.0, 0.0)*np.where(warm_patches == patch_no, 1.0, 0.0)):
            bp_no = patch_no
    
    my_sims[key]['lwp_sim'] = np.max(R*np.where(warm_patches == bp_no, 1.0, 0.0))/1000.
    my_sims[key]['lwp_est'] = lwp(float(key[1:]))/1000.

plt_colors = {}
my_cmap = mpl.cm.get_cmap('Reds')
for key in keys:
    plt_colors[key] = my_cmap((keys.index(key)+2)/float(len(keys)+3))

fig = plt.figure(figsize = (9.5, 5))
axa = fig.add_subplot(1, 2, 1)
my_x = np.array([1/(float(key[1:])**2) for key in keys])
my_y = np.array([my_sims[key]['vct_sim'] for key in keys])
[axa.plot(1/(float(key[1:])**2), my_sims[key]['vct_sim'], marker = 'o', color = plt_colors[key], markeredgecolor = 'none') for key in keys]
my_fit = np.polyfit(my_x, my_y, 1)
axa.plot([0,0.05], [np.sum([my_fit[i]*(x**(len(my_fit) - i - 1)) for i in range(len(my_fit))]) for x in [0,0.05]], 'k--')
axa.set_xlabel('$1/u^{2}$ of experiment (m$^{-2}$ s$^{2}$)')
axa.set_ylabel('$v_{CT}^{\prime}$ (m s$^{-1}$)')
axa.set_title('a) Scaling $v_{CT}^{\prime}$')
axa.set_xlim([0, 0.05])
axa.text(0.02, 1.75, '$y = $ ' + str(round(my_fit[0], 3)) + '$x$' + ['+' if my_fit[1] > 0 else ''][0] + str(round(my_fit[1], 3)))

axb = fig.add_subplot(1, 2, 2)
my_x = np.array([float(key[1:])**3 for key in keys])
my_y = [my_sims[key]['lwp_sim'] for key in keys]
[axb.plot(float(key[1:])**3, my_sims[key]['lwp_sim'], marker = 'o', color = plt_colors[key], markeredgecolor = 'none') for key in keys]
my_fit = np.polyfit(my_x, my_y, 1)
axb.plot([0,1100], [np.sum([my_fit[i]*(x**(len(my_fit) - i - 1)) for i in range(len(my_fit))]) for x in [0,1100]], 'k--')
axb.set_xlim([0, 1100])
axb.set_xlabel('$u^{3}$ (m$^{3}$ s$^{-3}$)')
axb.set_ylabel('$L_{WP}$ from simulation (km)')
axb.set_title('b) Scaling $L_{WP}$')
axb.text(400, 7, '$y = $ ' + str(round(my_fit[0], 3)) + '$x$' + ['+' if my_fit[1] > 0 else ''][0] + str(round(my_fit[1], 3)))
plt.savefig('../Ch5_Figure24.png', dpi = 250, bbox_inches = 'tight')
plt.show()

