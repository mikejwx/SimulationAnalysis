import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import u_key, v_key, w_key, zi_new_key, lcl_key, lwp_key
from multiprocessing import Pool


# Read data, u, v, w from winds_nc, lcl from zi_nc, lwp from lwp_nc
path = '/nerc/n02/n02/xb899100/CloudTrail/Control/'
hours = ['09']

for hour in hours:
    wind_nc = Dataset(path + 'wind_' + hour + '.nc', 'r')
    zi_nc   = Dataset(path + 'zi_' + hour + '.nc', 'r')
    lwp_nc  = Dataset(path + 'lwp_00.nc', 'r')
    
    # read the 4D fields
    print ' Reading data.'
    if hour == hours[0]:
        u_data   = wind_nc.variables[u_key][:]*1.
        v_data   = wind_nc.variables[v_key][:]*1.
        w_data   = wind_nc.variables[w_key][:]*1.
        lcl_data = zi_nc.variables[lcl_key][:]*1.
        zi_data  = zi_nc.variables[zi_new_key][:]*1.
        lwp_data = lwp_nc.variables[lwp_key][:]*1.
        t = zi_nc.variables['time'][:]*1.
        # create the coordinate system
        print ' Creating coordinate system.'
        z = wind_nc.variables['thlev_zsea_theta'][:]*1.
        x, y = np.meshgrid(np.arange(u_data.shape[3])*100.0, np.arange(u_data.shape[2])*100.0)
        
        time_key = [tkey for tkey in lwp_nc.variables.keys() if 'min' in tkey][0]
        t_lwp = lwp_nc.variables[time_key][:]*1.
    else:
        u_data   = np.concatenate((u_data, wind_nc.variables[u_key][:]*1.), axis = 0)
        v_data   = np.concatenate((v_data, wind_nc.variables[v_key][:]*1.), axis = 0)
        w_data   = np.concatenate((w_data, wind_nc.variables[w_key][:]*1.), axis = 0)
        lcl_data = np.concatenate((lcl_data, zi_nc.variables[lcl_key][:]*1.), axis = 0)
        zi_data  = np.concatenate((zi_data, zi_nc.variables[zi_new_key][:]*1.), axis = 0)
        t        = np.concatenate((t, zi_nc.variables['time'][:]*1.), axis = 0)
    wind_nc.close()
    zi_nc.close()
    lwp_nc.close()

# back track parcels from the lcl in the cloud trail band
print ' Identifying the cloud band region.'
t_0 = 540.0
t_start = 360.0
t_end   = 1080.0
it = np.where(np.abs(t - t_0) == np.min(np.abs(t - t_0)))[0][0]
it_lwp = np.where(np.abs(t_lwp - t_0) == np.min(np.abs(t_lwp - t_0)))[0][0]
it_start = np.where(np.abs(t_lwp - t_start) == np.min(np.abs(t_lwp - t_start)))[0][0]
it_end   = np.where(np.abs(t_lwp - t_end) == np.min(np.abs(t_lwp - t_end)))[0][0]+1

# Define cloud trail cloud band region
cloud_mask = np.where(lwp_data > 0, 1.0, 0.0)
cloud_band_region = np.nanmean(cloud_mask[it_start:it_end,:,:], axis = 0)
cloud_band_region_mask = np.where(cloud_band_region > 0.25, 1.0, 0.0)

# find the grid cells that contain lwp (1.0 if they have, np.nan if they don't)
# flatten from 2d - 1d array
# remove points that are np.nan
x_cld = np.array([x_val for x_val in np.where(cloud_mask[it_lwp,:,:]*cloud_band_region_mask, x, np.nan).flatten() if x_val == x_val])
y_cld = np.array([y_val for y_val in np.where(cloud_mask[it_lwp,:,:]*cloud_band_region_mask, y, np.nan).flatten() if y_val == y_val])

island_area = 50.0 # km2
island_radius = 1000.0*np.sqrt(island_area/np.pi) #m
x_c = 108000.0 - island_radius
y_c = 16000.0
range_from_island_centre = np.sqrt((x-x_c)**2 + (y-y_c)**2)
island_mask = np.where(range_from_island_centre <= island_radius, 1.0, 0.0)
x_cld = np.array([x_val for x_val in np.where(island_mask, x, np.nan).flatten() if x_val == x_val])
y_cld = np.array([y_val for y_val in np.where(island_mask, y, np.nan).flatten() if y_val == y_val])

# track the parcel forward in time from starting points over the island
print ' Tracking parcels.'
all_x = []
all_y = []
all_z = []
all_t = []
x_max = (x.max() + 100.0)
x_min = 0
y_max = (y.max() + 100.0)
y_min = 0
dt = 10.
"""
# Serial method
for parcel in range(0, len(x_cld), 25):
    print '  Parcel #' + str(parcel)
    iy_cld, ix_cld = np.where((x == x_cld[parcel])*(y == y_cld[parcel]))
    parcel_x, parcel_y, parcel_z = [[x_cld[parcel]],[y_cld[parcel]],[z[1]]] # start parcels at the lowest non-zero height
    iz_cld = np.where(np.abs(z - parcel_z[-1]) == np.min(np.abs(z - parcel_z[-1])))[0][0]
    parcel_t = [t_0*1]
    height_thresh = lcl_data[it,iy_cld[0],ix_cld[0]]
    while (parcel_z[-1] <= height_thresh) and (parcel_t[-1] <= t.max()):
        x_new = parcel_x[-1] + u_data[it,iz_cld,iy_cld[0],ix_cld[0]]*dt
        if (x_new >= x_max):
            x_new = x_new - x_max
        elif (x_new < x_min):
            x_new = x_new + x_max
        parcel_x.append(x_new)
        y_new = parcel_y[-1] + v_data[it,iz_cld,iy_cld[0],ix_cld[0]]*dt
        if (y_new >= y_max):
            y_new = y_new - y_max
        elif (y_new < y_min):
            y_new = y_new + y_max
        parcel_y.append(y_new)
        parcel_z.append(parcel_z[-1] + w_data[it,iz_cld,iy_cld[0],ix_cld[0]]*dt)
        parcel_t.append(parcel_t[-1] + dt/60.0)
        iy_cld, ix_cld = np.where((np.abs(x - parcel_x[-1]) == np.min(np.abs(x - parcel_x[-1])))*(np.abs(y - parcel_y[-1]) == np.min(np.abs(y - parcel_y[-1]))))
        iz_cld = np.where(np.abs(z - parcel_z[-1]) == np.min(np.abs(z - parcel_z[-1])))[0][0]
        it = np.where(np.abs(t - parcel_t[-1]) == np.min(np.abs(t - parcel_t[-1])))[0][0]
        height_thresh = lcl_data[it,iy_cld[0],ix_cld[0]]
    all_x.append(np.array(parcel_x))
    all_y.append(np.array(parcel_y))
    all_z.append(np.array(parcel_z))
    all_t.append(np.array(parcel_t))
"""
# Parallel method
def particle_trajectories(parcel):
    print '  Parcel #' + str(parcel)
    iy_cld, ix_cld = np.where((x == x_cld[parcel])*(y == y_cld[parcel]))
    parcel_x, parcel_y, parcel_z = [[x_cld[parcel]],[y_cld[parcel]],[z[1]]] # start parcels at the lowest non-zero height
    iz_cld = np.where(np.abs(z - parcel_z[-1]) == np.min(np.abs(z - parcel_z[-1])))[0][0]
    parcel_t = [t_0*1]
    height_thresh = lcl_data[it,iy_cld[0],ix_cld[0]]
    while (parcel_z[-1] <= height_thresh) and (parcel_t[-1] <= t.max()):
        x_new = parcel_x[-1] + u_data[it,iz_cld,iy_cld[0],ix_cld[0]]*dt
        if (x_new >= x_max):
            x_new = x_new - x_max
        elif (x_new < x_min):
            x_new = x_new + x_max
        parcel_x.append(x_new)
        y_new = parcel_y[-1] + v_data[it,iz_cld,iy_cld[0],ix_cld[0]]*dt
        if (y_new >= y_max):
            y_new = y_new - y_max
        elif (y_new < y_min):
            y_new = y_new + y_max
        parcel_y.append(y_new)
        parcel_z.append(parcel_z[-1] + w_data[it,iz_cld,iy_cld[0],ix_cld[0]]*dt)
        parcel_t.append(parcel_t[-1] + dt/60.0)
        iy_cld, ix_cld = np.where((np.abs(x - parcel_x[-1]) == np.min(np.abs(x - parcel_x[-1])))*(np.abs(y - parcel_y[-1]) == np.min(np.abs(y - parcel_y[-1]))))
        iz_cld = np.where(np.abs(z - parcel_z[-1]) == np.min(np.abs(z - parcel_z[-1])))[0][0]
        it = np.where(np.abs(t - parcel_t[-1]) == np.min(np.abs(t - parcel_t[-1])))[0][0]
        height_thresh = lcl_data[it,iy_cld[0],ix_cld[0]]
    return [np.array(parcel_x), np.array(parcel_y), np.array(parcel_z), np.array(parcel_t)]

p = Pool(24)
all_particles = p.map(particle_trajectories, range(0, len(x_cld), 25))
p.join()
p.close()
n_particles = len(all_particles)

# for parallel, need to sort all_particles out into coordinate arrays
all_x = [all_particles[particle][0] for particle in range(n_particles)]
all_y = [all_particles[particle][1] for particle in range(n_particles)]
all_z = [all_particles[particle][2] for particle in range(n_particles)]
all_t = [all_particles[particle][3] for particle in range(n_particles)]

# terminate parcels when they reach the lcl

### Fig 1
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, adjustable = 'box', aspect = 1)
ax.set_ylabel('y (km)')
ax.set_xlabel('x (km)')
ax.contourf(x/1000.0, y/1000.0, cloud_mask[it_lwp,:,:], levels = [0, 1-1e-16, 1], colors = ['w', 'k'])
ax.contour(x/1000.0, y/1000.0, cloud_band_region_mask, levels = [0.5], colors = ['r'], linewidths = [2])
for parcel in range(len(all_x)):
    ax.plot(np.array(all_x[parcel])/1000.0, np.array(all_y[parcel])/1000.0, 'b')

ax.set_title('Trajectories of all ' + str(len(all_x)) + ' particles from 2 m to the local LCL')
plt.savefig('./particle_trajectories.png', dpi = 150, bbox_inches = 'tight')
plt.show()

### Fig 2
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
for parcel in range(len(all_x)):
    ax.plot(all_t[parcel], all_z[parcel])

ax.set_title('Time evolution of particle height')
ax.set_ylabel('Height (m)')
ax.set_xlabel('Time (mins)')
ax.set_xticks(range(540, 721, 60))
ax.set_xlim([540, 720])
plt.savefig('./particle_heights.png', dpi = 150, bbox_inches = 'tight')
plt.show()

### Fig 3
endpoint_density = np.histogram2d([all_x[parcel][-1] for parcel in range(len(all_x))], [all_y[parcel][-1] for parcel in range(len(all_y))], bins = [np.arange(0, x.max()+100.0, 1000.0), np.arange(0, y.max()+100.0, 1000.0)])

from matplotlib.colors import BoundaryNorm

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, adjustable = 'box', aspect = 1)
ax.set_ylim([y.min()/1000.0, y.max()/1000.0])
ax.set_xlim([x.min()/1000.0, x.max()/1000.0])
X, Y = np.meshgrid(0.5*(endpoint_density[2][1:] + endpoint_density[2][:-1]), 0.5*(endpoint_density[1][1:] + endpoint_density[1][:-1]))
X = np.transpose(X/1000.0)
Y = np.transpose(Y/1000.0)
Z = np.transpose(endpoint_density[0])
my_cmap = plt.get_cmap('hot')
my_levels = range(int(Z.max())+1)
my_norm = BoundaryNorm(my_levels, ncolors=my_cmap.N, clip=True)
ax.contour(x/1000.0, y/1000.0, island_mask, levels = [0.5], colors = ['r'], linewidths = [1.5])
my_plt = ax.pcolormesh(Y, X, Z, cmap = my_cmap, norm = my_norm)
plt.colorbar(my_plt, ax = ax, orientation = 'horizontal', label = 'Number of particles')
ax.set_title('Density of where particles cross LCL')
ax.set_xlabel('x (km)')
ax.set_ylabel('y (km)')
plt.savefig('./particle_density.png', dpi = 150, bbox_inches = 'tight')
plt.show()

