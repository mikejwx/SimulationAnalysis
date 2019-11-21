import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import *
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from coarse_graining import coarse_grain
from scipy import ndimage, signal, spatial
import matplotlib.tri as mtri
from datetime import datetime as dt

print '[' + dt.now().strftime('%H:%M:%S') + '] Reading the data'
path = '/nerc/n02/n02/xb899100/CloudTrail/Control/'
iz_max = 30
with Dataset(path + 'bouy_09.nc', 'r') as bouy_nc:
    theta_data = bouy_nc.variables[theta_key][-1, :iz_max, :, :]
    z = bouy_nc.variables['thlev_zsea_theta'][:iz_max]*1.

# compute the theta anomalies
print '[' + dt.now().strftime('%H:%M:%S') + '] Computing theta anomalies'
theta_anom = np.array([theta_data[iz,:,:] - theta_data[iz,:,:].mean() for iz in range(theta_data.shape[0])])

# smooth the field to make reduce noise
print '[' + dt.now().strftime('%H:%M:%S') + '] Smoothing the anomalies'
smoothing_diameter = 6000.0
# Define a smoothing structure
n_r, n_c = np.array([smoothing_diameter, smoothing_diameter])/100
w_c = np.ones((n_r, n_c))
r, c = np.meshgrid(np.arange(n_r), np.arange(n_c))
w_c = np.where((np.sqrt((r - n_r/2)**2 + (c - n_c/2)**2) <= n_r/2), w_c, 0.)
w_c /= np.sum(w_c)
theta_anom = np.array([signal.convolve2d(theta_anom[iz,:,:], w_c, mode = 'same', boundary = 'wrap') for iz in range(theta_anom.shape[0])])

# coarse grain the theta anomalies to reduce computational cost
print '[' + dt.now().strftime('%H:%M:%S') + '] Coarse graining the anomalies'
scale = 10
if scale > 1:
    theta_anom = coarse_grain(theta_anom, scale, scale)

# manfacture a grid coordinate system
y, z, x = np.meshgrid(np.arange(theta_anom.shape[1])*100.0*scale, z, np.arange(theta_anom.shape[2])*100.0*scale)

# threshold to identify the warm spots
print '[' + dt.now().strftime('%H:%M:%S') + '] Thresholding the data'
threshold = 0.25
my_mask = np.where(theta_anom > threshold, 1.0, 0.0)

# only keep the contiguous warm spots...
my_warm_plume = np.zeros_like(my_mask)

# choose the largest spot at the surface
warm_spots, n_spots = ndimage.label(my_mask[0,:,:])
spot_sizes = np.array([np.nansum(np.where(warm_spots == spot, 1.0, 0.0)) for spot in range(1, n_spots+1)])
spot_idx = np.where(spot_sizes == np.max(spot_sizes))[0][0]+1
my_warm_plume[0,:,:] = np.where(warm_spots == spot_idx, 1.0, 0.0)

# only keep if the spot overlaps the surface spot we've just chosen
for iz in range(1, my_mask.shape[0]):
    warm_spots, n_spots = ndimage.label(my_mask[iz,:,:])
    for spot_idx in range(1, n_spots+1):
        if np.any(np.where(warm_spots == spot_idx, 1.0, np.nan) == my_warm_plume[iz-1,:,:]):
            my_warm_plume[iz,:,:] = np.max(np.array([my_warm_plume[iz,:,:], np.where(warm_spots == spot_idx, 1.0, 0.0)]), axis = 0)

# work our way back down and only keep spots if there is a match above as well
for iz in range(my_mask.shape[0]-2, 0, -1):
    warm_spots, n_spots = ndimage.label(my_mask[iz,:,:])
    for spot_idx in range(1, n_spots+1):
        if np.any(np.where(warm_spots == spot_idx, 1.0, np.nan) == my_warm_plume[iz+1,:,:]):
            my_warm_plume[iz,:,:] = np.max(np.array([my_warm_plume[iz,:,:], np.where(warm_spots == spot_idx, 1.0, 0.0)]), axis = 0)

# initialise empty lists to store the coordinates of contours
print '[' + dt.now().strftime('%H:%M:%S') + '] Computing plume outlines at each height'
contour_x, contour_y, contour_z = [[], [], []]
fig = plt.figure()
ax = fig.add_subplot(1,1,1,adjustable = 'box', aspect = 1)
for iz in range(my_mask.shape[0]):
    # use matplotlib.pyplot to get the coordinates of contours at each height
    plt.cla()
    c = ax.contour(x[0,:,:], y[0,:,:], my_mask[iz,:,:], levels = [0.5])
    ax.set_title(str(z[iz,0,0]) + ' m')
    plt.pause(0.01)
    for contour_path in c.collections[0].get_paths():
        # for each contour append all of the coordinates to the contour lists
        v = contour_path.vertices
        contour_x += list(v[:,0])
        contour_y += list(v[:,1])
        contour_z += list(np.zeros_like(v[:,0])+z[iz,0,0])

plt.close('all')

# convert the lists to arrays to do array manipulations/math
contour_x = np.array(contour_x)/1000.
contour_y = np.array(contour_y)/1000.
contour_z = np.array(contour_z)/1000.

# triangulate to determine the triangles
print '[' + dt.now().strftime('%H:%M:%S') + '] Performing delaunay triangulation'
triangle_idx = []
for idx in range(len(contour_x)):
    # find the two nearest points and include them in the triangle
    distance = np.sqrt((contour_x[idx] - contour_x)**2 + (contour_y[idx] - contour_y)**2)
    
    # remove the point itself
    distance[distance == 0.0] = np.nan
    
    # remove points that are more than one vertical level separated from idx
    # determine the height that we're at
    iz = np.where(z[:,0,0]/1000. == contour_z[idx])[0][0]
    our_hgts = [z[iz,0,0]/1000.]
    n_verts = 2
    
    if iz != 0:
        our_hgts.append(z[iz-1,0,0]/1000.)
        n_verts += 3
    
    if iz != (z.shape[0]-1):
        our_hgts.append(z[iz+1,0,0]/1000.)
        n_verts += 3
    
    distance = np.array([distance[idist] if contour_z[idist] in our_hgts else np.nan for idist in range(len(distance))])
    
    min_idxs = []
    while len(min_idxs) < n_verts:
        min_idxs.append(np.where(distance == np.nanmin(distance))[0][0])
        distance[min_idxs[-1]] = np.nan
    
    # combine the six adjacent points into triangles
    for i in range(n_verts):
        temp_triangle = [idx, min_idxs[i], min_idxs[(i+1)%n_verts]]
        triangle_idx.append(temp_triangle)


triang = mtri.Triangulation(contour_x, contour_y, triangles = triangle_idx)

print '[' + dt.now().strftime('%H:%M:%S') + '] plotting the data'
# use the contour coordinates on each height level to create a triangular surface
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection = '3d')
ax.plot_trisurf(triang, contour_z, edgecolor = 'none')
#ax.scatter(x2, y2, z2)
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.set_zlabel('Z (km)')

island_area = 50.0
island_radius = np.sqrt(island_area/np.pi)*1000.
x_c = 100000. + island_radius
y_c = 4*island_radius
R = np.sqrt((x[0,:,:] - x_c)**2 + (y[0,:,:] - y_c)**2)
ax.contour(x[0,:,:]/1000., y[0,:,:]/1000., R, levels = [0, island_radius], colors = ['r'], zdir = 'z', offset = 0)
ax.set_ylim([8, 24])
ax.set_xlim([100, 116])
ax.set_zlim([0, z.max()/1000.])
plt.show()

