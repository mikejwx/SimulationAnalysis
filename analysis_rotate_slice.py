"""
1. Generate real cartesian coordinates
2. Find the island center point
3. Find the domain mean wind direction at e.g. 750 m to estimate CT heaing
4. Find the distance and heading of each grid point from the island center (r, theta)
5. Find the difference between the wind direction at 750 m and 270 deg
    (So that we can rotate the domain so that flow is east-west.
6. Add the difference from the heading of each grid point
7. convert back from (r, theta) to (x,y)
"""
## Analysis of some netcdf output from island simulations
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

def hmean(aIN):
    """
    Takes in aIN as a 3D array (z, y, x) and computes the horizontal mean.
    """
    aOUT = np.zeros(aIN.shape[0])
    for k in range(aIN.shape[0]):
        aOUT[k] = np.mean(aIN[k,:,:])
    return aOUT

# Step 1
x = np.arange(0, 116000, 100.)
y = np.arange(0, 31900, 100.)
x, y = np.meshgrid(x, y)

# Step 2
"""
We define the domain size as (100 + 2D, 4D) (x, y) in km
Where D is the center point of the island
"""
A = 50.0 #Area of the island (km2)
D = 2.*np.sqrt(A/np.pi)
x_isl = (100+D)*1000.
y_isl = (2*D)*1000.

# Step 3
"""
Because the wind direction does change somewhat throughout the simulation, let's
assume that the wind direction is a constant 45deg (northeasterly).
"""
# My estimate for the cloud level wind direction
wind_r = np.sqrt((78726-x_isl)**2 + (9034 -y_isl)**2)
wind_dir = np.arcsin(-(78726-x_isl)/wind_r)
wind_dir = wind_dir*180./np.pi

# Step 4
def cartesian2polar(x, y, x_o = 0., y_o = 0.):
    """
    Convert from cartesian to polar coordinates
    x, y = the cartesian coordinates
    x_o, y_o = the coordinates of the origin
    """
    dx = x - x_o
    dy = y - y_o
    r = np.sqrt(dx**2 + dy**2)
    
    """
    Four cases:
    1. if dx and dy are positive
    2. if dx is positive and dy is negative
    3. if dx and dy are negative
    4. if dx is negative and dy is positive
    """
    theta1 = np.arcsin(dx/r)
    theta2 = np.arcsin(-dy/r) + np.pi/2.
    theta3 = np.arcsin(-dx/r) + np.pi
    theta4 = np.arcsin(dy/r)  + 3.*np.pi/2.
    
    theta1[dx < 0.] = 0.
    theta1[dy < 0.] = 0.
    
    theta2[dx < 0.] = 0.
    theta2[dy > 0.] = 0.

    theta3[dx > 0.] = 0.
    theta3[dy > 0.] = 0.
    
    theta4[dx > 0.] = 0.
    theta4[dy < 0.] = 0.
    
    # Add together and convert from radians to degrees
    theta = (theta1 + theta2 + theta3 + theta4)*180./np.pi
    
    return r, theta

r, theta = cartesian2polar(x, y, x_isl, y_isl)

# Step 5
"""
Wind direction is NE'ly (i.e. bearing = 45 degrees), heading is therefore SW'ly
and heading is 225 degrees.
"""

wind_heading = wind_dir + 180.
# check that the wind doesn't alias
if wind_heading > 360.:
    wind_heading -= 360.

# Find the offset to rotate the domain by
target_heading = 270.
heading_diff = target_heading - wind_heading

# Step 6
"""
Add the difference between the wind heading and the domain orientation
"""

theta_new = theta + heading_diff

# Step 7
"""
Convert back to (x,y) coordinates from the polar
"""

def polar2cartesian(r, theta, x_o, y_o):
    """
    Converts from polar coordinates to cartesian coordinates
    r = distance from the origin
    theta = angle from x = 0 in degrees
    x_o, y_o = the origin in cartesian space
    """
    # Convert from degrees to radians
    theta = theta*np.pi/180.
    x = r*np.sin(theta)
    y = r*np.cos(theta)
    
    # return the origin to 0,0
    x += x_o
    y += y_o
    
    return x, y

def cross_section_x(data, x_rot, x_pos, t_pos = -1):
    """
    data = data[time, z, y, x]
    x_rot = rotated x-coordinate
    x_pos = x_coordinate to do cross section at
    """
    data_xs = np.zeros_like(data[t_pos,:,:,x_pos/100])
    x_target = x_rot - x_pos
    
    my_ixs = []
    for iy in range(w.shape[2]):
        ix = np.where(abs(x_target[iy,:]) == np.min(abs(x_target[iy,:])))[0][0]
        my_ixs.append(ix)
        if x_target[iy, ix] > 0.:
            # interpolate ix and ix-1
            data_xs[:,iy] = (data[t_pos,:,iy,ix] - data[t_pos,:,iy,ix-1])*(x_pos - x_rot[iy,ix-1])/(x_rot[iy,ix] - x_rot[iy,ix-1]) + data[t_pos,:,iy,ix-1]
        elif x_target[iy, ix] > 0:
            # interpolate ix+1 and ix
            data_xs[:,iy] = (data[t_pos,:,iy,ix+1] - data[t_pos,:,iy,ix])*(x_pos - x_rot[iy,ix])/(x_rot[iy,ix+1] - x_rot[iy,ix]) + data[t_pos,:,iy,ix]
        else:
            # no need to interpolate as x_target == 0
            data_xs[:,iy] = data[t_pos,:,iy,ix]
    return data_xs, my_ixs

x_rot, y_rot = polar2cartesian(r, theta_new, x_isl, y_isl)
times = ["{0:02d}".format(t) for t in np.arange(0, 24, 3)]
xspos = 86000
count = 0
w_tmxs = np.zeros((141, 319))
print('Starting the plotting now')
for time in times:
    # Dip into some liquid water path data
    lwp_data = Dataset('lwp_'+time+'.nc', 'r')
    lwp_key = [key for key in lwp_data.variables.keys() if 'STASH' in key][0]
    lwp = lwp_data.variables[lwp_key][:]
    min1 = lwp_data.variables['min1'][:]
    lwp_data.close()

    # Grab a bit of vertical velocity data
    # Rotate a cross section out of it
    w_data = Dataset('bouy_'+time+'.nc', 'r')
    w_keys = w_data.variables.keys()

    w_key = w_keys[16]
    z_key = w_keys[9]

    w = w_data.variables[w_key][:]
    z = w_data.variables[z_key][:]
    min15 = w_data.variables['min15'][:]
    w_data.close()
    
    for it in range(len(min15)):
        w_xs, my_ixs = cross_section_x(w, x_rot, xspos, it)
        
        it1 = np.where(min1 == min15[it])[0][0]
        
        print('Making plots for ' + str(min15[it]))
        
        ### Make some plots ###
        plt.subplot(221)
        plt.contourf(x/1000., y/1000., lwp[it1,:,:], cmap = 'Greys_r', levels = np.arange(0, 1.01, 0.05), extend = 'max')
        cs1 = plt.contour(x/1000., y/1000., r/1000., levels = np.arange(10., 150.1, 10.), colors = 'r', lw = 2)
        plt.clabel(cs1, inline=1, fontsize=5, fmt = '%.0f')
        plt.ylabel('y (km)', fontsize=5)
        plt.xlabel('x (km)', fontsize=5)
        plt.title('Original Cartesian', fontsize = 7)
        
        plt.subplot(222)
        plt.contourf(x_rot/1000., y_rot/1000., lwp[it1,:,:], cmap = 'Greys_r', levels = np.arange(0, 1.01, 0.05), extend = 'max')
        cs2 = plt.contour(x_rot/1000., y_rot/1000., r/1000., levels = np.arange(10., 150.1, 10.), colors = 'r', lw = 2)
        plt.clabel(cs2,  inline=1, fontsize=5, fmt = '%.0f')
        plt.plot([x_rot[i,j]/1000. for i,j in enumerate(my_ixs)], [y_rot[i,j]/1000. for i,j in enumerate(my_ixs)], lw = 2)
        plt.ylabel('y-rotated (km)', fontsize=5)
        plt.xlabel('x-rotated (km)', fontsize=5)
        plt.title('Rotated into CT Heading Cartesian', fontsize = 7)
        
        ### Wrap around in the lwp?
        lwp_rot_wrap1 = lwp[it1,:,:]*1.
        lwp_rot_wrap1[y_rot > np.max(y)] = np.nan

        lwp_rot_wrap2 = lwp[it1,:,:]*1.
        lwp_rot_wrap2[y_rot < np.max(y)] = np.nan
        
        plt.subplot(223)
        plt.contourf(x_rot/1000., y_rot/1000., lwp_rot_wrap1, cmap = 'Greys_r', levels = np.arange(0, 1.01, 0.05), extend = 'max')
        plt.contourf(x_rot/1000., (y_rot-np.max(y))/1000., lwp_rot_wrap2, cmap = 'Greys_r', levels = np.arange(0, 1.01, 0.05), extend = 'max')
        plt.ylabel('y-rotated/wrapped (km)', fontsize=5)
        plt.xlabel('x-rotated (km)', fontsize=5)
        plt.title('Rotated into CT Heading and Wrapped Cartesian', fontsize = 7)
        
        # want to have a cross section in the y-direction at x = xspos m
        plt.subplot(224)
        plt.contourf([y_rot[i,j]/1000. for i,j in enumerate(my_ixs)], z/1000., w_xs, cmap = 'bwr', vmin = -8, vmax = 8, levels = np.arange(-8, 8.01, 0.5))
        plt.ylim([0, 10])
        plt.ylabel('Height (km)', fontsize=5)
        plt.xlabel('Rotated y (km)', fontsize=5)
        plt.colorbar()
        plt.title('Vertical Velocity (m/s) Cross section', fontsize = 7)
        
        plt.suptitle('T+' + "{0:04d}".format(int(min15[it])) + 'min', fontsize = 14)
        plt.tight_layout()
        plt.savefig('./slice_plots/rotate_wrap_around_lwp_'+"{0:04d}".format(int(min15[it]))+'.png', dpi = 300, bbox_inches = 'tight')
        if min15[it] == 1.:
            plt.show()
        plt.close('all')
        ### time mean w cross section
        w_tmxs += w_xs
        count += 1

w_tmxs /= count

plt.contourf([y_rot[i,j]/1000. for i,j in enumerate(my_ixs)], z/1000., w_tmxs, cmap = 'bwr', vmin = -1, vmax = 1)
plt.colorbar()
plt.savefig('./slice_plots/w_time_mean_cross_section_through_cloud_trail.png', dpi = 100, bbox_inches = 'tight')



