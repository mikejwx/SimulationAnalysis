"""
We are going to do a hovmoller diagram of different variables at different levels
along the length of the cloud trail
"""
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from netCDF4 import Dataset
from datetime import datetime as dt
from scipy import interpolate, integrate

def find_h(theta_v, u, v, z):
    """
    Iteratively solves equation 2 from Vogelezang and Holtslag (1996) to
    estimate the boundary layer depth.
    
    This method uses a modified Richardson number approach and vertical profiles
    of virtual potential temperature (theta_v), u- and v- wind components to
    estimate the boundary layer height.
    
    theta_v = f(z)
    u       = f(z)
    v       = f(z)
    """
    # Requires numpy arrays
    import numpy as np
    # Requires interpolation routines from scipy
    from scipy import interpolate
    
    # Constants
    g = 9.81
    
    # Need a first guess for the height of the boundary layer.
    h   = np.arange(1.0, 6000.1, 1.0) # m
    z_s = 0.1*h
    
    theta_v = interpolate.interp1d(x = z, y = theta_v, fill_value = 'extrapolate')
    u       = interpolate.interp1d(x = z, y = u, fill_value = 'extrapolate')
    v       = interpolate.interp1d(x = z, y = v, fill_value = 'extrapolate')
    
    Ri_g = (g/theta_v(z_s))*(theta_v(h) - theta_v(z_s))*(h - z_s)/((u(h) - u(z_s))**2. + (v(h) - v(z_s))**2.)
    
    z_i = h[np.where(np.min(np.abs(Ri_g - 0.30)) == np.abs(Ri_g - 0.30))[0][0]]
    
    return z_i

def bilinear_interpolation(x_in, y_in, z_in, x_out, y_out, kind = 0, d = 2000.0, p = 0.5, fv = np.nan):
    """
    Performs bilinear interpolation on a given 2D array, z_in.
    x_in = 2D array of x-coordinates, input
    y_in = 2D array of y-coordinates, input
    z_in = 3D array of data points, input z_in[nz, ny, nx]

    x_out = 1D array of x-coordinates to be interpolated onto
    y_out = 1D array of y-coordinates to be interpolated onto

    interpolation scheme finds the nearest points and uses scipy's bisplrep to 
    do the interpolation.
    """
    # Initialise our output array
    z_out = np.zeros((z_in.shape[0], len(x_out)))
    
    # For each point to be interpolated onto
    for i in xrange(len(x_out)):
        # Find the distance to the input data coordinates
        r = np.sqrt((x_in - x_out[i])**2 + (y_in - y_out[i])**2)
        
        if kind == 0:
            # Nearest neighbor approach
            iy, ix = np.where(r == np.min(r))
            z_out[:,i] = z_in[:,iy[0],ix[0]]
        elif kind == 1:
            # idw approach
            # Use all points within 'd' m of x_out, y_out for weighted mean
            iy, ix = np.where(r < d)
            w = 0
            
            for j in xrange(len(iy)):
                w += (1./r[iy[j], ix[j]]**p)
                z_out[:,i] += (1./r[iy[j], ix[j]]**p)*z_in[:,iy[j], ix[j]]
            
            z_out[:,i] /= w
        elif kind == 2:
            # Find the nearest point
            iy, ix = np.where(r == np.min(r))
            iy = iy[0]
            ix = ix[0]
            # Determine which quadrant (i.e. up and left, up and right, down and
            # right, or down and left) this point is with respect to the output 
            # point
            dx = x_in[iy, ix] - x_out[i] # if +ve, input point is to the right
            dy = y_in[iy, ix] - y_out[i] # if +ve, input point is up
            if np.min([j in xrange(z_in.shape[1]) for j in [iy-1, iy, iy+1]]) and np.min([j in xrange(z_in.shape[2]) for j in [ix-1, ix, ix+1]]):
                if dx > 0:
                    # nearest input point is to the right
                    if dy > 0:
                        # nearest input point is up
                        z_out_up = (z_in[:,iy,ix] - z_in[:,iy,ix-1])*(x_out[i] - x_in[iy,ix-1])/(x_in[iy,ix] - x_in[iy,ix-1]) + z_in[:,iy,ix-1]
                        z_out_down = (z_in[:,iy-1,ix] - z_in[:,iy-1,ix-1])*(x_out[i] - x_in[iy-1,ix-1])/(x_in[iy-1,ix] - x_in[iy-1,ix-1]) + z_in[:,iy-1,ix-1]
                        y_in_up = (y_in[iy,ix] - y_in[iy,ix-1])*(x_out[i] - x_in[iy,ix-1])/(x_in[iy,ix] - x_in[iy,ix-1]) + y_in[iy,ix-1]
                        y_in_down = (y_in[iy-1,ix] - y_in[iy-1,ix-1])*(x_out[i] - x_in[iy-1,ix-1])/(x_in[iy-1,ix] - x_in[iy-1,ix-1]) + y_in[iy-1,ix-1]
                        z_out[:,i] = (z_out_up - z_out_down)*(y_out[i] - y_in_down)/(y_in_up - y_in_down) + z_out_down
                    elif dy <= 0:
                        # nearest input point is down
                        z_out_up = (z_in[:,iy+1,ix] - z_in[:,iy+1,ix-1])*(x_out[i] - x_in[iy+1,ix-1])/(x_in[iy+1,ix] - x_in[iy+1,ix-1]) + z_in[:,iy+1,ix-1]
                        z_out_down = (z_in[:,iy,ix] - z_in[:,iy,ix-1])*(x_out[i] - x_in[iy,ix-1])/(x_in[iy,ix] - x_in[iy,ix-1]) + z_in[:,iy,ix-1]
                        y_in_up = (y_in[iy+1,ix] - y_in[iy+1,ix-1])*(x_out[i] - x_in[iy+1,ix-1])/(x_in[iy+1,ix] - x_in[iy+1,ix-1]) + y_in[iy+1,ix-1]
                        y_in_down = (y_in[iy,ix] - y_in[iy,ix-1])*(x_out[i] - x_in[iy,ix-1])/(x_in[iy,ix] - x_in[iy,ix-1]) + y_in[iy,ix-1]
                        z_out[:,i] = (z_out_up - z_out_down)*(y_out[i] - y_in_down)/(y_in_up - y_in_down) + z_out_down
                elif dx <= 0:
                    # nearest input point is to the left
                    if dy > 0:
                        # nearest input point is up
                        z_out_up = (z_in[:,iy,ix+1] - z_in[:,iy,ix])*(x_out[i] - x_in[iy,ix])/(x_in[iy,ix+1] - x_in[iy,ix]) + z_in[:,iy,ix]
                        z_out_down = (z_in[:,iy-1,ix+1] - z_in[:,iy-1,ix])*(x_out[i] - x_in[iy-1,ix])/(x_in[iy-1,ix+1] - x_in[iy-1,ix]) + z_in[:,iy-1,ix]
                        y_in_up = (y_in[iy,ix+1] - y_in[iy,ix])*(x_out[i] - x_in[iy,ix])/(x_in[iy,ix+1] - x_in[iy,ix]) + y_in[iy,ix]
                        y_in_down = (y_in[iy-1,ix+1] - y_in[iy-1,ix])*(x_out[i] - x_in[iy-1,ix])/(x_in[iy-1,ix+1] - x_in[iy-1,ix]) + y_in[iy-1,ix]
                        z_out[:,i] = (z_out_up - z_out_down)*(y_out[i] - y_in_down)/(y_in_up - y_in_down) + z_out_down
                    elif dy <= 0:
                        # nearest input point is down
                        z_out_up = (z_in[:,iy+1,ix+1] - z_in[:,iy+1,ix])*(x_out[i] - x_in[iy+1,ix])/(x_in[iy+1,ix+1] - x_in[iy+1,ix]) + z_in[:,iy+1,ix]
                        z_out_down = (z_in[:,iy,ix+1] - z_in[:,iy,ix])*(x_out[i] - x_in[iy,ix])/(x_in[iy,ix+1] - x_in[iy,ix]) + z_in[:,iy,ix]
                        y_in_up = (y_in[iy+1,ix+1] - y_in[iy+1,ix])*(x_out[i] - x_in[iy+1,ix])/(x_in[iy+1,ix+1] - x_in[iy+1,ix]) + y_in[iy+1,ix]
                        y_in_down = (y_in[iy,ix+1] - y_in[iy,ix])*(x_out[i] - x_in[iy,ix])/(x_in[iy,ix+1] - x_in[iy,ix]) + y_in[iy,ix]
                        z_out[:,i] = (z_out_up - z_out_down)*(y_out[i] - y_in_down)/(y_in_up - y_in_down) + z_out_down
            else:
                z_out[:,i] = fv
        elif kind == 3:
            iy, ix = np.where(r < d)
            z_out[:,i] = np.nanmean(z_in[:,iy,ix], axis = 1)
        
    return z_out

# Calculate the points along the cross section.
# Assume that the appropriate representative wind direction is equal to the 
# initial conditions for the wind

# Initial conditions taken directly from namelist
u_0 = np.array([-6.09,-7.02,-7.53,-7.89,-8.15,-8.36,-8.53,-8.68,-8.79,-8.89,
                -8.97,-9.02,-9.07,-9.1,-9.12,-9.13,-9.14,-9.15,-9.16,-9.16,
                -9.17,-9.19,-9.21,-9.24,-9.29,-9.35,-9.43,-9.54,-9.68,-9.84,
                -9.98,-10.09,-10.14,-10.15,-10.13,-10.1,-10.08,-10.06,-10.05,
                -10.04,-10.01,-10.0,-10.0,-10.0,-9.66,-0.09,0.0,0.0])
v_0 = np.array([-1.2,-1.38,-1.48,-1.55,-1.6,-1.64,-1.67,-1.7,-1.72,-1.74,
                -1.76,-1.78,-1.8,-1.81,-1.83,-1.84,-1.85,-1.86,-1.86,-1.85,
                -1.84,-1.82,-1.79,-1.76,-1.73,-1.69,-1.65,-1.61,-1.55,-1.45,
                -1.3,-1.08,-0.83,-0.62,-0.46,-0.34,-0.26,-0.2,-0.15,-0.11,-0.03,
                -0.01,0.0,0.0,0.0,0.0,0.0,0.0])
z_0 = np.array([1.0000004,3.6666676,7.666668,13.000004,19.666672,27.666672,
                37.000008,47.66668,59.66668,73.00004,87.66668,103.66668,
                121.00004,139.66672,159.66672,181.00004,203.66672,227.66672,
                337.00008,367.66676,399.66676,433.0,467.6668,503.6668,541.0,
                579.6668,619.6668,661.0,703.6668,747.6668,793.0004,839.6668,
                887.6668,937.0004,987.6668,1039.6668,1093.0004,1147.6672,
                1203.6668,1261.0004,1503.6672,1567.6668,1633.0004,8955.796,
                9205.932,14947.828,15802.464,15802.464])

# Calculate the mean wind direction in the boundary layer (lowest ~ 850 m)
z_1 = np.arange(0., 850.1, 1.)
u_1 = interpolate.interp1d(y = u_0, x = z_0, fill_value = 'extrapolate')(z_1)
v_1 = interpolate.interp1d(y = v_0, x = z_0, fill_value = 'extrapolate')(z_1)

U_0 = integrate.trapz(y = u_1, x = z_1)/850.
V_0 = integrate.trapz(y = v_1, x = z_1)/850.

wind_speed_0 = np.sqrt(U_0**2 + V_0**2)
# Wind speed should be ~ 9.5 m/s
wind_dir_0 = 360.*np.arctan(U_0/V_0)/(2.*np.pi)
# Wind direction should be ~ 80 deg.

landseamask = Dataset('/work/n02/n02/xb899100/island_masks/lsm50.nc', 'r')

lsm = landseamask.variables['lsm'][0,0,:,:]*1.
y = landseamask.variables['latitude'][:]*1.
x = landseamask.variables['longitude'][:]*1.

# create 2D coordinate mesh
X, Y = np.meshgrid(x, y)

landseamask.close()

# We know from our domain definition where the centre of the island should be
R_i = 1000.0*(50.0/np.pi)**0.5 # island radius
x_c = 100000.0 + R_i
y_c = 4*R_i

# From this center point, draw line in either direction parallel to wind_dir_0
# get the coordinates of this line at 100 m resolution.
h  = 100.0
dx = (h*np.sin(np.pi*wind_dir_0/180.0))
dy = (h*np.cos(np.pi*wind_dir_0/180.0))

# start lists of coordinates from the centre point of the island

x_cs = [x_c]
y_cs = [y_c]

# Populate my list of x and y coordinates along the cross section
while (x_cs[-1] < np.max(x))*(y_cs[-1] < np.max(y)):
    x_cs.append(x_cs[-1] + dx)
    y_cs.append(y_cs[-1] + dy)

x_cs = x_cs[::-1]
y_cs = y_cs[::-1]

while (x_cs[-1] > np.min(x))*(y_cs[-1] > np.min(y)):
    x_cs.append(x_cs[-1] - dx)
    y_cs.append(y_cs[-1] - dy)

# check that all the coordinates along the cross section are within the domain
i = 0
while i < len(x_cs):
    if (x_cs[i] > np.max(x)) or (x_cs[i] < np.min(x)) or (y_cs[i] > np.max(y)) or (y_cs[i] < np.min(y)):
        del x_cs[i]
        del y_cs[i]
    else:
        i += 1

# store to arrays
x_cs = np.array(x_cs)
y_cs = np.array(y_cs)

hours = ["{0:02d}".format(h) for h in xrange(0, 24, 3)]
theta_key = u'STASH_m01s00i004'
q_key = u'STASH_m01s00i010'

R_i /= 1000.
R = -np.sign(x_cs - x_c)*np.sqrt((x_cs - x_c)**2 + (y_cs - y_c)**2)/1000.
times = np.arange(10., 1440.1, 10.)/60.
for iz in xrange(126):
    for hour in hours:
        # read in the data
        bouy_nc = Dataset('../bouy_'+hour+'.nc', 'r')
        if hour == '00':
            # skip the first time step to avoid issue with moisture resetting
            it_start = 1
            hovmoller = np.zeros((18*len(hours), len(x_cs)))
            z = bouy_nc.variables['thlev_zsea_theta'][:]*1.
        else:
            it_start = 0
        hovmoller[18*hours.index(hour):(18*hours.index(hour) + 18),:] = bilinear_interpolation(X, Y, bouy_nc.variables[q_key][it_start:,iz,:,:], x_cs, y_cs, kind = 3)
        bouy_nc.close()
    print 'Making plot for height ' + str(int(z[iz])) + ' m'
    cbar_min = round(1000.*(np.nanmean(hovmoller) - 2*np.std(hovmoller)), 1)
    cbar_max = round(1000.*(np.nanmean(hovmoller) + 2*np.std(hovmoller)), 1)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    hov = ax.contourf(R, times, hovmoller*1000., cmap = 'BrBG', levels = np.linspace(cbar_min, cbar_max, 21), extend = 'both')
    plt.colorbar(hov, ax = ax, label = r'$q_v$ (K)')
    ax.plot([-R_i, -R_i], [0, 24], 'k--')
    ax.plot([R_i,R_i], [0, 24], 'k--')
    ax.set_yticks(np.arange(0., 24.1, 3.))
    ax.set_ylabel('Time (hrs)')
    ax.set_xlabel('Distance Downwind (km)')
    ax.set_title(r'Hovmoller of $q_v$ at ' + "{0:04d}".format(int(z[iz])) + 'm')
    plt.savefig('../Hovmoller_qv_'+"{0:04d}".format(int(z[iz]))+'m.png', dpi = 100)
    plt.show()


