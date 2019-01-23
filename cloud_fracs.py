"""
Analysis of the cloud fraction with height, collapsing the mcl dataset into a 
mask of cloud vs. no cloud f(t,z,y,x) then to a horizontal mean (a.k.a. cloud 
fraction) f(t,z)

Compare the time series of liquid water path to cloud fraction total, and at 
different heights.
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from analysis_tools import get_CTZ, bilinear_interpolation, get_cs_coords
import os
from scipy import interpolate, integrate

hours = ["{0:02d}".format(h) for h in xrange(0, 24, 3)]

################################################################################
#                                                                              #
# Use liquid water mixing ratio as where clouds exist.                         #
# Part 1: Cloud fractions                                                      #
# Part 2: Cloud top heights                                                    #
#                                                                              #
################################################################################

mcl_key = u'STASH_m01s00i392'
for hour in hours:
    print hour
    # Open the netCDF for mr
    mr_nc = Dataset('../mr_'+hour+'.nc', 'r')
    # Read the data for mcl from the mr netCDF
    mcl_data = mr_nc.variables[mcl_key][:]*1.
    # Read the height coordinate
    z = mr_nc.variables[u'thlev_zsea_theta'][:]*1.
    
    if hour == '00':
        # Mask out the mcl_data
        mcl_mask = np.where((mcl_data >= 1e-08), 1., 0.)
        CTZ      = get_CTZ(mcl_data, z) # returns CTZ = f(t,y,x)
        times    = mr_nc.variables[u'min10_0'][:]/60.
    else:
        mcl_mask = np.concatenate((mcl_mask, np.where((mcl_data >= 1e-08), 1., 0.)), axis = 0)
        CTZ      = np.concatenate((CTZ, get_CTZ(mcl_data, z)), axis = 0)
        times    = np.concatenate((times, mr_nc.variables[u'min10_0'][:]/60.), axis = 0)
    mr_nc.close()

print 'Finished getting the cloud mask and cloud top height fields.'
# Make a directory for these plots
if not os.path.isdir('../cloud_fractions/'):
    os.mkdir('../cloud_fractions/')

## Full Domain ##
# Collapse mcl_mask with horizontal mean to get cloud cover
print 'Starting averaging for the whole domain'
cc_tz  = np.nanmean(mcl_mask, axis = (2, 3))*100.                      # f(t,z,y,x) -> f(t,z)
cc_t   = np.nanmean(np.nanmax(mcl_mask, axis = 1), axis = (1, 2))*100. # f(t,z,y,x) -> f(t)
ctz_t = np.nanmean(CTZ, axis = (1, 2))                                 # f(t,y,x) -> f(t)

# heights
iz1 = np.where(np.abs(z - 750.0) == np.min(np.abs(z - 750.0)))[0][0]   # near cloud base
iz2 = np.where(np.abs(z - 1900.0) == np.min(np.abs(z - 1900.0)))[0][0] # mid-cloud
iz3 = np.where(np.abs(z - 3000.0) == np.min(np.abs(z - 3000.0)))[0][0] # near cloud top

print 'Making the cloud cover time series for the whole domain'
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(times, cc_t, 'k', lw = 2, label = 'Total Cloud Cover')
ax.plot(times, cc_tz[:,iz1], 'r', label = 'Cloud Cover ('+str(round(z[iz1], 2))+'m)')
ax.plot(times, cc_tz[:,iz2], 'purple', label = 'Cloud Cover ('+str(round(z[iz2], 2))+'m)')
ax.plot(times, cc_tz[:,iz3], 'b', label = 'Cloud Cover ('+str(round(z[iz3], 2))+'m)')
ax.legend(loc = 2)
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Cloud Fraction (%)')
fig.tight_layout()
plt.savefig('../cloud_fractions/cc_timeseries_whole.png', dpi = 100)
plt.close('all')

print 'Making the cloud top height time series for the whole domain'
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(times, ctz_t, 'k', lw = 2)
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Cloud Top Height (m)')
ax.set_title('Cloud Top Height time series')
fig.tight_layout()
plt.savefig('../cloud_fractions/ctz_timeseries_whole.png', dpi = 100)
plt.close('all')

print 'Making a cloud cover hovmoller for the whole domain'
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
H = ax.contourf(times, z, np.transpose(cc_tz), cmap = 'Greys_r', levels = np.linspace(0., 1., 21))
fig.colorbar(H, ax = ax, label = 'Cloud Fraction (%)')
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Height (m)')
plt.savefig('../cloud_fractions/cf_hovmoller_whole.png', dpi = 100)
plt.close('all')

## Just the CT region ##
# Collapse mcl_mask along the wind direction to get cloud cover along the cloud trail
### Step 1: determine the heading of the cloud trail ###
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

### Step 2: define a coordinate system along that heading ###
X = np.arange(0., 116000., 100.)
Y = np.arange(0., 31900., 100.)
X, Y = np.meshgrid(X, Y)

# We know from our domain definition where the centre of the island should be
R_i = 1000.0*(50.0/np.pi)**0.5 # island radius
x_c = 100000.0 + R_i
y_c = 4*R_i

# From this center point, draw line in either direction parallel to wind_dir_0
# get the coordinates of this line at 100 m resolution.
x_cs, y_cs = get_cs_coords(x_c, y_c, wind_dir_0, X, Y, h = 100.0)
R = - np.sign(x_c - x_cs)*np.sqrt((x_cs - x_c)**2 + (y_cs - y_c)**2)

### Step 3: collapse and average our arrays along our cloud trail in various ways ###
cc_tz_ct = np.zeros_like(cc_tz)
for it in xrange(cc_tz.shape[0]):
    # bilinear_interpolation should return an array f(z, R)
    # np.where should return only the bits that are downwind
    # np.nanmean should return the mean along the R axis i.e. f(z)
    # then store back to cc_tz_ct
    # do the above for each time step
    cc_tz_ct[it,:]  = np.nanmean(np.where((R >= R_i), bilinear_interpolation(X, Y, mcl_mask[it,:,:,:], x_cs, y_cs, kind = 3)*100., np.nan), axis = 1)

mcl_total_mask = np.nanmax(mcl_mask, axis = 1)*100. # f(t,z,y,x) -> f(t,y,x)
# bilinear interpolation should then return an array f(t, R)
# np.where should return only the bits that are downwind
# np.nanmean should collapse along the downwind direction f(t, R) -> f(t)
cc_t_ct = np.nanmean(np.where((R >= R_i), bilinear_interpolation(X, Y, mcl_total_mask, x_cs, y_cs, kind = 3), np.nan), axis = 1)

ctz_t_ct = np.nanmean(np.where((R >= R_i), bilinear_interpolation(X, Y, CTZ, x_cs, y_cs, kind = 3), np.nan), axis = 1)

### Step 4: plot it out ###
print 'Making the cloud cover time series for the cloud trail region'
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(times, cc_t_ct, 'k', lw = 2, label = 'Total Cloud Cover')
ax.plot(times, cc_tz_ct[:,iz1], 'r', label = 'Cloud Cover ('+str(round(z[iz1], 2))+'m)')
ax.plot(times, cc_tz_ct[:,iz2], 'purple', label = 'Cloud Cover ('+str(round(z[iz2], 2))+'m)')
ax.plot(times, cc_tz_ct[:,iz3], 'b', label = 'Cloud Cover ('+str(round(z[iz3], 2))+'m)')
ax.legend(loc = 2)
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Cloud Fraction (%)')
fig.tight_layout()
plt.savefig('../cloud_fractions/cc_timeseries_cloudtrail.png', dpi = 100)
plt.close('all')

print 'Making the cloud top height time series for the cloud trail region'
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(times, ctz_t_ct, 'k', lw = 2)
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Cloud Top Height (m)')
ax.set_title('Cloud Top Height time series')
fig.tight_layout()
plt.savefig('../cloud_fractions/ctz_timeseries_cloudtrail.png', dpi = 100)
plt.close('all')

print 'Making a cloud cover hovmoller for the cloud trail region'
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
H = ax.contourf(times, z, np.transpose(cc_tz_ct), cmap = 'Greys_r', levels = np.linspace(0., 1., 21))
fig.colorbar(H, ax = ax, label = 'Cloud Fraction (%)')
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Height (m)')
plt.savefig('../cloud_fractions/cf_hovmoller_cloudtrail.png', dpi = 100)
plt.close('all')
