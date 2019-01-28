#Plot skewT diagrams every 5 km downwind of the island.
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from scipy import interpolate, integrate
from netCDF4 import Dataset
from analysis_tools import bilinear_interpolation, get_cs_coords, send_email
from SkewT_archer import *
import os

def get_soundings(it):
    """
    Plots a swath of soundings for each time 'it' that is input
    """
    # For the skewT plots
    # Get a sounding every 5 km downwind
    h  = 100.0 # copied from above because for some reason, we cannot see this here.
    i_cs = np.where(R%5000. < h)[0]
    temp_cs_soundings = bilinear_interpolation(X, Y, bouy_nc.variables[temp_key][it,:,:,:]*1., x_cs[i_cs], y_cs[i_cs], kind = 3) # mean within a radius
    pres_cs_soundings = bilinear_interpolation(X, Y, fluxes_nc.variables[pres_key][it,:,:,:]*1., x_cs[i_cs], y_cs[i_cs], kind = 3)
    q_cs_soundings    = bilinear_interpolation(X, Y, bouy_nc.variables[q_key][it,:,:,:]*1., x_cs[i_cs], y_cs[i_cs], kind = 3)
    u_cs_soundings    = bilinear_interpolation(X, Y, wind_nc.variables[u_key][it,:,:,:]*1., x_cs[i_cs], y_cs[i_cs], kind = 3)
    v_cs_soundings    = bilinear_interpolation(X, Y, wind_nc.variables[v_key][it,:,:,:]*1., x_cs[i_cs], y_cs[i_cs], kind = 3)
    dew_cs_soundings  = getDew(q_cs_soundings, pres_cs_soundings/100.)

    # Ensure our soundings directory exists, if not create it
    if not os.path.isdir('../Soundings/'):
        os.mkdir('../Soundings/')

    #Plot the soundings as well
    for i in xrange(len(i_cs)):
        if not os.path.isdir('../Soundings/Radius_'+"{0:02d}".format(int(R[i_cs[i]]/1000.))+'/'):
            os.mkdir('../Soundings/Radius_'+"{0:02d}".format(int(R[i_cs[i]]/1000.))+'/')
        
        plotSkewT(temp_cs_soundings[:-1,i]-273.15, dew_cs_soundings[:-1,i]-273.15, pres_cs_soundings[:-1,i]/100., u_cs_soundings[:-1,i], v_cs_soundings[:-1,i], my_title = 'T+' + str(int(times[it])) + 'mins, ' + str(int(R[i_cs[i]]/1000.)) + 'km downwind', CAPE = True)
        plt.savefig('../Soundings/Radius_'+"{0:02d}".format(int(R[i_cs[i]]/1000.))+'/AlongWind_sounding_T'+"{0:04d}".format(int(times[it]))+'_R'+"{0:02d}".format(int(R[i_cs[i]]/1000.))+'.png', dpi = 100)
        plt.close('all')
        print "{0:02d}".format(int(R[i_cs[i]]/1000.)) + 'km Sounding Complete.'

print 'Getting the along-wind cloud trail coordinates'
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

print 'Taking a boundary layer mean wind'
# Calculate the mean wind direction in the boundary layer (lowest ~ 850 m)
z_1 = np.arange(0., 850.1, 1.)
u_1 = interpolate.interp1d(y = u_0, x = z_0, fill_value = 'extrapolate')(z_1)
v_1 = interpolate.interp1d(y = v_0, x = z_0, fill_value = 'extrapolate')(z_1)

U_0 = integrate.trapz(y = u_1, x = z_1)/850.
V_0 = integrate.trapz(y = v_1, x = z_1)/850.

print 'Calculating the wind speed and direction'
wind_speed_0 = np.sqrt(U_0**2 + V_0**2)
# Wind speed should be ~ 9.5 m/s
wind_dir_0 = 360.*np.arctan(U_0/V_0)/(2.*np.pi)
# Wind direction should be ~ 80 deg.

### Step 2: define a coordinate system along that heading ###
print 'Defining a coordinate system along the wind'
X = np.arange(0., 116000., 100.)
Y = np.arange(0., 31900., 100.)
X, Y = np.meshgrid(X, Y)

# We know from our domain definition where the centre of the island should be
print 'Determining where the island is'
R_i = 1000.0*(50.0/np.pi)**0.5 # island radius
x_c = 100000.0 + R_i
y_c = 4*R_i

# From this center point, draw line in either direction parallel to wind_dir_0
# get the coordinates of this line at 100 m resolution.
print 'Generating coordinates along our wind direction and across the island centre'
x_cs, y_cs = get_cs_coords(x_c, y_c, wind_dir_0, X, Y, h = 100.0)
R = np.sign(x_c - x_cs)*np.sqrt((x_cs - x_c)**2 + (y_cs - y_c)**2)

hours = ["{0:02d}".format(h) for h in xrange(0, 24, 3)]
for hour in hours:
    # read in the netCDFs
    bouy_nc   = Dataset('../bouy_' + hour + '.nc', 'r')
    fluxes_nc = Dataset('../fluxes_' + hour + '.nc', 'r')
    wind_nc   = Dataset('../wind_' + hour + '.nc', 'r')
    # keys
    temp_key  = u'STASH_m01s16i004' # bouy.nc
    pres_key  = u'STASH_m01s00i408' # fluxes.nc
    q_key     = u'STASH_m01s00i010' # bouy.nc
    u_key     = u'STASH_m01s00i002' # wind.nc
    v_key     = u'STASH_m01s00i003' # wind.nc
    
    z         = bouy_nc.variables[u'thlev_zsea_theta'][:]*1.
    times     = bouy_nc.variables[u'min10_0'][:]*1.

    for I in xrange(len(times)):
        if not ((hour == '00') and (I == 0)):
            print 'Hour = ' + hour + ', time = ' + str(int(times[I]))
            get_soundings(I)

    bouy_nc.close()
    fluxes_nc.close()
    wind_nc.close()

send_email(message = "Finished plotting downwind soundings!", subject = "downwind_soundings.py", attachments = ['../Soundings/Radius_10/AlongWind_sounding_T0720_R10.png', '../Soundings/Radius_20/AlongWind_sounding_T0720_R10.png'], isAttach = True)
