"""
What does the structure of a chunk of the cloud trail region look like?
Look at the along flow mean cloud liquid water mixing ratio, vertical velocity,
 potential temperature and in plane wind arrows/streamlines?

See Kirshbaum and Fairman (2015) figure 9(b) for inspiration
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from netCDF4 import Dataset
from analysis_tools import bilinear_interpolation, get_cs_coords, transform_winds, send_email
import os
from scipy import interpolate, integrate
from datetime import datetime as dt
from multiprocessing import Pool

f = open('../chunk_means/output.txt', 'w+')
# Calculate the points along the cross section.
# Assume that the appropriate representative wind direction is equal to the 
# initial conditions for the wind

f.write('[' + dt.now().strftime('%H:%M:%S') + '] Generating the cross section\n')
f.close()
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
z_1 = np.arange(0., 500.1, 1.)
u_1 = interpolate.interp1d(y = u_0, x = z_0, fill_value = 'extrapolate')(z_1)
v_1 = interpolate.interp1d(y = v_0, x = z_0, fill_value = 'extrapolate')(z_1)

U_0 = integrate.trapz(y = u_1, x = z_1)/500.
V_0 = integrate.trapz(y = v_1, x = z_1)/500.

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

# Get coordinates of the cross section along the flow:
h = 100. # resolution in chunk_length
x_cs, y_cs = get_cs_coords(x_c, y_c, wind_dir_0, X, Y, h = h)

f = open('../chunk_means/output.txt', 'a')
f.write('[' + dt.now().strftime('%H:%M:%S') + '] Along-wind cross section ready.\n')
f.close()
# convert island radius and along wind distance from m into km
R_i /= 1000.
R = -np.array([np.round(x, 0) for x in np.sign(x_cs - x_c)*np.sqrt((x_cs - x_c)**2 + (y_cs - y_c)**2)])

f = open('../chunk_means/output.txt', 'a')
f.write('[' + dt.now().strftime('%H:%M:%S') + '] Defining chunks downwind\n')
f.close()
# for each point in the cross section along the flow:
# only choose points that are within 5km of the exact cross section in 10 km chunks?
# makes 10 x 10 km squares
chunk_width  = 5000.
chunk_length = 20000.
res = 100. # the resolution in the chunk_width
chunks_x = {}
chunks_y = {}
chunk_Rs = {}
nchunks = int(np.max(R)/chunk_length - 1)

def get_chunk_coord(iC):
    iR = np.where(R == chunk_Rs['chunk_'+str(chunk)][iC])[0][0]
    chunks_x, chunks_y = get_cs_coords(x_cs[iR], y_cs[iR], wind_dir_0 + 90., X, Y, h = res, isPeriodic = True, max_r = chunk_width/2.)
    return chunks_x, chunks_y

for chunk in xrange(nchunks):
    f = open('../chunk_means/output.txt', 'a')
    f.write('[' + dt.now().strftime('%H:%M:%S') + '] Working on chunk ' + str(chunk) + '\n')
    f.close()
    print '[' + dt.now().strftime('%H:%M:%S') + '] Working on chunk ' + str(chunk)
    chunks_x['chunk_'+str(chunk)] = np.zeros((int(chunk_length/h + 1), int(chunk_width/res + 1))) #initialise an array for the chunk, f(along, across)
    chunks_y['chunk_'+str(chunk)] = np.zeros((int(chunk_length/h + 1), int(chunk_width/res + 1)))
    chunk_Rs['chunk_'+str(chunk)] = np.arange(chunk_length*chunk, (chunk_length + 0.1)*(chunk+1), h)
    
    p = Pool()
    tempChunk = p.map(get_chunk_coord, xrange(len(chunk_Rs['chunk_'+str(chunk)])))
    p.close()
    for iC in xrange(len(chunk_Rs['chunk_'+str(chunk)])):
        f = open('../chunk_means/output.txt', 'a')
        f.write('iC = ' + str(iC) + '\n')
        f.close()
        chunks_x['chunk_'+str(chunk)][iC,:], chunks_y['chunk_'+str(chunk)][iC,:] = tempChunk[iC] 

send_email(message = 'Finished defining the chunks', subject = 'along_wind_mean_plots.py', attachments = [''], isAttach = False)

# Above: chunks contains the coordinates for all the points needed.
f = open('../chunk_means/output.txt', 'a')
f.write('[' + dt.now().strftime('%H:%M:%S') + '] Defining keys for the netCDFs\n')
f.close()
w_key     = u'STASH_m01s00i150'
u_key     = u'STASH_m01s00i002'
v_key     = u'STASH_m01s00i003'
theta_key = u'STASH_m01s00i004'
mcl_key   = u'STASH_m01s00i392'
lwp_key   = u'STASH_m01s30i405'

#hours = ['06', '09']
hours = ['09']
#wanted_times = [480., 510., 540., 570., 600.]
wanted_times = [720.]
lwp_nc  = Dataset('../lwp_00.nc', 'r')
lwp_data = lwp_nc.variables[lwp_key][144,:,:]
lwp_nc.close()

my_cmap = mpl.cm.get_cmap('Greys')

plot_chunk = True
def interpolateMap(variable):
    return bilinear_interpolation(X, Y, variable, chunks_x[chunk_key].flatten(), chunks_y[chunk_key].flatten(), kind = 2).reshape((len(z), int(chunk_length/h + 1), int(chunk_width/res + 1)))

for hour in hours:
    send_email(message = 'Starting ' + hour, subject = 'along_wind_mean_plots.py', attachments = [''], isAttach = False)
    # read in some data 
    wind_nc = Dataset('../wind_' + hour + '.nc', 'r')
    bouy_nc = Dataset('../bouy_' + hour + '.nc', 'r')
    mr_nc   = Dataset('../mr_' + hour + '.nc', 'r')
    
    # determine which wanted_times are in this netCDF
    times = wind_nc.variables['min10_0'][:]
    its = [np.where(wanted_time == times)[0][0] for wanted_time in wanted_times if wanted_time in times]
    
    # read the variables
    f = open('../chunk_means/output.txt', 'a')
    f.write('[' + dt.now().strftime('%H:%M:%S') + '] Reading data from the netCDFs\n')
    f.close()
    z     = bouy_nc.variables['thlev_zsea_theta'][:]
    i_max = np.where(np.abs(z - 3000.) == np.min(np.abs(z - 3000.)))[0][0]
    z = z[:i_max]
    
    # do the winds first to be more time efficient
    u     = wind_nc.variables[u_key][its,:i_max,:,:] # need the time dimension for the transform winds function
    v     = wind_nc.variables[v_key][its,:i_max,:,:]
    
    # get along-flow (s) and across-flow (n) wind components
    f = open('../chunk_means/output.txt', 'a')
    f.write('[' + dt.now().strftime('%H:%M:%S') + '] Transforming the winds\n')
    f.close()
    s, n = transform_winds(u, v)
    for it in its:
        send_email(message = 'Starting ' + str(int(times[it])), subject = 'along_wind_mean_plots.py', attachments = [''], isAttach = False)
        theta = bouy_nc.variables[theta_key][it,:i_max,:,:]
        w     = wind_nc.variables[w_key][it,:i_max,:,:]
        mcl   = mr_nc.variables[mcl_key][it,:i_max,:,:]
        
        for chunk_key in chunks_x.keys():
            # interpolate the variables to our chunk coordinates
            f = open('../chunk_means/output.txt', 'a')
            f.write('[' + dt.now().strftime('%H:%M:%S') + '] Interpolating variables for ' + chunk_key + '\n')
            f.close()
            
            print '[' + dt.now().strftime('%H:%M:%S') + '] Interpolating variables for ' + chunk_key
            p = Pool()
            theta_chunk, w_chunk, n_chunk, mcl_chunk = p.map(interpolateMap, [theta, w, n[0,:,:,:], mcl])
            p.close()
            
            f = open('../chunk_means/output.txt', 'a')
            f.write('[' + dt.now().strftime('%H:%M:%S') + '] Collapsing variable chunks in the along-wind direction\n')
            f.close()
            theta_mean = np.mean(theta_chunk, axis = 1)
            w_mean = np.mean(w_chunk, axis = 1)
            n_mean = np.mean(n_chunk, axis = 1)
            mcl_mean = np.mean(mcl_chunk, axis = 1)
            
            if plot_chunk:
                f = open('../chunk_means/output.txt', 'a')
                f.write('[' + dt.now().strftime('%H:%M:%S') + '] Plot where the chunk is\n')
                f.close()
                fig = plt.figure()
                ax = fig.add_subplot(1, 1, 1, adjustable = 'box', aspect = 1)
                LWP_plt = ax.contourf(X, Y, lwp_data[:]*1000., colors = my_cmap(np.arange(0, 9.)/8.), levels = [0., 10., 20., 50., 100., 200., 300., 400., 500., 600.], extend = 'max')
                fig.colorbar(LWP_plt, ax = ax, label = r'LWP (g m$^{-2}$)')
                ax.contour(X, Y, lsm, colors = ['k'], linewidths = 2)
                ax.plot([chunks_x[chunk_key][0,0], chunks_x[chunk_key][0,-1], chunks_x[chunk_key][-1, -1], chunks_x[chunk_key][-1,0], chunks_x[chunk_key][0,0]], [chunks_y[chunk_key][0,0], chunks_y[chunk_key][0,-1], chunks_y[chunk_key][-1, -1], chunks_y[chunk_key][-1,0], chunks_y[chunk_key][0,0]], 'r')
                plt.savefig('../chunk_means/' + chunk_key + '.png', dpi = 100)
                plt.close('all')

            f = open('../chunk_means/output.txt', 'a')
            f.write('[' + dt.now().strftime('%H:%M:%S') + '] Plot along-wind chunk mean\n')
            f.close()
            y_p = np.arange(-chunk_width/2., chunk_width/2 + 0.1, res)/1000.
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            # potential temperature
            tc = ax.contour(y_p, z/1000., theta_mean, colors = ['k'], levels = range(298, 310), linewidths = [2])
            ax.clabel(tc, inline=1, fmt= '%1d')
            # vertical velocity
            ax.contourf(y_p, z/1000., w_mean, cmap = 'bwr', vmin = -2., vmax = 2., levels = [x for x in np.arange(-0.4, 2.01, 0.2) if x != 0], extend = 'max')
            ax.contour(y_p, z/1000., w_mean, colors = ['k'], levels = [x for x in np.arange(-0.4, 2.01, 0.2) if x != 0])
            # cloud liquid
            cld = ax.contourf(y_p, z/1000., mcl_mean*1000., cmap = 'Greys', levels = np.arange(0.05, 0.56, 0.05))
            fig.colorbar(cld, ax = ax, label = r'm$_{cl}$ (g kg$^{-1}$)')
            # winds
            q = ax.quiver(y_p, z[::3]/1000., n_mean[::3,:], w_mean[::3,:], scale = 25., width = 0.00125)
            ax.quiverkey(q, X = 1.1, Y = -0.1, U = 1, label = '1 m/s', labelpos = 'E')
            plt.xlabel("y' (km)")
            plt.ylabel('Height (km)')
            plt.title(str(chunk_length*chunks_x.keys().index(chunk_key)/1000) + ' to ' + str(chunk_length*(chunks_x.keys().index(chunk_key)+1)/1000) + ' km downwind of island, T+' + "{0:04d} mins".format(int(times[it])))
            plt.savefig('../chunk_means/along_wind_mean_' + chunk_key + "_T_{0:04d}".format(int(times[it])) + '.png', dpi = 100)
            plt.close('all')
            send_email(message = 'Finished along_wind_mean_' + chunk_key + "_T_{0:04d}".format(int(times[it])) + '.png', subject = 'along_wind_mean_plots.py', attachments = ['../chunk_means/along_wind_mean_' + chunk_key + "_T_{0:04d}".format(int(times[it])) + '.png'], isAttach = True)
            
        plot_chunk = False


