"""
Consider the wet-bulb potential temperature
"""
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from analysis_tools import get_cs_coords, bilinear_interpolation, getML_mean, lcl, get_theta_w
from SkewT_archer import getQ, PTtoTemp
from datetime import datetime as dt

print '[' + dt.now().strftime('%H:%M:%S') + '] Open some netCDFs'
bouy_nc = Dataset('../bouy_09.nc', 'r')
zi_nc = Dataset('../zi_09.nc', 'r')
fluxes_nc = Dataset('../fluxes_09.nc', 'r')

print '[' + dt.now().strftime('%H:%M:%S') + '] Define some keys to read the variables'
temp_key = u'STASH_m01s16i004'
theta_key = u'STASH_m01s00i004'
q_key = u'STASH_m01s00i010'
pres_key = u'STASH_m01s00i408'
zi_key = u'new boundary layer depth'
w_key = u'STASH_m01s00i150'

print '[' + dt.now().strftime('%H:%M:%S') + '] Read in array of heights and times'
z = bouy_nc.variables['thlev_zsea_theta'][:]
times = bouy_nc.variables['min10_0'][:]

it = len(times) - 1
print '[' + dt.now().strftime('%H:%M:%S') + '] Compute some theta_w and theta_w_s'
theta_w  = get_theta_w(bouy_nc.variables[temp_key][it,:,:,:]*1., bouy_nc.variables[q_key][it,:,:,:]*1., fluxes_nc.variables[pres_key][it,:,:,:]*1., t_units = 'K', q_units = 'kg/kg', p_units = 'Pa')
theta_ws = get_theta_w(bouy_nc.variables[temp_key][it,:,:,:]*1., getQ(bouy_nc.variables[temp_key][it,:,:,:]*1., 100., fluxes_nc.variables[pres_key][it,:,:,:]*1., t_units = 'K', p_units = 'Pa'), fluxes_nc.variables[pres_key][it,:,:,:]*1., t_units = 'K', q_units = 'kg/kg', p_units = 'Pa')

print '[' + dt.now().strftime('%H:%M:%S') + '] Compute the lcl from the mixed layer mean variables'
temperature_ML = PTtoTemp(getML_mean(bouy_nc.variables[theta_key][it,:,:,:]*1., z, zi_nc.variables[zi_key][it,:,:]*1.), fluxes_nc.variables[pres_key][it,0,:,:]/100.)
q_ML = getML_mean(bouy_nc.variables[q_key][it,:,:,:]*1., z, zi_nc.variables[zi_key][it,:,:]*1.)

print '[' + dt.now().strftime('%H:%M:%S') + '] Working on lcl for it = ' + str(it)
my_lcl = lcl(temperature_ML, q_ML, fluxes_nc.variables[pres_key][it,:,:,:], z)

print '[' + dt.now().strftime('%H:%M:%S') + '] Get the coordinates for our cross section'
x = np.arange(0., 116000., 100.)
y = np.arange(0., 31900., 100.)
X, Y = np.meshgrid(x,y)

x_c, y_c   = [108000., 15950.]
x_cs, y_cs = get_cs_coords(x_c, y_c, 80., X, Y)

theta_w_cs  = bilinear_interpolation(X, Y, theta_w, x_cs, y_cs, kind = 2)
theta_ws_cs = bilinear_interpolation(X, Y, theta_ws, x_cs, y_cs, kind = 2)
my_lcl_cs   = bilinear_interpolation(X, Y, my_lcl.reshape(1, 319, 1160), x_cs, y_cs, kind = 2)
w_cs        = bilinear_interpolation(X, Y, bouy_nc.variables[w_key][it,:,:,:], x_cs, y_cs, kind = 2)

print '[' + dt.now().strftime('%H:%M:%S') + '] Split the wet-bulb potential temperature; actual wet-bulb potential temperature below the lcl, and saturation wet-bulb potential temperature above the lcl'
theta_w_split = np.zeros_like(theta_w_cs)
theta_w_split_anom = np.zeros_like(theta_w_cs)
theta_ws_mean = np.nanmean(theta_ws, axis = (1, 2))
theta_w_mean = np.nanmean(theta_w, axis = (1, 2))

for iR in xrange(len(x_cs)):
    theta_w_split_anom[:,iR] = np.where((z > my_lcl_cs[:,iR]), theta_ws_cs[:,iR] - theta_ws_mean, theta_w_cs[:,iR] - theta_w_mean)
    theta_w_split[:,iR] = np.where((z > my_lcl_cs[:,iR]), theta_ws_cs[:,iR], theta_w_cs[:,iR])

R = - np.sign(x_cs - x_c)*np.sqrt((x_cs - x_c)**2. + (y_cs - y_c)**2.)

fig = plt.figure()
ax = fig.add_subplot(2, 1, 1)#, adjustable = 'box', aspect = 1)
tw = ax.contourf(R/1000., z/1000., theta_w_split_anom, cmap = 'bwr', levels = np.arange(-2., 2.01, 0.1), extend = 'both')
ax.plot(R/1000., my_lcl_cs[0]/1000., 'grey', lw = 2)
#ax.contour(R/1000., z/1000., w_cs, colors = ['k'], levels = np.arange(-2., 2.1, 0.1))
fig.colorbar(tw, ax = ax)
ax.set_ylim([0, 5])
#ax.set_xlim([10, 20])
ax.set_ylabel('Height (km)')

ax = fig.add_subplot(2, 1, 2)#, adjustable = 'box', aspect = 1)
tw = ax.contourf(R/1000., z/1000., theta_w_split, cmap = 'viridis', levels = np.arange(280., 300., 1.), extend = 'both')
ax.plot(R/1000., my_lcl_cs[0]/1000., 'grey', lw = 2)
#ax.contour(R/1000., z/1000., w_cs, colors = ['k'])
fig.colorbar(tw, ax = ax)
ax.set_xlabel('Distance downwind (km)')
ax.set_ylabel('Height (km)')
ax.set_ylim([0, 5])
#ax.set_xlim([10, 20])
plt.savefig('../theta_w_plot.png', dpi = 100)
plt.show()

