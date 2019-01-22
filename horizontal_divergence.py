# Code to plot horizontal divergence on top of vertical velocities.
import numpy as np
import matplotlib.pyplot as plt
from analysis_tools import ddx, ddy, circular_smoothing
from netCDF4 import Dataset

# Wind keys
u_key = u'STASH_m01s00i002'
v_key = u'STASH_m01s00i003'
w_key = u'STASH_m01s00i150'

x_u = np.arange(50., 116000., 100.) # x-coordinates for the u-grid
y_v = np.arange(50., 31900., 100.) # y-coordinates for the v-grid
x   = np.arange(0., 116000., 100.) # x-coordinates for the centre points
y   = np.arange(0., 31900., 100.) # y-coordinates for the centre points

X_u, Y_u = np.meshgrid(x_u, y)
X_v, Y_v = np.meshgrid(x, y_v)
X, Y     = np.meshgrid(x, y)

# read in the wind data
hours = ["{0:02d}".format(hour) for hour in xrange(0, 24, 3)]
#for hour in hours:
it = -1
iz = 35 # defined for theta levels, needs to be iz-1 for rho levels to get the 
        # nearest height below the theta levels

hour = '09'
u_nc = Dataset('../u_'+hour+'.nc', 'r')
v_nc = Dataset('../v_'+hour+'.nc', 'r')
w_nc = Dataset('../bouy_'+hour+'.nc', 'r')

z_theta = w_nc.variables['thlev_zsea_theta'][:]
z_rho   = u_nc.variables['rholev_zsea_rho'][:]

times = u_nc.variables['min10_0'][:]

# We use an Arakawa-C grid which means that the u wind is staggered in the 
# x-direction, and the v wind is staggered in the y-direction from the theta
# points grid.
# This means that the gradients we return from ddx() and ddy() will be on
# the same grid as the w points which are on the theta points grid.
DIV = -(ddx(u_nc.variables[u_key][it,iz-1,:,:], X_u) + ddy(v_nc.variables[v_key][it,iz-1,:,:], Y_v))
w   = w_nc.variables[w_key][it,iz,:,:]

u_nc.close()
v_nc.close()
w_nc.close()

w_vmin = 6*np.std(w)
w_vmax = -w_vmin
from datetime import datetime as dt
print '[' + dt.now().strftime("%H:%M:%S") + ']'
start = dt.now()
w_smoothed = circular_smoothing(X, Y, data = w, r = 1000.)
DIV_smoothed = circular_smoothing(X, Y, data = DIV, r = 1000.)
end = dt.now()
print '[' +  dt.now().strftime("%H:%M:%S") + '] Smoothing took ' + str(end - start)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
W = ax.contourf(X, Y, w_smoothed, cmap = 'bwr', vmin = w_vmin, vmax = w_vmax, extend = 'both', levels = np.linspace(-2., 2., 21))
ax.contour(X, Y, DIV, colors = ['k'], levels = np.linspace(-0.008, 0.008, 11))
ax.contour(X, Y, DIV, colors = ['k'], levels = [0], linewidths = [3])
fig.colorbar(W, ax = ax)
plt.show()



