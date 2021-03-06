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

print "Starting the wind analysis"
times = ["{0:02d}".format(x) for x in np.arange(0, 24, 3)]
for time in times:
    my_u = Dataset('u_'+time+'.nc', 'r')
    my_v = Dataset('v_'+time+'.nc', 'r')
    # Get u, and z_u
    variables_u = my_u.variables.keys()
    var_u = [var for var in variables_u if 'STASH' in var ][0]
    var_z = [var for var in variables_u if 'zsea' in var ][0]
    z_u = my_u.variables[var_z][:]
    # Get v and z_v
    variables_v = my_v.variables.keys()
    var_v = [var for var in variables_v if 'STASH' in var ][0]
    var_z = [var for var in variables_v if 'zsea' in var ][0]
    z_v = my_v.variables[var_z][:]
    for it in range(len(my_u.variables[var_u][:,0,0,0])):
        plt.subplot(121)
        lgdu, = plt.plot(hmean(my_u.variables[var_u][it,:,:,:]), z_u/1000., 
            label = 'u-wind', color = 'b', lw = 2)
        plt.plot([0,0], [0, 40], 'k--')
        plt.ylim([0, 10])
        plt.xlim([-7.5, 1])
        plt.ylabel('Height (km)')
        plt.subplot(122)
        lgdv, = plt.plot(hmean(my_v.variables[var_v][it,:,:,:]), z_v/1000., 
            label = 'v-wind', color = 'r', lw = 2)
        plt.plot([0,0], [0, 40], 'k--')
        plt.ylim([0, 10])
        plt.xlim([-7.5, 1])
        lgd = plt.legend([lgdu, lgdv], ["u-wind", "v-wind"], loc = 'upper center', 
            ncol = 2, bbox_to_anchor = (-0., -0.05))
        plt.suptitle('Horizontally Averaged Wind Profiles for T+' + str(my_u.variables['min15'][it]) + 'min', fontsize = 14)
        plt.savefig('./wind_plots/wind_'+"{0:03d}".format(it+int(time)*4)+'.png', 
            dpi = 100, bbox_to_inches = 'tight', bbox_extra_artists=(lgd,))
        plt.close('all')
