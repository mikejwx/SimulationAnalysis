## Analysis of some netcdf output from island simulations
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

times = ['00', '03', '06', '09', '12', '15', '18', '21']
for time in times:
    my_data = Dataset('lwp_'+time+'.nc', 'r')
    variables = my_data.variables.keys()
    for var in variables:
        if 'STASH_m01s30i405' in var:
            for it in range(len(my_data.variables[var][:,0,0])):
                plt.contourf(my_data.variables['longitude_t'][:], my_data['latitude_t'][:], my_data.variables[var][it, :, :], levels = np.arange(0, 1.01, 0.05), extend = 'max', cmap = 'Greys_r')
                plt.colorbar(label = 'LWP kg/m2')
                plt.title(str(my_data.variables['min1'][it]))
                plt.savefig('./lwp_plots/lwp_'+"{0:04d}".format(it+int(time)*60)+'.png', dpi = 100, bbox_inches = 'tight')
                plt.close('all')
                print it + int(time)*60
    my_data.close()

