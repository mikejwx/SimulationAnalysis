"""
What does the UM think equilibrium surface fluxes are?
"""
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

days = ["{0:02d}".format(day) for day in xrange(1, 11)]
ndays = len(days)
lhf_key = u'STASH_m01s03i234'
shf_key = u'STASH_m01s03i217'
E = np.zeros(1)
H = np.zeros(1)
for day in days:
    print 'STARTING: day ' + day
    fluxes = Dataset('fluxes_' + day + '.nc', 'r')
    E = np.concatenate((E, np.mean(fluxes.variables[lhf_key][:], axis = (1, 2))))
    H = np.concatenate((H, np.mean(fluxes.variables[shf_key][:], axis = (1, 2))))
    if day == '01':
        dt_i = fluxes.variables[lhf_key].shape[0]
    fluxes.close()

times = np.arange(1., 24*ndays*60., 60.)/1440.
ax1 = plt.subplot(211)
ax1.plot(times, H[1:], lw = 2, color = 'red')
plt.ylabel('Sensible Heat Flux (Wm$^{-2}$)')

ax2 = plt.subplot(212, sharex = ax1)
ax2.plot(times, E[1:], lw = 2, color = 'blue')
plt.ylabel('Latent Heat Flux (Wm$^{-2}$)')
plt.xlabel('Time (days)')

plt.savefig('Flux_equilibrium.png', dpi = 100)
plt.show()

with open('balanced_settings.txt', 'a') as my_file:
    my_file.write('\n')
    my_file.write('Initial Conditions and Forcing for a Balanced Simulation.\n')
    my_file.write('At least the last three days of a 10-day simulation were in approximate balance between the prescribed forcing and the environmental response.\n')
    my_file.write('\n')
    my_file.write('Surface Sensible Heat Flux = ' + str(round(np.mean(H[-4*dt_i:]), 4)) + ' W/m2\n')
    my_file.write('Surface Latent Heat Flux = ' + str(round(np.mean(E[-4*dt_i:]), 4)) + ' W/m2\n')

