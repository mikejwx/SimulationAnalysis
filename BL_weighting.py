"""
Code to estimate the profile of blending coefficient for the boundary layer
scheme used in our simulations.
"""
import numpy as np
import matplotlib.pyplot as plt

<<<<<<< HEAD
colors = ['blue', 'green', 'red', 'cyan', 'purple', 'gold']
dxs = [50., 100., 200., 400., 800., 1600.]
fig = plt.figure(figsize = (12,8))
ax = fig.add_subplot(1, 1, 1)
for par in [500., 700.]:
    z = np.arange(0., 40000.1, 1.) # z-coordinate
    zh = par*1. # boundary layer height
    beta_bl = 0.15
    beta_fa = 1.
    zfa = zh+1000. # free atmosphere height
    beta = np.where(z < zh, beta_bl, np.where(z > zfa, beta_fa, beta_bl*(zfa-z)/(zfa-zh) + beta_fa*(z-zh)/(zfa-zh)))
    l0 = 4.
    l1 = 0.25
    rf = 1./(l0 + l1)
    
    zsml = par*1. # depth of surface mixed layer
    zsc  = par*1. # depth of stratocumulus mixed layer (0 = no stratocumulus?)
    
    #zturb = np.array([np.min([np.max([z[i], zsml]), np.max([zsc, zh-z[i]])]) for i in range(len(z))])
    zturb = par*1.
    for dx in dxs:
        print 'Grid Spacing = ' + str(dx)
        W1d = 1. - np.tanh(beta*zturb/dx)*np.max([0., np.min([1., rf*(l0 - dx/zturb)])])
        print 'Max W1D = ' + str(W1d.max()) + '\nMin W1D = ' + str(W1d.min())
        
        ax.plot(W1d, z, label = 'dx = ' + str(dx) + u', z$_{turb}$ = ' + str(int(zturb)) + ' m', color = colors[dxs.index(dx)], ls = ['--' if zturb == 500. else '-'][0], lw = 2)

plt.legend(loc='center left', bbox_to_anchor=(1.025, 0.5), ncol=1)
ax.set_ylim([0, 3000])
ax.set_ylabel('height (m)', fontsize = 15)
ax.set_xlabel(u'W$_{1D}$', fontsize = 15)
plt.suptitle(u'W$_{1D}$ = 1 - tanh$(\\beta \\frac{z_{turb}}{\Delta x} )$max$[0,$min$[1,r_{f}(l_{0} - \\frac{\Delta x}{z_{turb}})]]$', fontsize = 15)
plt.subplots_adjust(right = 0.67)
plt.savefig('../BL_weighting.png', dpi = 150)
plt.show()

=======
z = np.arange(0., 40000.1, 1.) # z-coordinate
zh = 700. # boundary layer height
beta_bl = 0.15
beta_fa = 1.
zfa = zh+1000. # free atmosphere height
beta = np.where(z < zh, beta_bl, np.where(z > zfa, beta_fa, beta_bl*(zfa-z)/(zfa-zh) + beta_fa*(z-zh)/(zfa-zh)))
l0 = 4.
l1 = 0.25
rf = 1./(l0 + l1)

zsml = 800. # depth of surface mixed layer
zsc = 700. # depth of stratocumulus mixed layer (0 = no stratocumulus?)

zturb = np.where(np.where(z >= zsml, z, zsml) <= np.where(zsc >= zh - z, zsc, zh - z), np.where(z >= zsml, z, zsml), np.where(zsc >= zh - z, zsc, zh - z))

dx = 100. # model horizontal grid spacing
W1d_100 = 1. - np.tanh(beta*zturb/dx)*np.where(0. >= np.where(1. <= rf*(l0 - dx/zturb), 1., rf*(l0 - dx/zturb)), 0, np.where(1. <= rf*(l0 - dx/zturb), 1., rf*(l0 - dx/zturb)))
dx = 1000.
W1d_1000 = 1 - np.tanh(beta*zturb/dx)*np.where(0. >= np.where(1. <= rf*(l0 - dx/zturb), 1., rf*(l0 - dx/zturb)), 0, np.where(1. <= rf*(l0 - dx/zturb), 1., rf*(l0 - dx/zturb)))

plt.plot(W1d_100, z)
plt.plot(W1d_1000, z)
plt.ylim([0, 3000])
plt.show()
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696

