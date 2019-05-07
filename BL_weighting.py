"""
Code to estimate the profile of blending coefficient for the boundary layer
scheme used in our simulations.
"""
import numpy as np
import matplotlib.pyplot as plt

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

