import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from SkewT_archer import Lv, cpd, getQ
from STASH_keys import temp_key, pthe_key

def alphaL(T, p):
    """
    A scaled rate of change in the saturation specific humidity with temperature
    """
    alpha = dqsatdT(T, p)
    alpha_L = (1. + alpha*(Lv/cpd))**(-1.)
    return alpha_L

def dqsatdT(T, p):
    """
    The gradient ofthe saturation specific humidity curve with temperature, 
    assuming at a constant pressure.
    """
    dqsatdT = getQ(T+0.5, 100., p, t_units = 'K', p_units = 'Pa') - getQ(T-0.5, 100., p, t_units = 'K', p_units = 'Pa')
    return dqsatdT

def getbs(T, p, RHcrit):
    """
    half width of our assumed triangular pdf distribution
    """
    bs = alphaL(T, p)*getQ(T, 100., p, t_units = 'K', p_units = 'Pa')*(1. - RHcrit)
    return bs

def getRHcrit(dx):
    """
    Defines the mean RHcrit as a function of grid length.
    grid length is supplied as dx = metres, and is converted to kilometres.
    """
    RH_crit = (100. - 2.38*np.log(dx/1000.) - 4.09)/100.
    return RH_crit

# Get the horizontal mean temperature and pressure profiles for the Control simulation
with Dataset('/nerc/n02/n02/xb899100/CloudTrail/Control/bouy_00.nc', 'r') as bouy_nc:
    temperature = bouy_nc.variables[temp_key][0,:-1,:,:].mean(axis = (1, 2))
    z = bouy_nc.variables['thlev_zsea_theta'][:-1]

with Dataset('/nerc/n02/n02/xb899100/CloudTrail/Control/fluxes_00.nc', 'r') as fluxes_nc:
    pressure    = fluxes_nc.variables[pthe_key][0,:-1,:,:].mean(axis = (1, 2))

fig = plt.figure()
axa = fig.add_subplot(1, 1, 1)
for DX in [50., 100., 200., 400., 800., 1600.]:
    myRHcrit = getRHcrit(DX)
    print myRHcrit
    my_bs = [getbs(temperature[i], pressure[i], myRHcrit)*1000. for i in range(len(temperature))]
    axa.plot(my_bs, z, label = DX)

axa.set_ylim([0, 3500])
axa.legend(loc = 0)
plt.show()


