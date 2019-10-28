import numpy as np
import matplotlib.pyplot as plt
from SkewT_archer import *
from scipy import interpolate, integrate

# Read the initiatial conditions
exp_names = ['spinup_800m','Control','Spinup_Control_0800m_HRIC_INV']

for exp_name in exp_names:
    spinup_data = {}
    with open('../InitialFields_Temperature_' + exp_name + '.txt', 'r') as init_temp:
        for line in init_temp:
            if 'num_theta_init_heights' in line:
                spinup_data['num_theta_init_heights'] = int(line.strip('num_theta_init_heights: '))
            elif 'theta_init_height' in line:
                spinup_data['theta_init_height'] = np.array([float(hgt) for hgt in line.strip('theta_init_height: ').split(',')])
            elif 'theta_init_data' in line:
                spinup_data['theta_init_data'] = np.array([float(data) for data in line.strip('theta_init_data: ').split(',')])
    
    with open('../InitialFields_Moisture_' + exp_name + '.txt', 'r') as init_moist:
        for line in init_moist:
            if 'num_mv_init_heights' in line:
                spinup_data['num_mv_init_heights'] = int(line.strip('num_mv_init_heights: '))
            elif 'mv_init_height' in line:
                spinup_data['mv_init_height'] = np.array([float(hgt) for hgt in line.strip('mv_init_height: ').split(',')])
            elif 'mv_init_data' in line:
                spinup_data['mv_init_data'] = np.array([float(data) for data in line.strip('mv_init_data: ').split(',')])
    
    with open('../InitialFields_Wind_' + exp_name + '.txt', 'r') as init_wind:
        for line in init_wind:
            if 'num_uv_init_heights' in line:
                spinup_data['num_uv_init_heights'] = int(line.strip('num_uv_init_heights: '))
            elif 'uv_init_height' in line:
                spinup_data['uv_init_height'] = np.array([float(hgt) for hgt in line.strip('uv_init_height: ').split(',')])
            elif 'u_init_data' in line:
                spinup_data['u_init_data'] = np.array([float(data) for data in line.strip('u_init_data: ').split(',')])
            elif 'v_init_data' in line:
                spinup_data['v_init_data'] = np.array([float(data) for data in line.strip('v_init_data: ').split(',')])
    
    p_sfc = 101700. # the surface pressure used to initialise the simulations
    spinup_data['temperature'] = spinup_data['theta_init_data'] - g*spinup_data['theta_init_height']/cpd
    spinup_data['pressure']    = np.zeros_like(spinup_data['temperature'])
    spinup_data['pressure'][0] = p_sfc*1. # set the surface pressure
    
    dz = 1.
    z = 0
    temperature_interp = interpolate.interp1d(x = spinup_data['theta_init_height'], y = spinup_data['temperature'])
    for k in range(1, spinup_data['num_theta_init_heights']):
        rho_0 = spinup_data['pressure'][k-1]/(Rd*temperature_interp(z)) # get the air density at the level below
        p_temp = spinup_data['pressure'][k-1] - g*rho_0*dz
        z += dz
        
        # use that air density to compute the pressure slightly above and iterate until just below the next level
        while (z + dz) < spinup_data['theta_init_height'][k]:
            rho_0 = p_temp/(Rd*temperature_interp(z))
            p_temp = p_temp - g*rho_0*dz
            z += dz
        
        # do the remaining distance to get to the next level
        rho_0 = p_temp/(Rd*temperature_interp(z))
        spinup_data['pressure'][k] = p_temp - g*rho_0*(spinup_data['theta_init_height'][k] - z)
        z += (spinup_data['theta_init_height'][k] - z)
    
    # iterate to get a better temperature estimate
    spinup_data['temperature'] = PTtoTemp(spinup_data['theta_init_data'], spinup_data['pressure'], t_units = 'K', p_units = 'Pa')
    
    # iterate to get a better pressure estimate
    spinup_data['pressure'][0] = p_sfc*1. # set the surface pressure
    
    dz = 1.
    z = 0
    temperature_interp = interpolate.interp1d(x = spinup_data['theta_init_height'], y = spinup_data['temperature'])
    for k in range(1, spinup_data['num_theta_init_heights']):
        rho_0 = spinup_data['pressure'][k-1]/(Rd*temperature_interp(z)) # get the air density at the level below
        p_temp = spinup_data['pressure'][k-1] - g*rho_0*dz
        z += dz
        
        # use that air density to compute the pressure slightly above and iterate until just below the next level
        while (z + dz) < spinup_data['theta_init_height'][k]:
            rho_0 = p_temp/(Rd*temperature_interp(z))
            p_temp = p_temp - g*rho_0*dz
            z += dz
        
        # do the remaining distance to get to the next level
        rho_0 = p_temp/(Rd*temperature_interp(z))
        spinup_data['pressure'][k] = p_temp - g*rho_0*(spinup_data['theta_init_height'][k] - z)
        z += (spinup_data['theta_init_height'][k] - z)
    
    # get the dewpoint
    mv_interp = interpolate.interp1d(x = spinup_data['mv_init_height'], y = spinup_data['mv_init_data'], fill_value = 'extrapolate')
    spinup_data['q']        = getQ(spinup_data['temperature'], mv_interp(spinup_data['theta_init_height'])*100., spinup_data['pressure'], t_units = 'K', p_units = 'Pa')
    spinup_data['dewpoint'] = getDew(spinup_data['q'], spinup_data['pressure'], q_units = 'kg/kg', p_units = 'Pa')
    
    # get wind functions
    u_interp = interpolate.interp1d(x = spinup_data['uv_init_height'], y = spinup_data['u_init_data'], fill_value = 'extrapolate')
    v_interp = interpolate.interp1d(x = spinup_data['uv_init_height'], y = spinup_data['v_init_data'], fill_value = 'extrapolate')
    
    plotSkewT(spinup_data['temperature']-273.15, spinup_data['dewpoint']-273.15, spinup_data['pressure']/100., u_interp(spinup_data['theta_init_height']), v_interp(spinup_data['theta_init_height']), 
                  CAPE = True, my_title = exp_name + '\n')
    plt.savefig('../Spinup_SkewT_' + exp_name + '.png', dpi = 150)
    plt.show()



