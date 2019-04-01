"""
Analysis of a 10-day spin up simulation to arrive at balanced winds and surface
fluxes for our cloud trail simulations.

This script replaces 'temperature_profiles.py', 'surface_fluxes.py', and 
'hodographs.py' which are now depreciated. Code from those three scripts are 
combined using more efficient methods and code in general to improve performance
for the same results.

In addition, changes have been made to improve code readability.
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import integrate, interpolate
from analysis_tools import RDP, find_h, lcl, get_CC, get_CTZ, regrid
from SkewT_archer import getDew, PTtoTemp, getQ, Rd, cpd, g

# Toggle print statements
l_verbose = True

# Toggle thermodynamics
l_thermodynamics = 0

# Toggle dynamics
l_dynamics = 0

# Toggle surface fluxes
l_fluxes = 1

# file ID
#ID = 'Control'
#ID = 'U05'
ID = 'U05v2'

# Change the path to the data
#data_path = '/nerc/n02/n02/xb899100/CloudTrail/Control_Spinup/'
#data_path = '/nerc/n02/n02/xb899100/CloudTrail/U05_Spinup/'
data_path = '/work/n02/n02/xb899100/cylc-run/u-bg952/share/data/history/'

### Updated from 'temperature_profiles.py' ###
# Define a list of days to take advantage of regular expression in naming
days = ["{0:02d}".format(d) for d in xrange(1, 11)]
ndays = len(days)

# Define keys for the variables we want
theta_key    = u'STASH_m01s00i004'
pressure_key = u'STASH_m01s00i407'
q_key        = u'STASH_m01s00i010'
mv_key       = u'STASH_m01s00i391'
u_key        = u'STASH_m01s00i002'
v_key        = u'STASH_m01s00i003'

if l_thermodynamics:
    for day in days:
        if l_verbose: print 'Starting day: ' + day
        
        # Open the netCDF
        bouy_nc   = Dataset(data_path + 'bouy_' + day + '.nc', 'r')
        fluxes_nc = Dataset(data_path + 'fluxes_' + day + '.nc', 'r')
        mr_nc     = Dataset(data_path + 'mr_' + day + '.nc', 'r')
        u_nc      = Dataset(data_path + 'u_' + day + '.nc', 'r')
        v_nc      = Dataset(data_path + 'v_' + day + '.nc', 'r')
        
        theta_data      = bouy_nc.variables[theta_key][:]
        #regrid the pressures and specific humidity
        pressure_regrid = regrid(bouy_nc, fluxes_nc, pressure_key)
        q_regrid        = regrid(bouy_nc, mr_nc, q_key)
        mv_regrid       = regrid(bouy_nc, mr_nc, mv_key)
        u_regrid        = regrid(bouy_nc, u_nc, u_key)
        v_regrid        = regrid(bouy_nc, v_nc, v_key)
        
        if day == '01':
            if l_verbose: print 'Initialising arrays for data in time and height'
            temperature = np.zeros((theta_data.shape[0]*ndays, theta_data.shape[1]))
            pressure_rg = np.zeros((theta_data.shape[0]*ndays, theta_data.shape[1]))
            dewpoint_rg = np.zeros((theta_data.shape[0]*ndays, theta_data.shape[1]))
            u_rg        = np.zeros((theta_data.shape[0]*ndays, theta_data.shape[1]))
            v_rg        = np.zeros((theta_data.shape[0]*ndays, theta_data.shape[1]))
            q_rg        = np.zeros((theta_data.shape[0]*ndays, theta_data.shape[1]))
            mv_rg       = np.zeros((theta_data.shape[0]*ndays, theta_data.shape[1]))
            rh_rg       = np.zeros((theta_data.shape[0]*ndays, theta_data.shape[1]))
            theta_rg    = np.zeros((theta_data.shape[0]*ndays, theta_data.shape[1]))
            z           = bouy_nc.variables['thlev_zsea_theta'][:]
        
        # Close netCDF
        bouy_nc.close()
        fluxes_nc.close()
        mr_nc.close()
        u_nc.close()
        v_nc.close()
        
        # Convert from potential temperature to temperature
        temp = PTtoTemp(theta_data, pressure_regrid, t_units = 'K', p_units = 'Pa')
        
        # Horizontally average the temperatures
        temperature[theta_data.shape[0]*days.index(day):theta_data.shape[0]*(days.index(day)+1),:] = np.nanmean(temp, axis = (2, 3))
        pressure_rg[theta_data.shape[0]*days.index(day):theta_data.shape[0]*(days.index(day)+1),:] = np.nanmean(pressure_regrid, axis = (2, 3))
        dewpoint_rg[theta_data.shape[0]*days.index(day):theta_data.shape[0]*(days.index(day)+1),:] = np.nanmean(getDew(q_regrid, pressure_regrid, q_units = 'kg/kg', p_units = 'Pa'), axis = (2, 3))
        u_rg[theta_data.shape[0]*days.index(day):theta_data.shape[0]*(days.index(day)+1),:]        = np.nanmean(u_regrid, axis = (2, 3))
        v_rg[theta_data.shape[0]*days.index(day):theta_data.shape[0]*(days.index(day)+1),:]        = np.nanmean(v_regrid, axis = (2, 3))
        q_rg[theta_data.shape[0]*days.index(day):theta_data.shape[0]*(days.index(day)+1),:]        = np.nanmean(q_regrid, axis = (2, 3))
        mv_rg[theta_data.shape[0]*days.index(day):theta_data.shape[0]*(days.index(day)+1),:]       = np.nanmean(mv_regrid, axis = (2, 3))
        rh_rg[theta_data.shape[0]*days.index(day):theta_data.shape[0]*(days.index(day)+1),:]       = np.nanmean(q_regrid/getQ(temp, [100.], pressure_regrid, t_units = 'K', p_units = 'Pa'), axis = (2, 3))
        theta_rg[theta_data.shape[0]*days.index(day):theta_data.shape[0]*(days.index(day)+1),:]    = np.nanmean(theta_data, axis = (2, 3))

    # Time output every ten minutes (blindly manufactured)
    times = np.arange(1., 14400.*len(days), 10.)/60.

    dt_i = theta_data.shape[0] # Number of time steps per day
    if l_verbose: print 'Starting Temperature.'
    with open('../InitialFields_Temperature_' + ID + '.txt', 'w') as my_file:
        my_file.write('Specify initial temperature profiles\n')
        
        # Find minimum number of required levels to reproduce theta profile over the last four days
        theta_levels, theta_init = RDP(z, np.mean(theta_rg[-4*dt_i:,:], axis = 0), 0.1)
        n_thlev = len(theta_levels)
        
        # First namelist entry
        my_file.write('num_theta_init_heights: ')
        my_file.write(str(n_thlev) + '\n')
        
        # Second namelist entry
        my_file.write('theta_init_height: ')
        for k in xrange(n_thlev-1):
            my_file.write(str(theta_levels[k]) + ',')
        my_file.write(str(theta_levels[-1]) + '\n')
        
        # Third namelist entry
        my_file.write('theta_init_field_type: (10) dry potential temperature\n')
        
        # Fourth namelist entry
        my_file.write('theta_init_data: ')
        for k in xrange(n_thlev-1):
            my_file.write(str(theta_init[k]) + ',')
        my_file.write(str(theta_init[-1]) + '\n')

    if l_verbose: print 'Finished Temperature\n Starting Moisture.'
    with open('../InitialFields_Moisture_' + ID + '.txt', 'w') as my_file:
        my_file.write('Specify initial moisture profiles\n')
        
        # Find the minimum number of required levels to reproduce the mv profile
        mv_levels, mv_init = RDP(z, np.mean(rh_rg[-4*dt_i:, :], axis = 0), 0.01)
        n_mvlev = len(mv_levels)
        
        # First namelist entry
        my_file.write('num_mv_init_heights: ')
        my_file.write(str(n_mvlev) + '\n')
        
        # Second namelist entry
        my_file.write('mv_init_height: ')
        for k in xrange(n_mvlev-1):
            my_file.write(str(mv_levels[k]) + ',')
        my_file.write(str(mv_levels[-1]) + '\n')
        
        # Third namelist entry
        my_file.write('mv_init_field_type: (21) relative humidity\n')
        
        # Fourth namelist entry
        my_file.write('mv_init_data: ')
        for k in xrange(n_mvlev-1):
            my_file.write(str(mv_init[k]) + ',')
        my_file.write(str(mv_init[-1]) + '\n')

    if l_verbose: print 'Finished Moisture.'

### Updated from 'hodographs.py' ###
# Requires zi_lcl.py to be run on our simulation to create zi_DD.nc files
if l_dynamics:
    zi0_key = 'boundary layer depth'
    for day in days:
        # On the first day
        if l_verbose: print 'Starting day: ' + day
        
        # Open the netCDF for the horizontal winds
        u_nc = Dataset(data_path + 'u_' + day + '.nc', 'r')
        v_nc = Dataset(data_path + 'v_' + day + '.nc', 'r')
        
        if day == '01':
            # Read the height coordinate
            z_rho = u_nc.variables['rholev_zsea_rho'][:]*1.
        
            # Define the length of the time dimension for winds
            time_key = [key for key in u_nc.variables.keys() if 'min' in key][0]
            dt_i = u_nc.variables[u_key].shape[0]
            
            # populate our uwinds, vwinds, and time dimension with data from day 01
            zi = np.nanmean(Dataset(data_path + 'zi_' + day + '.nc', 'r').variables[zi0_key][:], axis = (1,2))
            u_mean = np.nanmean(u_nc.variables[u_key][:], axis = (2, 3))
            v_mean = np.nanmean(v_nc.variables[v_key][:], axis = (2, 3))
            times  = u_nc.variables[time_key][:]
        else:
            zi = np.concatenate((zi, np.nanmean(Dataset(data_path + 'zi_' + day + '.nc', 'r').variables[zi0_key][:], axis = (1,2))), axis = 0)
            u_mean = np.concatenate((u_mean, np.nanmean(u_nc.variables[u_key][:], axis = (2, 3))), axis = 0)
            v_mean = np.concatenate((v_mean, np.nanmean(v_nc.variables[v_key][:], axis = (2, 3))), axis = 0)
            times  = np.concatenate((times, u_nc.variables[time_key][:]), axis = 0)
            
        u_nc.close()
        v_nc.close()
        
    ### Get the mean wind speed profiles for the last 4 days ###
    # Focus alterations to the well-mixed layer

    #<------------------Geostrophic Forcing Needs to be Changed------------------->#
    u_g_z = np.array([0., 9000., 15000., 40000.])
    u_g = np.array([-10., -10., 0., 0.])
    u_g_interpolated = interpolate.interp1d(u_g_z, u_g, fill_value = 'extrapolate')(z_rho)
    v_g_interpolated = np.zeros_like(z_rho)
    #<---------------------------------------------------------------------------->#

    # Initialise the balanced wind arrays
    print len(zi)
    print u_mean.shape[0]
    target_zi  = np.mean(zi[-4*dt_i:])
    u_balanced = np.nanmean([interpolate.interp1d(z_rho/zi[timestep], u_mean[timestep+1,:], fill_value = 'extrapolate')(z_rho/target_zi) for timestep in xrange(-4*dt_i, u_mean.shape[0]-1)], axis = 0)
    v_balanced = np.nanmean([interpolate.interp1d(z_rho/zi[timestep], v_mean[timestep+1,:], fill_value = 'extrapolate')(z_rho/target_zi) for timestep in xrange(-4*dt_i, v_mean.shape[0]-1)], axis = 0)

    if l_verbose: print 'Weighting winds by boundary layer depth'
    z_balanced = z_rho*1.
    
    ### Blended transition from mean in the boundary layer to geostrophic wind above ###
    factor = z_rho/target_zi
    u_balanced_new = np.array([u_balanced[k] if factor[k] < 1.0 else (factor[k])*u_g_interpolated[k] + (1. - factor[k])*u_balanced[k] if factor[k] < 2.0 else u_g_interpolated[k] for k in xrange(len(z_rho))])
    v_balanced_new = np.array([v_balanced[k] if factor[k] < 1.0 else (factor[k])*v_g_interpolated[k] + (1. - factor[k])*v_balanced[k] if factor[k] < 2.0 else v_g_interpolated[k] for k in xrange(len(z_rho))])

    ### Minimize required points to reproduce the profiles using the RDP function ###
    # Create a well-reproduces u-profile
    z_u, u_u = RDP(z_balanced, u_balanced_new, 0.001)

    # Create a v-profile at the same heights as the u-profile
    v_u = interpolate.interp1d(z_balanced, v_balanced_new, fill_value = 'extrapolate')(z_u)

    if l_verbose: print 'Making the wind initial conditions.'
    n_uvlev = len(z_u)
    with open('../InitialFields_Wind_' + ID + '.txt', 'w') as my_new_file:
        my_new_file.write('Specify initial wind profiles.\n')
        
        # First namelist entry
        my_new_file.write('num_uv_init_heights: ')
        my_new_file.write(str(n_uvlev) + '\n')
        
        # Second namelist entry
        my_new_file.write('uv_init_height: ')
        for k in xrange(n_uvlev-1):
            my_new_file.write(str(z_u[k])+',')
        my_new_file.write(str(z_u[-1])+'\n')
        
        # Third namelist entry
        my_new_file.write('u_init_data: ')
        for k in xrange(n_uvlev-1):
            my_new_file.write(str(round(u_u[k],2))+',')
        my_new_file.write(str(round(u_u[-1],2))+'\n')
        
        # Fourth namelist entry
        my_new_file.write('v_init_data: ')
        for k in xrange(n_uvlev-1):
            my_new_file.write(str(round(v_u[k],2))+',')
        my_new_file.write(str(round(v_u[-1],2)))

### Updates from 'surface_fluxes.py' ###

if l_fluxes:
    lhf_key = u'STASH_m01s03i234'
    shf_key = u'STASH_m01s03i217'
    for day in days:
        if l_verbose: print 'STARTING: day ' + day
        fluxes_nc = Dataset(data_path + 'fluxes_' + day + '.nc', 'r')
        
        if day == '01':
            dt_i = fluxes_nc.variables[lhf_key].shape[0]
            
            LHF   = np.nanmean(fluxes_nc.variables[lhf_key][:], axis = (1, 2))
            SHF   = np.nanmean(fluxes_nc.variables[shf_key][:], axis = (1, 2))
            times = fluxes_nc.variables['hr1'][:]*1.
        else:
            LHF   = np.concatenate((LHF, np.mean(fluxes_nc.variables[lhf_key][:], axis = (1, 2))))
            SHF   = np.concatenate((SHF, np.mean(fluxes_nc.variables[shf_key][:], axis = (1, 2))))
            times = np.concatenate((times, fluxes_nc.variables['hr1'][:]*1.), axis = 0)
        fluxes_nc.close()

    SHF_mean = round(np.nanmean(SHF[-4*dt_i:]), 4)
    LHF_mean = round(np.nanmean(LHF[-4*dt_i:]), 4)
    with open('../SurfaceFluxes_balanced_' + ID + '.txt', 'w') as my_file:
        my_file.write('Specify or calculate surface fluxes of latent and sensible heat.\n')
        
        # First namelist entry
        my_file.write('idlsurffluxseaoption: (3) Constant\n')
        
        # Second namelist entry
        my_file.write('idlsurffluxseaparams: (1) ')
        my_file.write(str(SHF_mean))
        my_file.write(' (2) ')
        my_file.write(str(LHF_mean))
        my_file.write(' (3) 0 (4) 0\n')
        
        # Third namelist entry
        my_file.write('l_specz0: True\n')
        
        # Fourth namelist entry
        my_file.write('roughlen_z0m: 0.0002\n')
        
        # Fifth namelist entry
        my_file.write('roughlen_z0h: 0.00002')
        

