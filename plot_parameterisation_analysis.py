import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import mcl_key, mv_key, mci_key, q_key, w_key, theta_key, 
                       temp_key, u_key, v_key, rho_key, pthe_key, lcl_key, 
                       shf_key, lhf_key
import SkewT_archer as skewt
import sys
import os

################################################################################
############################ read in the data ##################################
################################################################################
base_path = '/nerc/n02/n02/xb899100/CloudTrail/'

# create a default output path
output_path = '../output_from_param_analysis/'
# if that directory doesn't already exist, make it
if not os.path.exists(output_path):
    os.mkdir(output_path)

sub_levels = ['Control/', 'Control_0800m/', 'Control_1600m/']
sub_level_colors = {'Control/' : 'darkred',
                    'Control_0800m/' : 'red',
                    'Control_1600m/' : 'orange'}

my_keys = [mcl_key, mv_key, mci_key, q_key, w_key, theta_key, temp_key, u_key, v_key, rho_key, pthe_key, lcl_key, shf_key, lhf_key]
dim_keys = ['min', 'zsea', 'latitude', 'longitude']

# initialise a dictionary for the data and a dictionary for the dimensions
data_dict = {}
dims_dict = {}

for sub_level in sub_levels:
    path = base_path + sub_level
    print path
    
    print 'Output path set as = ' + output_path
    
    label = path.split('/')[-2]
    print 'Experiment label = ' + label
    
    # give each experiment a subdict entry
    data_dict[sub_level] = {}
    dims_dict[sub_level] = {}
    
    nc_files = [filename for filename in os.listdir(path) if '.nc' in filename]
    
    # Don't care about the raw u and v
    nc_files = [filename for filename in nc_files if ('u_' not in filename) and ('v_' not in filename)]
    
    # Only consider the period leading up to the peak heating
    if sub_level == 'Control/':
        l_short = False
    else:
        l_short = True
        
    if l_short:
        nc_files = [filename for filename in nc_files if '_04' in filename]
    else:
        nc_files = [filename for filename in nc_files if '_09' in filename]
    
    print 'Reading in the data'
    for nc_file in nc_files:
        print 'Opening ' + nc_file
        # open the netcdf as read-only
        data_nc = Dataset(path + nc_file, 'r')
        
        for my_key in my_keys:
            # loop through all the keys we care about
            if (my_key in data_nc.variables.keys()):
                if ('swath' in nc_file) and (my_key + '_swath' not in data_dict[sub_level].keys()):
                    print 'Reading swath variable: ' + my_key
                    # treat the swath variables slightly differently
                    data_dict[sub_level][my_key + '_swath'] = data_nc.variables[my_key][:]
                elif (my_key not in data_dict[sub_level].keys()):
                    print 'Reading variable: ' + my_key
                    # if a given key is in the open netcdf and not already read in from elsewhere, read it in
                    data_dict[sub_level][my_key] = data_nc.variables[my_key][:]
                
        for dim_key in dim_keys:
            # loop through all the parts of the dimension keys, if there are keys in the netcdf which contain those parts, include them in the list 'matches'
            matches = [dimension_key for dimension_key in data_nc.variables.keys() if (dim_key in dimension_key) and ('rotated' not in dimension_key)]
            if (len(matches) >= 1):
                # if there are matches
                for matching_key in matches:
                    if matching_key not in dims_dict[sub_level].keys():
                        print 'Reading dimension: ' + matching_key
                        # and the matches aren't already in the dimensions dictionary, add them
                        dims_dict[sub_level][matching_key] = data_nc.variables[matching_key][:]
                        if ('min' in matching_key) and l_short:
                            # Add the offset for short nights so that the peak heating is at 'noon'
                            dims_dict[sub_level][matching_key] += 240.0
        print 'Closing ' + nc_file + '\n'
        data_nc.close() # close the netcdf

    missing_keys = [my_key for my_key in my_keys if my_key not in data_dict[sub_level].keys()]

    print 'Unable to locate data for the following keys:'
    print missing_keys

    print 'Creating the faux bl fluxes height coordinate...'
    dims_dict[sub_level]['rholev_zsea_rho_fluxes'] = dims_dict[sub_level]['rholev_zsea_rho'][:]*1.
    dims_dict[sub_level]['rholev_zsea_rho_fluxes'][0] = 0.0

    print 'Finished reading in the data.\n'

################################################################################
######################## do mass flux computations #############################
################################################################################

print 'Starting computation routine: MASS FLUX'

# timeseries of LCL mass flux
fig1 = plt.figure(tight_layout = True, figsize = (9, 7))
axa = fig1.add_subplot(1, 1, 1)
axa.set_ylabel(u'mass flux (m s$^{-1}$ kg m$^{-3}$)')
axa.set_xlabel('time (hours)')
axa.set_title('Time-series of mean LCL mass flux in swath downwind of island')
axa.set_ylim([0.0, 0.5])
axa.set_xlim([9, 12])

for sub_level in sub_levels:
    # target the lcl
    target_height = np.nanmean(data_dict[sub_level][lcl_key], axis = (1, 2))
    test = [np.abs(dims_dict[sub_level]['thlev_zsea_theta'] - height) for height in target_height]
    iz = [np.where(test1 == np.min(test1))[0][0] for test1 in test]
    
    # compute the time evolution of mass flux
    print '---> Calculating the time evolution of mass flux through this period...'
    mass_flux = data_dict[sub_level][w_key+'_swath']*data_dict[sub_level][rho_key+'_swath']
    time_key = [tkey for tkey in dims_dict[sub_level].keys() if dims_dict[sub_level][tkey].shape[0] == mass_flux.shape[0]][0]
    
    # create a time series of cloud base mass flux in the cloud trail swath
    mean_cb_mass_flux_ts = np.array([np.nanmean(mass_flux[it,iz[it],:,:]) for it in range(mass_flux.shape[0])])
    
    # create a time series of cloud base mass flux where there are clouds in the cloud trail swath
    mean_cb_cond_mass_flux_ts = np.array([np.nanmean(np.where(np.nanmax(data_dict[sub_level][mcl_key+'_swath'][it,:,:,:], axis = 0) > 0., mass_flux[it,iz[it],:,:], np.nan)) for it in range(mass_flux.shape[0])])
    
    ############################################################################
    ####################### plot the mass flux data ############################
    ############################################################################
    
    print 'Starting plotting routine: MASS FLUX'
    # plot the time evolution of the mass flux
    print '---> Plotting the time evolution of mass fluxes...'
    axa.plot(dims_dict[sub_level][time_key]/60., mean_cb_mass_flux_ts, lw = 2, color = sub_level_colors[sub_level], label = sub_level[:-1])
    axa.plot(dims_dict[sub_level][time_key]/60., mean_cb_cond_mass_flux_ts, lw = 2, ls = ':', color = sub_level_colors[sub_level], label = sub_level[:-1] + ' (cloudy)')

axa.legend(loc = 0, ncol = 2)
plt.savefig(output_path + label + '_swath_mean_massflux_at_lcl_timeseries.png', dpi = 150)
plt.show()

# profiles of mass flux + cloudy mass flux
fig2 = plt.figure(tight_layout = True, figsize = (9, 13.5))
axb = fig2.add_subplot(1, 1, 1)
axb.set_xlabel(u'mass flux (m s$^{-1}$ kg m$^{-3}$)')
axb.set_ylabel('height (km)')
axb.set_title('Mass flux profiles averaged from t=660 to 720 mins over swath downwind of island')
axb.set_xlim([-0.01, 2.5])
axb.set_ylim([0, 5])

for sub_level in sub_levels:
    # compute the mass flux
    print '---> Calculating the mass flux through this period...'
    mass_flux = data_dict[sub_level][w_key+'_swath']*data_dict[sub_level][rho_key+'_swath']
    time_key = [tkey for tkey in dims_dict[sub_level].keys() if dims_dict[sub_level][tkey].shape[0] == mass_flux.shape[0]][0]
    
    # Find times between t=660 and t=720 mins
    its = [IT for IT in range(dims_dict[sub_level][time_key].shape[0]) if 660. <= dims_dict[sub_level][time_key][IT] <= 720.]
    
    # create a time series of cloud base mass flux in the cloud trail swath
    mean_cb_mass_flux_prof = np.nanmean(mass_flux[its,:,:,:], axis = (0, 2, 3))
    
    # create a time series of cloud base mass flux where there are clouds in the cloud trail swath
    mean_cb_cond_mass_flux_prof = np.nanmean(np.where(data_dict[sub_level][mcl_key+'_swath'][its,:,:,:] > 0., mass_flux[its,:,:,:], np.nan), axis = (0, 2, 3))
    
    ############################################################################
    ####################### plot the mass flux data ############################
    ############################################################################
    
    print 'Starting plotting routine: MASS FLUX'
    # plot the time evolution of the mass flux
    print '---> Plotting the time evolution of mass fluxes...'
    axb.plot(mean_cb_mass_flux_prof, dims_dict[sub_level]['thlev_zsea_theta']/1000., lw = 2, color = sub_level_colors[sub_level], label = sub_level[:-1])
    axb.plot(mean_cb_cond_mass_flux_prof, dims_dict[sub_level]['thlev_zsea_theta']/1000., lw = 2, ls = ':', color = sub_level_colors[sub_level], label = sub_level[:-1] + ' (cloudy)')

axb.legend(loc = 0)
plt.savefig(output_path + label + '_swath_mean_massflux_at_lcl_profile.png', dpi = 150)
plt.show()

################################################################################
##################### compute the mean thermodynamic state #####################
################################################################################

print '---> Starting computation routine: SKEWT + CLD'

# make sure output directories exist
output_mcl_path = output_path + 'mcl/'
output_skewt_path = output_path + 'skewt/'

# if that directory doesn't already exist, make it
if not os.path.exists(output_mcl_path):
    os.mkdir(output_mcl_path)
if not os.path.exists(output_skewt_path):
    os.mkdir(output_skewt_path)

# do earlier
fig = plt.figure(tight_layout = True, figsize = (9, 9))
axa = fig.add_subplot(1, 1, 1)
axa.set_xlabel(u'm$_{cl}$ or m$_{cf}$ (g kg$^{-1}$)')
axa.set_ylabel(u' height (km)')
axa.set_xlim([0, 0.05])
axa.set_title('Domain-mean cloud water (m$_{cl}$ + m$_{cf}$) \nT = 9 to 10 hrs')

for sub_level in sub_levels:
    time_key = [tkey for tkey in dims_dict[sub_level].keys() if dims_dict[sub_level][tkey].shape[0] == data_dict[sub_level][mcl_key].shape[0]][0]
    # Choose 'its' to be 9-10am
    test09 = np.abs(dims_dict[sub_level][time_key] - 9*60)
    test10 = np.abs(dims_dict[sub_level][time_key] - 10*60)
    its09 = range(np.where(test09 == np.min(test09))[0][0], np.where(test10 == np.min(test10))[0] + 1)
    # compute the domain mean liquid water profiles
    print '---> Calculating the domain mean liquid water path at the end of this period...'
    mcl_domain_mean = np.nanmean(data_dict[sub_level][mcl_key][its09,:,:,:], axis = (0, 2, 3))
    print '---> Plotting the domain mean liquid water profile at the end of this period...'
    axa.plot(mcl_domain_mean*1000., dims_dict[sub_level]['thlev_zsea_theta']/1000., lw = 2, color = sub_level_colors[sub_level], label = sub_level[:-1] + ' m$_{cl}$')

for sub_level in sub_levels:
    time_key = [tkey for tkey in dims_dict[sub_level].keys() if dims_dict[sub_level][tkey].shape[0] == data_dict[sub_level][mcl_key].shape[0]][0]
    # Choose 'its' to be 9-10am
    test09 = np.abs(dims_dict[sub_level][time_key] - 9*60)
    test10 = np.abs(dims_dict[sub_level][time_key] - 10*60)
    its09 = range(np.where(test09 == np.min(test09))[0][0], np.where(test10 == np.min(test10))[0] + 1)
    # compute the domain mean liquid water profiles
    print '---> Calculating the domain mean liquid water path at the end of this period...'
    mcf_domain_mean = np.nanmean(data_dict[sub_level][mci_key][its09,:,:,:], axis = (0, 2, 3))
    print '---> Plotting the domain mean liquid water profile at the end of this period...'
    axa.plot(mcf_domain_mean*1000., dims_dict[sub_level]['thlev_zsea_theta']/1000., lw = 2, color = sub_level_colors[sub_level], ls = ':', label = sub_level[:-1] + ' m$_{cf}$')

axa.legend(loc = 0, ncol = 2)
plt.savefig(output_mcl_path + '_domain_mean_mcl_' + "{0:04d}".format(int(dims_dict[sub_level][time_key][its09[0]])) + 'mins.png', dpi = 150)
plt.show()

# do middle
fig = plt.figure(tight_layout = True, figsize = (9, 9))
axa = fig.add_subplot(1, 1, 1)
axa.set_xlabel(u'm$_{cl}$ or m$_{cf}$ (g kg$^{-1}$)')
axa.set_ylabel(u' height (km)')
axa.set_xlim([0, 0.05])
axa.set_title('Domain-mean cloud water (m$_{cl}$ + m$_{cf}$) \nT = 9 to 10 hrs')

for sub_level in sub_levels:
    time_key = [tkey for tkey in dims_dict[sub_level].keys() if dims_dict[sub_level][tkey].shape[0] == data_dict[sub_level][mcl_key].shape[0]][0]
    # Choose 'its' to be 10-11am
    test10 = np.abs(dims_dict[sub_level][time_key] - 10*60)
    test11 = np.abs(dims_dict[sub_level][time_key] - 11*60)
    its10 = range(np.where(test10 == np.min(test10))[0][0], np.where(test11 == np.min(test11))[0] + 1)
    # compute the domain mean liquid water profiles
    print '---> Calculating the domain mean liquid water path at the end of this period...'
    mcl_domain_mean = np.nanmean(data_dict[sub_level][mcl_key][its10,:,:,:], axis = (0, 2, 3))
    print '---> Plotting the domain mean liquid water profile at the end of this period...'
    axa.plot(mcl_domain_mean*1000., dims_dict[sub_level]['thlev_zsea_theta']/1000., lw = 2, color = sub_level_colors[sub_level], label = sub_level[:-1])

for sub_level in sub_levels:
    time_key = [tkey for tkey in dims_dict[sub_level].keys() if dims_dict[sub_level][tkey].shape[0] == data_dict[sub_level][mcl_key].shape[0]][0]
    # Choose 'its' to be 9-10am
    test10 = np.abs(dims_dict[sub_level][time_key] - 10*60)
    test11 = np.abs(dims_dict[sub_level][time_key] - 11*60)
    its10 = range(np.where(test10 == np.min(test10))[0][0], np.where(test11 == np.min(test11))[0] + 1)
    # compute the domain mean liquid water profiles
    print '---> Calculating the domain mean liquid water path at the end of this period...'
    mcf_domain_mean = np.nanmean(data_dict[sub_level][mci_key][its10,:,:,:], axis = (0, 2, 3))
    print '---> Plotting the domain mean liquid water profile at the end of this period...'
    axa.plot(mcf_domain_mean*1000., dims_dict[sub_level]['thlev_zsea_theta']/1000., lw = 2, color = sub_level_colors[sub_level], ls = ':', label = sub_level[:-1] + ' m$_{cf}$')

axa.legend(loc = 0, ncol = 2)
plt.savefig(output_mcl_path + '_domain_mean_mcl_' + "{0:04d}".format(int(dims_dict[sub_level][time_key][its10[0]])) + 'mins.png', dpi = 150)
plt.show()

# do later
fig = plt.figure(tight_layout = True, figsize = (9, 9))
axa = fig.add_subplot(1, 1, 1)
axa.set_xlabel(u'm$_{cl}$ or m$_{cf}$ (g kg$^{-1}$)')
axa.set_ylabel(u' height (km)')
axa.set_xlim([0, 0.05])
axa.set_title('Domain-mean cloud water (m$_{cl}$ + m$_{cf}$) \nT = 9 to 10 hrs')

for sub_level in sub_levels:
    time_key = [tkey for tkey in dims_dict[sub_level].keys() if dims_dict[sub_level][tkey].shape[0] == data_dict[sub_level][mcl_key].shape[0]][0]
    # Choose 'its' to be 11-12pm
    test11 = np.abs(dims_dict[sub_level][time_key] - 11*60)
    test12 = np.abs(dims_dict[sub_level][time_key] - 12*60)
    its12 = range(np.where(test11 == np.min(test11))[0][0], np.where(test12 == np.min(test12))[0][0] + 1)
    # compute the domain mean liquid water profiles
    print '---> Calculating the domain mean liquid water path at the end of this period...'
    mcl_domain_mean = np.nanmean(data_dict[sub_level][mcl_key][its12,:,:,:], axis = (0, 2, 3))
    print '---> Plotting the domain mean liquid water profile at the end of this period...'
    axa.plot(mcl_domain_mean*1000., dims_dict[sub_level]['thlev_zsea_theta']/1000., lw = 2, color = sub_level_colors[sub_level], label = sub_level[:-1])

for sub_level in sub_levels:
    time_key = [tkey for tkey in dims_dict[sub_level].keys() if dims_dict[sub_level][tkey].shape[0] == data_dict[sub_level][mcl_key].shape[0]][0]
    # Choose 'its' to be 9-10am
    test11 = np.abs(dims_dict[sub_level][time_key] - 11*60)
    test12 = np.abs(dims_dict[sub_level][time_key] - 12*60)
    its12 = range(np.where(test11 == np.min(test11))[0][0], np.where(test12 == np.min(test12))[0] + 1)
    # compute the domain mean liquid water profiles
    print '---> Calculating the domain mean liquid water path at the end of this period...'
    mcf_domain_mean = np.nanmean(data_dict[sub_level][mci_key][its12,:,:,:], axis = (0, 2, 3))
    print '---> Plotting the domain mean liquid water profile at the end of this period...'
    axa.plot(mcf_domain_mean*1000., dims_dict[sub_level]['thlev_zsea_theta']/1000., lw = 2, color = sub_level_colors[sub_level], ls = ':', label = sub_level[:-1] + ' m$_{cf}$')

axa.legend(loc = 0, ncol = 2)
plt.savefig(output_mcl_path + '_domain_mean_mcl_' + "{0:04d}".format(int(dims_dict[sub_level][time_key][its12[-1]])) + 'mins.png', dpi = 150)
plt.show()

for sub_level in sub_levels:
    for it in range(dims_dict[sub_level][time_key].shape[0]):
        # plot the mean thermodynamic state
        print '---> Plotting the domain mean thermodynamic state...'
        skewT_dict = {}
        skewT_dict['temperature'] = np.nanmean(data_dict[sub_level][temp_key][it,:,:,:], axis = (1, 2))
        skewT_dict['pressure']    = np.nanmean(data_dict[sub_level][pthe_key][it,:,:,:], axis = (1, 2))
        skewT_dict['q']           = np.nanmean(data_dict[sub_level][q_key][it,:,:,:], axis = (1, 2))
        skewT_dict['dewpoint']    = skewt.getDew(QIN = skewT_dict['q'][:]*1., PIN1 = skewT_dict['pressure'][:]*1., q_units = 'kg/kg', p_units = 'Pa')
        skewT_dict['u']           = np.nanmean(data_dict[sub_level][u_key][it,:,:,:], axis = (1, 2))
        skewT_dict['v']           = np.nanmean(data_dict[sub_level][v_key][it,:,:,:], axis = (1,2))
        
        skewt.plotSkewT(temp = skewT_dict['temperature'][:]-273.15, t_dew = skewT_dict['dewpoint'][:]*1., p = skewT_dict['pressure'][:]/100., u = skewT_dict['u'][:]/0.5144, v = skewT_dict['v'][:]/0.5144, CAPE = False, date = -999, my_title = label + ' domain-mean skewT-logP \nT = ' + "{0:04d}".format(int(dims_dict[sub_level][time_key][it])) + ' min', temp_col = ['r'], dew_col = ['b'])
        plt.savefig(output_skewt_path + sub_level[:-1] + '_domain_mean_skewt_' + "{0:04d}".format(int(dims_dict[sub_level][time_key][it])) + 'mins.png', dpi = 150)
        plt.close('all')

################################################################################
############## compute the turbulent fluxes (heat and moisture) ################
################################################################################

print '---> Calculating the turbulent heat and moisture fluxes...'

for hour in [9, 10, 11]:
    fig = plt.figure(tight_layout = True, figsize = (12, 7.5))
    axa = fig.add_subplot(1, 2, 1)
    axa.set_xlabel(u'($\overline{w^{\prime}\\theta^{\prime}_{l}}$) (W m$^{-2}$)')
    axa.set_ylabel('Height (km)')
    axa.set_ylim([0, 9])
    
    axb = fig.add_subplot(1, 2, 2, sharey = axa)
    axb.set_xlabel(u'($\overline{w^{\prime}q^{\prime}_{T}}$) (W m$^{-2}$)')
    axb.set_ylim([0, 9])
    for sub_level in sub_levels:
        print '------> Computing thermodynamic variables...'
        thetal = data_dict[sub_level][theta_key] - (data_dict[sub_level][theta_key]/data_dict[sub_level][temp_key])*(skewt.Lv/skewt.cpd)*data_dict[sub_level][mcl_key]
        # use itsHH from previous routine
        
        time_key = [tkey for tkey in dims_dict[sub_level].keys() if dims_dict[sub_level][tkey].shape[0] == data_dict[sub_level][w_key].shape[0]][0]
        # Choose 'its' to be 11-12pm
        test1 = np.abs(dims_dict[sub_level][time_key] - hour*60)
        test2 = np.abs(dims_dict[sub_level][time_key] - (hour+1)*60)
        its = range(np.where(test1 == np.min(test1))[0][0], np.where(test2 == np.min(test2))[0][0] + 1)
        print '------> Calculating perturbations...'
        wp      = np.array([np.transpose(np.transpose(data_dict[sub_level][w_key][it,:,:,:]) - np.nanmean(data_dict[sub_level][w_key][it,:,:,:], axis = (1, 2))) for it in its])
        thetalp = np.array([np.transpose(np.transpose(thetal[it,:,:,:]) - np.nanmean(thetal[it,:,:,:], axis = (1, 2))) for it in its])
        qtp     = np.array([np.transpose(np.transpose(data_dict[sub_level][q_key][it,:,:,:]) - np.nanmean(data_dict[sub_level][q_key][it,:,:,:], axis = (1, 2))) for it in its])
        
        print '------> Interpolating stuff from boundary layer scheme from rho to theta levels...'
        from scipy import interpolate
        shf_bl = np.array([interpolate.interp1d(y = data_dict[sub_level][shf_key][it,:,:,:], x = dims_dict[sub_level]['rholev_zsea_rho_fluxes'], fill_value = 'extrapolate', axis = 0)(dims_dict[sub_level]['thlev_zsea_theta']) for it in its])
        lhf_bl = np.array([interpolate.interp1d(y = data_dict[sub_level][lhf_key][it,:,:,:], x = dims_dict[sub_level]['rholev_zsea_rho_fluxes'], fill_value = 'extrapolate', axis = 0)(dims_dict[sub_level]['thlev_zsea_theta']) for it in its])
        rho_bl = np.array([interpolate.interp1d(y = data_dict[sub_level][rho_key][it,:,:,:], x = dims_dict[sub_level]['rholev_zsea_rho'], fill_value = 'extrapolate', axis = 0)(dims_dict[sub_level]['thlev_zsea_theta']) for it in its])
        
        print '------> Combining resolved and subgrid scale fluxes...'
        wpthetalp = np.nanmean(skewt.cpd*rho_bl*wp*thetalp + shf_bl, axis = (0, 2, 3))
        wpqtp     = np.nanmean(skewt.Lv*(rho_bl*wp*qtp + lhf_bl), axis = (0, 2, 3))
        
        print 'Finished computation routine.\n'
        
        print '---> Plotting the vertical profiles of turbulent heat and moisture fluxes...'
        
        axa.plot(wpthetalp, dims_dict[sub_level]['thlev_zsea_theta']/1000., lw = 2, color = sub_level_colors[sub_level], label = sub_level[:-1])
        axb.plot(wpqtp, dims_dict[sub_level]['thlev_zsea_theta']/1000., lw = 2, color = sub_level_colors[sub_level], label = sub_level[:-1])
    
    axa.set_title(u'a) time-mean (w$^{\prime} \\theta_{l}^{\prime}$) ' + "{0:02d}".format(hour) + ' to ' + "{0:02d}".format(hour+1) + ' hrs')
    axb.set_title(u'b) time-mean (w$^{\prime} q_{T}^{\prime}$) ' + "{0:02d}".format(hour) + ' to ' + "{0:02d}".format(hour+1) + ' hrs')
    axa.legend(loc = 0)
    plt.savefig(output_path + 'horizontal_mean_flux_profiles_' + "{0:04d}".format(int(dims_dict[sub_level][time_key][its[0]])) + 'mins.png', dpi = 150)
    plt.show()

print 'Finished plotting routine.\n'

