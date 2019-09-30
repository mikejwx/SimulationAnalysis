import numpy as np
import matplotlib.pyplot as plt
from SkewT_archer import *
from analysis_tools import toComponents, round5
"""
Plots the radiosonde from a case CT day.
Plots the simplified initial conditions based on that radiosonde in a comparison
    three-panel plot.
Plots the spin-up-based initial conditions in a comparison with the simplified
    initial conditions in a three-panel plot.
"""

def three_panel(theta_1, theta_2, RH_1, RH_2, U_1, U_2, height, label_1, label_2, figname):
    """
    Style for the three-panel plot to compare radiosonde and simplified initial
    profiles. Or the simplified initial profiles to the spin-up.
    
    These are done in height coordinates for the lowest few km to highlight
    low-level differences.
    theta_N = the potential temperature in Kelvin
    RH_N    = the relative humidity as a percentage
    U_N     = the wind speed in m/s
    height  = the height of the observations in km
    """
    fig = plt.figure(figsize = (18, 6))
    axa = fig.add_subplot(1, 3, 1)
    axa.set_ylim([0, 5])
    axb = fig.add_subplot(1, 3, 2)
    axb.set_ylim([0, 5])
    axc = fig.add_subplot(1, 3, 3)
    axc.set_ylim([0, 5])
    
    # find the index of the maximum height
    iz_max = np.where(np.abs(5 - height) == np.min(np.abs(5 - height)))[0][0]+2
    axa.plot(theta_1, height, 'k')
    axa.plot(theta_2, height, 'k', lw = 2)
    axa.set_ylabel('Height (km)')
    axa.set_xlabel(u'$\\theta$ (K)')
    x_min = np.min([theta_1[:iz_max].min(), theta_2[:iz_max].min()])
    x_max = np.max([theta_1[:iz_max].max(), theta_2[:iz_max].max()])
    axa.set_xlim([round5(x_min) - [5 if round5(x_min) > x_min else 0][0], round5(x_max) + [5 if round5(x_max) < x_max else 0][0]])
    
    axb.plot(RH_1, height, 'k', label = label_1)
    axb.plot(RH_2, height, 'k', lw = 2, label = label_2)
    axb.set_yticklabels([''])
    axb.set_xlabel('RH (%)')
    axb.legend(loc = 0)
    x_min = np.min([RH_1[:iz_max].min(), RH_2[:iz_max].min()])
    x_max = np.max([RH_1[:iz_max].max(), RH_2[:iz_max].max()])
    axb.set_xlim([round5(x_min) - [5 if round5(x_min) > x_min else 0][0], round5(x_max) + [5 if round5(x_max) < x_max else 0][0]])
    
    axc.plot(U_1, height, 'k')
    axc.plot(U_2, height, 'k', lw = 2)
    axc.set_yticklabels([''])
    axc.set_xlabel(u'Wind Speed (m s$^{-1}$)')
    x_min = np.min([U_1[:iz_max].min(), U_2[:iz_max].min()])
    x_max = np.max([U_1[:iz_max].max(), U_2[:iz_max].max()])
    axc.set_xlim([round5(x_min) - [5 if round5(x_min) > x_min else 0][0], round5(x_max) + [5 if round5(x_max) < x_max else 0][0]])
    
    plt.savefig('../' + figname, dpi = 150, bbox_inches = 'tight')
    plt.show()

### Read the radiosonde data ###
radiosonde_data = {}
with open('../IC_radiosonde.txt', 'r') as radiosonde_txt:
    radiosonde_read = radiosonde_txt.readlines()
    for line in radiosonde_read:
        line_split = line.split()
        # 1st column is the pressure in hPa
        if 'pressure' in radiosonde_data.keys():
            radiosonde_data['pressure'].append(float(line_split[0]))
        else:
            radiosonde_data['pressure'] = [float(line_split[0])]
        
        # 2nd column is height in metres
        if 'height' in radiosonde_data.keys():
            radiosonde_data['height'].append(float(line_split[1]))
        else:
            radiosonde_data['height'] = [float(line_split[1])]
        # 3rd column is temperature in C
        if 'temperature' in radiosonde_data.keys():
            radiosonde_data['temperature'].append(float(line_split[2]))
        else:
            radiosonde_data['temperature'] = [float(line_split[2])]
        # 4th column in dew point in C
        if 'dewpoint' in radiosonde_data.keys():
            radiosonde_data['dewpoint'].append(float(line_split[3]))
        else:
            radiosonde_data['dewpoint'] = [float(line_split[3])]
        # 5th column is RH in %
        if 'RH' in radiosonde_data.keys():
            radiosonde_data['RH'].append(float(line_split[4]))
        else:
            radiosonde_data['RH'] = [float(line_split[4])]
        # 6th column is mixing ratio in g/kg
        # 7th column is wind direction in degrees
        if 'wind_dir' in radiosonde_data.keys():
            radiosonde_data['wind_dir'].append(float(line_split[6]))
        else:
            radiosonde_data['wind_dir'] = [float(line_split[6])]
        # 8th column is wind speed in kts
        if 'wind_spd' in radiosonde_data.keys():
            radiosonde_data['wind_spd'].append(float(line_split[7]))
        else:
            radiosonde_data['wind_spd'] = [float(line_split[7])]
        # 9th column is potential temperature in Kelvin
        if 'theta' in radiosonde_data.keys():
            radiosonde_data['theta'].append(float(line_split[8]))
        else:
            radiosonde_data['theta'] = [float(line_split[8])]
        #10th column is equivalent potential temperature in Kelvin
        #11th column is virtual potential temperature in Kelvin


# Convert everything to arrays
for key in radiosonde_data.keys():
    radiosonde_data[key] = np.array(radiosonde_data[key])

radiosonde_data['u'], radiosonde_data['v'] = toComponents(radiosonde_data['wind_spd'], radiosonde_data['wind_dir'])

# Plot the radiosonde
plotSkewT(radiosonde_data['temperature'], radiosonde_data['dewpoint'], radiosonde_data['pressure'], u = radiosonde_data['u'], v = radiosonde_data['v'], date = '16th July 2015', CAPE = True)
plt.savefig('../ic_radiosonde_skewT.png', dpi = 150, bbox_inches = 'tight')
plt.show()

# Plot the three panel for radiosonde and simplified
# Interpolate onto the same height levels
my_heights = np.linspace(0.0, 40000.0, 1000)
from scipy import interpolate
simple_theta = np.array([299., 299., 299.01, 302.47, 311.38, 318.05, 323.86, 328.95, 335.98, 339.1, 342.42, 354.57, 402.01, 1200.0])
simple_theta_z = np.array([0, 148, 509, 857, 3217, 3840, 5931, 7635, 9714, 10963, 12408, 14189, 16649, 40000])
simple_theta = interpolate.interp1d(x = simple_theta_z, y = simple_theta, fill_value = 'extrapolate')(my_heights)
radiosonde_data['theta_int'] = interpolate.interp1d(x = radiosonde_data['height'], y = radiosonde_data['theta'], fill_value = 'extrapolate')(my_heights)

simple_RH = np.array([77.5, 81.4, 89.0, 55.4, 38.6, 20., 20.])
simple_RH_z = np.array([0, 148, 509, 857, 3217, 3509, 40000])
simple_RH = interpolate.interp1d(x = simple_RH_z, y = simple_RH, fill_value = 'extrapolate')(my_heights)
radiosonde_data['RH_int'] = interpolate.interp1d(x = radiosonde_data['height'], y = radiosonde_data['RH'], fill_value = 'extrapolate')(my_heights)

simple_U = np.array([10, 10, 0, 0])
simple_U_z = np.array([0, 9000, 15000, 40000])
simple_U = interpolate.interp1d(x = simple_U_z, y = simple_U, fill_value = 'extrapolate')(my_heights)
radiosonde_data['wind_spd_int'] = interpolate.interp1d(x = radiosonde_data['height'], y = radiosonde_data['wind_spd'], fill_value = 'extrapolate')(my_heights)

three_panel(theta_1 = radiosonde_data['theta_int'], theta_2 = simple_theta, RH_1 = radiosonde_data['RH_int'], RH_2 = simple_RH, U_1 = radiosonde_data['wind_spd_int']*0.5144, U_2 = simple_U, height = my_heights/1000.0, label_1 = 'Radiosonde', label_2 = 'Simplified', figname = 'radiosonde_simplified_comp.png')

# Plot the three panel for simplified and spin-up
### Read the spinup data ###
spinup_data = {}
with open('../InitialFields_Temperature_Control.txt', 'r') as init_txt:
    init_read = init_txt.readlines()
    for line in init_read:
        if 'theta_init_height: ' in line:
            spinup_data['theta_init_height'] = np.array([float(point) for point in line.strip('theta_init_height: ').split(',')])
        elif 'theta_init_data: ' in line:
            spinup_data['theta_init_data'] = np.array([float(point) for point in line.strip('theta_init_data: ').split(',')])

with open('../InitialFields_Moisture_Control.txt', 'r') as init_txt:
    init_read = init_txt.readlines()
    for line in init_read:
        if 'mv_init_height: ' in line:
            spinup_data['mv_init_height'] = np.array([float(point) for point in line.strip('mv_init_height: ').split(',')])
        elif 'mv_init_data: ' in line:
            spinup_data['mv_init_data'] = np.array([float(point) for point in line.strip('mv_init_data: ').split(',')])

with open('../InitialFields_Wind_Control.txt', 'r') as init_txt:
    init_read = init_txt.readlines()
    for line in init_read:
        if 'uv_init_height: ' in line:
            spinup_data['uv_init_height'] = np.array([float(point) for point in line.strip('uv_init_height: ').split(',')])
        elif 'u_init_data: ' in line:
            spinup_data['u_init_data'] = np.array([float(point) for point in line.strip('u_init_data: ').split(',')])
        elif 'v_init_data: ' in line:
            spinup_data['v_init_data'] = np.array([float(point) for point in line.strip('v_init_data: ').split(',')])

# Interpolate everything onto the same height levels for the three-panel
simple_theta = np.array([299., 299., 299.01, 302.47, 311.38, 318.05, 323.86, 328.95, 335.98, 339.1, 342.42, 354.57, 402.01, 1200.0])
simple_theta_z = np.array([0, 148, 509, 857, 3217, 3840, 5931, 7635, 9714, 10963, 12408, 14189, 16649, 40000])
simple_theta = interpolate.interp1d(x = simple_theta_z, y = simple_theta, fill_value = 'extrapolate')(my_heights)
spinup_data['theta'] = interpolate.interp1d(x = spinup_data['theta_init_height'], y = spinup_data['theta_init_data'], fill_value = 'extrapolate')(my_heights)

simple_RH = np.array([77.5, 81.4, 89.0, 55.4, 38.6, 20., 20.])
simple_RH_z = np.array([0, 148, 509, 857, 3217, 3509, 40000])
simple_RH = interpolate.interp1d(x = simple_RH_z, y = simple_RH, fill_value = 'extrapolate')(my_heights)
spinup_data['RH'] = 100*interpolate.interp1d(x = spinup_data['mv_init_height'], y = spinup_data['mv_init_data'], fill_value = 'extrapolate')(my_heights)

simple_U = np.array([10, 10, 0, 0])
simple_U_z = np.array([0, 9000, 15000, 40000])
simple_U = interpolate.interp1d(x = simple_U_z, y = simple_U, fill_value = 'extrapolate')(my_heights)
spinup_data['U'] = np.sqrt(interpolate.interp1d(x = spinup_data['uv_init_height'], y = spinup_data['u_init_data'], fill_value = 'extrapolate')(my_heights)**2 + interpolate.interp1d(x = spinup_data['uv_init_height'], y = spinup_data['v_init_data'], fill_value = 'extrapolate')(my_heights)**2)

three_panel(theta_1 = simple_theta, theta_2 = spinup_data['theta'], RH_1 = simple_RH, RH_2 = spinup_data['RH'], U_1 = simple_U, U_2 = spinup_data['U'], height = my_heights/1000.0, label_1 = 'Simplified', label_2 = 'Spin-up', figname = 'simplified_spinup_comp.png')



