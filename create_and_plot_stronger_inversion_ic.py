"""
Create and plot the new initial conditions with a stronger inversion at 3.5 km.
"""

import numpy as np
import matplotlib.pyplot as plt
from SkewT_archer import *
from analysis_tools import toComponents, round5, RDP, send_email
from scipy import interpolate

# Start with the high resolutino spin-up data and then follow a warmer moist 
# adaibat from the height of the base of the original inversion.

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

# Find the lowest height at whidh RH < 0.2
idx_20 = np.where(spinup_data['mv_init_data'] <= 0.2)[0]

# Then reset all RH above that level to 0.2
spinup_data['mv_init_data'][idx_20[0]:] = 0.2

# Then remove all points in between lowest 0.2 and highest 0.2
spinup_data['mv_init_height'] = np.array([spinup_data['mv_init_height'][i] for i in range(len(spinup_data['mv_init_height'])) if i not in idx_20[1:-1]])
spinup_data['mv_init_data'] = np.array([spinup_data['mv_init_data'][i] for i in range(len(spinup_data['mv_init_data'])) if i not in idx_20[1:-1]])

with open('../InitialFields_Wind_Control.txt', 'r') as init_txt:
    init_read = init_txt.readlines()
    for line in init_read:
        if 'uv_init_height: ' in line:
            spinup_data['uv_init_height'] = np.array([float(point) for point in line.strip('uv_init_height: ').split(',')])
        elif 'u_init_data: ' in line:
            spinup_data['u_init_data'] = np.array([float(point) for point in line.strip('u_init_data: ').split(',')])
        elif 'v_init_data: ' in line:
            spinup_data['v_init_data'] = np.array([float(point) for point in line.strip('v_init_data: ').split(',')])

my_heights = np.linspace(0.0, 40000.0, 1000)
spinup_data['height'] = my_heights*1.
spinup_data['theta'] = interpolate.interp1d(x = spinup_data['theta_init_height'], y = spinup_data['theta_init_data'], fill_value = 'extrapolate')(my_heights)
spinup_data['RH'] = 100*interpolate.interp1d(x = spinup_data['mv_init_height'], y = spinup_data['mv_init_data'], fill_value = 'extrapolate')(my_heights)
spinup_data['U'] = np.sqrt(interpolate.interp1d(x = spinup_data['uv_init_height'], y = spinup_data['u_init_data'], fill_value = 'extrapolate')(my_heights)**2 + interpolate.interp1d(x = spinup_data['uv_init_height'], y = spinup_data['v_init_data'], fill_value = 'extrapolate')(my_heights)**2)

# skew t of the spinup data
spinup_data['pressure'] = interpolate.interp1d(radiosonde_data['height'], radiosonde_data['pressure'], fill_value = 'extrapolate')(my_heights) #assume the pressure profile hasn't changed from the radiosonde
spinup_data['temperature'] = PTtoTemp(spinup_data['theta'], spinup_data['pressure'], t_units = 'K', p_units = 'hPa')
spinup_data['q'] = getQ(spinup_data['temperature'], spinup_data['RH'], spinup_data['pressure'], t_units = 'K', p_units = 'hPa')
spinup_data['dewpoint'] = getDew(spinup_data['q'], spinup_data['pressure'], q_units = 'kg/kg', p_units = 'hPa')
spinup_data['u'] = interpolate.interp1d(x = spinup_data['uv_init_height'], y = spinup_data['u_init_data'], fill_value = 'extrapolate')(my_heights)
spinup_data['v'] = interpolate.interp1d(x = spinup_data['uv_init_height'], y = spinup_data['v_init_data'], fill_value = 'extrapolate')(my_heights)

##### Plot the origina spin-up data #####
print 'spin-up skew t'
plotSkewT(spinup_data['temperature'][:415]-273.15, spinup_data['dewpoint'][:415]-273.15, spinup_data['pressure'][:415], u = spinup_data['u'][:415], v = spinup_data['v'][:415], CAPE = True)
plt.savefig('../original_spinup_skewT.png', dpi = 150)
plt.close('all')

##### Make the inversion stronger #####
my_CAPE, my_CIN, my_ParcelT, my_ParcelP, LCLp, LFCp = getCAPE(spinup_data['temperature'][:415]+7.5, spinup_data['q'][:415], spinup_data['pressure'][:415], parcel_type = 1, t_units = 'K', q_units = 'kg/kg', p_units = 'hPa')

spinup_data['temperature_int'] = interpolate.interp1d(spinup_data['pressure'], spinup_data['temperature'], fill_value = 'extrapolate')(my_ParcelP)
lnb = np.where(spinup_data['temperature_int'] > my_ParcelT)[0]
i = lnb.max()
while i in lnb:
    i -= 1

lnb = i + 1
p_inv = interpolate.interp1d(radiosonde_data['height'], radiosonde_data['pressure'], fill_value = 'extrapolate')(3500)
iz = np.where(np.abs(my_ParcelP - p_inv) == np.min(np.abs(my_ParcelP - p_inv)))[0][0]

spinup_data['temperature_new'] = spinup_data['temperature_int']*1.
spinup_data['temperature_new'][iz:lnb] = my_ParcelT[iz:lnb]
spinup_data['temperature_new'] = interpolate.interp1d(my_ParcelP, spinup_data['temperature_new'], fill_value = 'extrapolate')(spinup_data['pressure'])

# recompute the dewpoint....
spinup_data['q'] = getQ(spinup_data['temperature_new'], spinup_data['RH'], spinup_data['pressure'], t_units = 'K', p_units = 'hPa')
spinup_data['dewpoint'] = getDew(spinup_data['q'], spinup_data['pressure'], q_units = 'kg/kg', p_units = 'hPa')

print 'stronger inversion skew t'
plotSkewT(spinup_data['temperature_new'][:415]-273.15, spinup_data['dewpoint'][:415]-273.15, spinup_data['pressure'][:415], u = spinup_data['u'][:415], v = spinup_data['v'][:415], CAPE = True)
plt.savefig('../strong_inversion_skewT.png', dpi = 150)
plt.close('all')

spinup_data['theta_new'] = spinup_data['temperature_new']*(1000.0/spinup_data['pressure'])**(Rd/cpd)
spinup_data['theta_new'] = np.array([spinup_data['theta_new'][idx] for idx in range(len(spinup_data['theta_new'])) if spinup_data['height'][idx] < 16649])
spinup_data['theta_new_z'] = np.array([spinup_data['height'][idx] for idx in range(len(spinup_data['height'])) if spinup_data['height'][idx] < 16649])
spinup_data['theta_new'] = np.concatenate((spinup_data['theta_new'], np.array([1200.0])), axis = 0)
spinup_data['theta_new_z'] = np.concatenate((spinup_data['theta_new_z'], np.array([40000.0])), axis = 0)
# RDP to a few heights
heights, spinup_data['theta_new'] = RDP(spinup_data['theta_new_z'], spinup_data['theta_new'], 0.1)
print [int(h) for h in heights]
print [round(t, 2) for t in spinup_data['theta_new']]
print len(heights)

my_message = 'heights:'
for h in heights:
    my_message = my_message +str(int(h)) + ','

# remove the last comma and add a new line
my_message = my_message[:-1] + '\n'

for t in spinup_data['theta_new']:
    my_message = my_message + str(round(t, 2)) + ','

# remove the last comma and add a new line
my_message = my_message[:-1] + '\n'
my_message = my_message + str(len(heights))

send_email(message = my_message, subject = 'Stronger inversion', attachments = ['../original_spinup_skewT.png', '../strong_inversion_skewT.png'])



