import numpy as np
import matplotlib.pyplot as plt
from SkewT_archer import *
from analysis_tools import toComponents, round5, RDP
"""
Plots the radiosonde from a case CT day.
Plots the simplified initial conditions based on that radiosonde in a comparison
    three-panel plot.
Plots the spin-up-based initial conditions in a comparison with the simplified
    initial conditions in a three-panel plot.
"""

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
print 'original skew t'
plotSkewT(radiosonde_data['temperature'], radiosonde_data['dewpoint'], radiosonde_data['pressure'], u = radiosonde_data['u'], v = radiosonde_data['v'], date = '1200 UTC 16th July 2015', CAPE = True)
#plt.savefig('../ic_radiosonde_skewT.png', dpi = 150, bbox_inches = 'tight')
plt.show()

# Plot the three panel for radiosonde and simplified
# Interpolate onto the same height levels
my_heights = np.linspace(0.0, 40000.0, 1000)
from scipy import interpolate
simple_data = {}
simple_data['height'] = my_heights*1.

simple_theta = np.array([299., 299., 299.01, 302.47, 311.38, 318.05, 323.86, 328.95, 335.98, 339.1, 342.42, 354.57, 402.01, 1200.0])
simple_theta_z = np.array([0, 148, 509, 857, 3217, 3840, 5931, 7635, 9714, 10963, 12408, 14189, 16649, 40000])
simple_theta = interpolate.interp1d(x = simple_theta_z, y = simple_theta, fill_value = 'extrapolate')(my_heights)
radiosonde_data['theta_int'] = interpolate.interp1d(x = radiosonde_data['height'], y = radiosonde_data['theta'], fill_value = 'extrapolate')(my_heights)
simple_data['theta'] = simple_theta*1.

simple_RH = np.array([77.5, 81.4, 89.0, 55.4, 38.6, 20., 20.])
simple_RH_z = np.array([0, 148, 509, 857, 3217, 3509, 40000])
simple_RH = interpolate.interp1d(x = simple_RH_z, y = simple_RH, fill_value = 'extrapolate')(my_heights)
radiosonde_data['RH_int'] = interpolate.interp1d(x = radiosonde_data['height'], y = radiosonde_data['RH'], fill_value = 'extrapolate')(my_heights)
simple_data['RH'] = simple_RH*1.

simple_U = np.array([10, 10, 0, 0])
simple_U_z = np.array([0, 9000, 15000, 40000])
simple_U = interpolate.interp1d(x = simple_U_z, y = simple_U, fill_value = 'extrapolate')(my_heights)
radiosonde_data['wind_spd_int'] = interpolate.interp1d(x = radiosonde_data['height'], y = radiosonde_data['wind_spd'], fill_value = 'extrapolate')(my_heights)
simple_data['u'] = -simple_U*1.
simple_data['v'] = np.zeros_like(simple_U)

three_panel(theta_1 = radiosonde_data['theta_int'], theta_2 = simple_theta, RH_1 = radiosonde_data['RH_int'], RH_2 = simple_RH, U_1 = radiosonde_data['wind_spd_int']*0.5144, U_2 = simple_U, height = my_heights/1000.0, label_1 = 'Radiosonde', label_2 = 'Simplified', figname = 'radiosonde_simplified_comp.png')

# skew t of the simple data
simple_data['pressure'] = interpolate.interp1d(radiosonde_data['height'], radiosonde_data['pressure'], fill_value = 'extrapolate')(simple_data['height']) #assume the pressure profile hasn't changed from the radiosonde
simple_data['temperature'] = PTtoTemp(simple_data['theta'], simple_data['pressure'], t_units = 'K', p_units = 'hPa')
simple_data['q'] = getQ(simple_data['temperature'], simple_data['RH'], simple_data['pressure'], t_units = 'K', p_units = 'hPa')
simple_data['dewpoint'] = getDew(simple_data['q'], simple_data['pressure'], q_units = 'kg/kg', p_units = 'hPa')

print 'simplified skew t'
plotSkewT(simple_data['temperature'][:415]-273.15, simple_data['dewpoint'][:415]-273.15, simple_data['pressure'][:415], u = simple_data['u'][:415], v = simple_data['v'][:415], CAPE = True)
plt.show()

# make the inversion at ~3km deeper by 500 m
# follow a moist adiabat from the top of the inversion until it intersects the original simple
new_simple_data = {}
new_simple_data['height'] = my_heights*1.

new_simple_theta   = np.array([299., 299., 299.01, 302.47, 311.38, 318.05, 323.86, 328.95, 335.98, 339.1, 342.42, 354.57, 402.01, 1200.0])
new_simple_theta_z = np.array([0,    148,  509,    857,    3217,   3840,   5931, 7635, 9714, 10963, 12408, 14189, 16649, 40000])
new_simple_theta = interpolate.interp1d(x = new_simple_theta_z, y = new_simple_theta, fill_value = 'extrapolate')(my_heights)
new_simple_data['theta'] = new_simple_theta*1.

new_simple_RH = np.array([77.5, 81.4, 89.0, 55.4, 38.6, 20., 20.])
new_simple_RH_z = np.array([0, 148, 509, 857, 3217, 3509, 40000])
new_simple_RH = interpolate.interp1d(x = new_simple_RH_z, y = new_simple_RH, fill_value = 'extrapolate')(my_heights)
new_simple_data['RH'] = new_simple_RH*1.

new_simple_U = np.array([10, 10, 0, 0])
new_simple_U_z = np.array([0, 9000, 15000, 40000])
new_simple_U = interpolate.interp1d(x = new_simple_U_z, y = new_simple_U, fill_value = 'extrapolate')(my_heights)
new_simple_data['u'] = -new_simple_U*1.
new_simple_data['v'] = np.zeros_like(simple_U)

# skew t of the simple data
new_simple_data['pressure'] = interpolate.interp1d(radiosonde_data['height'], radiosonde_data['pressure'], fill_value = 'extrapolate')(new_simple_data['height']) #assume the pressure profile hasn't changed from the radiosonde
new_simple_data['temperature'] = PTtoTemp(new_simple_data['theta'], new_simple_data['pressure'], t_units = 'K', p_units = 'hPa')
new_simple_data['q'] = getQ(new_simple_data['temperature'], new_simple_data['RH'], new_simple_data['pressure'], t_units = 'K', p_units = 'hPa')
new_simple_data['dewpoint'] = getDew(new_simple_data['q'], new_simple_data['pressure'], q_units = 'kg/kg', p_units = 'hPa')

my_CAPE, my_CIN, my_ParcelT, my_ParcelP, LCLp, LFCp = getCAPE(new_simple_data['temperature'][:415]+7.5, new_simple_data['q'][:415], new_simple_data['pressure'][:415], parcel_type = 1, t_units = 'K', q_units = 'kg/kg', p_units = 'hPa')

new_simple_data['temperature_int'] = interpolate.interp1d(new_simple_data['pressure'], new_simple_data['temperature'], fill_value = 'extrapolate')(my_ParcelP)
lnb = np.where(new_simple_data['temperature_int'] > my_ParcelT)[0]
i = lnb.max()
while i in lnb:
    i -= 1

lnb = i + 1
p_inv = interpolate.interp1d(radiosonde_data['height'], radiosonde_data['pressure'], fill_value = 'extrapolate')(3500)
iz = np.where(np.abs(my_ParcelP - p_inv) == np.min(np.abs(my_ParcelP - p_inv)))[0][0]

new_simple_data['temperature_new'] = new_simple_data['temperature_int']*1.
new_simple_data['temperature_new'][iz:lnb] = my_ParcelT[iz:lnb]
new_simple_data['temperature_new'] = interpolate.interp1d(my_ParcelP, new_simple_data['temperature_new'], fill_value = 'extrapolate')(new_simple_data['pressure'])

print 'stronger inversion skew t'
plotSkewT(new_simple_data['temperature_new'][:415]-273.15, new_simple_data['dewpoint'][:415]-273.15, new_simple_data['pressure'][:415], u = new_simple_data['u'][:415], v = new_simple_data['v'][:415], CAPE = True)
plt.show()

new_simple_data['theta_new'] = new_simple_data['temperature_new']*(1000.0/new_simple_data['pressure'])**(Rd/cpd)
new_simple_data['theta_new'] = np.array([new_simple_data['theta_new'][idx] for idx in range(len(new_simple_data['theta_new'])) if new_simple_data['height'][idx] < 16649])
new_simple_data['theta_new_z'] = np.array([new_simple_data['height'][idx] for idx in range(len(new_simple_data['height'])) if new_simple_data['height'][idx] < 16649])
new_simple_data['theta_new'] = np.concatenate((new_simple_data['theta_new'], np.array([1200.0])), axis = 0)
new_simple_data['theta_new_z'] = np.concatenate((new_simple_data['theta_new_z'], np.array([40000.0])), axis = 0)
# RDP to a few heights
heights, new_simple_data['theta_new'] = RDP(new_simple_data['theta_new_z'], new_simple_data['theta_new'], 0.1)
print [int(h) for h in heights]
print [round(t, 2) for t in new_simple_data['theta_new']]
print len(heights)

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

three_panel(theta_1 = radiosonde_data['theta_int'], theta_2 = spinup_data['theta'], RH_1 = radiosonde_data['RH_int'], RH_2 = spinup_data['RH'], U_1 = radiosonde_data['wind_spd_int']*0.5144, U_2 = spinup_data['U'], height = my_heights/1000.0, label_1 = 'Simplified', label_2 = 'Spin-up', figname = 'radiosonde_spinup_comp.png')

# skew t of the spinup data
spinup_data['pressure'] = interpolate.interp1d(radiosonde_data['height'], radiosonde_data['pressure'], fill_value = 'extrapolate')(my_heights) #assume the pressure profile hasn't changed from the radiosonde
spinup_data['temperature'] = PTtoTemp(spinup_data['theta'], spinup_data['pressure'], t_units = 'K', p_units = 'hPa')
spinup_data['q'] = getQ(spinup_data['temperature'], spinup_data['RH'], spinup_data['pressure'], t_units = 'K', p_units = 'hPa')
spinup_data['dewpoint'] = getDew(spinup_data['q'], spinup_data['pressure'], q_units = 'kg/kg', p_units = 'hPa')
spinup_data['u'] = interpolate.interp1d(x = spinup_data['uv_init_height'], y = spinup_data['u_init_data'], fill_value = 'extrapolate')(my_heights)
spinup_data['v'] = interpolate.interp1d(x = spinup_data['uv_init_height'], y = spinup_data['v_init_data'], fill_value = 'extrapolate')(my_heights)

print 'spin-up skew t'
plotSkewT(spinup_data['temperature'][:415]-273.15, spinup_data['dewpoint'][:415]-273.15, spinup_data['pressure'][:415], u = spinup_data['u'][:415], v = spinup_data['v'][:415], CAPE = True)
plt.show()


