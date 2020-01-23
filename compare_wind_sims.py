import numpy as np
import matplotlib.pyplot as plt
from SkewT_archer import getQ, g, cpd, p0, Rd, PTtoTemp
from scipy import interpolate, integrate

def get_IC_from_txt(experiment):
    """
    Opens and reads the initial conditions text files.
    """
    ### Open the initial conditions text files ###
    with open('../InitialFields_Moisture/InitialFields_Moisture_' + experiment + '.txt') as moisture:
        moisture = moisture.read()
    with open('../InitialFields_Temperature/InitialFields_Temperature_' + experiment + '.txt') as temperature:
        temperature = temperature.read()
    with open('../InitialFields_Wind/InitialFields_Wind_' + experiment + '.txt') as wind:
        wind = wind.read()
    
    ### Read in the data ###
    # The heights for the moisture data, km
    mv_z = np.array([float(x) for x in moisture.split('\n')[2].split(' ')[1].split(',')])/1000.
    # The moisture data
    mv_data = np.array([float(x) for x in moisture.split('\n')[4].split(' ')[1].split(',')])*100.

    # The heights for the theta data, km
    theta_z = np.array([float(x) for x in temperature.split('\n')[2].split(' ')[1].split(',')])
    # The theta data, K
    theta_data = np.array([float(x) for x in temperature.split('\n')[4].split(' ')[1].split(',')])
    
    # The heights for the wind data, km
    u_z = np.array([float(x) for x in wind.split('\n')[2].split(' ')[1].split(',')])/1000.
    # The wind data, m/s
    u_data = np.array([float(x) for x in wind.split('\n')[3].split(' ')[1].split(',')])
    v_data = np.array([float(x) for x in wind.split('\n')[4].split(' ')[1].split(',')])
    
    ### Estimate the pressure so that we can compute the specific humidity ###
    p_sfc = 101700. # the surface pressure used to initialise the simulations
    temperature_data = theta_data - g*theta_z/cpd
    
    pressure_data    = np.zeros_like(temperature_data)
    pressure_data[0] = p_sfc*1. # set the surface pressure
    
    dz = 1.
    z = 0
    temperature_interp = interpolate.interp1d(x = theta_z, y = temperature_data)
    for k in range(1, len(theta_z)):
        rho_0 = pressure_data[k-1]/(Rd*temperature_interp(z)) # get the air density at the level below
        p_temp = pressure_data[k-1] - g*rho_0*dz
        z += dz
        
        # use that air density to compute the pressure slightly above and iterate until just below the next level
        while (z + dz) < theta_z[k]:
            rho_0 = p_temp/(Rd*temperature_interp(z))
            p_temp = p_temp - g*rho_0*dz
            z += dz
        
        # do the remaining distance to get to the next level
        rho_0 = p_temp/(Rd*temperature_interp(z))
        pressure_data[k] = p_temp - g*rho_0*(theta_z[k] - z)
        z += (theta_z[k] - z)
    
    # iterate to get a better temperature estimate
    temperature_data = PTtoTemp(theta_data, pressure_data, t_units = 'K', p_units = 'Pa')
    
    # iterate to get a better pressure estimate
    pressure_data[0] = p_sfc*1. # set the surface pressure
    
    dz = 1.
    z = 0
    temperature_interp = interpolate.interp1d(x = theta_z, y = temperature_data)
    for k in range(1, len(theta_z)):
        rho_0 = pressure_data[k-1]/(Rd*temperature_interp(z)) # get the air density at the level below
        p_temp = pressure_data[k-1] - g*rho_0*dz
        z += dz
        
        # use that air density to compute the pressure slightly above and iterate until just below the next level
        while (z + dz) < theta_z[k]:
            rho_0 = p_temp/(Rd*temperature_interp(z))
            p_temp = p_temp - g*rho_0*dz
            z += dz
        
        # do the remaining distance to get to the next level
        rho_0 = p_temp/(Rd*temperature_interp(z))
        pressure_data[k] = p_temp - g*rho_0*(theta_z[k] - z)
        z += (theta_z[k] - z)
    
    ### Store the data to a dictionary and return ###
    data_dict = {'mv_z' : mv_z, 'theta_z' : theta_z/1000., 'u_z' : u_z,
                 'RH_data' : mv_data, 'theta_data' : theta_data, 'u_data' : u_data, 'v_data' : v_data, 'mv_data' : getQ(interpolate.interp1d(x = theta_z/1000., y = temperature_data, fill_value = 'extrapolate')(mv_z), mv_data, interpolate.interp1d(x = theta_z/1000., y = pressure_data, fill_value = 'extrapolate')(mv_z), t_units = 'K', p_units = 'Pa')*1000.0}
    
    return data_dict

def plot_wind(ax, dict_in, plt_color):
    """
    Creates a panel of the initial conditions plot from an input dictionary
    of the form returned by get_IC_from_txt.
    
    This panel contains the u- and v- wind components.
    
    The axis to plot the data on, ax, and the dictionary containing the data,
    dict_in, must be provided.
    """
    ax.plot(dict_in['u_data'], dict_in['u_z'], color = plt_color)
    ax.plot(dict_in['v_data'], dict_in['u_z'], ls = '--', color = plt_color)
    

def plot_theta(ax, dict_in, plt_color, data_label):
    """
    Creates a panel of the initial conditions plot from an input dictionary
    of the form returned by get_IC_from_txt.
    
    This panel contains the theta data.
    
    The axis to plot the data on, ax, and the dictionary containing the data,
    dict_in, must be provided.
    """
    ax.plot(dict_in['theta_data'], dict_in['theta_z'], color = plt_color, label = data_label)
    

def plot_moisture(ax, dict_in, plt_color):
    """
    Creates a panel of the initial conditions plot from an input dictionary
    of the form returned by get_IC_from_txt.
    
    This panel contains the moisture data.
    
    The axis to plot the data on, ax, and the dictionary containing the data,
    dict_in, must be provided.
    """
    ax.plot(dict_in['RH_data'], dict_in['mv_z'], color = plt_color)


def plot_q(ax, dict_in, plt_color):
    """
    Creates a panel of the initial conditions plot from an input dictionary of
    the form returned by get_IC_from_txt.
    
    This panel contains the specific humidity data.
    
    The axis to plot the data on, ax, and the dictionary containinf the data,
    dict_in, must be provided.
    """
    ax.plot(dict_in['mv_data'], dict_in['mv_z'], color = plt_color)


### Grab the data ###
# Write a list of expected experiment names
experiments = ['Spinup_U05', 'Spinup_U06', 'Spinup_U07', 'Spinup_U08', 'Spinup_U09', 'Control']
# Initialise a dictionary to store the data read in for each experiment
all_data = {}
# Initialise a dictionary to store colors for each experiment
plt_colors = {}
import matplotlib
my_cmap = matplotlib.cm.get_cmap('Reds')
for experiment in experiments:
    # Create a dictionary entry for each experiment
    # each entry contains a dictionary with the initial conditions
    all_data[experiment] = get_IC_from_txt(experiment)
    plt_colors[experiment] = my_cmap((experiments.index(experiment)+2)/float(len(experiments)+3))


### Make the plot Comparison ###
my_ylim = [0, 4]
fig = plt.figure(figsize=(8,8))
ax1 = fig.add_subplot(2, 2, 1)
ax1.set_xlim([-12, 2])
ax1.set_ylim(my_ylim)
ax1.set_xlabel(u'Wind (m s$^{-1}$)')
ax1.set_ylabel('Height (km)')
ax1.text(1, 3.0, 'V')
ax1.text(-4, 3.0, 'U')
ax1.set_title('a) Wind')

ax2 = fig.add_subplot(2, 2, 2)
ax2.set_xlabel('RH (%)')
ax2.set_ylim(my_ylim)
ax2.set_yticklabels([''])
ax2.set_title('b) Relative Humidity')

ax3 = fig.add_subplot(2, 2, 3)
ax3.set_xlabel('$\\theta$ (K)')
ax3.set_ylim(my_ylim)
ax3.set_xlim([298, 320])
ax3.set_ylabel('Height (km)')
ax3.set_title('c) Potential Temperature')

ax4 = fig.add_subplot(2, 2, 4)
ax4.set_xlabel('q$_{v}$ (g kg$^{-1}$)')
ax4.set_ylim(my_ylim)
ax4.set_yticklabels([''])
ax4.set_xlim([0, 20])
ax4.set_title('d) Specific Humidity')

# Do the plotting
for experiment in experiments:
    plot_wind(ax1, all_data[experiment], plt_color = plt_colors[experiment])
    plot_moisture(ax2, all_data[experiment], plt_color = plt_colors[experiment])
    plot_theta(ax3, all_data[experiment], plt_color = plt_colors[experiment], data_label = [experiment[-3:] if experiment != 'Control' else 'U10'][0])
    plot_q(ax4, all_data[experiment], plt_color = plt_colors[experiment])

ax3.legend(loc = 'lower right', frameon = False)
plt.subplots_adjust(hspace = 0.3)
plt.savefig('../spinup_plots/wind_IC.png', dpi = 500, bbox_inches = 'tight')
plt.show()

