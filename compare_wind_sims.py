import numpy as np
import matplotlib.pyplot as plt

def get_IC_from_txt(experiment):
    """
    Opens and reads the initial conditions text files.
    """
    ### Open the initial conditions text files ###
    with open('../InitialFields_Moisture_' + experiment + '.txt') as moisture:
        moisture = moisture.read()
    with open('../InitialFields_Temperature_' + experiment + '.txt') as temperature:
        temperature = temperature.read()
    with open('../InitialFields_Wind_' + experiment + '.txt') as wind:
        wind = wind.read()
    
    ### Read in the data ###
    # The heights for the moisture data, km
    mv_z = np.array([float(x) for x in moisture.split('\n')[2].split(' ')[1].split(',')])/1000.
    # The moisture data
    mv_data = np.array([float(x) for x in moisture.split('\n')[4].split(' ')[1].split(',')])*100.

    # The heights for the theta data, km
    theta_z = np.array([float(x) for x in temperature.split('\n')[2].split(' ')[1].split(',')])/1000.
    # The theta data, K
    theta_data = np.array([float(x) for x in temperature.split('\n')[4].split(' ')[1].split(',')])
    
    # The heights for the wind data, km
    u_z = np.array([float(x) for x in wind.split('\n')[2].split(' ')[1].split(',')])/1000.
    # The wind data, m/s
    u_data = np.array([float(x) for x in wind.split('\n')[3].split(' ')[1].split(',')])
    v_data = np.array([float(x) for x in wind.split('\n')[4].split(' ')[1].split(',')])

    ### Store the data to a dictionary and return ###
    data_dict = {'mv_z' : mv_z, 'theta_z' : theta_z, 'u_z' : u_z,
                 'mv_data' : mv_data, 'theta_data' : theta_data, 'u_data' : u_data, 'v_data' : v_data}
    
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
fig = plt.figure(figsize = (15, 6))
ax1 = fig.add_subplot(1, 3, 1)
ax1.set_xlim([-12, 2])
ax1.set_ylim(my_ylim)
ax1.set_xlabel(u'u (m s$^{-1}$) and v (m s$^{-1}$)')
ax1.set_ylabel('height (km)')
ax1.text(1, 3.0, 'V')
ax1.text(-4, 3.0, 'U')

ax2 = fig.add_subplot(1, 3, 2)
ax2.set_xlabel('RH (%)')
ax2.set_ylim(my_ylim)
ax2.set_yticklabels([''])

ax3 = fig.add_subplot(1, 3, 3)
ax3.set_xlabel('$\\theta$ (K)')
ax3.set_ylim(my_ylim)
ax3.set_xlim([298, 320])
ax3.set_yticklabels([''])

# Do the plotting
for experiment in experiments:
    plot_wind(ax1, all_data[experiment], plt_color = plt_colors[experiment])
    plot_moisture(ax2, all_data[experiment], plt_color = plt_colors[experiment])
    plot_theta(ax3, all_data[experiment], plt_color = plt_colors[experiment], data_label = [experiment[-3:] if experiment != 'Control' else 'U10'][0])

ax3.legend(loc = 'lower right')
plt.savefig('../spinup_plots/winds_comparison.png', dpi = 150, bbox_inches = 'tight')
plt.show()

