from regrid_winds import main as winds
from zi_lcl import main as zi_lcl
from create_balanced_initial_conditions import main as create_IC

"""
Overarching code structure to run 'regrid_winds.py' and 'zi_lcl.py'

'regrid_winds.py'
Takes the UM output wind components (u, v, and w) and regrids them from the 
native c-grid so that they are all on w grid points.

'zi_lcl.py'
Takes thermodynamic data and the regridded wind data to compute the boundary 
layer height, the lifting condensation level, and the cloud top height (if 
cloud liquid water content is available). This is then stored in netCDF files 
following the time splits of the original data. Finally, a time-series of the 
evolution in these newly computed parameters is plotted and saved.
"""

path          = raw_input('Input path:')
ID            = raw_input('Input experiment ID:')
l_spinup      = bool(int(raw_input('Is this a spinup experiment? (yes = 1, no = 0):')))
if not l_spinup:
    l_short   = bool(int(raw_input('Is this a short experiment? (yes = 1, no = 0):')))
else:
    l_short   = False
create_netCDF = bool(int(raw_input('Do you want to create zi_lcl netCDF? (yes = 1, no = 0):')))

print 'path         : ' + path
print 'ID           : ' + ID
print 'l_spinup     : ' + str(l_spinup)
print 'l_short      : ' + str(l_short)
print 'create_netCDF: ' + str(create_netCDF)

winds(path, l_spinup, l_short)
zi_lcl(path, ID, l_spinup, l_short, create_netCDF)
if l_spinup:
    UG = float(raw_input('Input zonal geostrophic wind:'))
    create_IC(path, ID, UG)
