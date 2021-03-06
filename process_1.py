from regrid_winds import main as winds
from zi_lcl import main as zi_lcl
from create_balanced_initial_conditions import main as create_IC
<<<<<<< HEAD
from analysis_tools import send_email
=======

>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
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
<<<<<<< HEAD
paths = ['/nerc/n02/n02/xb899100/CloudTrail/Control2/']
for path in paths:
    ID            = path.split('/')[-2] # raw_input('Input experiment ID:')
    l_spinup      = False # bool(int(raw_input('Is this a spinup experiment? (yes = 1, no = 0):')))
    l_short       = False
    #if not l_spinup:
    #    l_short   = bool(int(raw_input('Is this a short experiment? (yes = 1, no = 0):')))
    #else:
    #    l_short   = False
    create_WIND   = True # bool(int(raw_input('Do you want to regrid the winds? (yes = 1, no = 0):')))
    create_netCDF = True # bool(int(raw_input('Do you want to create zi_lcl netCDF? (yes = 1, no = 0):')))
    l_diagnostics = False
    
    print 'path         : ' + path
    print 'ID           : ' + ID
    print 'l_spinup     : ' + str(l_spinup)
    print 'l_short      : ' + str(l_short)
    print 'create_netCDF: ' + str(create_netCDF)
    
    print '\nRegridding winds...'
    if create_WIND:
        winds(path, l_spinup, l_short, l_diagnostics)
    
    print '\nComputing boundary layer and cloud characteristics...'
    zi_lcl(path, ID, l_spinup, l_short, create_netCDF, l_diagnostics)
    
    if l_spinup:
        UG = -10.0#float(raw_input('Input zonal geostrophic wind:'))
        create_IC(path, ID, UG, l_diagnostics)
    
    send_email(message = 'Finished Processing ' + ID, subject='Process 1', attachments=['../zi_lcl_cc_lwp_'+ID+'.png'], isAttach = True)
=======

path          = '/nerc/n02/n02/xb899100/CloudTrail/H150E250/' # raw_input('Input path:')
ID            = 'H150E250' # raw_input('Input experiment ID:')
l_spinup      = False # bool(int(raw_input('Is this a spinup experiment? (yes = 1, no = 0):')))
l_short       = False
#if not l_spinup:
#    l_short   = bool(int(raw_input('Is this a short experiment? (yes = 1, no = 0):')))
#else:
#    l_short   = False
create_WIND   = False # bool(int(raw_input('Do you want to regrid the winds? (yes = 1, no = 0):')))
create_netCDF = True # bool(int(raw_input('Do you want to create zi_lcl netCDF? (yes = 1, no = 0):')))
l_diagnostics = False

print 'path         : ' + path
print 'ID           : ' + ID
print 'l_spinup     : ' + str(l_spinup)
print 'l_short      : ' + str(l_short)
print 'create_netCDF: ' + str(create_netCDF)

print '\nRegridding winds...'
if create_WIND:
    winds(path, l_spinup, l_short)

print '\nComputing boundary layer and cloud characteristics...'
zi_lcl(path, ID, l_spinup, l_short, create_netCDF)

if l_spinup:
    UG = float(raw_input('Input zonal geostrophic wind:'))
    create_IC(path, ID, UG)
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
