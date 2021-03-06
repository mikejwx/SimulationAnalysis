import os
import logging
#from configparser import ConfigParser
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
<<<<<<< HEAD
#import iris # no module named iris
=======
import iris
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
from utils import label_clds, cloud_center, cloud_edge, cloud_radius, cloud_geometry_center
from thermodynamics import Epsilon, Cp_da, Cp_v, Cp_lw, Lv, Pref, Rs_da, Ls
import json
import pickle
from scipy.interpolate import interp1d
import math
import netCDF4
<<<<<<< HEAD
from STASH_keys import *
from analysis_tools import send_email
=======
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696

##################################################################################
############   Set some parameters   #############################################
############   Set file path         #############################################
############   Set the condition for object identification    ####################
##################################################################################

<<<<<<< HEAD
dx  = 1600.0
=======
dx  = 100
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696

ud_thres = 0.0
dd_thres = -0.0

starttime  = 0  # start time of analysis in seconds
endtime    = 0  # end time of analysis in seconds

nsig = 1.0
step = 0.1
<<<<<<< HEAD

=======
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
################## Percentile Setup ###################
############################################################
percent    = np.arange(0,100+step,step)

pickup_ud  = [990]  # set the percentile threshold to pick up updraft
pickup_dd  = [10]   # set the percentile threshold to pick up downdraft
pickup_env = [0]

#### different conditions
nameoption = ['plume_only', 'clw_only', 'plume_pls_clw', 'clw_pls_plume', 'plume_no_clw', 'clw_no_plume']
nanum = 1
<<<<<<< HEAD
print 'name option = ' + nameoption[nanum]

## set the path for the files to be processed ###
filepath   = '/nerc/n02/n02/xb899100/CloudTrail/Control_1600m_HRIC_INV/'

## set path for outputing the results
resultpath_base = filepath + 'CT_results/'

if len(pickup_ud)==3:
    resultpath = resultpath_base + nameoption[nanum] + "/test_ud_" + str(pickup_ud[0]) + "_"+str(pickup_ud[1]) + "_" + str(pickup_ud[2]) + "/"
elif len(pickup_ud)==2:
    resultpath = resultpath_base + str(len(pickup_ud)) + "/nothres/m0400_g330/new_time/cld_slice/new_buoy/1000percent/highli_thres/" + nameoption[nanum] + "/test_ud_" + str(pickup_ud[0]) + "_" + str(pickup_ud[1]) + "/"
elif len(pickup_ud)==1:
=======


## set the path for the files to be processed ###
#filepath   = "/gws/nopw/j04/paracon_rdg/users/toddj/MONC_RCE_1.5/MONC_RCE_1.5_m0400_g0330/diagnostic_files/"
filepath   = '/nerc/n02/n02/xb899100/CloudTrail/Control/'

## set path for outputing the results
resultpath_base = '/nerc/n02/n02/xb899100/CloudTrail/Control/results/'
if len(pickup_ud)==3:
    #resultpath = "/home/users/jfgu/Analysis/PLUMEDECOM/results/MONC/TODD/Decom/Updraft"+str(len(pickup_ud))+"/nothres/m0400_g330/new_time/cld_slice/new_buoy/1000percent/highli_thres/"+nameoption[nanum]+"/test_ud_"+str(pickup_ud[0])+"_"+str(pickup_ud[1])+"_"+str(pickup_ud[2])+"/"
    resultpath = resultpath_base + nameoption[nanum] + "/test_ud_" + str(pickup_ud[0]) + "_"+str(pickup_ud[1]) + "_" + str(pickup_ud[2]) + "/"
elif len(pickup_ud)==2:
    #resultpath = "/home/users/jfgu/Analysis/PLUMEDECOM/results/MONC/TODD/Decom/Updraft"+str(len(pickup_ud))+"/nothres/m0400_g330/new_time/cld_slice/new_buoy/1000percent/highli_thres/"+nameoption[nanum]+"/test_ud_"+str(pickup_ud[0])+"_"+str(pickup_ud[1])+"/"
    resultpath = resultpath_base + str(len(pickup_ud)) + "/nothres/m0400_g330/new_time/cld_slice/new_buoy/1000percent/highli_thres/" + nameoption[nanum] + "/test_ud_" + str(pickup_ud[0]) + "_" + str(pickup_ud[1]) + "/"
elif len(pickup_ud)==1:
    #resultpath = "/home/users/jfgu/Analysis/PLUMEDECOM/results/MONC/TODD/Decom/Updraft"+str(len(pickup_ud))+"/nothres/m0400_g330/new_time/cld_slice/new_buoy/1000percent/highli_thres/"+nameoption[nanum]+"/test_ud_"+str(pickup_ud[0])+"/"
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
    resultpath = resultpath_base + str(len(pickup_ud)) + "/nothres/m0400_g330/new_time/cld_slice/new_buoy/1000percent/highli_thres/" + nameoption[nanum] + "/test_ud_" + str(pickup_ud[0]) + "/"

if not os.path.exists(resultpath):
    os.makedirs(resultpath)

<<<<<<< HEAD
############# Get lists of all of the files    #####################
all_files = os.listdir(filepath)

# remove 'swath' files
all_files = [my_file for my_file in all_files if 'swath' not in my_file]

#### begin to get time dimensions for each file and also x,y,z dimensions  ####
# just go through one data type to get this information, e.g. the lwp file
lwp_fn = [my_file for my_file in all_files if 'lwp' in my_file][0]

## get time dimensions for each file ##
data_nc = netCDF4.Dataset(filepath + lwp_fn, 'r')
nt = data_nc.variables[lwp_key].shape[0]
data_nc.close()
bouy_fn = [my_file for my_file in all_files if 'bouy' in my_file][0]
data_nc = netCDF4.Dataset(filepath + bouy_fn, 'r')
numz = data_nc.variables[theta_key].shape[1]
numx = data_nc.variables[theta_key].shape[3]
numy = data_nc.variables[theta_key].shape[2]

###############################################
##########       Part II           ############
##########    now read the data    ############
###############################################

#hours = [my_file.split('_')[-1].strip('.nc') for my_file in all_files if ('results' not in my_file) and ('output' not in my_file)]
#hours = list(set(hours))
#hours.sort()
hours = ['09','21']
print 'Reading the data'
## basic variables
for hour in hours:
    with netCDF4.Dataset(filepath + 'fluxes_' + hour + '.nc', 'r') as fluxes_nc:
        if hour == hours[0]:
            prs = fluxes_nc.variables[pthe_key][1:,:,:,:]*1.
            rho = fluxes_nc.variables[rho_key][1:,:,:,:]*1.
            hgt_rho = fluxes_nc.variables['rholev_zsea_rho'][:]*1.
        else:
            prs = np.concatenate((prs, fluxes_nc.variables[pthe_key][:]), axis = 0)
            rho = np.concatenate((rho, fluxes_nc.variables[rho_key][:]), axis = 0)
    
    with netCDF4.Dataset(filepath + 'bouy_' + hour + '.nc', 'r') as bouy_nc:
        if hour == hours[0]:
            theta = bouy_nc.variables[theta_key][1:,:,:,:]*1.
            qv    = bouy_nc.variables[q_key][1:,:,:,:]*1.
        else:
            theta = np.concatenate((theta, bouy_nc.variables[theta_key][:]), axis = 0)
            qv    = np.concatenate((qv, bouy_nc.variables[q_key][:]), axis = 0)
    
    with netCDF4.Dataset(filepath + 'wind_' + hour + '.nc', 'r') as wind_nc:
        if hour == hours[0]:
            WW  = wind_nc.variables[w_key][1:,:,:,:]*1.
            hgt = wind_nc.variables['thlev_zsea_theta'][:]*1.
        else:
            WW  = np.concatenate((WW, wind_nc.variables[w_key][:]), axis = 0)
    
    with netCDF4.Dataset(filepath + 'mr_' + hour + '.nc', 'r') as mr_nc:
        if hour == hours[0]:
            mv    = mr_nc.variables[mv_key][1:,:,:,:]*1.
            mcl   = mr_nc.variables[mcl_key][1:,:,:,:]*1.
            mrain = mr_nc.variables[mr_key][1:,:,:,:]*1.
            mci   = mr_nc.variables[mci_key][1:,:,:,:]*1.
            mg    = mr_nc.variables[mg_key][1:,:,:,:]*1.
        else:
            mv    = np.concatenate((mv, mr_nc.variables[mv_key][:]), axis = 0)
            mcl   = np.concatenate((mcl, mr_nc.variables[mcl_key][:]), axis = 0)
            mrain = np.concatenate((mrain, mr_nc.variables[mr_key][:]), axis = 0)
            mci   = np.concatenate((mci, mr_nc.variables[mci_key][:]), axis = 0)
            mg    = np.concatenate((mg, mr_nc.variables[mg_key][:]), axis = 0)


thetav = theta*(1.0 + mv/Epsilon)/(1.0 + mv + mcl + mrain + mci + mg)
thetal = theta - Lv/Cp_da*mcl*(Pref/prs)**(Rs_da/Cp_da) - Ls/Cp_da*mci*(Pref/prs)**(Rs_da/Cp_da)
mt     = mcl + mv + mrain + mci + mg
qcld   = mcl + mci   # for deep convection; if shallow convection, qcld=ql

### Get perturbations
thetav_prime = np.array([[thetav[it,iz,:,:] - thetav[it,iz,:,:].mean() for iz in range(thetav.shape[1])] for it in range(thetav.shape[0])])
thetal_prime = np.array([[thetal[it,iz,:,:] - thetal[it,iz,:,:].mean() for iz in range(thetal.shape[1])] for it in range(thetal.shape[0])])
theta_prime = np.array([[theta[it,iz,:,:] - theta[it,iz,:,:].mean() for iz in range(theta.shape[1])] for it in range(theta.shape[0])])
ww_prime = np.array([[WW[it,iz,:,:] - WW[it,iz,:,:].mean() for iz in range(WW.shape[1])] for it in range(WW.shape[0])])
mt_prime = np.array([[mt[it,iz,:,:] - mt[it,iz,:,:].mean() for iz in range(mt.shape[1])] for it in range(mt.shape[0])])

## cloud field
ud_field          = np.zeros((nt,len(pickup_ud),numz+1,numy,numx))
dd_field          = np.zeros((nt,len(pickup_dd),numz+1,numy,numx))
env_field         = np.zeros((nt,len(pickup_env),numz+1,numy,numx))
ud_sizes          = np.zeros((nt,len(pickup_ud),numz+1))
ud_field_flg      = np.zeros((nt,len(pickup_ud),numz+1,numy,numx))

# flags
dd_flg            = np.zeros((nt,len(pickup_dd),numz+1,numy,numx))  # downdraft flag
env_flg           = np.zeros((nt,len(pickup_env),numz+1,numy,numx))  # environment flag

## cloud statistics
thetav_dd_stat       = {}
thetav_p_dd_stat     = {}
thetal_dd_stat       = {}
thetal_p_dd_stat     = {}
theta_dd_stat        = {}
theta_p_dd_stat      = {}
ww_dd_stat           = {}
ww_p_dd_stat         = {}
qv_dd_stat           = {}
ql_dd_stat           = {}
qt_dd_stat           = {}
qt_p_dd_stat         = {}

###################################### Lets's Begin ##############################################
####################################### Get Data #################################################
percenval_ud = np.array([[np.nanpercentile(WW[it,iz,WW[it,iz,:,:] <= dd_thres], percent) for iz in range(WW.shape[1])] for it in range(WW.shape[0])])
percenval_dd = np.array([[np.nanpercentile(WW[it,iz,WW[it,iz,:,:] <= dd_thres], percent) for iz in range(WW.shape[1])] for it in range(WW.shape[0])])

# Initialise lists for the updraught stats
# variables
thetav_ud_stat        = []
thetav_p_ud_stat      = []
thetal_ud_stat        = []
thetal_p_ud_stat      = []
theta_ud_stat         = []
theta_p_ud_stat       = []
ww_ud_stat            = []
ww_p_ud_stat          = []
qv_ud_stat            = []
mcl_ud_stat           = []
qcld_ud_stat          = []
mt_ud_stat            = []
mt_p_ud_stat          = []

# other stats
ud_sizes_stat         = []
ud_cld_pts_stat       = []

ud_pos_stat           = []
ud_wcenter_stat       = []
ud_mcenter_stat       = []
ud_gcenter_stat       = []
ud_radius_stat        = []
ud_radius_var_stat    = []
ud_radius_m_stat      = []
ud_edge_stat          = []
ud_in_cld_stat        = []

# create an updraught flag
ud_flg  = np.array([[np.where(WW[it,iz,:,:] < percenval_ud[it,iz,pickup_ud], 0.0, 1.0) for iz in range(WW.shape[1])] for it in range(WW.shape[0])])

# create a cloud flag
cld_flg = np.where(mcl > 1e-4, 1.0, 0.0)

if nanum == 1:
    # set the updraught flag to be wherever there are clouds
    ud_flg = cld_flg[:]*1.
elif nanum == 2:
    # set the updraught to be wherever there is ascent and clouds
    ud_flg = np.where(cld_flg > 0.0, ud_flg, 0.0)
elif nanum == 3:
    # set the cloud flag to be wherever there is ascent and cloud
    cld_flg = np.where(ud_flg > 0.0, cld_flg, 0.0)     
    # then reset the updraught to be where there is ascent and cloud
    ud_flg  = cld_flg[:]*1.
elif nanum == 4:
    # only ascending non-cloudy regions
    ud_flg  = np.where(cld_flg > 0.0, 0.0, ud_flg)     
elif nanum == 5:
    # only cloud regions that are not ascending
    cld_flg = np.where(ud_flg > 0.0, 0.0, cld_flg)
    ud_flg  = cld_flg

print 'Gathering the statistics'
for it in range(WW.shape[0]):
    # append a time list
    thetav_ud_stat.append([])
    thetav_p_ud_stat.append([])
    thetal_ud_stat.append([])
    thetal_p_ud_stat.append([])
    theta_ud_stat.append([])
    theta_p_ud_stat.append([])
    ww_ud_stat.append([])
    ww_p_ud_stat.append([])
    qv_ud_stat.append([])
    mcl_ud_stat.append([])
    qcld_ud_stat.append([])
    mt_ud_stat.append([])
    mt_p_ud_stat.append([])
    
    # other stats
    ud_sizes_stat.append([])
    ud_cld_pts_stat.append([])
    
    ud_pos_stat.append([])
    ud_wcenter_stat.append([])
    ud_mcenter_stat.append([])
    ud_gcenter_stat.append([])
    ud_radius_stat.append([])
    ud_radius_var_stat.append([])
    ud_radius_m_stat.append([])
    ud_edge_stat.append([])
    ud_in_cld_stat.append([])
    
    for iz in range(WW.shape[1]):
        ### append a heights list ##
        thetav_ud_stat[it].append([])
        thetav_p_ud_stat[it].append([])
        thetal_ud_stat[it].append([])
        thetal_p_ud_stat[it].append([])
        theta_ud_stat[it].append([])
        theta_p_ud_stat[it].append([])
        ww_ud_stat[it].append([])
        ww_p_ud_stat[it].append([])
        qv_ud_stat[it].append([])
        mcl_ud_stat[it].append([])
        qcld_ud_stat[it].append([])
        mt_ud_stat[it].append([])
        mt_p_ud_stat[it].append([])
        
        # other stats
        ud_cld_pts_stat[it].append([])
        
        ud_pos_stat[it].append([])
        ud_wcenter_stat[it].append([])
        ud_mcenter_stat[it].append([])
        ud_gcenter_stat[it].append([])
        ud_radius_stat[it].append([])
        ud_radius_var_stat[it].append([])
        ud_radius_m_stat[it].append([])
        ud_edge_stat[it].append([])
        ud_in_cld_stat[it].append([])
        #################################
        
        # to label connected updrafts
        ud_field     = label_clds(ud_flg[it,iz,:,:], diagonal=True, min_cells=5)[1]
        ud_field_flg = np.where(ud_field < 1.0, ud_field, 1.0)
        
        # colect statistics for updrafts sizes
        max_label      = int(np.max(ud_field))
        ud_sizes_stat[it].append(np.histogram(ud_field, range(1, max_label + 2))[0]) # number of grid cells in each cloud at iz
        
        # get cloud properties for each cloud 
        for num_cld in range(1, max_label+1):
            cld_num = ud_sizes_stat[it][iz][num_cld-1]
            ## select grid points for each plume
            ww_ud        = WW[it,iz,ud_field == num_cld]
            ww_p_ud      = ww_prime[it,iz,ud_field == num_cld]
            rho_ud       = rho[it,iz,ud_field == num_cld]
            theta_ud     = theta[it,iz,ud_field == num_cld]
            theta_p_ud   = theta_prime[it,iz,ud_field == num_cld]
            thetav_ud    = thetav[it,iz,ud_field == num_cld]
            thetav_p_ud  = thetav_prime[it,iz,ud_field == num_cld]
            thetal_ud    = thetal[it,iz,ud_field == num_cld]
            thetal_p_ud  = thetal_prime[it,iz,ud_field == num_cld]
            qv_ud        = qv[it,iz,ud_field == num_cld]
            mcl_ud       = mcl[it,iz,ud_field == num_cld]
            qcld_ud      = qcld[it,iz,ud_field == num_cld]
            mt_ud        = mt[it,iz,ud_field == num_cld]
            mt_p_ud      = mt_prime[it,iz,ud_field == num_cld]
            
            pos_ud        = np.where(ud_field == num_cld)
            # check to see if any of the updraught region intersects the lateral boundaries
            if (((numy-1) in pos_ud[0]) and (0 in pos_ud[0])) or ((numx-1) in pos_ud[1]) and (0 in pos_ud[1]):
                continue
            else:
                if nanum == 0 or nanum == 2 or nanum == 4:
                    pos_cx, pos_cy        = cloud_center(pos_ud, ww_ud)
                elif nanum == 1 or nanum == 3 or nanum == 5:
                    pos_cx, pos_cy        = cloud_center(pos_ud, qcld_ud)
                
                cld_pts = ud_sizes_stat[it][iz][num_cld-1]
                
                pos_mcx, pos_mcy          = cloud_center(pos_ud, ww_ud/ww_ud)
                pos_cld_edge, pos_cld_in  = cloud_edge(pos_ud, ud_field, wrap=True, diagonal=True)
                cld_radius                = cloud_radius([pos_cx, pos_cy], pos_cld_edge, dx)[0]
                pos_gcx, pos_gcy          = cloud_geometry_center(pos_cld_edge, pos_cld_in)
                
                ## get statistics
                ud_pos_stat[it][iz].append(pos_ud)
                ud_wcenter_stat[it][iz].append([pos_cx, pos_cy])
                ud_mcenter_stat[it][iz].append([pos_mcx, pos_mcy])
                ud_radius_stat[it][iz].append(cld_radius)
                ud_radius_var_stat[it][iz].append(np.std(cld_radius)/np.max(cld_radius))
                ud_radius_m_stat[it][iz].append(np.mean(cld_radius))
                ud_edge_stat[it][iz].append(pos_cld_edge)
                ud_in_cld_stat[it][iz].append(pos_cld_in)
                ud_cld_pts_stat[it][iz].append(cld_pts)
                
                if pos_gcx==0 and pos_gcy==0:
                    ud_gcenter_stat[it][iz].append([])
                else:
                    ud_gcenter_stat[it][iz].append([pos_gcx, pos_gcy])
                
                thetav_ud_stat[it][iz].append(thetav_ud)
                thetav_p_ud_stat[it][iz].append(thetav_p_ud)
                theta_ud_stat[it][iz].append(theta_ud)
                theta_p_ud_stat[it][iz].append(theta_p_ud)
                thetal_ud_stat[it][iz].append(thetal_ud)
                thetal_p_ud_stat[it][iz].append(thetal_p_ud)
                ww_ud_stat[it][iz].append(ww_ud)
                ww_p_ud_stat[it][iz].append(ww_p_ud)
                qv_ud_stat[it][iz].append(qv_ud)
                mcl_ud_stat[it][iz].append(mcl_ud)
                qcld_ud_stat[it][iz].append(qcld_ud)
                mt_ud_stat[it][iz].append(mt_ud)
                mt_p_ud_stat[it][iz].append(mt_p_ud)
=======

############# Following part is sorting files  #####################
## get the times of files and sorted them in an ascending order ####
all_file = []
for filename in sorted(os.listdir(filepath)):
    all_file.append(filename)

all_files = []
filenum   = []
filename1 = []
filename2 = []
for num in range(0, len(all_file)):
    strtmp = all_file[num].split('.')[0]
    for ii in range(0, len(strtmp)):
        if strtmp[ii].isdigit() and (not strtmp[ii-1].isdigit()):
           filenum.append(int(strtmp[ii:]))
           filename1.append(strtmp[:ii])

# sort the times
filenum_sort = sorted(filenum)
# get the sorted names
sort_ind     = sorted(range(len(filenum)), key=lambda k: filenum[k])
for num in range(0, len(sort_ind)):
    filename2.append(filename1[sort_ind[num]])

filenum_array = np.array(filenum_sort)
num_start     = np.where(filenum_array==starttime)
num_end       = np.where(filenum_array==endtime)
for num in range(num_start[0][0], num_end[0][0]+1):
    all_files.append(filename2[num]+str(filenum_sort[num])+'.0.nc')

#### begin to get time dimensions for each file and also x,y,z dimensions  ####
num = 0
nt  = []
for filename in all_files:

    ## get time dimensions for each file ##
    data = netCDF4.Dataset(filepath+filename)
    nt.append(data.variables['u'].shape[0])
    numz = data.variables['u'].shape[3]
    numx = data.variables['u'].shape[1]
    numy = data.variables['u'].shape[2]
    num += 1

### get the total number of times and x,y,z dimensions ###
ntimes = np.sum(nt)
nx     = numx
ny     = numy
nz     = numz - 1

print ntimes

###############################################
##########       Part II           ############
##########  now allocate the data  ############
###############################################

## basic variables
prs               = np.zeros((ntimes,nz,ny,nx))
rho               = np.zeros((ntimes,nz,ny,nx))
prs_w             = np.zeros((ntimes,nz+1,ny,nx))
rho_w             = np.zeros((ntimes,nz+1,ny,nx))
WW_z              = np.zeros((ntimes,nz+1,ny,nx))
WW                = np.zeros((ntimes,nz+1,ny,nx))
pref              = np.zeros((ntimes,nz+1))
prs2d             = np.zeros((ntimes,nz+1))
thref             = np.zeros((ntimes,nz+1))
theta             = np.zeros((ntimes,nz+1,ny,nx))
thetav            = np.zeros((ntimes,nz+1,ny,nx))
thetal            = np.zeros((ntimes,nz+1,ny,nx))
qcld              = np.zeros((ntimes,nz+1,ny,nx))
qt                = np.zeros((ntimes,nz+1,ny,nx))
qv                = np.zeros((ntimes,nz+1,ny,nx))
rv                = np.zeros((ntimes,nz+1,ny,nx))
ql                = np.zeros((ntimes,nz+1,ny,nx))
rl                = np.zeros((ntimes,nz+1,ny,nx))
qrain             = np.zeros((ntimes,nz+1,ny,nx))
rrain             = np.zeros((ntimes,nz+1,ny,nx))
qice              = np.zeros((ntimes,nz+1,ny,nx))
rice              = np.zeros((ntimes,nz+1,ny,nx))
qsnow             = np.zeros((ntimes,nz+1,ny,nx))
rsnow             = np.zeros((ntimes,nz+1,ny,nx))
qgrap             = np.zeros((ntimes,nz+1,ny,nx))
rgrap             = np.zeros((ntimes,nz+1,ny,nx))
hgt               = np.zeros((ntimes,nz+1))

## cloud field
ud_field          = np.zeros((ntimes,len(pickup_ud),nz+1,ny,nx)) 
dd_field          = np.zeros((ntimes,len(pickup_dd),nz+1,ny,nx)) 
env_field         = np.zeros((ntimes,len(pickup_env),nz+1,ny,nx)) 
ud_sizes          = np.zeros((ntimes,len(pickup_ud),nz+1)) 
ud_field_flg      = np.zeros((ntimes,len(pickup_ud),nz+1,ny,nx)) 

# flags
ud_flg            = np.zeros((ntimes,len(pickup_ud),nz+1,ny,nx))  # updraft flag
dd_flg            = np.zeros((ntimes,len(pickup_dd),nz+1,ny,nx))  # downdraft flag
env_flg           = np.ones((ntimes,len(pickup_env),nz+1,ny,nx))  # environment flag
cld_flg           = np.zeros((ntimes,len(pickup_ud),nz+1,ny,nx))  # updraft flag

# perturbations and fluxes over the whole domain
thetav_prime      = np.zeros((ntimes,nz+1,ny,nx))  # thetav perturbation with respect to domain averaged value
thetal_prime      = np.zeros((ntimes,nz+1,ny,nx))  # thetal perturbation with respect to domain averaged value
theta_prime       = np.zeros((ntimes,nz+1,ny,nx))  # theta perturbation with respect to domain averaged value
ww_prime          = np.zeros((ntimes,nz+1,ny,nx))  # vertical velocity perturbation with respect to domain averged value
qt_prime          = np.zeros((ntimes,nz+1,ny,nx))  # total water perturbation with respect to domain averged value

## values for splitting the plumes
percenval_ud      = np.zeros((ntimes,nz+1,len(percent)))
percenval_dd      = np.zeros((ntimes,nz+1,len(percent)))


## cloud statistics
thetav_ud_stat       = {}
thetav_dd_stat       = {}
thetav_p_ud_stat     = {}
thetav_p_dd_stat     = {}
thetal_ud_stat       = {}
thetal_dd_stat       = {}
thetal_p_ud_stat     = {}
thetal_p_dd_stat     = {}
theta_ud_stat        = {}
theta_dd_stat        = {}
theta_p_ud_stat      = {}
theta_p_dd_stat      = {}
ww_ud_stat           = {}
ww_dd_stat           = {}
ww_p_ud_stat         = {}
ww_p_dd_stat         = {}
qv_ud_stat           = {}
qv_dd_stat           = {}
ql_ud_stat           = {}
qcld_ud_stat         = {}
ql_dd_stat           = {}
qt_ud_stat           = {}
qt_dd_stat           = {}
qt_p_ud_stat         = {}
qt_p_dd_stat         = {}

ud_sizes_stat        = {}
ud_cld_pts_stat      = {}

ud_pos_stat          = {}
ud_wcenter_stat      = {}
ud_mcenter_stat      = {}
ud_gcenter_stat      = {}
ud_radius_stat       = {}
ud_radius_var_stat   = {}
ud_radius_m_stat     = {}
ud_edge_stat         = {}
ud_in_cld_stat       = {}

###################################### Lets's Begin ##############################################
####################################### Get Data #################################################
num = 0
nfile = 0
#for filename in all_files:
for file_num in range(0, 1):
    filename = all_files[0]
    print filename

    ## get data from each file ##
    data = netCDF4.Dataset(filepath+filename)
    for ntt in range(0, nt[nfile]):
        # get variables
        print 'get variables'
        temp             = data.variables['thref'][ntt,:]
        thref[num,:]     = temp
        temp             = data.variables['prefn'][ntt,:]
        pref[num,:]      = temp
        temp             = data.variables['w'][ntt,:,:,:]
        WW_z[num,:,:,:]  = np.transpose(temp, (2,1,0))
        temp             = data.variables['th'][ntt,:,:,:]
        theta[num,:,:,:] = np.transpose(temp, (2,1,0))
        temp             = data.variables['q_vapour'][ntt,:,:,:]
        rv[num,:,:,:]    = np.transpose(temp, (2,1,0))
        temp             = data.variables['q_cloud_liquid_mass'][ntt,:,:,:]
        rl[num,:,:,:]    = np.transpose(temp, (2,1,0))
        temp             = data.variables['q_rain_mass'][ntt,:,:,:]
        rrain[num,:,:,:] = np.transpose(temp, (2,1,0))
        temp             = data.variables['q_ice_mass'][ntt,:,:,:]
        rice[num,:,:,:]  = np.transpose(temp, (2,1,0))
        temp             = data.variables['q_snow_mass'][ntt,:,:,:]
        rsnow[num,:,:,:] = np.transpose(temp, (2,1,0))
        temp             = data.variables['q_graupel_mass'][ntt,:,:,:]
        rgrap[num,:,:,:] = np.transpose(temp, (2,1,0))
        hgt[num,:]       = data.variables['zn'][:]

        prs2d[num,:]     = pref[num,:]
        ## interpolate the vertical velocity to theta level
        WW[num,1:nz+1,:,:] = 0.5*(WW_z[num,0:nz,:,:] + WW_z[num,1:nz+1,:,:])
        WW[num,0,:,:]      = WW_z[num,1,:,:]
        
        ## convert specific humidity to mixing ratio  
        ## for MONC, no need to convert
        qv[num,:,:,:]     = rv[num,:,:,:]
        ql[num,:,:,:]     = rl[num,:,:,:]
        qrain[num,:,:,:]  = rrain[num,:,:,:]
        qice[num,:,:,:]   = rice[num,:,:,:]
        qsnow[num,:,:,:]  = rsnow[num,:,:,:]
        qgrap[num,:,:,:]  = rgrap[num,:,:,:]

        #qv[num,:,:,:]     = rv[num,:,:,:]/(1.0 - rv[num,:,:,:] - rl[num,:,:,:] - rrain[num,:,:,:] - rice[num,:,:,:] - rsnow[num,:,:,:] - rgrap[num,:,:,:])
        #ql[num,:,:,:]     = rl[num,:,:,:]/(1.0 - rv[num,:,:,:] - rl[num,:,:,:] - rrain[num,:,:,:] - rice[num,:,:,:] - rsnow[num,:,:,:] - rgrap[num,:,:,:])
        #qrain[num,:,:,:]  = rrain[num,:,:,:]/(1.0 - rv[num,:,:,:] - rl[num,:,:,:] - rrain[num,:,:,:] - rice[num,:,:,:] - rsnow[num,:,:,:] - rgrap[num,:,:,:])
        #qice[num,:,:,:]   = rice[num,:,:,:]/(1.0 - rv[num,:,:,:] - rl[num,:,:,:] - rrain[num,:,:,:] - rice[num,:,:,:] - rsnow[num,:,:,:] - rgrap[num,:,:,:])
        #qsnow[num,:,:,:]  = rsnow[num,:,:,:]/(1.0 - rv[num,:,:,:] - rl[num,:,:,:] - rrain[num,:,:,:] - rice[num,:,:,:] - rsnow[num,:,:,:] - rgrap[num,:,:,:])
        #qgrap[num,:,:,:]  = rgrap[num,:,:,:]/(1.0 - rv[num,:,:,:] - rl[num,:,:,:] - rrain[num,:,:,:] - rice[num,:,:,:] - rsnow[num,:,:,:] - rgrap[num,:,:,:])

        ## get the total value of theta from the reference value and perturbation  
        ## for UM, don't need to add the reference conference
        for kk in range(0, nz+1):
            theta[num,kk,:,:]     = theta[num,kk,:,:] + thref[num,kk]
            thetal[num,kk,:,:]    = theta[num,kk,:,:] - Lv/Cp_da*ql[num,kk,:,:]*(Pref/prs2d[num,kk])**(Rs_da/Cp_da) - Ls/Cp_da*qice[num,kk,:,:]*(Pref/prs2d[num,kk])**(Rs_da/Cp_da)

        ## theta_v
        thetav[num,:,:,:]     = theta[num,:,:,:]*(1.0+qv[num,:,:,:]/Epsilon)/(1.0+qv[num,:,:,:]+ql[num,:,:,:]+qrain[num,:,:,:]+qice[num,:,:,:]+qsnow[num,:,:,:]+qgrap[num,:,:,:])

        qt[num,:,:,:]         = ql[num,:,:,:] + qv[num,:,:,:] + qrain[num,:,:,:] + qice[num,:,:,:] + qsnow[num,:,:,:] + qgrap[num,:,:,:]
        qcld[num,:,:,:]       = ql[num,:,:,:] + qice[num,:,:,:]   # for deep convection; if shallow convection, qcld=ql

        ### Get perturbations
        for kk in range(0,nz+1):
            thetav_prime[num,kk,:,:]   = thetav[num,kk,:,:] - np.mean(thetav[num,kk,:,:])
            thetal_prime[num,kk,:,:]   = thetal[num,kk,:,:] - np.mean(thetal[num,kk,:,:])
            theta_prime[num,kk,:,:]    = theta[num,kk,:,:] - np.mean(theta[num,kk,:,:])
            ww_prime[num,kk,:,:]       = WW[num,kk,:,:] - np.mean(WW[num,kk,:,:])
            qt_prime[num,kk,:,:]       = qt[num,kk,:,:] - np.mean(qt[num,kk,:,:])

        ## Begin to mask clouds
        thetav_ud_stat[num]       = {}
        thetav_dd_stat[num]       = {}
        thetav_p_ud_stat[num]     = {}
        thetav_p_dd_stat[num]     = {}
        thetal_ud_stat[num]       = {}
        thetal_dd_stat[num]       = {}
        thetal_p_ud_stat[num]     = {}
        thetal_p_dd_stat[num]     = {}
        theta_ud_stat[num]        = {}
        theta_dd_stat[num]        = {}
        theta_p_ud_stat[num]      = {}
        theta_p_dd_stat[num]      = {}
        ww_ud_stat[num]           = {}
        ww_dd_stat[num]           = {}
        ww_p_ud_stat[num]         = {}
        ww_p_dd_stat[num]         = {}
        qv_ud_stat[num]           = {}
        qv_dd_stat[num]           = {}
        ql_ud_stat[num]           = {}
        qcld_ud_stat[num]         = {}
        ql_dd_stat[num]           = {}
        qt_ud_stat[num]           = {}
        qt_dd_stat[num]           = {}
        qt_p_ud_stat[num]         = {}
        qt_p_dd_stat[num]         = {}

        ud_sizes_stat[num]         = {}
        ud_cld_pts_stat[num]       = {}

        ud_pos_stat[num]           = {}
        ud_wcenter_stat[num]       = {}
        ud_mcenter_stat[num]       = {}
        ud_gcenter_stat[num]       = {}
        ud_radius_stat[num]        = {}
        ud_radius_var_stat[num]    = {}
        ud_radius_m_stat[num]      = {}
        ud_edge_stat[num]          = {}
        ud_in_cld_stat[num]        = {}
        for kk in range(0,nz+1):

            thetav_ud_stat[num][kk]        = {}
            thetav_dd_stat[num][kk]        = {}
            thetav_p_ud_stat[num][kk]      = {}
            thetav_p_dd_stat[num][kk]      = {}
            thetal_ud_stat[num][kk]        = {}
            thetal_dd_stat[num][kk]        = {}
            thetal_p_ud_stat[num][kk]      = {}
            thetal_p_dd_stat[num][kk]      = {}
            theta_ud_stat[num][kk]         = {}
            theta_dd_stat[num][kk]         = {}
            theta_p_ud_stat[num][kk]       = {}
            theta_p_dd_stat[num][kk]       = {}
            ww_ud_stat[num][kk]            = {}
            ww_dd_stat[num][kk]            = {}
            ww_p_ud_stat[num][kk]          = {}
            ww_p_dd_stat[num][kk]          = {}
            qv_ud_stat[num][kk]            = {}
            qv_dd_stat[num][kk]            = {}
            ql_ud_stat[num][kk]            = {}
            qcld_ud_stat[num][kk]          = {}
            ql_dd_stat[num][kk]            = {}
            qt_ud_stat[num][kk]            = {}
            qt_dd_stat[num][kk]            = {}
            qt_p_ud_stat[num][kk]          = {}
            qt_p_dd_stat[num][kk]          = {}        

            ud_sizes_stat[num][kk]         = {}
            ud_cld_pts_stat[num][kk]       = {}
            ud_pos_stat[num][kk]           = {}
            ud_wcenter_stat[num][kk]       = {}
            ud_mcenter_stat[num][kk]       = {}
            ud_gcenter_stat[num][kk]       = {}
            ud_radius_stat[num][kk]        = {}
            ud_radius_var_stat[num][kk]    = {}
            ud_radius_m_stat[num][kk]      = {}
            ud_edge_stat[num][kk]          = {}
            ud_in_cld_stat[num][kk]        = {}
            if len(WW[num,kk,WW[num,kk,:,:]>=ud_thres])!=0 and len(WW[num,kk,WW[num,kk,:,:]<=dd_thres])!=0:

                percenval_ud[num,kk,:]   = np.percentile(WW[num,kk,WW[num,kk,:,:]>=ud_thres], percent)
                percenval_dd[num,kk,:]   = np.percentile(WW[num,kk,WW[num,kk,:,:]<=dd_thres], percent)

                # updraft mask
                for ndd in range(0, len(pickup_ud)):
               
                    thetav_ud_stat[num][kk][ndd]        = []
                    thetav_p_ud_stat[num][kk][ndd]      = []
                    thetal_ud_stat[num][kk][ndd]        = []
                    thetal_p_ud_stat[num][kk][ndd]      = []
                    theta_ud_stat[num][kk][ndd]         = []
                    theta_p_ud_stat[num][kk][ndd]       = []
                    ww_ud_stat[num][kk][ndd]            = []
                    ww_p_ud_stat[num][kk][ndd]          = []
                    qv_ud_stat[num][kk][ndd]            = []
                    ql_ud_stat[num][kk][ndd]            = []
                    qcld_ud_stat[num][kk][ndd]          = []
                    qt_ud_stat[num][kk][ndd]            = []
                    qt_p_ud_stat[num][kk][ndd]          = []

                    ud_sizes_stat[num][kk][ndd]         = []
                    ud_cld_pts_stat[num][kk][ndd]       = []

                    ud_pos_stat[num][kk][ndd]           = []
                    ud_wcenter_stat[num][kk][ndd]       = []
                    ud_mcenter_stat[num][kk][ndd]       = []
                    ud_gcenter_stat[num][kk][ndd]       = []
                    ud_radius_stat[num][kk][ndd]        = []
                    ud_radius_var_stat[num][kk][ndd]    = []
                    ud_radius_m_stat[num][kk][ndd]      = []
                    ud_edge_stat[num][kk][ndd]          = []
                    ud_in_cld_stat[num][kk][ndd]        = []

                    if ndd < len(pickup_ud)-1:
                        ud_flg[num,ndd,kk,:,:]   = np.where(WW[num,kk,:,:]<percenval_ud[num,kk,pickup_ud[ndd]], ud_flg[num,ndd,kk,:,:], 1.0)
                        ud_flg[num,ndd,kk,:,:]   = np.where(WW[num,kk,:,:]>=percenval_ud[num,kk,pickup_ud[ndd+1]], 0.0, ud_flg[num,ndd,kk,:,:])
                    elif ndd == len(pickup_ud)-1:
                        ud_flg[num,ndd,kk,:,:]   = np.where(WW[num,kk,:,:]<percenval_ud[num,kk,pickup_ud[ndd]], ud_flg[num,ndd,kk,:,:], 1.0)
 
                    cld_flg[num,ndd,kk,:,:]      = np.where(qcld[num,kk,:,:]>1e-4, 1.0, cld_flg[num,ndd,kk,:,:])

                    if nanum == 1:
                        ud_flg[num,ndd,kk,:,:]    = cld_flg[num,ndd,kk,:,:]
                    elif nanum == 2:
                        ud_flg[num,ndd,kk,:,:]    = np.where(cld_flg[num,ndd,kk,:,:]>0.0, ud_flg[num,ndd,kk,:,:], 0.0)
                    elif nanum == 3:
                        cld_flg[num,ndd,kk,:,:]   = np.where(ud_flg[num,ndd,kk,:,:]>0.0, cld_flg[num,ndd,kk,:,:], 0.0)     
                        ud_flg[num,ndd,kk,:,:]    = cld_flg[num,ndd,kk,:,:]
                    elif nanum == 4:     
                        ud_flg[num,ndd,kk,:,:]    = np.where(cld_flg[num,ndd,kk,:,:]>0.0, 0.0, ud_flg[num,ndd,kk,:,:])     
                    elif nanum == 5:
                        cld_flg[num,ndd,kk,:,:]   = np.where(ud_flg[num,ndd,kk,:,:]>0.0, 0.0, cld_flg[num,ndd,kk,:,:])
                        ud_flg[num,ndd,kk,:,:]    = cld_flg[num,ndd,kk,:,:]
 
                    # to label connected updrafts
                    ud_field[num,ndd,kk,:,:] = label_clds(ud_flg[num,ndd,kk,:,:], diagonal=True, min_cells=5)[1]
                    ud_field_flg[num,ndd,kk,:,:] = np.where(ud_field[num,ndd,kk,:,:]<1.0, ud_field[num,ndd,kk,:,:], 1.0)
                    # colect statistics for updrafts sizes
                    max_label = int(np.max(ud_field[num,ndd,kk,:,:]))
                    ud_sizes_stat[num][kk][ndd] = np.histogram(ud_field[num,ndd,kk,:,:], range(1, max_label + 2))[0]
                    # get cloud properties for each cloud 
                    for num_cld in range(1, max_label+1):
                        cld_num = ud_sizes_stat[num][kk][ndd][num_cld-1]
                        ## select grid points for each plume
                        ww_ud         = WW[num,kk,ud_field[num,ndd,kk,:,:]==num_cld]
                        ww_p_ud       = ww_prime[num,kk,ud_field[num,ndd,kk,:,:]==num_cld]
                        rho_ud        = rho_w[num,kk,ud_field[num,ndd,kk,:,:]==num_cld]
                        theta_ud      = theta[num,kk,ud_field[num,ndd,kk,:,:]==num_cld]
                        theta_p_ud    = theta_prime[num,kk,ud_field[num,ndd,kk,:,:]==num_cld]
                        thetav_ud     = thetav[num,kk,ud_field[num,ndd,kk,:,:]==num_cld]
                        thetav_p_ud   = thetav_prime[num,kk,ud_field[num,ndd,kk,:,:]==num_cld]
                        thetal_ud     = thetal[num,kk,ud_field[num,ndd,kk,:,:]==num_cld]
                        thetal_p_ud   = thetal_prime[num,kk,ud_field[num,ndd,kk,:,:]==num_cld]
                        qv_ud         = qv[num,kk,ud_field[num,ndd,kk,:,:]==num_cld]
                        ql_ud         = ql[num,kk,ud_field[num,ndd,kk,:,:]==num_cld]
                        qcld_ud       = qcld[num,kk,ud_field[num,ndd,kk,:,:]==num_cld]
                        qt_ud         = qt[num,kk,ud_field[num,ndd,kk,:,:]==num_cld]
                        qt_p_ud       = qt_prime[num,kk,ud_field[num,ndd,kk,:,:]==num_cld]

                        pos_ud        = np.where(ud_field[num,ndd,kk,:,:]==num_cld)
                        if (((ny-1) in pos_ud[0]) and (0 in pos_ud[0])) or ((nx-1) in pos_ud[1]) and (0 in pos_ud[1]):
 
                            continue
 
                        else:
                            if nanum == 0 or nanum == 2 or nanum == 4:
                                pos_cx, pos_cy        = cloud_center(pos_ud, ww_ud)
                            elif nanum == 1 or nanum ==3 or nanum == 5:
                                pos_cx, pos_cy        = cloud_center(pos_ud, qcld_ud)

                            cld_pts = ud_sizes_stat[num][kk][ndd][num_cld-1] 

                            pos_mcx, pos_mcy          = cloud_center(pos_ud, ww_ud/ww_ud)
                            pos_cld_edge, pos_cld_in  = cloud_edge(pos_ud, ud_field[num,ndd,kk,:,:], wrap=True, diagonal=True)
                            cld_radius                = cloud_radius([pos_cx, pos_cy], pos_cld_edge, dx)[0]
                            pos_gcx, pos_gcy          = cloud_geometry_center(pos_cld_edge, pos_cld_in)
 
                            ## get statistics
                            ud_pos_stat[num][kk][ndd].append(pos_ud)
                            ud_wcenter_stat[num][kk][ndd].append([pos_cx, pos_cy])
                            ud_mcenter_stat[num][kk][ndd].append([pos_mcx, pos_mcy])
                            ud_radius_stat[num][kk][ndd].append(cld_radius)
                            ud_radius_var_stat[num][kk][ndd].append(np.std(cld_radius)/np.max(cld_radius))
                            ud_radius_m_stat[num][kk][ndd].append(np.mean(cld_radius))
                            ud_edge_stat[num][kk][ndd].append(pos_cld_edge)
                            ud_in_cld_stat[num][kk][ndd].append(pos_cld_in)
                            ud_cld_pts_stat[num][kk][ndd].append(cld_pts)
                            if pos_gcx==0 and pos_gcy==0:
                                ud_gcenter_stat[num][kk][ndd].append([])
                            else:
                                ud_gcenter_stat[num][kk][ndd].append([pos_gcx, pos_gcy])
 
                            thetav_ud_stat[num][kk][ndd].append(thetav_ud)
                            thetav_p_ud_stat[num][kk][ndd].append(thetav_p_ud)
                            theta_ud_stat[num][kk][ndd].append(theta_ud)
                            theta_p_ud_stat[num][kk][ndd].append(theta_p_ud)
                            thetal_ud_stat[num][kk][ndd].append(thetal_ud)
                            thetal_p_ud_stat[num][kk][ndd].append(thetal_p_ud)
                            ww_ud_stat[num][kk][ndd].append(ww_ud)
                            ww_p_ud_stat[num][kk][ndd].append(ww_p_ud)
                            qv_ud_stat[num][kk][ndd].append(qv_ud)
                            ql_ud_stat[num][kk][ndd].append(ql_ud)
                            qcld_ud_stat[num][kk][ndd].append(qcld_ud)
                            qt_ud_stat[num][kk][ndd].append(qt_ud)
                            qt_p_ud_stat[num][kk][ndd].append(qt_p_ud)
 

        num += 1

    nfile += 1
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696

### write to the file ###

### write the 3D field into file ###
<<<<<<< HEAD
np.savez(resultpath+'primary_field', ud_field_flg=ud_field_flg, theta=theta, thetav=thetav, thetal=thetal, theta_p=theta_prime, thetav_p=thetav_prime, thetal_p=thetal_prime, ww=WW, ww_prime=ww_prime, qv=qv, ql=mcl, qt=mt, qt_prime=mt_prime, qcld=qcld, hgt=hgt)
=======
np.savez(resultpath+'primary_field', ud_field_flg=ud_field_flg, theta=theta, thetav=thetav, thetal=thetal, theta_p=theta_prime, thetav_p=thetav_prime, thetal_p=thetal_prime, ww=WW, ww_prime=ww_prime, qv=qv, ql=ql, qt=qt, qt_prime=qt_prime, qcld=qcld, hgt=hgt)
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696

### write plume statistics ####
## updraft size ##
f = open(resultpath+'/ud_sizes_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
pickle.dump(ud_sizes_stat,f)
f.close()

## updraft cloud points ##
f = open(resultpath+'/ud_cld_pts_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
pickle.dump(ud_cld_pts_stat,f)
f.close()

## updraft position ##
f = open(resultpath+'/ud_pos_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
pickle.dump(ud_pos_stat,f)
f.close()

## updraft weighted center ##
f = open(resultpath+'/ud_wcenter_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
pickle.dump(ud_wcenter_stat,f)
f.close()

## updraft mean center ##
f = open(resultpath+'/ud_mcenter_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
pickle.dump(ud_mcenter_stat,f)
f.close()

## updraft geometric center ##
f = open(resultpath+'/ud_gcenter_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
pickle.dump(ud_gcenter_stat,f)
f.close()

## updraft radius ##
f = open(resultpath+'/ud_radius_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
pickle.dump(ud_radius_stat,f)
f.close()

## updraft radius variance ##
f = open(resultpath+'/ud_radius_var_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
pickle.dump(ud_radius_var_stat,f)
f.close()

## updraft mean radius ##
f = open(resultpath+'/ud_radius_m_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
pickle.dump(ud_radius_m_stat,f)
f.close()

## updraft edge position ##
f = open(resultpath+'/ud_edge_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
pickle.dump(ud_edge_stat,f)
f.close()

## in updraft position ##
f = open(resultpath+'/ud_in_cld_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
pickle.dump(ud_in_cld_stat,f)
f.close()

## thetav in updraft ##
f = open(resultpath+'/thetav_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
pickle.dump(thetav_ud_stat,f)
f.close()

## thetav prime in updraft ##
f = open(resultpath+'/thetav_p_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
pickle.dump(thetav_p_ud_stat,f)
f.close()

## theta in updraft ##
f = open(resultpath+'/theta_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
pickle.dump(theta_ud_stat,f)
f.close()

## theta prime in updraft ##
f = open(resultpath+'/theta_p_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
pickle.dump(theta_p_ud_stat,f)
f.close()

## thetal in updraft ##
f = open(resultpath+'/thetal_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
pickle.dump(thetal_ud_stat,f)
f.close()

## thetal prime in updraft ##
f = open(resultpath+'/thetal_p_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
pickle.dump(thetal_p_ud_stat,f)
f.close()

## vertical motion in updraft ##
f = open(resultpath+'/ww_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
pickle.dump(ww_ud_stat,f)
f.close()

## vertical motion perturbation in updraft ##
f = open(resultpath+'/ww_p_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
pickle.dump(ww_ud_stat,f)
f.close()

## qv in updraft ##
f = open(resultpath+'/qv_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
pickle.dump(qv_ud_stat,f)
f.close()

## ql in updraft ##
f = open(resultpath+'/ql_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
<<<<<<< HEAD
pickle.dump(mcl_ud_stat,f)
=======
pickle.dump(ql_ud_stat,f)
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
f.close()

## qcld in updraft ##
f = open(resultpath+'/qcld_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
pickle.dump(qcld_ud_stat,f)
f.close()

## qt in updraft ##
f = open(resultpath+'/qt_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
<<<<<<< HEAD
pickle.dump(mt_ud_stat,f)
=======
pickle.dump(qt_ud_stat,f)
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
f.close()

## qt perturbation in updraft ##
f = open(resultpath+'/qt_p_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
<<<<<<< HEAD
pickle.dump(mt_p_ud_stat,f)
f.close()

send_email('Gathering cloud statistics complete.\n\nRegards,\ncld_composite.py', 'cld_composite.py --> complete!', attachments = [], isAttach = False)
=======
pickle.dump(qt_p_ud_stat,f)
f.close()

>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696

