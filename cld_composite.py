import os
import logging
#from configparser import ConfigParser
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#import iris # no module named iris
from utils import label_clds, cloud_center, cloud_edge, cloud_radius, cloud_geometry_center
from thermodynamics import Epsilon, Cp_da, Cp_v, Cp_lw, Lv, Pref, Rs_da, Ls
import json
import pickle
from scipy.interpolate import interp1d
import math
import netCDF4
from STASH_keys import *
from analysis_tools import send_email

##################################################################################
############   Set some parameters   #############################################
############   Set file path         #############################################
############   Set the condition for object identification    ####################
##################################################################################

dx  = 1600.0

ud_thres = 0.0
dd_thres = -0.0

starttime  = 0  # start time of analysis in seconds
endtime    = 0  # end time of analysis in seconds

nsig = 1.0
step = 0.1

################## Percentile Setup ###################
############################################################
percent    = np.arange(0,100+step,step)

pickup_ud  = [990]  # set the percentile threshold to pick up updraft
pickup_dd  = [10]   # set the percentile threshold to pick up downdraft
pickup_env = [0]

#### different conditions
nameoption = ['plume_only', 'clw_only', 'plume_pls_clw', 'clw_pls_plume', 'plume_no_clw', 'clw_no_plume']
nanum = 1
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
    resultpath = resultpath_base + str(len(pickup_ud)) + "/nothres/m0400_g330/new_time/cld_slice/new_buoy/1000percent/highli_thres/" + nameoption[nanum] + "/test_ud_" + str(pickup_ud[0]) + "/"

if not os.path.exists(resultpath):
    os.makedirs(resultpath)

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

### write to the file ###

### write the 3D field into file ###
np.savez(resultpath+'primary_field', ud_field_flg=ud_field_flg, theta=theta, thetav=thetav, thetal=thetal, theta_p=theta_prime, thetav_p=thetav_prime, thetal_p=thetal_prime, ww=WW, ww_prime=ww_prime, qv=qv, ql=mcl, qt=mt, qt_prime=mt_prime, qcld=qcld, hgt=hgt)

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
pickle.dump(mcl_ud_stat,f)
f.close()

## qcld in updraft ##
f = open(resultpath+'/qcld_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
pickle.dump(qcld_ud_stat,f)
f.close()

## qt in updraft ##
f = open(resultpath+'/qt_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
pickle.dump(mt_ud_stat,f)
f.close()

## qt perturbation in updraft ##
f = open(resultpath+'/qt_p_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','wb')
pickle.dump(mt_p_ud_stat,f)
f.close()

send_email('Gathering cloud statistics complete.\n\nRegards,\ncld_composite.py', 'cld_composite.py --> complete!', attachments = [], isAttach = False)

