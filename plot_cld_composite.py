import os
import logging
#from configparser import ConfigParser
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#import iris
from utils import label_clds, cloud_center, cloud_edge, cloud_radius
from thermodynamics import Epsilon, Cp_da, Cp_v, Cp_lw, Lv, Pref, Rs_da
import json
import pickle
from scipy.interpolate import interp1d
import math
from analysis_tools import send_email

##############################################################
############   To read the model result   ####################
############   and pick up the sounding   ####################
############   at some specific points    ####################
##############################################################

ntimes = 144
nx  = 1160
ny  = 319
nz  = 127
dtt = 2
dt  = 3
dx  = 100
dz  = 40

ud_thres = 0.0
dd_thres = -0.0

nsig = 1.0
step = 0.1
################## Percentile Setup ###################
############################################################
percent    = np.arange(0,100+step,step)
############################################################

pickup_ud  = [990]
pickup_dd  = [10]
pickup_env = [0]

pickup_cld = [1,1]

####
nameoption = ['plume_only', 'clw_only', 'plume_pls_clw', 'clw_pls_plume', 'plume_no_clw', 'clw_no_plume']
nanum = 1

## set path for outputing the results
resultpath_base = '/nerc/n02/n02/xb899100/CloudTrail/Control/CT_results/'

if len(pickup_ud)==3:
    resultpath = resultpath_base + nameoption[nanum] + "/test_ud_" + str(pickup_ud[0]) + "_"+str(pickup_ud[1]) + "_" + str(pickup_ud[2]) + "/"
elif len(pickup_ud)==2:
    resultpath = resultpath_base + str(len(pickup_ud)) + "/nothres/m0400_g330/new_time/cld_slice/new_buoy/1000percent/highli_thres/" + nameoption[nanum] + "/test_ud_" + str(pickup_ud[0]) + "_" + str(pickup_ud[1]) + "/"
elif len(pickup_ud)==1:
    resultpath = resultpath_base + str(len(pickup_ud)) + "/nothres/m0400_g330/new_time/cld_slice/new_buoy/1000percent/highli_thres/" + nameoption[nanum] + "/test_ud_" + str(pickup_ud[0]) + "/"

### path for saving the figures
figurepath = resultpath + "Figs/"

if not os.path.exists(figurepath):
    os.makedirs(figurepath)

print figurepath

### allocate
ww_prime = np.zeros((ntimes, nz+1, ny, nx))
mt_prime = np.zeros((ntimes, nz+1, ny, nx))

### load primary variables

data               = np.load(resultpath+'primary_field.npz')

ud_field_flg       = data['ud_field_flg']
thetav             = data['thetav']
thetav_prime       = data['thetav_p']
theta              = data['theta']
theta_prime        = data['theta_p']
thetal             = data['thetal']
thetal_prime       = data['thetal_p']
WW                 = data['ww']
#ww_prime           = data['ww_prime']
qv                 = data['qv']
mcl                = data['ql']
qcld               = data['qcld']
mt                 = data['qt']
hgt                = data['hgt']
#qt_prime           = data['qt_prime']

for num in range(0, mt.shape[0]):
    for kk in range(0, mt.shape[1]):
        ww_prime[num,kk,:,:] = WW[num,kk,:,:] - np.mean(WW[num,kk,:,:])
        mt_prime[num,kk,:,:] = mt[num,kk,:,:] - np.mean(mt[num,kk,:,:])


### load cloud statistics
ud_pos_stat        = pickle.load(open(resultpath+'/ud_pos_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

ud_wcenter_stat    = pickle.load(open(resultpath+'/ud_wcenter_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

ud_mcenter_stat    = pickle.load(open(resultpath+'/ud_mcenter_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

ud_gcenter_stat    = pickle.load(open(resultpath+'/ud_gcenter_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

ud_radius_stat     = pickle.load(open(resultpath+'/ud_radius_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

ud_radius_var_stat = pickle.load(open(resultpath+'/ud_radius_var_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

ud_radius_m_stat   = pickle.load(open(resultpath+'/ud_radius_m_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

ud_edge_stat       = pickle.load(open(resultpath+'/ud_edge_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

ud_in_cld_stat     = pickle.load(open(resultpath+'/ud_in_cld_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

thetav_ud_stat     = pickle.load(open(resultpath+'/thetav_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

thetav_p_ud_stat   = pickle.load(open(resultpath+'/thetav_p_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

theta_ud_stat      = pickle.load(open(resultpath+'/theta_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

theta_p_ud_stat    = pickle.load(open(resultpath+'/theta_p_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

thetal_ud_stat     = pickle.load(open(resultpath+'/thetal_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

thetal_p_ud_stat   = pickle.load(open(resultpath+'/thetal_p_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

ww_ud_stat         = pickle.load(open(resultpath+'/ww_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

ww_p_ud_stat       = pickle.load(open(resultpath+'/ww_p_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

qv_ud_stat         = pickle.load(open(resultpath+'/qv_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

mcl_ud_stat        = pickle.load(open(resultpath+'/ql_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

qcld_ud_stat       = pickle.load(open(resultpath+'/qcld_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

mt_ud_stat         = pickle.load(open(resultpath+'/qt_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

mt_p_ud_stat       = pickle.load(open(resultpath+'/qt_p_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

## plot
my_hgts = [1500, 2000, 3000] 
kk_list = [np.where(np.abs(hgt-my_hgt) == np.min(np.abs(hgt-my_hgt)))[0][0] for my_hgt in my_hgts]

## basic setup for composite
nstep = 12                            # number of points on each side of the center for plot
ncld  = 8                             # number of points we want to interpolate within cloud
xplt  = np.arange(-nstep, nstep+1)
ndd   = 0
sx    = 14                            # figure size x
sy    = 16                            # figure size y
v1  = [1.0e-5, 1.0e-4, 5.0e-4, 1.0e-3, 5.0e-3, 1.0e-2]       # contour levels for ql plot

###################

for kk in kk_list:
    # kk = a select height index
    ### set slices list at each height to be empty
    thetav_p_merge = []
    thetal_p_merge = []
    theta_p_merge  = []
    ww_merge       = []
    mt_merge       = []
    mcl_merge      = []
    qcld_merge     = []
    
    ### number of slices for composite
    num_slice = 0
    for num in range(0, thetav.shape[0]):
        # num is a time index
        if (len(ud_in_cld_stat[num][kk]) > 0):
            # if there are any clouds at this time and height
            for nn in range(0, len(ud_in_cld_stat[num][kk])):
                # for each cloud, nn, at this select height and at this time
                if (len(ud_in_cld_stat[num][kk][nn]) > 0) and (len(ud_gcenter_stat[num][kk][nn]) > 0):
                
                ## find out maximum value of variables within the plume
                    thetav_p_max  = np.max(np.abs(thetav_p_ud_stat[num][kk][nn]))
                    thetal_p_max  = np.max(np.abs(thetal_p_ud_stat[num][kk][nn]))
                    theta_p_max   = np.max(np.abs(theta_p_ud_stat[num][kk][nn]))
                    ww_p_max      = np.max(np.abs(ww_p_ud_stat[num][kk][nn]))
                    mt_p_max      = np.max(np.abs(mt_p_ud_stat[num][kk][nn]))
                    mcl_max       = np.max(np.abs(mcl_ud_stat[num][kk][nn]))
                    qcld_max      = np.max(np.abs(qcld_ud_stat[num][kk][nn]))
                    
                ## weighted plume center
                    xcenter = ud_gcenter_stat[num][kk][nn][1]
                    ycenter = ud_gcenter_stat[num][kk][nn][0]
                    
                    cycle = 0
                    while cycle < 4: 
                        
                        if cycle == 0:
                            cx    = int(xcenter)
                            cy    = int(ycenter) 
                        elif cycle == 1:
                            cx    = int(xcenter) + 1
                            cy    = int(ycenter)
                        elif cycle == 2:
                            cx    = int(xcenter)
                            cy    = int(ycenter) + 1
                        else:
                            cx    = int(xcenter) + 1
                            cy    = int(ycenter) + 1
                        
                        ## x, y positions for plume
                        xp = ud_pos_stat[num][kk][nn][1]
                        yp = ud_pos_stat[num][kk][nn][0]
                        
                        ## x,y slice position at center
                        if (cx in xp) and (cy in yp):
                            ys_pos = yp[xp==cx]
                            xs_pos = xp[yp==cy]
                            
                            ## work out the length of x, y slice and intervals
                            ys_length = np.max(ys_pos) - np.min(ys_pos)
                            ys_interval = ys_length/float(ncld)
                            
                            xs_length = np.max(xs_pos) - np.min(xs_pos)
                            xs_interval = xs_length/float(ncld)
                            
                            ysign = (cy - np.max(ys_pos))*(cy - np.min(ys_pos))
                            xsign = (cx - np.max(xs_pos))*(cx - np.min(xs_pos))
                            
                            # make sure that the slice has enough points and doesn't break, and also the center is not outside the plume
                            if ys_length >=2 and xs_length >= 2 and (xs_length+1)==len(xs_pos) and (ys_length+1) == len(ys_pos) and ysign < 0 and xsign < 0:      
                                ## points to interploate
                                ys_interp_pts = np.arange(cy-ys_interval*nstep, cy+ys_interval*(nstep+1), ys_interval)
                                xs_interp_pts = np.arange(cx-xs_interval*nstep, cx+xs_interval*(nstep+1), xs_interval)
                                
                                ## points for interpolation
                                ys_pts = np.arange(np.min([cy-nstep, np.min(ys_interp_pts)]), np.max([cy+nstep+1, np.max(ys_interp_pts+1)]), 1)
                                xs_pts = np.arange(np.min([cx-nstep, np.min(xs_interp_pts)]), np.max([cx+nstep+1, np.max(xs_interp_pts+1)]), 1)
                                
                                ys_temp_pts = np.zeros((len(ys_pts)))
                                xs_temp_pts = np.zeros((len(xs_pts)))
                                
                                ## in case the points for interpolation is near the edge of the domain
                                for pts_num in range(0, len(ys_pts)):
                                    ys_temp_pts[pts_num] = ys_pts[pts_num]
                                    if ys_temp_pts[pts_num] >(ny-1) or ys_temp_pts[pts_num] <0:
                                        ys_temp_pts[pts_num] %= ny
                                
                                for pts_num in range(0, len(xs_pts)):
                                    xs_temp_pts[pts_num] = xs_pts[pts_num]
                                    if xs_temp_pts[pts_num] >(nx-1) or xs_temp_pts[pts_num] <0:
                                        xs_temp_pts[pts_num] %= nx
                                
                                ys_orig_pts = ys_temp_pts.astype(int)
                                xs_orig_pts = xs_temp_pts.astype(int)
                                
                                ## variables at the points for interpolation
                                thetav_p_orig_ys = thetav_prime[num,kk,ys_orig_pts,cx]/thetav_p_max
                                thetal_p_orig_ys = thetal_prime[num,kk,ys_orig_pts,cx]/thetal_p_max
                                theta_p_orig_ys  = theta_prime[num,kk,ys_orig_pts,cx]/theta_p_max
                                ww_orig_ys       = ww_prime[num,kk,ys_orig_pts,cx]/ww_p_max
                                mt_orig_ys       = mt_prime[num,kk,ys_orig_pts,cx]/mt_p_max
                                mcl_orig_ys      = mcl[num,kk,ys_orig_pts,cx]/mcl_max
                                qcld_orig_ys     = qcld[num,kk,ys_orig_pts,cx]/qcld_max
                                thetav_p_orig_xs = thetav_prime[num,kk,cy,xs_orig_pts]/thetav_p_max
                                thetal_p_orig_xs = thetal_prime[num,kk,cy,xs_orig_pts]/thetal_p_max
                                theta_p_orig_xs  = theta_prime[num,kk,cy,xs_orig_pts]/theta_p_max
                                ww_orig_xs       = ww_prime[num,kk,cy,xs_orig_pts]/ww_p_max
                                mt_orig_xs       = mt_prime[num,kk,cy,xs_orig_pts]/mt_p_max
                                mcl_orig_xs      = mcl[num,kk,cy,xs_orig_pts]/mcl_max
                                qcld_orig_xs     = qcld[num,kk,cy,xs_orig_pts]/qcld_max
                                
                                ## interpolation
                                func_thetav_p_ys = interp1d(ys_pts, thetav_p_orig_ys)
                                func_thetal_p_ys = interp1d(ys_pts, thetal_p_orig_ys)
                                func_theta_p_ys  = interp1d(ys_pts, theta_p_orig_ys)
                                func_ww_ys       = interp1d(ys_pts, ww_orig_ys)
                                func_mt_ys       = interp1d(ys_pts, mt_orig_ys)
                                func_mcl_ys      = interp1d(ys_pts, mcl_orig_ys)
                                func_qcld_ys     = interp1d(ys_pts, qcld_orig_ys)
                                func_thetav_p_xs = interp1d(xs_pts, thetav_p_orig_xs)
                                func_thetal_p_xs = interp1d(xs_pts, thetal_p_orig_xs)
                                func_theta_p_xs  = interp1d(xs_pts, theta_p_orig_xs)
                                func_ww_xs       = interp1d(xs_pts, ww_orig_xs)
                                func_mt_xs       = interp1d(xs_pts, mt_orig_xs)
                                func_mcl_xs      = interp1d(xs_pts, mcl_orig_xs)
                                func_qcld_xs     = interp1d(xs_pts, qcld_orig_xs)
                                
                                thetav_p_ys      = func_thetav_p_ys(ys_interp_pts)
                                thetal_p_ys      = func_thetal_p_ys(ys_interp_pts)
                                theta_p_ys       = func_theta_p_ys(ys_interp_pts)
                                ww_ys            = func_ww_ys(ys_interp_pts)
                                mt_ys            = func_mt_ys(ys_interp_pts)
                                mcl_ys           = func_mcl_ys(ys_interp_pts)
                                qcld_ys          = func_qcld_ys(ys_interp_pts)
                                thetav_p_xs      = func_thetav_p_xs(xs_interp_pts)
                                thetal_p_xs      = func_thetal_p_xs(xs_interp_pts)
                                theta_p_xs       = func_theta_p_xs(xs_interp_pts)
                                ww_xs            = func_ww_xs(xs_interp_pts)
                                mt_xs            = func_mt_xs(xs_interp_pts)
                                mcl_xs           = func_mcl_xs(xs_interp_pts)
                                qcld_xs          = func_qcld_xs(xs_interp_pts)
                                
                                thetav_p_merge.append(thetav_p_ys)
                                thetav_p_merge.append(thetav_p_xs)
                                thetal_p_merge.append(thetal_p_ys)
                                thetal_p_merge.append(thetal_p_xs)
                                theta_p_merge.append(theta_p_ys)
                                theta_p_merge.append(theta_p_xs)
                                ww_merge.append(ww_ys)
                                ww_merge.append(ww_xs)
                                mt_merge.append(mt_ys)
                                mt_merge.append(mt_xs)
                                mcl_merge.append(mcl_ys)
                                mcl_merge.append(mcl_xs)
                                qcld_merge.append(qcld_ys)
                                qcld_merge.append(qcld_xs)
                                
                                num_slice += 1
                        
                        cycle += 1
        
    print 'number of slice', num_slice
    
    ######################
    ### plot composite ###
    ######################
    
    # convert lists to np arrays
    thetav_p_array = np.asarray(thetav_p_merge)
    thetal_p_array = np.asarray(thetal_p_merge)
    theta_p_array  = np.asarray(theta_p_merge)
    ww_array       = np.asarray(ww_merge)
    mt_array       = np.asarray(mt_merge)
    #mcl_array      = np.asarray(mcl_merge)
    qcld_array     = np.asarray(qcld_merge)
    
    ## ql might have nan values
    ## remove the rows that have nan values
    #mask           = np.any(np.isnan(mcl_array), axis=1)
    #mcl_array       = mcl_array[~mask]
    mask           = np.any(np.isnan(qcld_array), axis=1)
    qcld_array     = qcld_array[~mask]
    
    ## composite slices
    thetav_p_comp  = np.mean(thetav_p_array, axis=0)
    thetal_p_comp  = np.mean(thetal_p_array, axis=0)
    theta_p_comp   = np.mean(theta_p_array, axis=0)
    ww_comp        = np.mean(ww_array, axis=0)
    mt_comp        = np.mean(mt_array, axis=0)
    #mcl_comp        = np.mean(mcl_array, axis=0)
    qcld_comp      = np.mean(qcld_array, axis=0)
    
    fig = plt.figure(figsize=[sx,sy])
    plt.plot(xplt, thetav_p_comp, '-b', label='$\\theta_v^\prime$')
    plt.plot(xplt, ww_comp, '-r', label='$w$')
    plt.plot(xplt, mt_comp, '-g', label='$m_t$')
    plt.plot(xplt, np.zeros((len(xplt))), '--k')
    plt.legend(loc=0)
    plt.ylabel('normalized vars')
    plt.xlabel('distance from center')
    plt.title('cloud slice ('+str(hgt[kk])+' m, '+str(num_slice*2)+' slices)')
    axes = plt.gca()
    axes.set_xlim([-12, 12])
    axes.set_xticks([-12, -8, -4, 0, 4, 8, 12])
    axes.set_xticklabels(['-1.5 L', '-L', '-0.5 L', 'C', '0.5 L', 'L', '1.5 L'])
    plt.savefig(figurepath + 'cloud_slice_at' + str(int(round(hgt[kk],0))) + '.png', dpi = 150, bbox_inches = 'tight')
    plt.show()

send_email('Plotting cloud slices complete.', 'plot_cld_composite.py --> complete!', attachments = [figurepath + 'cloud_slice_at' + str(int(round(hgt[kk],0))) + '.png' for kk in kk_list])

