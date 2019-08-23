import os
import logging
#from configparser import ConfigParser
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import iris
from utils import label_clds, cloud_center, cloud_edge, cloud_radius
from thermodynamics import Epsilon, Cp_da, Cp_v, Cp_lw, Lv, Pref, Rs_da
import json
import pickle
from scipy.interpolate import interp1d
import math

##############################################################
############   To read the model result   ####################
############   and pick up the sounding   ####################
############   at some specific points    ####################
##############################################################

ntimes = 6
nx  = 330
ny  = 330
nz  = 98
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

### result path
if len(pickup_ud)==3:
    resultpath = "/home/users/jfgu/Analysis/PLUMEDECOM/results/MONC/TODD/Decom/Updraft"+str(len(pickup_ud))+"/nothres/m0400_g330/cld_slice/new_buoy/1000percent/highli_thres/"+nameoption[nanum]+"/test_ud_"+str(pickup_ud[0])+"_"+str(pickup_ud[1])+"_"+str(pickup_ud[2])+"/"
elif len(pickup_ud)==2:
    resultpath = "/home/users/jfgu/Analysis/PLUMEDECOM/results/MONC/TODD/Decom/Updraft"+str(len(pickup_ud))+"/nothres/m0400_g330/cld_slice/new_buoy/1000percent/highli_thres/"+nameoption[nanum]+"/test_ud_"+str(pickup_ud[0])+"_"+str(pickup_ud[1])+"/"
elif len(pickup_ud)==1:
    resultpath = "/home/users/jfgu/Analysis/PLUMEDECOM/results/MONC/TODD/Decom/Updraft"+str(len(pickup_ud))+"/nothres/m0400_g330/cld_slice/new_buoy/1000percent/highli_thres/"+nameoption[nanum]+"/test_ud_"+str(pickup_ud[0])+"/"

### path for saving the figures
figurepath = resultpath + "Figs/"

if not os.path.exists(figurepath):
    os.makedirs(figurepath)

### allocate
ww_prime = np.zeros((ntimes, nz+1, ny, nx))
qt_prime = np.zeros((ntimes, nz+1, ny, nx))

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
ql                 = data['ql']
qcld               = data['qcld']
qt                 = data['qt']
hgt                = data['hgt']
#qt_prime           = data['qt_prime']

for num in range(0, qt.shape[0]):
    for kk in range(0, qt.shape[1]):
        ww_prime[num,kk,:,:] = WW[num,kk,:,:] - np.mean(WW[num,kk,:,:])
        qt_prime[num,kk,:,:] = qt[num,kk,:,:] - np.mean(qt[num,kk,:,:])


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

ww_p_ud_stat         = pickle.load(open(resultpath+'/ww_p_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

qv_ud_stat         = pickle.load(open(resultpath+'/qv_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

ql_ud_stat         = pickle.load(open(resultpath+'/ql_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

qcld_ud_stat       = pickle.load(open(resultpath+'/qcld_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

qt_ud_stat         = pickle.load(open(resultpath+'/qt_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

qt_p_ud_stat         = pickle.load(open(resultpath+'/qt_p_ud_udthres_'+str(ud_thres)+'_'+str(pickup_ud[0])+'.pkl','rb'))

## plot
kk_list = [12, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70]
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
    ### set slices list at each height to be empty
    thetav_p_merge = []
    thetal_p_merge = []
    theta_p_merge  = []
    ww_merge       = []
    qt_merge       = []
    ql_merge       = []
    qcld_merge     = []

    ### number of slices for composite
    num_slice = 0
    for num in range(0, thetav.shape[0]):
        for nn in range(0, len(ud_in_cld_stat[num][kk][ndd])):
            if len(ud_in_cld_stat[num][kk][ndd][nn])>1 and len(ud_gcenter_stat[num][kk][ndd][nn])>0 :

            ## find out maximum value of variables within the plume
                thetav_p_max  = np.max(np.abs(thetav_p_ud_stat[num][kk][ndd][nn]))
                thetal_p_max  = np.max(np.abs(thetal_p_ud_stat[num][kk][ndd][nn]))
                theta_p_max   = np.max(np.abs(theta_p_ud_stat[num][kk][ndd][nn]))
                ww_p_max      = np.max(np.abs(ww_p_ud_stat[num][kk][ndd][nn]))
                qt_p_max      = np.max(np.abs(qt_p_ud_stat[num][kk][ndd][nn]))
                ql_max        = np.max(np.abs(ql_ud_stat[num][kk][ndd][nn]))
                qcld_max      = np.max(np.abs(qcld_ud_stat[num][kk][ndd][nn]))

            ## weighted plume center
                xcenter = ud_gcenter_stat[num][kk][ndd][nn][1]
                ycenter = ud_gcenter_stat[num][kk][ndd][nn][0]

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
                    xp = ud_pos_stat[num][kk][ndd][nn][1]
                    yp = ud_pos_stat[num][kk][ndd][nn][0]

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
                            qt_orig_ys       = qt_prime[num,kk,ys_orig_pts,cx]/qt_p_max
                            ql_orig_ys       = ql[num,kk,ys_orig_pts,cx]/ql_max
                            qcld_orig_ys     = qcld[num,kk,ys_orig_pts,cx]/qcld_max
                            thetav_p_orig_xs = thetav_prime[num,kk,cy,xs_orig_pts]/thetav_p_max
                            thetal_p_orig_xs = thetal_prime[num,kk,cy,xs_orig_pts]/thetal_p_max
                            theta_p_orig_xs  = theta_prime[num,kk,cy,xs_orig_pts]/theta_p_max
                            ww_orig_xs       = ww_prime[num,kk,cy,xs_orig_pts]/ww_p_max
                            qt_orig_xs       = qt_prime[num,kk,cy,xs_orig_pts]/qt_p_max
                            ql_orig_xs       = ql[num,kk,cy,xs_orig_pts]/ql_max
                            qcld_orig_xs     = qcld[num,kk,cy,xs_orig_pts]/qcld_max

                            ## interpolation
                            func_thetav_p_ys = interp1d(ys_pts, thetav_p_orig_ys)
                            func_thetal_p_ys = interp1d(ys_pts, thetal_p_orig_ys)
                            func_theta_p_ys  = interp1d(ys_pts, theta_p_orig_ys)
                            func_ww_ys       = interp1d(ys_pts, ww_orig_ys)
                            func_qt_ys       = interp1d(ys_pts, qt_orig_ys)
                            func_ql_ys       = interp1d(ys_pts, ql_orig_ys)
                            func_qcld_ys     = interp1d(ys_pts, qcld_orig_ys)
                            func_thetav_p_xs = interp1d(xs_pts, thetav_p_orig_xs)
                            func_thetal_p_xs = interp1d(xs_pts, thetal_p_orig_xs)
                            func_theta_p_xs  = interp1d(xs_pts, theta_p_orig_xs)
                            func_ww_xs       = interp1d(xs_pts, ww_orig_xs)
                            func_qt_xs       = interp1d(xs_pts, qt_orig_xs)
                            func_ql_xs       = interp1d(xs_pts, ql_orig_xs)
                            func_qcld_xs     = interp1d(xs_pts, qcld_orig_xs)

                            thetav_p_ys      = func_thetav_p_ys(ys_interp_pts)
                            thetal_p_ys      = func_thetal_p_ys(ys_interp_pts)
                            theta_p_ys       = func_theta_p_ys(ys_interp_pts)
                            ww_ys            = func_ww_ys(ys_interp_pts)
                            qt_ys            = func_qt_ys(ys_interp_pts)
                            ql_ys            = func_ql_ys(ys_interp_pts)
                            qcld_ys          = func_qcld_ys(ys_interp_pts)
                            thetav_p_xs      = func_thetav_p_xs(xs_interp_pts)
                            thetal_p_xs      = func_thetal_p_xs(xs_interp_pts)
                            theta_p_xs       = func_theta_p_xs(xs_interp_pts)
                            ww_xs            = func_ww_xs(xs_interp_pts)
                            qt_xs            = func_qt_xs(xs_interp_pts)
                            ql_xs            = func_ql_xs(xs_interp_pts)
                            qcld_xs          = func_qcld_xs(xs_interp_pts)

                            thetav_p_merge.append(thetav_p_ys)
                            thetav_p_merge.append(thetav_p_xs)
                            thetal_p_merge.append(thetal_p_ys)
                            thetal_p_merge.append(thetal_p_xs)
                            theta_p_merge.append(theta_p_ys)
                            theta_p_merge.append(theta_p_xs)
                            ww_merge.append(ww_ys)
                            ww_merge.append(ww_xs)
                            qt_merge.append(qt_ys)
                            qt_merge.append(qt_xs)
                            ql_merge.append(ql_ys)
                            ql_merge.append(ql_xs)
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
    qt_array       = np.asarray(qt_merge)
    #ql_array       = np.asarray(ql_merge)
    qcld_array     = np.asarray(qcld_merge)

    ## ql might have nan values
    ## remove the rows that have nan values
    #mask           = np.any(np.isnan(ql_array), axis=1)
    #ql_array       = ql_array[~mask]
    mask           = np.any(np.isnan(qcld_array), axis=1)
    qcld_array     = qcld_array[~mask]

    ## composite slices
    thetav_p_comp  = np.mean(thetav_p_array, axis=0)
    thetal_p_comp  = np.mean(thetal_p_array, axis=0)
    theta_p_comp   = np.mean(theta_p_array, axis=0)
    ww_comp        = np.mean(ww_array, axis=0)
    qt_comp        = np.mean(qt_array, axis=0)
    #ql_comp        = np.mean(ql_array, axis=0)
    qcld_comp      = np.mean(qcld_array, axis=0)

    fig = plt.figure(figsize=[sx,sy])
    plt.plot(xplt, thetav_p_comp, '-b', label='$\\theta_v^\prime$')
    plt.plot(xplt, ww_comp, '-r', label='$w$')
    plt.plot(xplt, qt_comp, '-g', label='$q_t$')
    plt.plot(xplt, np.zeros((len(xplt))), '--k')
    plt.legend(loc=0)
    plt.ylabel('normalized vars')
    plt.xlabel('distance from center')
    plt.title('cloud slice ('+str(hgt[num-1,kk])+' m, '+str(num_slice*2)+' slices)')
    axes = plt.gca()
    axes.set_xlim([-12, 12])
    axes.set_xticks([-12, -8, -4, 0, 4, 8, 12])
    axes.set_xticklabels(['-1.5 L', '-L', '-0.5 L', 'C', '0.5 L', 'L', '1.5 L'])

plt.show()



