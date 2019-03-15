import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate, integrate
from netCDF4 import Dataset
from multiprocessing import Pool

################################################################################
#                                                                              #
# Boundary Layer diagnostics                                                   #
#                                                                              #
################################################################################

def zi(theta_v, z):
    """
    Function to calculate the mixed layer depth as a function of the virtual
    potential temeprature.
    
    Input:
    theta_v = virtual potential temperature (K), theta_v = theta_v(z, y, x)
    z       = the heights (m) at which theta_v are computed
    
    Output:
    zi      = the mixed layer depth
    
    Assumes that the mixed layer depth is defined as the level to which a
    surface based parcel of air (i.e. theta_v_parcel = theta_v(z = 0)) must be 
    lifted in order for the parcel to become neutrally buoyant (i.e. 
    theta_v_parcel = theta_v_env(z = zi).
    """
    
    # Use a surface based parcel
    theta_v_parcel = theta_v[0,:,:]
    # Find the height of the level of neutral buoyancy
    z_i = np.zeros_like(theta_v_parcel)
    for j in xrange(theta_v.shape[1]):
        for i in xrange(theta_v.shape[2]):
            ik = np.max(np.where(np.abs(theta_v_parcel[j,i] - theta_v[2:,j,i]) == np.min(np.abs(theta_v_parcel[j,i] - theta_v[2:,j,i])))[0])
            if theta_v[ik,j,i] > theta_v_parcel[j,i]:
                # Means the nearest point is higher height than the boundary layer depth
                # Linearly interpolate to get z_i
                m = 1
                while theta_v[ik-m,j,i] == theta_v[ik,j,i]:
                    m += 1
                n = m - 1
                z_i[j,i] = (z[ik-n] - z[ik-m])*(theta_v_parcel[j,i] - theta_v[ik-m,j,i])/(theta_v[ik-n,j,i] - theta_v[ik-m,j,i]) + z[ik-m]
            elif theta_v[ik,j,i] < theta_v_parcel[j,i]:
                # Means the nearest point is at a lower height than the boundary layer depth
                m = 1
                while theta_v[ik+m,j,i] == theta_v[ik,j,i]:
                    m += 1
                n = m - 1
                z_i[j,i] = (z[ik+m] - z[ik+n])*(theta_v_parcel[j,i] - theta_v[ik+n,j,i])/(theta_v[ik+m,j,i] - theta_v[ik+n,j,i]) + z[ik+n]
            else:
                z_i[j,i] = z[ik]
            #z_i[j, i] = interpolate.interp1d(theta_v[t2:, j, i], z[2:], fill_value = "extrapolate")(theta_v_parcel[j,i])
            #if z_i[j, i] < 0.:
            #    z_i[j, i] = 0.
            
            #k = 2
            #while theta_v_parcel[j,i] > theta_v[k,j,i]:
            #    k += 1
            #z_i[j,i] = z[k]
    return z_i

def getML_mean(var, z, zi, axis = 0):
    """
    Find and return the mean of 'var' in height over the mixed layer 'zi'.
    ----------------------------------------------------------------------------
    INPUT:
    var  = variable that is some function of height along the given axis
    var  = f([t], z, y, x)
    z    = height observations
    zi   = height of the mixed layer determined from elsewhere
    zi   = f([t], y, x)
    axis = the axis in var that contains the height variation
    
    i.e. if var = f(t, z, y, x) axis should be provided as axis = 1
    N.B. if just a single profile, this should be able to pick out the correct 
    ML mean value
    OUTPUT:
    var_ML = the mean of the given variable over the mixed layer depth, zi
    ----------------------------------------------------------------------------
    """
    from scipy import interpolate, integrate
    if len(var.shape) <= 2:
        # This is a single profile, either f(t, z) or f(z)
        var_f  = interpolate.interp1d(z, var, axis = axis)
        my_z   = np.arange(np.min(z), zi)
        var_ML = integrate.trapz(var_f(my_z), my_z, axis = axis)/zi
    elif len(var.shape) == 3:
        # This is f(z, y, x)
        var_ML = np.zeros_like(zi)
        for j in xrange(var.shape[1]):
            for i in xrange(var.shape[2]):
                var_f       = interpolate.interp1d(z, var[:,j,i])
                my_z        = np.arange(np.min(z), zi[j,i])
                var_ML[j,i] = integrate.trapz(var_f(my_z), my_z)/zi[j,i]
    elif len(var.shape) == 4:
        # This is f(t, z, y, x)
        var_ML = np.zeros_like(zi)
        for it in xrange(var.shape[0]):
            for j in xrange(var.shape[2]):
                for i in xrange(var.shape[3]):
                    var_f          = interpolate.interp1d(z, var[it,:,j,i])
                    my_z           = np.arange(np.min(z), zi[it,j,i])
                    var_ML[it,j,i] = integrate.trapz(var_f(my_z), my_z)/zi[it,j,i]
    return var_ML

def lcl(temp, q, pres, z):
    """
    lcl calculates the surface-based lifting condensation level as a function of
    the temperature (temp), specific humidity (q), and pressure (pres).
    These input are used to find the temperature and dew point which are used in
    Bolton (1980)'s equation for the LCL.
    
    Input:
    temp = air temperature (K), temp = temp[z=0, y, x]
    q    = specific humidity (kg/kg), q = q[z=0, y, x]
    pres = pressure (Pa), pres = pres[z, y, x]
    z    = heights (m)
    
    Output:
    lcl  = the lifting condensation level (m)
    """
    # Import some constants
    from SkewT_archer import Lv, Rv, Rd, EPS, cpd
    T0 = 273.15 # freezing point of water
    e0 = 611.2  # vapour pressure of water vapour at T0
    
    # calculate the dew point temperature
    e = q*pres[0,:,:]/(EPS - q*EPS + q)
    dewpoint = (1./T0 - (Rv/Lv)*np.log(e/e0))**(-1.)
    
    term_A = 1./(dewpoint - 56)
    term_B = np.log(temp/dewpoint)/800.
    T_LCL = 1/(term_A + term_B) + 56
    
    p_LCL = pres[0,:,:]*(T_LCL/temp)**(cpd/Rd)
    z_LCL = np.zeros_like(p_LCL)
    for j in xrange(pres.shape[1]):
        for i in xrange(pres.shape[2]):
            z_LCL[j,i] = interpolate.interp1d(pres[:,j,i], z)(p_LCL[j,i])
    
    return z_LCL

def find_h(theta_v, u, v, z):
    """
    Iteratively solves equation 2 from Vogelezang and Holtslag (1996) to
    estimate the boundary layer depth.
    
    This method uses a modified Richardson number approach and vertical profiles
    of virtual potential temperature (theta_v), u- and v- wind components to
    estimate the boundary layer height.
    ----------------------------------------------------------------------------
    INPUT:
    theta_v = f(z,y,x)
    u       = f(z,y,x) -> regridded to theta points
    v       = f(z,y,x) -> regridded to theta points
    z       = f(z)     -> heights in theta levels
    
    OUTPUT:
    z_i     = f(y,x)
    """
    
    # Requires numpy arrays
    import numpy as np
    # Requires interpolation routines from scipy
    from scipy import interpolate
    
    # Constants
    g    = 9.81
    Ri_c = 0.25 # Critical Richardson Number distinguishing between laminar and turbulent flow
    
    # Need a first guess for the height of the boundary layer.
    z_s = 0.1*z # approximation to the surface layer depth
    
    Ri_g = np.zeros((len(z),theta_v.shape[1], theta_v.shape[2]))
    
    theta_v = interpolate.interp1d(x = z, y = theta_v, fill_value = 'extrapolate', axis = 0)
    u       = interpolate.interp1d(x = z, y = u, fill_value = 'extrapolate', axis = 0)
    v       = interpolate.interp1d(x = z, y = v, fill_value = 'extrapolate', axis = 0)
    
    for k in xrange(1, len(z)):
        Ri_g[k,:,:] = (g/theta_v(z_s[k]))*(theta_v(z[k]) - theta_v(z_s[k]))*(z[k] - z_s[k])/((u(z[k]) - u(z_s[k]))**2. + (v(z[k]) - v(z_s[k]))**2.)
    
    z_i = np.zeros((Ri_g.shape[1], Ri_g.shape[2]))
    for j in xrange(Ri_g.shape[1]):
        for i in xrange(Ri_g.shape[2]):
            # Assume that Ri_g is monotonically increasing in the boundary layer
            k = 0
            while Ri_g[k,j,i] <= Ri_c:
                k += 1
            z_i[j,i] = (z[k] - z[k-1])*(Ri_c - Ri_g[k-1,j,i])/(Ri_g[k,j,i] - Ri_g[k-1,j,i]) + z[k-1]
        
    return z_i

def get_TKE(U, V, W, rho = 1., start = 0, end = None):
    """
    Function to compute the vertical profile of turbulent kinetic energy from 
    given 4D U, V, and W fields.
    
    First, calculates the perturbations from horizontal mean profiles of U, V, 
    and W. Then, calculates the turbulent kinetic energy, and then the 
    horizontal mean so that there is a time series of profiles.
    
    U, V, and W must all be on the same grid.
    """
    U_mean = horizontalMean(U)
    V_mean = horizontalMean(V)
    W_mean = horizontalMean(W)
    U_p = np.zeros_like(U)
    V_p = np.zeros_like(V)
    W_p = np.zeros_like(W)
    for t in xrange(U.shape[0]):
        for k in xrange(U.shape[1]):
            U_p[t,k,:,:] = U[t,k,:,:] - U_mean[t,k]
            V_p[t,k,:,:] = V[t,k,:,:] - V_mean[t,k]
            W_p[t,k,:,:] = W[t,k,:,:] - W_mean[t,k]
    if type(rho) != float:
        TKE_profile = 0.5*horizontalMean(rho)*(horizontalMean(U_p**2)+ horizontalMean(V_p**2) + horizontalMean(W_p**2))
    else:
        TKE_profile = 0.5*rho*(horizontalMean(U_p**2)+ horizontalMean(V_p**2) + horizontalMean(W_p**2))
    
    return TKE_profile

################################################################################
#                                                                              #
# Cloud diagnostics                                                            #
#                                                                              #
################################################################################

def get_CC(m_cl, start = 0, end = None, threshold = 0.):
    """
    Function to compute the vertical profile of cloud cover from given 4D m_cl 
    field.
    
    First, finds the number of cloudy pixels, i.e. those that contain non-zero 
    m_cl. This is a binary, 0 = no cloud, 1 = cloud. Then calculates the cloud 
    cover on each vertical level as the mean of the binary field at each level. 
    Finally, does a time mean over a given period. Returns a vertical profile of
    this cloud cover estimate.
    """
    cloud_mask = m_cl[:]
    cloud_mask[m_cl > threshold] = 1.
    cloud_mask[m_cl <= threshold] = 0.
    cc_profile = np.mean(horizontalMean(cloud_mask)[start:end,:], axis = 0)
    
    return cc_profile

def get_CTZ(mc, z, threshold = 1e-16):
    """
    takes a 3- or 4-dimensional array of total cloud water mixing ratio (i.e. mcl
    + mcf) and a 1-dimensional array of heights.
    Then finds the highest level at which the total cloud water mixing ratio is
    non-zero. Options exist to interpolate to this level using scipy functions, 
    which might be useful for deep convection, but is more expensive.
    ----------------------------------------------------------------------------
    Input:
    mc = a 3- or 4-dimensional array of total cloud water mixing ratio 
    mc = mcl + mcf = f([t],z,y,x)
    z = heights corresponding to the z coordinate of mc
    z = f(z)
    threshold = some value below which total water mixing ratio is not considered
        to be cloud
    
    Output:
    Z_top = 2- or 3-dimensional array of cloud top heights, if there is no cloud
    in the column, a np.nan should fill in.
    Z_top = f([t],y,x)
    ----------------------------------------------------------------------------
    """
    if len(mc.shape) == 3:
        Z = np.repeat(z, mc.shape[1]*mc.shape[2], axis = 0).reshape(mc.shape)
        Z_top = np.nanmax(np.where((mc >= threshold), Z, np.nan), axis = 0)
    else:
        Z = np.repeat(z, mc.shape[2]*mc.shape[3], axis = 0).reshape(mc.shape[1:])
        Z_top = np.zeros((mc.shape[0], mc.shape[2], mc.shape[3]))
        for it in xrange(mc.shape[0]):
            Z_top[it,:,:] = np.nanmax(np.where((mc[it,:,:,:] >= threshold), Z, np.nan), axis = 0)
    
    return Z_top

def get_theta_w(temperature, q, pressure, t_units = 'K', q_units = 'kg/kg', p_units = 'hPa', ndp = 20):
    """
    Calculates the wet bulb potential temperature.
    First, the parcel must be lifted dry adiabatically to the LCL, then the 
    parcel must be brought pseudoadiabatically to the reference pressure 
    (1000 hPa) to arrive at the wet bulb potential temperature.
    
    If the saturated wet-bulb potential temperature is sought, replace q with 
    getQ(temperature, RH, pressure) from SkewT_archer and use RH = 100%.
    """
    # Do unit conversions
    if t_units == 'C':
        # convert to K
        temperature += 273.15
    
    if q_units == 'g/kg':
        # convert to kg/kg
        q /= 1000.
    
    if p_units == 'Pa':
        # Convert to hPa
        pressure /= 100.
    
    from SkewT_archer import getDew, cpd, Rd, p0, g, getGM
    # First, find the LCL temperature (requires dew point from getDew)
    dewpoint = getDew(q*1., pressure*1., q_units = 'kg/kg', p_units = 'hPa')
    A = 1./(dewpoint - 56.)
    B = np.log(temperature/dewpoint)/800.
    LCL_temperature  = (1./(A + B)) + 56.
    
    # Second, find the LCL pressure by inverting Poisson's relation
    LCL_pressure = pressure/((temperature/LCL_temperature)**(cpd/Rd)) # hPa
    
    # Third, bring the parcel pseudoadiabatically to the reference pressure
    # using getGM, an integration with height and an estimate for dz.
    """
    # cheap and nasty one step integration
    
    dp = p0 - LCL_pressure

    rho = 0.5*(LCL_pressure + p0)*100./(Rd*LCL_temperature)
    dz = - dp*100./(rho*g)
    
    theta_w = LCL_temperature + getGM(LCL_temperature*1., 0.5*(LCL_pressure+p0), t_units = 'K', p_units = 'hPa')*dz
    """
    
    # multi-step integration
    dp = (p0 - LCL_pressure)/ndp
    p1 = LCL_pressure*1.
    theta_w = LCL_temperature*1.
    for i in xrange(ndp):
        # integrate down until we are at p0
        rho = (p1 + dp/2.)*100./(Rd*theta_w)
        dz = - dp*100./(rho*g)
        theta_w += getGM(theta_w*1., (p1 + dp/2.), t_units = 'K', p_units = 'hPa')*dz
        p1 += dp
    
    return theta_w

################################################################################
#                                                                              #
# Regridding, interpolation and spatial smoothing routines                     #
#                                                                              #
################################################################################

def cartesian2polar(x, y, x_o = 0., y_o = 0.):
    """
    Convert from cartesian to polar coordinates
    x, y = the cartesian coordinates
    x_o, y_o = the coordinates of the origin
    """
    dx = x - x_o
    dy = y - y_o
    r = np.sqrt(dx**2 + dy**2)
    
    """
    Four cases:
    1. if dx and dy are positive
    2. if dx is positive and dy is negative
    3. if dx and dy are negative
    4. if dx is negative and dy is positive
    """
    theta1 = np.arcsin(dx/r)
    theta2 = np.arcsin(-dy/r) + np.pi/2.
    theta3 = np.arcsin(-dx/r) + np.pi
    theta4 = np.arcsin(dy/r)  + 3.*np.pi/2.
    
    theta1[dx < 0.] = 0.
    theta1[dy < 0.] = 0.
    
    theta2[dx < 0.] = 0.
    theta2[dy > 0.] = 0.

    theta3[dx > 0.] = 0.
    theta3[dy > 0.] = 0.
    
    theta4[dx > 0.] = 0.
    theta4[dy < 0.] = 0.
    
    # Add together and convert from radians to degrees
    theta = (theta1 + theta2 + theta3 + theta4)*180./np.pi
    
    return r, theta

def polar2cartesian(r, theta, x_o, y_o):
    """
    Converts from polar coordinates to cartesian coordinates
    r = distance from the origin
    theta = angle from x = 0 in degrees
    x_o, y_o = the origin in cartesian space
    """
    # Convert from degrees to radians
    theta = theta*np.pi/180.
    x = r*np.sin(theta)
    y = r*np.cos(theta)
    
    # return the origin to 0,0
    x += x_o
    y += y_o
    
    return x, y

def bilinear_interpolation(x_in, y_in, z_in0, x_out, y_out, kind = 0, d = 2000.0, p = 0.5, operation = np.nanmean):
    """
    Performs bilinear interpolation on a given 2D array, z_in.
    ----------------------------------------------------------------------------
    INPUT:
    x_in = 2D array of x-coordinates, input
    y_in = 2D array of y-coordinates, input
    z_in = 3D array of data points, input z_in[nz, ny, nx]

    x_out = 1D array of x-coordinates to be interpolated onto
    y_out = 1D array of y-coordinates to be interpolated onto
    
    kind = option for how to interpolate
    option 0: nearest neighbor
    option 1: inverse distance weighted within radius d to power p
    option 2: bilinear interpolation
    option 3: mean within radius d
    
    OUTPUT:
    z_out = 2D array, output z_out[nz, nr] where nr is the length of x_out
    e.g. nz could be the number of height levels, and nr could be the number of 
    downwind points
    
    interpolation scheme finds the nearest points and uses scipy's bisplrep to 
    do the interpolation.
    """
    import numpy as np
    # make use of the bi-periodic boundary conditions to not truncate at boundaries
    dx = x_in[0,1] - x_in[0,0]
    dy = y_in[1,0] - y_in[0,0]
    x_in_p = np.concatenate((x_in-np.max(x_in)-dx, x_in, x_in+np.max(x_in)+dx), axis = 1)
    x_in_p = np.concatenate((x_in_p, x_in_p, x_in_p), axis = 0)
    y_in_p = np.concatenate((y_in-np.max(y_in)-dy, y_in, y_in+np.max(y_in)+dy), axis = 0)
    y_in_p = np.concatenate((y_in_p, y_in_p, y_in_p), axis = 1)
    # Initialise our output array
    z_out = np.zeros((z_in0.shape[0], len(x_out)))
    
    # For each point to be interpolated onto
    for i in xrange(len(x_out)):
        #print 'bilinearInterpolation [' + str(i) +']'
        # Find the distance to the input data coordinates
        r = np.sqrt((x_in_p - x_out[i])**2 + (y_in_p - y_out[i])**2)
        
        if kind == 0:
            z_in = z_in0*1.
            # Nearest neighbor approach
            iy, ix = np.where(r == np.min(r))
            z_out[:,i] = z_in[:,iy[0]%z_in.shape[1],ix[0]%z_in.shape[2]]
        elif kind == 1:
            z_in = z_in0*1.
            # idw approach
            # Use all points within 'd' m of x_out, y_out for weighted mean
            iy, ix = np.where(r < d)
            w = 0
            
            for j in xrange(len(iy)):
                w += (1./r[iy[j], ix[j]]**p)
                z_out[:,i] += (1./r[iy[j], ix[j]]**p)*z_in[:,iy[j]%z_in.shape[1], ix[j]%z_in.shape[2]]
            z_out[:,i] /= w
        elif kind == 2:
            # Find the nearest point
            iy0, ix0 = np.where(r == np.min(r))
            Dx = x_out[i] - x_in_p[iy, ix] # if +ve, input point is to the right
            Dy = y_out[i] - y_in_p[iy, ix] # if +ve, input point is up
            iy = iy0[0]%z_in0.shape[1]
            ix = ix0[0]%z_in0.shape[2]
            iy0 = iy0[0]
            ix0 = ix0[0]
            if Dx < 0:
                # nearest point is to the right, move one to the left
                ix -= 1
            
            if Dy < 0:
                # nearest point is above, move one down
                iy -= 1
            
            z_in = np.zeros((z_in0.shape[0], 3, 3))
            for J in xrange(-1,2):
                for I in xrange(-1,2):
                    z_in[:,J,I] = z_in0[:,(iy+J)%z_in0.shape[1],(ix+I)%z_in0.shape[2]]
            # bilinearly interpolate
            z_out_up = (z_in[:,1,1] - z_in[:,1,0])*(x_out[i] - x_in_p[iy0,ix0-1])/(x_in_p[iy0,ix0] - x_in_p[iy0,ix0-1]) + z_in[:,1,0]
            z_out_down = (z_in[:,0,1] - z_in[:,0,0])*(x_out[i] - x_in_p[iy0-1,ix0-1])/(x_in_p[iy0-1,ix0] - x_in_p[iy0-1,ix0-1]) + z_in[:,0,0]
            y_in_up = (y_in_p[iy0,ix0] - y_in_p[iy0,ix0-1])*(x_out[i] - x_in_p[iy0,ix0-1])/(x_in_p[iy0,ix0] - x_in_p[iy0,ix0-1]) + y_in_p[iy0,ix0-1]
            y_in_down = (y_in_p[iy0-1,ix0] - y_in_p[iy0-1,ix0-1])*(x_out[i] - x_in_p[iy0-1,ix0-1])/(x_in_p[iy0-1,ix0] - x_in_p[iy0-1,ix0-1]) + y_in_p[iy0-1,ix0-1]
            z_out[:,i] = (z_out_up - z_out_down)*(y_out[i] - y_in_down)/(y_in_up - y_in_down) + z_out_down
        elif kind == 3:
            z_in = z_in0*1.
            iy, ix = np.where(r < d)
            z_out[:,i] = operation(z_in[:,iy%z_in.shape[1],ix%z_in.shape[2]], axis = 1)
    
    return z_out

def get_cs_coords(x_c, y_c, direction, x, y, h = 100., isPeriodic = False, max_r = 80000.):
    """
    Uses an input coordinate (x_c, y_c) and a wind direction (direction) to 
    calculate the coordinates of a line that goes through the input coordinate
    and would be parallel to the wind direction. These coordinates are given
    at a horizontal resolution (h) in metres.
    ---------------------------------------------------------------------------
    x_c: the x-coordinate of the input point (e.g. island centre)
    y_c: the y-coordinate of the input point
    direction: the direction of the cross section
    x: the 2D array of X coordinates
    y: the 2D array of Y coordinates
    h: the resolution at which to take points along the cross section
    """
    if isPeriodic:
        # make use of the bi-periodic boundary conditions to not truncate at boundaries
        x = np.concatenate((x-np.max(x), x, x+np.max(x)), axis = 1)
        x = np.concatenate((x, x, x), axis = 0)
        y = np.concatenate((y-np.max(y), y, y+np.max(y)), axis = 0)
        y = np.concatenate((y, y, y), axis = 1)
    
    dx = (h*np.sin(np.pi*direction/180.0))
    dy = (h*np.cos(np.pi*direction/180.0))

    x_cs = [x_c]
    y_cs = [y_c]

    # Populate my list of x and y coordinates along the cross section
    while (x_cs[-1] < np.max(x))*(y_cs[-1] < np.max(y))*(round(np.sqrt((x_cs[-1] - x_c)**2 + (y_cs[-1] - y_c)**2),2) < max_r):
        # checks that the x and y cross section coordinate are within the range
        # on the right hand side, and top of the domain, also checks that the 
        # cross section isn't longer than max_r.
        x_cs.append(x_cs[-1] + dx)
        y_cs.append(y_cs[-1] + dy)

    x_cs = x_cs[::-1]
    y_cs = y_cs[::-1]

    while (x_cs[-1] > np.min(x))*(y_cs[-1] > np.min(y))*(round(np.sqrt((x_cs[-1] - x_c)**2 + (y_cs[-1] - y_c)**2),2) < max_r):
        # checks that the x and y cross section coordinate are within the range 
        # of the domain on the left hand side and bottom of the domain, also
        # also checks that the cross section isn't longer than max_r
        x_cs.append(x_cs[-1] - dx)
        y_cs.append(y_cs[-1] - dy)

    # double check that all the coordinates along the cross section are within the domain
    i = 0
    while i < len(x_cs):
        if (round(np.sqrt((x_cs[i] - x_c)**2 + (y_cs[i] - y_c)**2),2) > max_r) or (x_cs[i] > np.max(x)) or (x_cs[i] < np.min(x)) or (y_cs[i] > np.max(y)) or (y_cs[i] < np.min(y)):
            del x_cs[i]
            del y_cs[i]
        else:
            i += 1

    # store to arrays
    x_cs = np.array(x_cs)
    y_cs = np.array(y_cs)
    
    return x_cs, y_cs

def regrid(target_grid, current_grid, current_data_key):
    """
    Function to attempt to regrid an already read-in netCDF B (current_grid) to 
    the same grid as already read-in netCDF A (target_grid).
    Should only work if there is some overlap between the range of the functions
    otherwise the interpolation won't make sense.
    
    Regrids in all four dimensions, although going from e.g. 1 minute data to
    10 or 15 minute data shouldn't really require regridding. But going from
    10 minute to 1 minute data definitely requires regridding and the linear 
    interpolation approach used here is probably not appropriate.
    """
    # Read the data that contains the target grid
    target_x_key     = [key for key in target_grid.variables.keys() if ('longitude' in key) and ('_l' not in key)][0]
    target_x = target_grid.variables[target_x_key][:]
    target_y_key     = [key for key in target_grid.variables.keys() if ('latitude' in key) and ('_l' not in key)][0]
    target_y = target_grid.variables[target_y_key][:]
    target_z_key     = [key for key in target_grid.variables.keys() if ('zsea' in key) and ('bounds' not in key)][0]
    target_z = target_grid.variables[target_z_key][:]
    target_t_key     = [key for key in target_grid.variables.keys() if 'min' in key][0]
    target_t = target_grid.variables[target_t_key][:]
    
    # Read the netCDF that needs to be transformed
    current_x_key     = [key for key in current_grid.variables.keys() if ('longitude' in key) and ('_l' not in key)][0]
    current_x = current_grid.variables[current_x_key][:]
    current_y_key     = [key for key in current_grid.variables.keys() if ('latitude' in key) and ('_l' not in key)][0]
    current_y = current_grid.variables[current_y_key][:]
    current_z_key     = [key for key in current_grid.variables.keys() if ('zsea' in key) and ('bounds' not in key)][0]
    current_z = current_grid.variables[current_z_key][:]
    current_t_key     = [key for key in current_grid.variables.keys() if 'min' in key][0]
    current_t = current_grid.variables[current_t_key][:]
    current_data = current_grid.variables[current_data_key][:]
     
    ### Procedure ###
    # Will usually only need to interpolate in z, sometimes in x or y, probably
    # never required to interpolate in more than two dimensions
    # 1. interpolate in time
    if target_t_key != current_t_key:
        print 'Regridding in time...'
        # Interpolate in time
        current_data = interpolate.interp1d(current_t, current_data, axis = 0, fill_value = 'extrapolate')(target_t)
        print 'Complete.'
    # 2. interpolate in z
    if target_z_key != current_z_key:
        print 'Regridding in height...'
        # Interpolate in z
        current_data = interpolate.interp1d(current_z, current_data, axis = 1, fill_value = 'extrapolate')(target_z)
        print 'Complete.'
    # 3. interpolate in y
    if target_y_key != current_y_key:
        print 'Regridding in latitude...'
        # Interpolate in y
        for i in xrange(current_data.shape[3]):
            current_data[:,:,:,i] = interpolate.interp1d(current_y[:,i], current_data[:,:,:,i], axis = 2, fill_value = 'extrapolate')(target_y[:,i])
        print 'Complete.'
    # 4. interpolate in x
    if target_x_key != current_x_key:
        print 'Regridding in longitude...'
        # Interpolate in x
        for j in xrange(current_data.shape[2]):
            current_data[:,:,j,:] = interpolate.interp1d(current_x[j,:], current_data[:,:,j,:], axis = 2, fill_value = 'extrapolate')(target_x[j,:])
        print 'Complete.'
    
    return current_data

def cross_section_x(data, x_rot, x_pos, t_pos = -1):
    """
    data = data[time, z, y, x]
    x_rot = rotated x-coordinate
    x_pos = x_coordinate to do cross section at
    """
    data_xs = np.zeros_like(data[t_pos,:,:,x_pos/100])
    x_target = x_rot - x_pos
    
    my_ixs = []
    for iy in range(w.shape[2]):
        ix = np.where(abs(x_target[iy,:]) == np.min(abs(x_target[iy,:])))[0][0]
        my_ixs.append(ix)
        if x_target[iy, ix] > 0.:
            # interpolate ix and ix-1
            data_xs[:,iy] = (data[t_pos,:,iy,ix] - data[t_pos,:,iy,ix-1])*(x_pos - x_rot[iy,ix-1])/(x_rot[iy,ix] - x_rot[iy,ix-1]) + data[t_pos,:,iy,ix-1]
        elif x_target[iy, ix] > 0:
            # interpolate ix+1 and ix
            data_xs[:,iy] = (data[t_pos,:,iy,ix+1] - data[t_pos,:,iy,ix])*(x_pos - x_rot[iy,ix])/(x_rot[iy,ix+1] - x_rot[iy,ix]) + data[t_pos,:,iy,ix]
        else:
            # no need to interpolate as x_target == 0
            data_xs[:,iy] = data[t_pos,:,iy,ix]
    return data_xs, my_ixs

def ddx(data, X, periodic = True):
    """
    Calculate the horizontal gradient in the x-direction.

    data     = 2, 3, or 4 dimensional array
    data     = f([t],[z],y,x)
    X        = 2D array of x-coordinates
    X        = f(y,x)
    periodic = Boolean, 
               .True. means periodic boundary conditions
               .False. means zero gradient boundary conditions.
    """
    print 'WARNING: gradients are returned at staggered coordinates'
    d_dx = np.zeros_like(data)
    if len(data.shape) == 2:
        d_dx[:,:-1] = (data[:,1:] - data[:,:-1])/(X[:,1:] - X[:,:-1])
        d_dx[:,-1] = (data[:,0] - data[:,-1])/(X[:,0] - X[:,-1])
    elif len(data.shape) == 3:
        d_dx[:,:,:-1] = (data[:,:,1:] - data[:,:,:-1])/(X[:,1:] - X[:,:-1])
        d_dx[:,:,-1] = (data[:,:,0] - data[:,:,-1])/(X[:,0] - X[:,-1])
    elif len(data.shape) == 4:
        d_dx[:,:,:,:-1] = (data[:,:,:,1:] - data[:,:,:,:-1])/(X[:,1:] - X[:,:-1])
        d_dx[:,:,:,-1] = (data[:,:,:,0] - data[:,:,:,-1])/(X[:,0] - X[:,-1])
    
    return d_dx

def ddy(data, Y, periodic = True):
    """
    Calculate the horizontal gradient in the x-direction.

    data     = 2, 3, or 4 dimensional array
    data     = f([t],[z],y,x)
    Y        = 2D array of y-coordinates
    Y        = f(y,x)
    periodic = Boolean, 
               .True. means periodic boundary conditions
               .False. means zero gradient boundary conditions.
    """
    print 'WARNING: gradients are returned at staggered coordinates'
    d_dy = np.zeros_like(data)
    if len(data.shape) == 2:
        d_dy[:-1,:] = (data[1:,:] - data[:-1,:])/(Y[1:,:] - Y[:-1,:])
        d_dy[-1,:] = (data[0,:] - data[-1,:])/(Y[0,:] - Y[-1,:])
    elif len(data.shape) == 3:
        d_dy[:,:-1,:] = (data[:,1:,:] - data[:,:-1,:])/(Y[1:,:] - Y[:-1,:])
        d_dy[:,-1,:] = (data[:,0,:] - data[:,-1,:])/(Y[0,:] - Y[-1,:])
    elif len(data.shape) == 4:
        d_dy[:,:,:-1,:] = (data[:,:,1:,:] - data[:,:,:-1,:])/(Y[1:,:] - Y[:-1,:])
        d_dy[:,:,-1,:] = (data[:,:,0,:] - data[:,:,-1,:])/(Y[0,:] - Y[-1,:])
    
    return d_dy

def circular_smoothing(X, Y, data, r, scale = 1):
    """
    Smooth the data by taking the mean of all points within radius r.
    Requires the following input:
    X, a 2D array f(y,x) of x-coordinates
    Y, a 2D array f(y,x) of y-coordinates
    data a 2, 3, or 4D array f([t],[z],y,x) of data to be smoothed
    r, a single float in the same units as the X and Y coordinates to define a 
    radius with which to smooth the data.
    scale = a parameter to reduce the size of the arrays to process quicker,
        must be an integer greater than or equal to 1.
    
    N.B. a future version of this function may attempt to do some form of
    weighting to make the points closest to the actual point weigh more than
    the points further away.
    """
    x_periodic = np.concatenate((X-np.max(X), X, X+np.max(X)), axis = 1)
    x_periodic = np.concatenate((x_periodic, x_periodic, x_periodic), axis = 0)
    y_periodic = np.concatenate((Y-np.max(Y), Y, Y+np.max(Y)), axis = 0)
    y_periodic = np.concatenate((y_periodic, y_periodic, y_periodic), axis = 1)
    smoothed_data = np.zeros_like(data)
    if len(data.shape) == 2:
        for i in xrange(0, data.shape[0], scale):
            for j in xrange(0, data.shape[1], scale):
                iy, ix = np.where(np.sqrt((X[i,j] - x_periodic)**2 + (Y[i,j] - y_periodic)**2) <= r)
                smoothed_data[i,j] = np.nanmean(data[iy%X.shape[0],ix%X.shape[1]])
        smoothed_data = smoothed_data[::scale,::scale]
    elif len(data.shape) == 3:
        for tz0 in xrange(data.shape[0]):
            for i in xrange(0, data.shape[1], scale):
                for j in xrange(0, data.shape[2], scale):
                    iy, ix = np.where(np.sqrt((X[i,j] - x_periodic)**2 + (Y[i,j] - y_periodic)**2) <= r)
                    smoothed_data[tz0,i,j] = np.nanmean(data[tz0,iy%X.shape[0],ix%X.shape[1]])
        smoothed_data = smoothed_data[:,::scale,::scale]
    elif len(data.shape) == 4:
        for tz0 in xrange(data.shape[0]):
            for tz1 in xrange(data.shape[1]):
                for i in xrange(0, data.shape[2], scale):
                    for j in xrange(0, data.shape[3], scale):
                        iy, ix = np.where(np.sqrt((X[i,j] - x_periodic)**2 + (Y[i,j] - y_periodic)**2) <= r)
                        smoothed_data[tz0,tz1,i,j] = np.nanmean(data[tz0,tz1,iy%X.shape[0],ix%X.shape[1]])
        smoothed_data = smoothed_data[:,:,::scale,::scale]
    
    return smoothed_data

def get_rectangle(x_c, y_c, X, Y, direction, half_width, min_distance, max_distance):#mask, sector, length, my_res = 10.0):
    """
    x_c, y_c are the coordinates of the point of interest (e.g. the centre of the island)
    X = a 2d array of the x-corrdinates
    Y = a 2d array of the y-coordinates
    direction = the direction in which to look (i.e. if the wind is from 90 deg, we look positive in the 270 deg direction)
    half_width = the width of the rectangle divided by two (user inputs in the same units are the coordinate system (i.e. degrees lat/lon or metres)
    min_distance = the minimum distance from the x_c,y_c coordinate
    max_distance = the maximum distance from the x_c,y_c coordinate
    """
    # Convert from the wind direction to the heading
    direction = (direction + 180.)%360.
    # Find the distance of each point from the x_c, y_c coordinate of interest
    R = np.sqrt((x_c - X)**2 + (y_c - Y)**2)
    # Find the bearing of each point from the x_c, y_c coordinate of interest
    T = np.where(((y_c - Y) >= 0), (np.arctan((x_c - X)/(y_c - Y))*180./np.pi-180.)%360., (np.arctan((x_c - X)/(y_c - Y))*180./np.pi)%360.)
    # Find the index for the point with R nearest to min_distance and T nearest to direction
    iy, ix = np.where(np.abs(R - min_distance)*np.abs(T - direction) == np.min(np.abs(R - min_distance)*np.abs(T - direction)))
    iy = iy[0]
    ix = ix[0]
    
    # draw line perpendicular to direction
    if direction <= 90.:
        x_fact = 1.
        y_fact = 1.
    elif direction <= 180.:
        x_fact = 1.
        y_fact = -1.
    elif direction <= 270.:
        x_fact = -1.
        y_fact = -1.
    else:
        x_fact = -1.
        y_fact = 1.
    
    dx1 = half_width*np.cos(direction*np.pi/180.)
    dy1 = half_width*np.sin(direction*np.pi/180.)
    
    left_X = X[iy, ix]-dx1
    left_Y = Y[iy, ix]+dy1
    right_X = X[iy, ix]+dx1
    right_Y = Y[iy, ix]-dy1
    
    ## draw lines at 90deg to our first line
    opposite = abs(left_X - right_X)
    adjacent = abs(left_Y - right_Y)
    hypotenuse = np.sqrt(opposite**2 + adjacent**2)
    angle = np.arcsin(opposite/hypotenuse)*360./(2.*np.pi)
    ## use angle_b to find the coordinates for the other two corners of the rectangle
    length = max_distance - min_distance
    dx = length*np.cos(angle*np.pi/180.)
    dy = length*np.sin(angle*np.pi/180.)
    new_left_X = left_X + dx*x_fact
    new_left_Y = left_Y + dy*y_fact
    new_right_X = right_X + dx*x_fact
    new_right_Y = right_Y + dy*y_fact
    
    ## get the indices enclosed by that rectangle
    #A. left to right
    ma = (right_Y - left_Y)/(right_X - left_X)
    ca = left_Y
    #B. right to new_right
    mb = (new_right_Y - right_Y)/(new_right_X - right_X)
    cb = right_Y
    #C. new_right to new_left
    mc = (new_left_Y - new_right_Y)/(new_left_X - new_right_X)
    cc = new_right_Y
    #D. new_left to left
    md = (left_Y - new_left_Y)/(left_X - new_left_X)
    cd = new_left_Y
    
    if direction == 0.:
        i_rect_xa, i_rect_ya = np.where(Y > ma*(X-left_X)+ca)
        i_rect_xb, i_rect_yb = np.where(X < (Y-cb)/mb + right_X)
        i_rect_xc, i_rect_yc = np.where(Y < mc*(X-new_right_X)+cc)
        i_rect_xd, i_rect_yd = np.where(X > (Y-cd)/md + new_left_X)
    elif direction < 90.:
        i_rect_xa, i_rect_ya = np.where(Y > ma*(X-left_X)+ca)
        i_rect_xb, i_rect_yb = np.where(Y > mb*(X-right_X)+cb)
        i_rect_xc, i_rect_yc = np.where(Y < mc*(X-new_right_X)+cc)
        i_rect_xd, i_rect_yd = np.where(Y < md*(X-new_left_X)+cd)
    elif direction == 90.:
        i_rect_xa, i_rect_ya = np.where(X > (Y-ca)/ma + left_X)
        i_rect_xb, i_rect_yb = np.where(Y > mb*(X-right_X)+cb)
        i_rect_xc, i_rect_yc = np.where(X < (Y-cc)/mc + new_right_X)
        i_rect_xd, i_rect_yd = np.where(Y < md*(X-new_left_X)+cd)
    elif direction < 180.:
        i_rect_xa, i_rect_ya = np.where(Y < ma*(X-left_X)+ca)
        i_rect_xb, i_rect_yb = np.where(Y > mb*(X-right_X)+cb)
        i_rect_xc, i_rect_yc = np.where(Y > mc*(X-new_right_X)+cc)
        i_rect_xd, i_rect_yd = np.where(Y < md*(X-new_left_X)+cd)
    elif direction == 180.:
        i_rect_xa, i_rect_ya = np.where(Y < ma*(X-left_X)+ca)
        i_rect_xb, i_rect_yb = np.where(X > (Y-cb)/mb + right_X)
        i_rect_xc, i_rect_yc = np.where(Y > mc*(X-new_right_X)+cc)
        i_rect_xd, i_rect_yd = np.where(X < (Y-cd)/md + new_left_X)
    elif direction < 270.:
        i_rect_xa, i_rect_ya = np.where(Y < ma*(X-left_X)+ca)
        i_rect_xb, i_rect_yb = np.where(Y < mb*(X-right_X)+cb)
        i_rect_xc, i_rect_yc = np.where(Y > mc*(X-new_right_X)+cc)
        i_rect_xd, i_rect_yd = np.where(Y > md*(X-new_left_X)+cd)
    elif direction == 270:
        i_rect_xa, i_rect_ya = np.where(X < (Y-ca)/ma + left_X)
        i_rect_xb, i_rect_yb = np.where(Y < mb*(X-right_X)+cb)
        i_rect_xc, i_rect_yc = np.where(X > (Y-cc)/mc + new_right_X)
        i_rect_xd, i_rect_yd = np.where(Y > md*(X-new_left_X)+cd)
    else:
        i_rect_xa, i_rect_ya = np.where(Y > ma*(X-left_X)+ca)
        i_rect_xb, i_rect_yb = np.where(Y < mb*(X-right_X)+cb)
        i_rect_xc, i_rect_yc = np.where(Y < mc*(X-new_right_X)+cc)
        i_rect_xd, i_rect_yd = np.where(Y > md*(X-new_left_X)+cd)
    
    my_mask_a = np.zeros_like(X)
    my_mask_b = np.zeros_like(X)
    my_mask_c = np.zeros_like(X)
    my_mask_d = np.zeros_like(X)
    for i in range(len(i_rect_xa)):
        my_mask_a[i_rect_xa[i],i_rect_ya[i]] = 1.
    for i in range(len(i_rect_xb)):
        my_mask_b[i_rect_xb[i],i_rect_yb[i]] = 1.
    for i in range(len(i_rect_xc)):
        my_mask_c[i_rect_xc[i],i_rect_yc[i]] = 1.
    for i in range(len(i_rect_xd)):
        my_mask_d[i_rect_xd[i],i_rect_yd[i]] = 1.
    
    my_mask = my_mask_a*my_mask_b*my_mask_c*my_mask_d
    
    return my_mask
################################################################################
#                                                                              #
# Minimise number of points in a profile                                       #
#                                                                              #
################################################################################

def RDP(x, y, e, l_testing = 0):
    """
    Algorithm to choose the minimum number of points that retains the shape of 
    y with maximum error, "error". Returns the values of y, and the coordinates 
    x for those values.
    """
    # find the end of the data
    end = len(x)
    index0 = 0
    index1 = end -1
    # The first point of our data is the same as the first point in our simplification
    x_simplified = [x[index0]]
    y_simplified = [y[index0]]
    if l_testing: fig = plt.figure()
    if l_testing: ax = fig.add_subplot(1, 1, 1)
    if l_testing: ax.plot(x, y, 'k-', marker = 'o', lw = 2)
    if l_testing: plt.pause(0.1)
    while (index0 < (end - 1)):
        #first estimate is a straight line between the first and last point of the data
        y_est = (y[index1] - y[index0])*(x - x[index0])/(x[index1] - x[index0]) + y[index0]
        d = np.abs(y - y_est)[index0:index1]
        dmax = np.max(d)
        while (dmax > e)*(index0 < index1):
            if l_testing: ax.cla()
            if l_testing: ax.plot(x, y, 'k-', marker = 'o', lw = 2)
            if l_testing: ax.plot(x[index0:(index1+1)], y_est[index0:(index1+1)], 'b--')
            index1 = index0 + np.where(d == dmax)[0][0] # find the index of distance_max > error
            if l_testing: ax.plot(x[index1], y[index1], 'ro')
            if l_testing: plt.pause(0.1)
            y_est = (y[index1] - y[index0])*(x - x[index0])/(x[index1] - x[index0]) + y[index0] # draw straight line to that point
            d = np.abs(y[index0:index1] - y_est[index0:index1]) # recalculate the differences
            dmax = np.max(d) # find the maximum differences
        x_simplified.append(x[index1])
        y_simplified.append(y[index1])
        index0 = index1
        index1 = end - 1
    return np.array(x_simplified), np.array(y_simplified)

################################################################################
#                                                                              #
# Convert winds from u, v, to speed, direction (and reverse), and along flow   #
#                                                                              #
################################################################################

def toComponents(speed, direction):
    "takes the wind speed and wind direction and breaks the wind speed into"
    "u- and v- components to be returned."
    "speed = wind speed"
    "direction = wind direction"
    u = - speed*np.sin(direction*np.pi/180)
    v = - speed*np.cos(direction*np.pi/180)
    
    return u, v

def fromComponents(u, v, isList = False):
    "Takes the two horizontal wind components and combines them into a wind"
    "speed and direction to be returned."
    "u = the zonal wind component"
    "v = the meridional wind component"
    speed = np.sqrt(np.array(u)**2 + np.array(v)**2)
    if isList:
        direction = np.arcsin(v/speed)*180/np.pi
        for iV in range(len(u)):
            if speed[iV] > 0:
                if (u[iV] >= 0) and (v[iV] > 0):
                    #if u > 0 and v > 0 SW quadrant
                    direction[iV] = 270.0 - direction[iV]
                elif (u[iV] > 0) and (v[iV] <= 0):
                    #if u > 0 and v < 0 NW quadrant
                    direction[iV] = 270.0 - direction[iV]
                elif (u[iV] <= 0) and (v[iV] > 0):
                    #if u < 0 and v > 0 SE quadrant
                    direction[iV] = 90.0 + direction[iV]
                elif (u[iV] < 0) and (v[iV] <= 0):
                    #if u < 0 and v < 0 NE quadrant
                    direction[iV] = 90.0 + direction[iV]
            else:
                direction[iV] = 0.
    else:
        direction = np.arcsin(v/speed)*180/np.pi
        if speed > 0:
            if (u >= 0) and (v > 0):
                #if u > 0 and v > 0 SW quadrant
                direction = 270.0 - direction
            elif (u > 0) and (v <= 0):
                #if u > 0 and v < 0 NW quadrant
                direction = 270.0 - direction
            elif (u <= 0) and (v > 0):
                #if u < 0 and v > 0 SE quadrant
                direction = 90.0 + direction
            elif (u < 0) and (v <= 0):
                #if u < 0 and v < 0 NE quadrant
                direction = 90.0 + direction
        else:
            direction = 0.
    return speed, direction

def transform_winds(u, v):
    """
    Computes the flow-relative wind anomalies.
    At each vertical level, computes the along and across flow anomalies (s & n)
    
    These anomalies are  with respect to the domain total horizontal-mean winds
    at that level.
    
    N.B.
        u and v must be regridded to the same horizontal grid points BEFORE they
        are passed to this function!
    
    IN:
        u = u[nt, nz, ny, nx] -> wind component in x-direction
        v = v[nt, nz, ny, nx] -> wind component in y-direction
    OUT:
        s = s[nt, nz, ny, nx] -> wind component in direction of mean flow
        n = n[nt, nz, ny, nx] -> wind component across direction of mean flow
    """
    # Init arrays for s and n
    s = np.zeros_like(u) # along flow
    n = np.zeros_like(u) # cross flow
    
    # For every time step
    for it in xrange(u.shape[0]):
        # For every height
        for k in xrange(u.shape[1]):
            u_slice = u[it,k,:,:]
            v_slice = v[it,k,:,:]
            # Find the mean wind direction
            wind_speed = np.sqrt(u_slice**2 + v_slice**2)
            wind_direction = np.where((v_slice >= 0), (np.arctan(u_slice/v_slice)*180./np.pi-180.)%360., (np.arctan(u_slice/v_slice)*180./np.pi)%360.)
            mean_dir = np.where((v_slice.mean() >= 0), (np.arctan(u_slice.mean()/v_slice.mean())*180./np.pi-180.)%360., (np.arctan(u_slice.mean()/v_slice.mean())*180./np.pi)%360.)
            # Rotate the coordinates by subtracting mean_dir
            rotated_dir = wind_direction - mean_dir
            n[it,k,:,:] = wind_speed*np.sin(rotated_dir*np.pi/180.) # cross-wind
            s[it,k,:,:] = wind_speed*np.cos(rotated_dir*np.pi/180.) # along-wind
            
    return s, n

################################################################################
#                                                                              #
# Other functions                                                              #
#                                                                              #
################################################################################
def send_email(message, subject, attachments, isAttach = True):
    import smtplib
    from email.mime.text import MIMEText
    from email.mime.image import MIMEImage
    from email.mime.multipart import MIMEMultipart
    
    # Create the container (outer) email message.
    msg = MIMEMultipart()
    msg['Subject'] = subject
    
    me = 'm.c.johnston@pgr.reading.ac.uk'
    
    msg['From'] = me
    msg['To'] = me
    msg.preamble = '\n'
    
    if isAttach:
        # Assume we know that the image files are all in PNG format
        if type(attachments) != list:
            attachments = [attachments]
        for my_file in attachments:
            # Open the files in binary mode.  Let the MIMEImage class automatically
            # guess the specific image type.
            fp = open(my_file, 'rb')
            img = MIMEImage(fp.read())
            fp.close()
            msg.attach(img)
    
    # Create a text/plain message
    body = MIMEText(message) # convert the body to a MIME compatible string
    msg.attach(body) # attach it to your main message
    
    # Send the message via our own SMTP server, but don't include the
    # envelope header.
    s = smtplib.SMTP('localhost')
    s.sendmail(me, me, msg.as_string())
    s.quit()
    
    return None

def summary(aIN):
    """
    Function to print some basic statistics for an input array.
    Requires numpy, prints max, 95th, 75th, 50th, mean 25th, 5th percentiles, and min.
    """
    
    print "\n===================="
    print "\n  MAX: "+str(np.nanmax(aIN))
    print "\n  95th: " + str(np.percentile(aIN, 95))
    print "\n  75th: " + str(np.percentile(aIN, 75))
    print "\n  50th: " + str(np.percentile(aIN, 50))
    print "\n  MEAN: " + str(np.nanmean(aIN))
    print "\n  25th: " + str(np.percentile(aIN, 25))
    print "\n  5th: " + str(np.percentile(aIN, 5))
    print "\n  MIN: " + str(np.nanmin(aIN))
    print "\n===================="

