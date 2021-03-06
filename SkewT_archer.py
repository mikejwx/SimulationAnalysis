import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy import interpolate

def getThetaE(theta, temperature, pressure, t_units = 'K', p_units = 'hPa'):
    """
    Calculates the simplified version of equivalent potential temperature from 
    Holton and Hakim (5th edition) An Introduction to Dynamic Meteorology.
    ----------------------------------------------------------------------------
    theta_e ~ theta*exp(Lv*qs/(cp*T))
    ----------------------------------------------------------------------------
    """
    # do unit conversions if neccessary
    if t_units == 'C':
        temperature = temperature + 273.15
    elif t_units == 'F':
        temperature = (temperature - 32.)*1.8 + 273.15
    
    if p_units == 'Pa':
        pressure = pressure/100.
    
    # calculate the saturated specific humidity from temperature and pressure
    qs = getQ(temperature*1., [100.], pressure*1., t_units = 'K', p_units = 'hPa')
    
    # calculate the equivalent potential temperature and return it
    theta_e = theta*np.exp(Lv*qs/(cpd*temperature))
    return theta_e

def getQ(TIN, RHIN, PIN, t_units = 'C', p_units = 'hPa'):
    """
    Uses equations from Bolton (1980) to calculate the specific humidity from the temperature (TIN), 
    the relative humidity (RHIN), and the pressure (PIN).
    
    TIN: input temperatures (deg C)
    RHIN: input relative humidities (%)
    PIN: input pressures (hPa)
    
    This method is not appropriate for temperatures lower than -80 C.
    Specific humidity is set to zero where temperature is lower than -80 C.
    """
    
    if (type(TIN) == int) or (type(TIN) in [float, np.float16, np.float32, np.float64]):
        TIN = np.array([TIN])
    if (type(RHIN) == int) or (type(RHIN) in [float, np.float16, np.float32, np.float64]):
        RHIN = np.array([RHIN])
    if (type(PIN) == int) or (type(PIN) in [float, np.float16, np.float32, np.float64]):
        PIN = np.array([PIN])
    
    #convert to degC
    if t_units == 'K':
        TIN = TIN - 273.15
    #convert to Pa
    if p_units == 'Pa':
        p = PIN
    elif p_units == 'hPa':
        p = PIN*100.
    
    # Step one - calculate saturation vapour pressure (es)    
    # Define some coefficients
    c0=0.6105851e+03
    c1=0.4440316e+02
    c2=0.1430341e+01
    c3=0.2641412e-01
    c4=0.2995057e-03
    c5=0.2031998e-05
    c6=0.6936113e-08
    c7=0.2564861e-11
    c8=-.3704404e-13
    # Calculate the saturation vapor pressure
    es = c0 + TIN*(c1 + TIN*(c2 + TIN*(c3 + TIN*(c4 + TIN*(c5 + TIN*(c6 + TIN*(c7 + TIN*c8)))))))
    
    # Step two - calculate the actual vapor pressure
    ev = RHIN*es/100. #in hPa
    
    # Step three - calculate the specific humidity
    q = ev*EPS/(p-ev*(1-EPS)) # kg kg-1
    
    q = np.where((TIN < -80.), 0., q)
    
    return q

def getGM(TIN, PIN, t_units = 'K', p_units = 'hPa'):
    """
    Calculate the moist adiabatic lapse rate, GM, for a given temperature and pressure.
    Assumes that the air is saturated at the given temperature and that all moisture is
    condensed and removed immediately from the parcel.
    
    Equation from Holton and Hakim, An Introduction to Dynamic Meteorology 5th 
    edition. Pages 504-505.
    
    TIN: Input temperature (K)
    PIN: Input pressure (hPa)
    """
    # Convert to Kelvin
    if t_units == 'C':
        # Convert from C to K
        TIN += 273.15
    elif t_units == 'K':
        # Do nothing
        TIN *= 1.
    
    # Convert to hPa
    if p_units == 'Pa':
        # Convert from Pa to hPa
        PIN /= 100.
    elif p_units == 'hPa':
        # Do nothing
        PIN *= 1.
    
    # Calculate the specific humidity at saturation
    #QIN  = getQ(TIN*1., 100., PIN*1., t_units = 'K', p_units = 'hPa')[0]
    es = e_0*np.exp((Lv/Rv)*((1./T_0) - (1./TIN)))
    QIN = EPS*es/(PIN*100.)
    # Calculate the moist adiabatic lapse rate
    A = Lv*QIN/(Rd*TIN)
    B = EPS*QIN*Lv**2./(cpd*Rd*TIN**2.)
    G_m = (g/cpd)*(1. + A)/(1. + B)
    
    return G_m

def getDew(QIN, PIN1, q_units = 'kg/kg', p_units = 'hPa'):
    """
    Convert the specific humidity into a dew point temperature given a pressure.
    Convert the specific humidity to mixing ratio, then mixing ratio to vapour pressure,
    Finally, inverts the Clausius-Clapeyron Equation to get the dew point temperature.
    
    QIN: Specific humidity (kg kg^-1)
    PIN1: Pressure (hPa)
    """
    
    # Convert pressure from hPa to Pa
    if p_units == 'hPa':
        PIN = PIN1*100.
    elif p_units == 'Pa':
        PIN = PIN1*1.
    
    if q_units == 'g/kg':
        QIN = QIN/1000.
    
    # Convert the specific humidity to a mixing ratio
    W = QIN/(1. - QIN)
    
    # Convert the mixing ratio to a vapour pressure
    e = W*PIN/(EPS + W)
    
    # Calculate the dew point temperature from the vapour pressure
    Td = ((1./273.15) - (Rv/Lv)*np.log(e/e_0))**(-1.)
    
    return Td

def getCAPE(temperature, q, pressure, parcel_type = 0, t_units = 'K', q_units = 'kg/kg', p_units = 'hPa'):
    """
    Calculate a parcel ascent, first following the dry adiabat until the 
    Lifiting Condensation Level (LCL), then following a pseudo-adiabat through 
    to 100 hPa or the highest input pressure level.
    
    Once a parcel ascent has been computed, the Convective Available Potential
    Energy (CAPE) is calculated between the lowest level of free convection 
    (LFC) and the highest level of neutral buoyancy (LNB). Similarly, Convective 
    Inhibition (CIN) is calculated between the LCL and the lowest LFC.
    ----------------------------------------------------------------------------
    INPUT:
    temperature = An environmental temperature profile in Kelvin, ie. f(z)
    q           = An environmental specific humidity profile in kg/kg, ie. f(z)
    pressure    = An environmental pressure profile in hPa, i.e. f(z)
    
    parcel_type = type of parcel used to compute ascent, either 0 = surface, or 
                  1 = mean of lowest 500 m ~ 50 hPa
    x_units     = the units of the provided temperature, q, or pressure
    """
    # First convert everything to the correct units
    if t_units == 'C':
        temperature += 273.15
    
    if q_units == 'g/kg':
        q /= 1000.
    
    if p_units == 'Pa':
        pressure /= 100.
    
    # Next, calculate the parcel
    from scipy import interpolate, integrate
    if parcel_type == 0:
        # Surface Parcel
        T0 = temperature[0]
        Q0 = q[0]
        P0 = pressure[0]
    elif parcel_type == 1:
        # Lowest ~500 m
        DP = g*500./100. # hPa change in pressure between the surface and 500 m for rho = 1 kg/m3
        
        theta_fun = interpolate.interp1d(x = pressure, y = temp2theta(temperature, pressure), fill_value = 'extrapolate')
        q_fun     = interpolate.interp1d(x = pressure, y = q, fill_value = 'extrapolate')
        p_layer = pressure[0] - np.arange(DP)
        T0 = PTtoTemp(integrate.trapz(x = p_layer, y = theta_fun(p_layer)*p_layer)/integrate.trapz(p_layer, p_layer), pressure[0])
        Q0 = integrate.trapz(x = p_layer, y = q_fun(p_layer)*p_layer)/integrate.trapz(p_layer, p_layer)
        P0 = pressure[0]

    # Then start parcel ascent:
    # Need to know the LCL_temperature and pressure
    Td0 = getDew(Q0, P0, q_units = 'kg/kg', p_units = 'hPa')
    # Use Bolton (1980)
    A = 1./(Td0 - 56.)
    B = np.log(T0/Td0)/800.
    LCL_temperature  = (1./(A + B)) + 56.
    LCL_pressure     = pressure[0]*(LCL_temperature/temperature[0])**(cpd/Rd)

    n_steps = 1000
    T_parcel = np.zeros(n_steps)
    P_parcel = np.zeros(n_steps)
    T_parcel[0] = T0
    P_parcel[0] = P0
    dp = (np.min([pressure[-1], 100.]) - pressure[0])/n_steps
    for i in xrange(1, n_steps):
        # do 1000 step integration between the surface and highest (or 100. hPa) pressure
        rho = 100.*P_parcel[i-1]/(Rd*T_parcel[i-1]) # Equation of state
        dz = - 100.*dp/(rho*g) # Hydrostatic balance
        if P_parcel[i-1] > LCL_pressure:
            # Parcel ascends dry adiabatically to the LCL
            T_parcel[i] = T_parcel[i-1] - g*dz/cpd # Dry Adiabatic Lapse Rate
            P_parcel[i] = P_parcel[i-1] + dp # Pressure at new level
        elif P_parcel[i-1] <= LCL_pressure:
            # Parcel ascends pseudoadiabatically through the rest of the atmosphere
            T_parcel[i] = T_parcel[i-1] - getGM(T_parcel[i-1]*1., P_parcel[i-1]*1., t_units = 'K', p_units = 'hPa')*dz # Pseudoadiabatic lapse rate
            P_parcel[i] = P_parcel[i-1] + dp # pressure at new level

    # Find LFC, LNB, Calculate CAPE and CIN
    T_environment = interpolate.interp1d(x = pressure, y = temperature, fill_value = 'extrapolate')(P_parcel)
    rho = 100.*P_parcel/(Rd*T_parcel)
    dz = -100.*dp/(rho*g)
    # CAPE above LCL, above lowest LFC, below highest LNB
    CAPE = np.sum(np.where(((T_parcel - T_environment) > 0.)*(P_parcel < LCL_pressure), g*(T_parcel - T_environment)*dz/T_environment, 0))

    # lowest LFC, the lowest level above or equal to the LCL that is positively buoyant
    LFC_pressure = np.max(np.where(((T_parcel - T_environment) > 0.)*(P_parcel <= LCL_pressure), P_parcel, 0.))
    # CIN above LCL, below lowest LFC
    CIN = np.sum(np.where(((T_parcel - T_environment) < 0.)*(P_parcel < LCL_pressure)*(P_parcel > LFC_pressure), g*(T_parcel - T_environment)*dz/T_environment, 0))

    return CAPE, CIN, T_parcel, LCL_pressure, LFC_pressure

def getCAPEv1(TIN1, QIN, PIN, parcel_type = 0):
    """
    Function to calculate surface based CAPE.
    Tp = User provided parcel (e.g. surface based virtual temperature)
    TIN1 = User provided profile of temperatures in K
    QIN  = User provided profile of specific humidity in kg/kg
    PIN  = User provided profile of pressure in hPa
    """
    from scipy import interpolate
    
    # Convert to the correct units
    if (TIN1[0] - 273.15) < 0:
        # then the given temperatures are in celcius and need to be converted
        Te = TIN1 + 273.15 # converts the environment temperature to Kelvin
    else:
        Te = TIN1*1.
    
    # Locally define some constants
    G_d = g/cpd # dry adiabatic lapse rate
    
    ### Parcel Type
    if parcel_type == 0:
        # surface based parcel
        Tp = Te[0]
        Q_0  = QIN[0] # get the surface specific humidity
        T_d0 = getDew(Q_0*1., PIN[0]) # get the surface dew point
    elif parcel_type == 1:
        # mean over lowest ~ 500 m
        rho = PIN[0]/(Rd*Te[0])
        dp500 = 500.*rho*g/100.
        p_dummy = PIN[0] - np.arange(0., dp500+0.1)
        Tp = PTtoTemp(np.mean(interpolate.interp1d(PIN, temp2theta(Te, PIN))(p_dummy)), PIN[0])
        Q_0  = np.mean(interpolate.interp1d(PIN, QIN)(p_dummy)) # get the surface specific humidity
        T_d0 = getDew(Q_0*1., PIN[0]) # get the surface dew point
    else:
        print 'WARNING: Unknown parcel type!!!'
    
    # Calculate the LCL temperature using the University of Wyoming's method (Bolton, 1980)
    A = 1./(T_d0 - 56.)
    B = np.log(Tp/T_d0)/800.
    LCL  = (1./(A + B)) + 56.
    
    ## Find the LFC ##
    LFC  = getLFC(Tp*1., Te*1., QIN*1., PIN*1.) # parcel temperature at level of free convection
    
    ## Initialise the parcel ascent list and the CAPE to accumulate ##
    Tps  = [Tp] # parcel temperature
    CAPE = 0 # convective available potential energy
    CIN  = 0 # Convective Inhibition
    
    ### Start accumulating CAPE ###
    if LFC != None:
        # an LFC has been found, there should be some CAPE at some level
        for level in range(1,len(Te)):
            p_level_half = 0.5*(PIN[level] + PIN[level-1])
            # convert the thickness of the layer from pressure to height units
            rho = (p_level_half*100.)/(Rd*Tp)
            dz = - (PIN[level] - PIN[level-1])*100./(rho*g) # dp/rho*g, rho = p/RT
            #print 'Tp:' + str(int(Tp)) + '\nLCL:' + str(int(LCL)) + '\nLFC:' + str(int(LFC))
            if PIN[level] != 0.0:
                #print round(Tp, 1)
                if Tp >= LCL:
                    #print 'Original parcel below LCL'
                    # if the parcel is warmer than the LCL it is below the LCL and should cool at the dry adiabatic lapse rate
                    Tp -= G_d*dz # dry adiabatic lapse rate
                    if Tp >= LCL:
                        #print 'Lifted parcel is still below LCL'
                        # if the new parcel temperature is still below the LCL append it and move on
                        Tps.append(Tp)
                    else:
                        ## Step 1: undo the ascent
                        Tp += G_d*dz
                        ## Step 2: find the difference between the parcel temperature and the LCL temperature
                        dT = Tp - LCL
                        # find the height change needed for the parcel to get to the LCL temperature
                        delta_z = dT/G_d
                        ## Step 3: set parcel temperature to the LCL and then lift the parcel the remaining distance but moist adiabatically
                        Tp = LCL
                        p_level_middle = 0.5*(PIN[level] + (PIN[level] - PIN[level-1])*delta_z/dz + PIN[level-1])
                        Tp -= getGM(LCL, p_level_middle)*(dz - delta_z)
                        # append the correct parcel temperature
                        Tps.append(Tp)
                        if Tps[-1] >= LFC:
                            #print 'Lifted parcel is now above LCL, but below the LFC'
                            # If we're still below the LFC, only accumulate some CIN
                            Tp_level_half = 0.5*(Tps[-1] + Tps[-2])
                            Te_level_half = 0.5*(Te[level] + Te[level-1])
                            if Tp_level_half < Te_level_half:
                                # Only accumulate CIN if negatively buoyant!
                                CIN += g*(Tp_level_half - Te_level_half)*dz/Te_level_half
                            #print 'Accumulated some CIN between LCL and Tps[-1]'
                        else:
                            #print 'Lifted parcel is now above the LCL, and above the LFC'
                            # If we're above the LFC, accumulate CIN between LCL and LFC
                            # Also accumulate CAPE from LFC to Tp
                            
                            # Find the distance between the LCL and LFC
                            dz0 = (LFC - LCL)/getGM(0.5*(LCL + LFC), p_level_middle)
                            
                            # accumulate CIN in the LCL -> LFC layer
                            Tp_level_half = 0.5*(LCL + LFC)
                            Te_level_half = 0.5*(Te[level] + Te[level-1]) # An easy approximation
                            if Tp_level_half < Te_level_half:
                                # Only accumulate CIn if negatively buoyant
                                CIN += g*(Tp_level_half - Te_level_half)*dz0/Te_level_half
                            
                            # find the distance between the LFC and Tp
                            dz1 = (Tps[-1] - LFC)/getGM(0.5*(LFC + Tps[-1]), p_level_middle)
                            
                            # accumulate CAPE in the LFC -> Tp layer
                            Tp_level_half = 0.5*(Tps[-1] + LFC)
                            Te_level_half = 0.5*(Te[level] + Te[level-1]) # An easy approximation
                            if Tp_level_half > Te_level_half:
                                # if the parcel is warmer than the environment, accumulate CAPE
                                CAPE += g*(Tp_level_half - Te_level_half)*dz1/Te_level_half
                            #print 'Accumulated some CIN between LCL and LFC, then some CAPE between LFC and Tps[-1]'
                elif Tp >= LFC:
                    #print 'Original parcel between LCL and LFC'
                    # if the parcel is cooler than the LCL but warmer than the LFC it should cool at the moist adiabatic lapse rate, but not accumulate CAPE
                    Tp -= getGM(Tps[-1], p_level_half)*dz
                    # append the correct parcel temperature
                    Tps.append(Tp)
                    if Tp >= LFC:
                        #print 'Lifted parcel is still below the LFC'
                        # If we're still below the LFC, only accumulate CIN
                        Tp_level_half = 0.5*(Tps[-1] + Tps[-2])
                        Te_level_half = 0.5*(Te[level] + Te[level-1])
                        if Tp_level_half < Te_level_half:
                            # Only accumulate CIN if negatively buoyant
                            CIN += g*(Tp_level_half - Te_level_half)*dz/Te_level_half
                            #print 'Accumulated some CIN'
                    else:
                        #print 'Lifted parcel is now above the LFC'
                        # If we're above the LFC, accumulate CIN between Tps[-2] and LFC
                        # Also accumulate CAPE from LFC to Tps[-1]
                        
                        # Find the distance between the Tps[-2] and LFC
                        dz0 = (LFC - Tps[-2])/getGM(0.5*(LCL + LFC), p_level_half) # An easy approximation for pressure
                        
                        # accumulate CIN in the Tps[-2] -> LFC layer
                        Tp_level_half = 0.5*(LFC - Tps[-2])
                        Te_level_half = 0.5*(Te[level] + Te[level-1]) # An easy approximation
                        if Tp_level_half < Te_level_half:
                            # Only accumulate CIN if negatively buoyant
                            CIN += g*(Tp_level_half - Te_level_half)*dz0/Te_level_half
                        
                        # find the distance between the LFC and Tps[-1]
                        dz1 = (Tps[-1] - LFC)/getGM(0.5*(Tps[-1] + LFC), p_level_half)
                        
                        # accumulate CAPE in the LFC -> Tp layer
                        Tp_level_half = 0.5*(Tps[-1] + LFC)
                        Te_level_half = 0.5*(Te[level] + Te[level-1]) # An easy approximation
                        if Tp_level_half > Te_level_half:
                            # if the parcel is warmer than the environment, accumulate CAPE
                            CAPE += g*(Tp_level_half - Te_level_half)*dz1/Te_level_half
                        #print 'Accumulated some CIN between Tps[-2] and LFC, then some CAPE between LFC and Tps[-1]'
                else:
                    # Parcel is above the LFC, cool at moist adiabatic lapse rate and accumulate CAPE
                    Tp -= getGM(Tps[-1], p_level_half)*dz
                    Tps.append(Tp)
                    
                    # get the mean temperature between levels for the parcel and environment
                    Tp_level_half = 0.5*(Tps[-1] + Tps[-2])
                    Te_level_half = 0.5*(Te[level] + Te[level-1])
                    if Tp_level_half > Te_level_half:
                        # if the parcel is warmer than the environment, accumulate CAPE
                        CAPE += g*(Tp_level_half - Te_level_half)*dz/Te_level_half
                    # Don't accumulate CIN because we've already passed the LFC
            else:
                Tp -= G_d*dz
                Tps.append(Tp)
    else:
        # there is no LFC
        for level in range(1,len(Te)):
            p_level_half = 0.5*(PIN[level] + PIN[level-1])
            # convert the thickness of the layer from pressure to height units
            rho = (p_level_half*100.)/(Rd*Tp)
            dz = - (PIN[level] - PIN[level-1])*100./(rho*g) # dp/rho*g, rho = p/RT
            
            if PIN[level] != 0.0:
                if Tp >= LCL:
                    # if the parcel is warmer than the LCL it is below the LCL and should cool at the dry adiabatic lapse rate
                    Tp -= G_d*dz # dry adiabatic lapse rate
                    if Tp >= LCL:
                        # if the new parcel temperature is still below the LCL append it and move on
                        Tps.append(Tp)
                    else:
                        # the new parcel is now above the LCL
                        # undo the ascent
                        Tp += G_d*dz
                        # find the difference between the parcel temperature and the LCL temperature
                        dT = Tp - LCL
                        # find the height change needed for the parcel to get to the LCL temperature
                        delta_z = dT/G_d
                        # set parcel temperature to the LCL and then lift the parcel the remaining distance but moist adiabatically
                        Tp = LCL
                        p_level_middle = 0.5*(PIN[level] + (PIN[level] - PIN[level-1])*delta_z/dz + PIN[level-1])
                        Tp -= getGM(LCL, p_level_middle)*(dz - delta_z)
                        # append the correct parcel temperature
                        Tps.append(Tp)
                else:
                    # above the LCL and cool moist adiabatcally
                    Tp -= getGM(Tps[-1], p_level_half)*dz
                    Tps.append(Tp)
                    
                    # Accumulate some CIN
                    Tp_level_half = 0.5*(Tps[-1] + Tps[-2])
                    Te_level_half = 0.5*(Te[level] + Te[level-1])
                    CIN += g*(Tp_level_half - Te_level_half)*dz/Te_level_half
            else:
                Tp -= G_d*dz
                Tps.append(Tp)
    
    # calculate the pressure of the LCL and the LFC
    LCLp = interpolate.interp1d(Tps, PIN, fill_value = 'extrapolate')(LCL)#PIN[0]*(LCL/Te[0])**(cpd/R)
    if LFC != None:
        LFCp = interpolate.interp1d(Tps, PIN, fill_value = 'extrapolate')(LFC)#PIN[0]*(LFC/Te[0])**(cpd/R)
    else:
        LFCp = None
    
    return CAPE, CIN, Tps, LCLp, LFCp
    

def PTtoTemp(theta, PIN, t_units = 'K', p_units = 'hPa'):
    """
    Converts potential temperature to temperature.
    theta: potential temperature (K)
    PIN: pressure (hPa)
    
    Other recognised units include, C - Celsius and Pa - Pascals.
    """
    
    # Perform any required unit conversions
    if t_units == 'C':
        theta = theta + 273.15
    
    if p_units == 'Pa':
        PIN = PIN/100.
    
    Temperature = theta/((1000./PIN)**(Rd/cpd))
    
    return Temperature

def temp2theta(TIN, PIN, P0 = 1000.):
    """
    Converts temperature to potential temperature
    TIN: input temperature (K)
    PIN: input pressures (hPa)
    """
    
    # check that a temperature is passed in Kelvin.
    if type(TIN) == np.ndarray:
        if TIN[0] - 273.15 < 0:
            TIN += 273.15
    elif TIN - 273.15 < -100:
        TIN += 273.15
    
    theta = TIN*(P0/PIN)**(Rd/cpd)
    
    return theta

def getLFC(Tp, Te, QIN, PIN):
    """
    Function to calculate the level of free convection LFC.
    Tp = the temperature of a parcel, K
    Te = the temperature of the environment, K
    QIN = the specific humidity profile, kg/kg
    PIN = the pressure profile, hPa
    """
    # Locally define dry adiabatic lapse rate
    G_d  = g/cpd
    
    # Lets get some unit conversions going
    Q_0  = QIN[0] # surface specific humidity
    T_d0 = getDew(Q_0, PIN[0]) # surface dew point
    
    # Use the University of Wyoming's formula (Bolton, 1980) to find LCL temperature, K
    A    = 1./(T_d0 - 56.)
    B    = np.log(Tp/T_d0)/800.
    LCL  = (1./(A + B)) + 56.
    
    # Begin lifting the parcel to find where it is positively buoyant
    for level in range(len(Te)-1):
        # for each observed vertical level
        p_level_half = 0.5*(PIN[level+1] + PIN[level]) # pressure in between current level and the level above
        rho = (p_level_half*100.)/(Rd*Tp) # air density in the layer between the current level and the level above
        dz = -(PIN[level+1] - PIN[level])*100./(rho*g) # hydrostatic balance to get the height between the two levels e.g. dp/rho*g, rho = p/RT
        if Tp > LCL:
            #print 'Parcel is below the LCL'
            # if the parcel is warmer than the LCL it is below the LCL and should cool at the dry adiabatic lapse rate
            Tp -= G_d*dz # dry adiabatic lapse rate
            if Tp < LCL:
                # the new parcel is now above the LCL
                # undo the ascent
                Tp += G_d*dz
                # find the difference between the parcel temperature and the LCL temperature
                dT = Tp - LCL
                # find the height change needed for the parcel to get to the LCL temperature
                delta_z = dT/G_d
                # set parcel temperature to the LCL and then lift the parcel the remaining distance but moist adiabatically
                Tp = LCL
                p_level_middle = 0.5*(PIN[level] + (PIN[level] - PIN[level-1])*delta_z/dz + PIN[level-1])
                Tp -= getGM(LCL, p_level_middle)*(dz - delta_z)
        else:
            # above LCL cool at moist lapse rate
            Tp -= getGM(Tp, p_level_half)*dz
            # print [PIN[level], Tp, Te[level], level]
            if Tp >= Te[level+1]:
                # if the parcel is warmer than environment it is free
                LFCT = Tp # set LFC temperature
                break
    if level == (len(Te)-2):
        LFCT = None
    
    return LFCT
################################################################################
# Methods to plot a skew T when passed temperature and specific humidity       #
# values varying with height and pressure                                      #
################################################################################

### Constants
g = 9.81 #gravity
cpd = 1005. #specific heat capacity of dry air
cl = 4200. #specific heat capacity of liquid water
R = 8.314 #gas constant
Mrd = 0.0289 #molar mass of dry air
Mrv = 0.018 #molar mass of water vapour
EPS = Mrv/Mrd #ratio of molar masses of water vapour and dry air
Rd = 287.05 #gas constant for dry air
Rv = R/Mrv #gas constant for water vapour
Lv = 2.501e06 #latent heat of vapourisation
Lm = 3.3e05 #latent heat of melting
Ls = Lv + Lm #latent heat of sublimation
gamma_d = -g/cpd #dry adiabatic lapse rate
p0 = 1000. #standard pressure level
adiabat_p = np.arange(1050., 99., -5.) #pressure levels to calculate adiabats
q_p = np.arange(1050., 500., -5.) # pressure levels on which to plot lines of constant specific humidity (q)
target_q = np.array([0.1, 1, 2, 4, 7, 10, 15, 20, 30])*1e-3 # constant specific humidities to plot lines
T0 = np.arange(-40, 300., 20.) # p0 level temperatures to calculate dry adiabats
T1 = np.arange(-40., 90.1, 5)
e_0 = 611.2 # saturation vapour pressure for T_0 = 273.15
T_0 = 273.15

def dry_adiabats():
    """
	Uses T0 to draw dry adiabats. Uses pressure levels every 5hPa and 
    poisson equation for potential temperature
	"""
    adiabats = np.zeros((len(adiabat_p), len(T0)))
    for T in range(len(T0)):
        adiabats[:,T] = (T0[T]+273.15)/(p0/adiabat_p)**(Rd/cpd) - 273.15
    	
    return adiabats

def get_saturation_mixing_ratio(temp):
    """
	Calculate the saturation mixing ratio for a given temperature (temp).
    temp: input temperature (deg C)
    """
    es = 6.112*np.exp(17.67*temp/(temp+243.5))
    ws = EPS*es/(adiabat_p - es)

    return ws

def moist_adiabats():
    """
    Draws the moist adiabats for the skew T.
    Uses T1, temperatures in deg C, at 1000 hPa to draw these adiabats.
    Also uses pressure every 5 hPa. 
    """
    moist_adiabats = np.zeros((len(adiabat_p), len(T1)))
    i_1000 = np.where(adiabat_p == 1000)[0][0]
    moist_adiabats[i_1000,:] = T1
    for T in range(i_1000+1,len(adiabat_p)):
        # loop over the pressure levels
        for t in range(len(T1)):
            rho = (0.5*(adiabat_p[T-1]+adiabat_p[T])*100.)/(Rd*(moist_adiabats[T-1,t] + 273.15)*(1 + 0.608*getQ(moist_adiabats[T-1,t], 100., 0.5*(adiabat_p[T-1]+adiabat_p[T]))[0]))
            dz = (500./(rho*g)) # dp/rho*g, rho = p/RT
            moist_adiabats[T,t] = moist_adiabats[T-1, t] - \
                getGM(moist_adiabats[T-1,t]+273.15, 0.5*(adiabat_p[T-1] + adiabat_p[T]))*dz
    for T in range(i_1000, 0, -1):
        # loop over the pressure levels
        for t in range(len(T1)):
            rho = (0.5*(adiabat_p[T-1]+adiabat_p[T])*100.)/(Rd*(moist_adiabats[T,t] + 273.15)*(1 + 0.608*getQ(moist_adiabats[T,t], 100., 0.5*(adiabat_p[T-1]+adiabat_p[T]))[0]))
            dz = (500./(rho*g)) # dp/rho*g, rho = p/RT
            moist_adiabats[T-1,t] = moist_adiabats[T, t] + \
                getGM(moist_adiabats[T,t]+273.15, 0.5*(adiabat_p[T-1] + adiabat_p[T]))*dz

    return moist_adiabats

def const_qs():
    """
    Creates curves for lines of constant specific humidity
    """
    const_qs = np.zeros((len(q_p), len(target_q)))
    for iq in range(len(target_q)):
        const_qs[:,iq] = getDew(target_q[iq]/(1. + target_q[iq]), q_p)-273.15
    
    return const_qs

wmo_p = [1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100]

def skew(temperature, pressure):
    """
    Creates a skew for an input temperature.
    The temperature could represent, actual air temperature, dry adiabats, moist adiabats,
    lines of constant specific humidity, actual specific humidity...
    The amount of skew is designed to roughly match the amount of skew on UWyoming's skewT plots
    """
    
    c = 0.01219513#np.tan(30*np.pi/180.)
    skewedT = np.zeros_like(temperature) + temperature
    for p_i in range(len(pressure)):
        dp = np.log10(pressure[p_i]/1000.)
        skewedT[p_i] -= dp/c
    
    return skewedT

def paper():
    """
    Plots the paper used to plot skew T diagrams on without requiring input data.
    Additional plots of skewed temperatures can be overlaid.
    It might be possible to include this is a subplot and wind barbs in a 
    separate subplot using windBarbs().
    """
    theta = dry_adiabats()
    thetae = moist_adiabats()
    my_q = const_qs()
    for sfc in range(len(T0)):
        theta[:,sfc] = skew(theta[:,sfc], adiabat_p)
    for sfc in range(len(T1)):
        thetae[:,sfc] = skew(thetae[:,sfc], adiabat_p)
    for sfc in range(len(target_w)):
        my_q[:,sfc] = skew(my_q[:,sfc], q_p)
    
    therms = np.arange(-120., 40.1, 5.)
    isotherms = np.zeros((len(adiabat_p), len(therms)))
    for level in range(len(adiabat_p)):
        isotherms[level,:] = therms
    for sfc in range(len(therms)):
        isotherms[:,sfc] = skew(isotherms[:,sfc], adiabat_p)
    
    fig = plt.gcf()
    fig.set_size_inches(6, 8)
    plt.semilogy(theta, adiabat_p, color = 'g', lw = 0.25)
    plt.semilogy(thetae, adiabat_p, color = 'b', lw = 0.25)
    plt.semilogy(isotherms, adiabat_p, color = 'gray', lw = 0.25, ls = '--')
    plt.semilogy(my_q, q_p, color = 'purple', lw = 0.5, ls = '--')
    q_labels = [q*1000. if q*1000 < 1 else int(q*1000) for q in target_q]
    for iq in range(len(target_q)):
        plt.text(my_q[np.where(q_p == 850.)[0][0], iq], 850., str(q_labels[iq]), color = 'purple')
    plt.gca().invert_yaxis()
    plt.ylim([1050, 100])
    plt.gca().yaxis.grid()
    plt.yticks(wmo_p, wmo_p)
    plt.xlim([-40, 40])

def windBarbs(u_in, v_in, p_in, gs):
    """
    Plots wind barbs for given u and v wind components at wmo pressure levels.
    """
    # Only use winds on the wmo levels
    u_wmo = []
    v_wmo = []
    for i in range(len(p_in)):
        if p_in[i] in wmo_p:
            u_wmo.append(u_in[i])
            v_wmo.append(v_in[i])
    u_wmo = np.array(u_wmo)
    v_wmo = np.array(v_wmo)
    # Adds wind barbs to the SkewT
    ax2 = plt.subplot(gs[1])
    ax2.barbs(np.zeros_like(wmo_p), wmo_p, u_wmo, v_wmo, barbcolor = 'k', clip_on=False, zorder=100)
    plt.xlim([-1, 1])
    plt.ylim([1005., 100.])
    plt.yticks([])
    plt.xticks([])
    plt.axis('off')

def plotSkewT(temp, t_dew, p, u = np.array([-999]), v = np.array([-999]), 
              CAPE = False, date = -999, my_title = '', temp_col = ['r'], 
              dew_col = ['b']):
    """
    Plots a skewT-logP diagram on the paper I've defined above.
    Requires an input temperature, dewpoint, pressure, and optional u- and v- 
    wind components. If CAPE is set to .True. then a call is made to calculate
    the surface-based parcel ascent. This provides CAPE, CIN computed using
    temperature.
    
    Here, CAPE is defined as the positively buoyant region between the level
    of free convection and the highest level of neutral buoyancy.
    
    CIN, however, is defined as the negatively buoyant region between the 
    lifting condensation level and the level of free convection.
    ---------------------------------------------------------------------------
    Optimal input units:
    temp    = deg C
    t_dew   = deg C
    p       = hPa
    u       = kts
    v       = kts
    """
    ### get the correct units for temp and t_dew ###
    # should be in deg C
    from scipy import interpolate 
    
    gs = gridspec.GridSpec(1, 2, width_ratios = [4, 1])
    gs.update(wspace=0.0, hspace=0.0)
    #paper()
    
    ## plot the skew T paper ##
    theta = dry_adiabats()
    thetae = moist_adiabats()
    my_q = const_qs()
    # skew the adiabats and the mixing ratio lines
    for sfc in range(len(T0)):
        theta[:,sfc] = skew(theta[:,sfc], adiabat_p)
    for sfc in range(len(T1)):
        thetae[:,sfc] = skew(thetae[:,sfc], adiabat_p)
    for sfc in range(len(target_q)):
        my_q[:,sfc] = skew(my_q[:,sfc], q_p)
    # get some isotherms
    therms = np.arange(-120., 40.1, 5.)
    isotherms = np.zeros((len(adiabat_p), len(therms)))
    for level in range(len(adiabat_p)):
        isotherms[level,:] = therms
    for sfc in range(len(therms)):
        isotherms[:,sfc] = skew(isotherms[:,sfc], adiabat_p)
    
    #start the plot for the paper
    fig = plt.gcf()
    fig.set_size_inches(10, 8)
    ax = fig.add_subplot(gs[0])
    # dry adiabats
    ax.semilogy(theta, adiabat_p, color = 'g', lw = 0.25)
    # pseudo adiabats
    ax.semilogy(thetae, adiabat_p, color = 'b', lw = 0.25)
    # isotherms
    ax.semilogy(isotherms, adiabat_p, color = 'gray', lw = 0.25, ls = '--')
    # specific humidities
    ax.semilogy(my_q, q_p, color = 'purple', lw = 0.5, ls = '--')
    # label the specific humidities
    q_labels = [q*1000. if q*1000 < 1 else int(q*1000) for q in target_q]
    for iq in range(len(target_q)):
        ax.text(my_q[np.where(q_p == 850.)[0][0], iq], 850., str(q_labels[iq]), color = 'purple')
    
    ## plot the sounding ##
    # make sure our input are lists
    if type(temp) != list:
        #print( 'WARNING! converting to a list of temperature profiles.')
        temp = [temp]
    if type(t_dew) != list:
        #print( 'WARNING! Converting to a list of dewpoint profiles.')
        t_dew = [t_dew]
    if type(p) != list:
        #print( 'WARNING! Converting to a list of pressure profiles.')
        p = [p]
    if type(u) != list:
        #print( 'WARNING! Converting to a list of uwind profiles.')
        u = [u]
    if type(v) != list:
        #print( 'WARNING! Converting to a list of vwind profiles.')
        v = [v]
    
    # some error checking to ensure that there is a matching temperature and dewpoint list
    if len(temp) != len(t_dew):
        print( 'ERROR! You have an unequal number of temperature and dewpoint profiles')
    if len(temp_col) != len(temp):
        #print( 'WARNING! The number of temperature colors provided does not match the number of temperature profiles. Taking the same dewpoint color for each profile.')
        while len(temp_col) != len(temp):
            temp_col.append(temp_col[-1])
    if len(dew_col) != len(t_dew):
        #print( 'WARNING! The number of dewpoint colors provided does not match the number of dewpoint profiles. Taking the same dewpoint color for each profile.')
        while len(dew_col) != len(t_dew):
            dew_col.append(dew_col[-1])
    
    for i in range(len(temp)):
        #each i in the list is a profile of t, td, and p
        # make sure they are the right units
        if (np.min(temp[i]) > 0) and (np.max(temp[i] > 100.)):
            temp[i] -= 273.15
        if (np.min(t_dew[i]) > 0) and (np.max(t_dew[i] > 100.)):
            t_dew[i] -= 273.15
        ax.semilogy(skew(temp[i], p[i]), p[i], color = temp_col[i], lw = 2)
        ax.semilogy(skew(t_dew[i], p[i]), p[i], color = dew_col[i], lw = 2)
        if CAPE:
            my_q = getQ(t_dew[i]+273.15, np.zeros_like(t_dew[i])+100., p[i])
            my_CAPE, my_CIN, my_Parcel, LCLp, LFCp = getCAPE(temp[i], my_q, p[i], parcel_type = 1)
            my_title = my_title + '\nCAPE = ' + str(round(my_CAPE, 0)) + 'Jkg$^{-1}$, CIN = ' + str(round(my_CIN, 0)) + 'Jkg$^{-1}$'
            if np.nanmin(my_Parcel) > 0:
                my_Parcel = np.array([obs - 273.15 for obs in my_Parcel])
            ax.semilogy(skew(my_Parcel, p[i]), p[i], color = 'gray', lw = 2)
            ax.semilogy([-100, 100], [LCLp, LCLp], color = 'k')
            ax.semilogy([-100, 100], [LFCp, LFCp], color = 'k')
    ax.set_ylim([np.max([np.max([np.max(p_i) for p_i in p]), 1025.]), 100])
    ax.grid(axis = 'y')
    ax.set_yticks(wmo_p)
    ax.set_yticklabels(wmo_p)
    ax.set_xlim([-40, 40])
    ax.set_xlabel('Temperature ($^{\circ}$C)')
    ax.set_ylabel('Pressure (hPa)')
    
    if date != -999:
        my_title = date + '\n' + my_title
    ax.set_title(my_title)    
    
    ## plot the wind barbs ##
    for i in range(len(u)):
        if u[i][0] != -999:
            # windBarbs(u, v, p, gs)
            # only use wind obs on wmo pressure levels so we can see all the wind barbs
            #print 'Length p[i] = ' + str(len(p[i]))
            #print 'Length u[i] = ' + str(len(u[i]))
            u_wmo = interpolate.interp1d(p[i], u[i], fill_value = 'extrapolate')(wmo_p)
            v_wmo = interpolate.interp1d(p[i], v[i], fill_value = 'extrapolate')(wmo_p)
            
            # add wind barbs to the SkewT
            ax2 = fig.add_subplot(gs[1], sharey=ax)
            ax2.barbs(np.zeros_like(wmo_p), wmo_p, u_wmo, v_wmo, barbcolor = 'k', \
                clip_on = False, zorder=100)
            ax2.set_xlim([-1, 1])
            ax2.set_ylim([np.max([np.max([np.max(p_i) for p_i in p]), 1025.]), 100])
            #print 'Plotted wind barbs'
            ax2.set_xticks([])
            ax2.axis('off')


    
###real data to plot
#path = '/home/xb899100/bin/'
#import scipy.io as io
#data = io.readsav('/home/xb899100/bin/profilearmsonde_save_20010401-20060816.dat')
#sounding = 1130
#temp = data['tdry'][sounding,3:-1]
#Q = data['qsat'][sounding,3:-1]*data['rh'][sounding,3:-1]/100.
#q_kgkg = Q/1000.
#temp_v = temp #*(1 + 0.608*q_kgkg)
#Z = data['alt'][sounding,3:-1]
#p = data['pres'][3:-1]
#t_dew = getDew(q_kgkg, p) - 273.15
#Tvp = np.nanmean(theta(temp_v[0:11], p[0:11]))
#myCAPE = getCAPE(Tvp, temp_v-273.15, q_kgkg, p, Z)
# test_data = []
# with open(path + 'test_sounding.csv', 'r') as test_file:
    # for line in test_file:
        # test_data.append(line)

# p = []
# temp = []
# t_dew = []
# Z = []
# for line in test_data:
	# myline = line.split(',')
	# p.append(float(myline[0]))
	# temp.append(float(myline[2]))
	# t_dew.append(float(myline[3]))
	# Z.append(float(myline[1]))

# execfile('/glusterfs/greybls/users/xb899100/SatData/readCloudFrequencies.py')
# ###Key inputs
# lsmQ = True
# innerR = 0.16
# outerR = 0.25
# myres = 22.5
# lowerP = 950
# upperP = 850
# execfile('/glusterfs/greybls/users/xb899100/SatData/soundingParse.py')
# i = 0
# i += 10
# #p = np.array(p)
# p = interpolated['p'][:,i]
# #temp = np.array(temp)
# temp = interpolated['Temp'][:,i] -273.15
# #t_dew = np.array(t_dew)
# RH = interpolated['RH'][:,i]
# #Z = np.array(Z)


# #get parcel
# #q_kgkg = getQ(t_dew, 100., p)
# q_kgkg = getQ(temp, RH, p)
# t_dew = getDew(q_kgkg, p) - 273.15
# temp_v = (temp+273.15)*(1 + 0.608*q_kgkg)


# CAPE = getCAPE(temp_v[0], temp_v -273.15, q_kgkg, p, Z)

# #CAPE = getCAPE(temp_v[0], temp, q_kgkg, p)

# parcel_temp = np.array(CAPE[1])-273.15
#paper()
#plt.plot(skew(np.array(temp)-273.15, p), p, color = 'r')
#plt.plot(skew(t_dew, p), p, color = 'g')
#plt.plot(skew(np.array(myCAPE[1])-273.15, p), p, color = 'gray')
#plt.title('CAPE using $T$ = ' + str(round(myCAPE[0], 1)) + '$Jkg^{-1}$')
#plt.xlabel('Temperature ($^o$C)')
#plt.ylabel('Pressure (hPa)')
#plt.show()
# plt.savefig('/home/xb899100/bin/test_skewt.png', bbox_tight = True)

