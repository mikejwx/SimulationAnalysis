"""
Script containing all of the nc keys for variables
"""
# UM Model Keys
u_key           = u'STASH_m01s00i002' # U COMPNT OF WIND AFTER TIMESTEP
v_key           = u'STASH_m01s00i003' # V COMPNT OF WIND AFTER TIMESTEP
theta_key       = u'STASH_m01s00i004' # THETA AFTER TIMESTEP
q_key           = u'STASH_m01s00i010' # SPECIFIC HUMIDITY AFTER TIMESTEP
w_key           = u'STASH_m01s00i150' # W COMPNT OF WIND AFTER TIMESTEP
rho_key         = u'STASH_m01s00i389' # DRY RHO AFTER TIMESTEP
mv_key          = u'STASH_m01s00i391' # VAPOUR MIXING RATIO (mv) AFTER TS
mcl_key         = u'STASH_m01s00i392' # CLD LIQ MIXING RATIO (mcl) AFTER TS
mci_key         = u'STASH_m01s00i393' # CLD ICE MIXING RATIO (mcf) AFTER TS
mr_key          = u'STASH_m01s00i394' # RAIN MIXING RATIO (mr) AFTER TS
mg_key          = u'STASH_m01s00i395' # GRAUPEL MIXING RATIO (mg) AFTER TS
prho_key        = u'STASH_m01s00i407' # PRESSURE AT RHO LEVELS AFTER TS
pthe_key        = u'STASH_m01s00i408' # PRESSURE AT THETA LEVELS AFTER TS
tinc_bl_key     = u'STASH_m01s03i181' # 
qinc_bl_key     = u'STASH_m01s03i182' # 
shf_key         = u'STASH_m01s03i216' # BOUNDARY LAYER HEAT FLUXES W/M2
sshf_key        = u'STASH_m01s03i217' # SURFACE SENSIBLE HEAT FLUX W/M2
lhf_key         = u'STASH_m01s03i222' # B.LAYER TOTAL MOISTURE FLUXS KG/M2/S
slhf_key        = u'STASH_m01s03i234' # SURFACE LATENT HEAT FLUX W/M2
tinc_rain_key   = u'STASH_m01s04i181' # TEMPERATURE INCR: ls rain
qinc_rain_key   = u'STASH_m01s04i182' # SPEC HUMIDITY INC: ls rain
ls_rain_amt_key = u'STASH_m01s04i201' # LARGE SCALE RAIN AMOUNT KG/M2/TS
tinc_bl_cld_key = u'STASH_m01s09i181' # TEMPERATURE INC: bdy layer + ls cld
qinc_bl_cld_key = u'STASH_m01s09i182' # SPEC HUMIDITYINC: bdy layer + ls cld
tinc_adv_key    = u'STASH_m01s12i181' # TEMPERATURE INCR: advect K/Timestep
qinc_adv_key    = u'STASH_m01s12i182' # Q INCR: advect kg/kg/timestep
tinc_diff_key   = u'STASH_m01s13i181' # dT DIFFUSION INC ON MODEL LEVELS
qinc_diff_key   = u'STASH_m01s13i182' # dQ DIFFUSION INC ON MODEL LEVELS
smag_m_key      = u'STASH_m01s13i190' # 
smag_h_key      = u'STASH_m01s13i191' # 
smag_s_key      = u'STASH_m01s13i192' # 
tinc_qt_key     = u'STASH_m01s15i181' # dT INC FROM QT_BAL_CLD
qinc_qt_key     = u'STASH_m01s15i182' # DQ INC FROM QT_BAL_CLD
temp_key        = u'STASH_m01s16i004' # TEMPERATURE ON THETA LEVELS
tinc_total_key  = u'STASH_m01s30i181' # T TOTAL INCREMENT ON MODEL LEVELS
qinc_total_key  = u'STASH_m01s30i182' # Q TOTAL INCREMENT ON MODEL LEVELS
lwp_key         = u'STASH_m01s30i405' # TOTAL COLUMN QCL RHO GRID
iwp_key         = u'STASH_m01s30i406' # TOTAL COLUMN QCF RHO GRID
tinc_ideal_key  = u'STASH_m01s53i181' # Temperatue Inc: Idealised (K/step)
qinc_ideal_key  = u'STASH_m01s53i182' # water vapour Inc: Ideal (kg/kg/step)

# My keys
s_key          = u'Flow-Parallel'              # Component of wind parallel to horizontal mean 0-500 m wind direction
n_key          = u'Flow-Perpendicular'         # Component of wind perpendicular to horizontal mean 0-500 m wind direction
np_key         = u'n-prime'                    # Component of wind in the plane of the 0-500 m mean wind direction
zi_new_key     = u'new boundary layer depth'   # Parcel method mixed layer depth
zi_old_key     = u'boundary layer depth'       # Richardson number method mixed layer depth
lcl_key        = u'lifting condensation level' # Lifting condensation level from surface parcel
ctz_key        = u'cloud top height'           # Cloud top height as highest level containing non-zero mcl
theta_anom_key = u'STASH_m01s00i004_anomaly'   # Theta anomaly from horizontal mean at each time

