###################################################################################################

from mpl_toolkits.basemap import Basemap
import matplotlib
import math
from scipy import *
from scipy import spatial
import pylab as P
import numpy as np
import sys
import os
import time
import netCDF4
from optparse import OptionParser
from news_e_post_cbook import *

####################################### File Variables: ######################################################

parser = OptionParser()
parser.add_option("-d", dest="summary_dir", type="string", default= None, help="Directory of summary files (both input and output dir)")
parser.add_option("-t", dest="t", type="int", help = "Forecast timestep being processed")
parser.add_option("-n", dest="fcst_nt", type="int", help = "Total number of forecast timesteps")

(options, args) = parser.parse_args()

if ((options.summary_dir == None) or (options.t == None) or (options.fcst_nt == None)):
   print
   parser.print_help()
   print
   sys.exit(1)
else:
   summary_dir = options.summary_dir
   t = options.t
   fcst_nt = options.fcst_nt

print 'SWATH STARTED: ', t

############################ Set Thresholds for convolution and neighborhood probabilities: #################################

radius_max_9km             = 3			#grid point radius for maximum value filter (3x3 square neighborhood)
radius_max_15km            = 5			#grid point radius for maximum value filter (9x9 square neighborhood)
radius_max_27km            = 9			#grid point radius for maximum value filter (14x14 square neighborhood)

radius_gauss               = 2			#grid point radius of convolution operator

kernel                     = gauss_kern(radius_gauss)   #convolution operator kernel

neighborhood               = 15                 #15 gridpoint radius for probability matched mean 

comp_dz_thresh             = 40.		#40 dBZ
rain_thresh_low            = 2.	        	#2 in.
rain_thresh_high           = 3.                 #3 in.
ws_80_thresh_ts            = 34.		#58 kts
ws_80_thresh_h1            = 64.		#58 kts
ws_80_thresh_h3            = 96.		#58 kts
wz_0to2_thresh_low         = 0.002		#0.003 s^-1
wz_0to2_thresh_high        = 0.003		#0.003 s^-1


############################ Find WRFOUT files to process: #################################

### Find ENS Summary files ### 

ne = 18

############hack to try and fix realtime for the last timestep: 
#if (t == fcst_nt): 
#   print 'ENTERING WHILE LOOP (NEWS-E-POST-SWATH): ', t, fcst_nt
#   lastfile = 0
#   while (lastfile == 0):  
#      ens_files = []
#      summary_files_temp = os.listdir(summary_dir)
#      for f, file in enumerate(summary_files_temp):
#         if (file[-28:-25] == 'ENS'):                               #assumes filename format of: news-e_ENS_20170516_2200_0000.nc
#            ens_files.append(file)
#      ens_files.sort()
#      tempfile = ens_files[-1]
#      if (tempfile[-24:-22] == '36'): 
#         print 'swathswathswathswathswath', tempfile
#         lastfile == 1
#         time.sleep(10) #give it 10 s to finish writing ... just in case
#      else: 
#         print 'Not Yet, SWATH'
#         time.sleep(10) #try waiting 30 s for final ens. file
#############

ens_files = []
summary_files_temp = os.listdir(summary_dir)

for f, file in enumerate(summary_files_temp): 
   if (file[-28:-25] == 'ENS'):                               #assumes filename format of: news-e_ENS_20170516_2200_0000.nc
      ens_files.append(file) 
      
ens_files.sort() 

for tt in range(0, t+1):
   ens_file = ens_files[tt]
   infile = os.path.join(summary_dir, ens_file)

   if (tt == t): 
   ### Set output path ###
      outname = ens_file[0:7] + 'LTC' + ens_file[10:]
      output_path = summary_dir + outname

   try:                                                 #open WRFOUT file
      fin = netCDF4.Dataset(infile, "r")
      print "Opening %s \n" % infile
   except:
      print "%s does not exist! \n" %infile
      sys.exit(1)

   if (tt == 0):
      dx = fin.DX                                             #east-west grid spacing (m)
      dy = fin.DY                                             #north-south grid spacing (m)
      cen_lat = fin.CEN_LAT                                   #center of domain latitude (dec deg)
      cen_lon = fin.CEN_LON                                   #center of domain longitude (dec deg)
      stand_lon = fin.STAND_LON                               #center lon of Lambert conformal projection
      true_lat_1 = fin.TRUE_LAT1                               #true lat value 1 for Lambert conformal conversion (dec deg)
      true_lat_2 = fin.TRUE_LAT2                               #true lat value 2 for Lambert conformal conversion (dec deg)
      projection = fin.PROJECTION                             #domain projection
      init_time_seconds = fin.INIT_TIME_SECONDS               #initialization time in seconds from 0000 UTC of case
      valid_time_seconds = fin.VALID_TIME_SECONDS             #valid time in seconds from 0000 UTC of case
      forecast_timestep = fin.FORECAST_TIME_STEP              #index of forecast time step

      xlat = fin.variables["xlat"][:,:]                     #latitude (dec deg; Lambert conformal)
      xlon = fin.variables["xlon"][:,:]                     #longitude (dec deg; Lambert conformal)
      hgt = fin.variables["hgt"][:,:]                       #terrain height above MSL (m)

      ny = xlat.shape[0]
      nx = xlat.shape[1]

      wz_0to2                            = fin.variables["wz_0to2"][:,:,:]
      rain                               = fin.variables["rain"][:,:,:]
      comp_dz                            = fin.variables["comp_dz"][:,:,:]
      ws_80                              = fin.variables["ws_80"][:,:,:]

      wz_0to2_hourly                     = fin.variables["wz_0to2"][:,:,:]
      rain_hourly                        = fin.variables["rain"][:,:,:]
      comp_dz_hourly                     = fin.variables["comp_dz"][:,:,:]
      ws_80_hourly                       = fin.variables["ws_80"][:,:,:]

      ### Calc probability matched mean for reflectivity only ###
      temp_dz = np.where(comp_dz > 100000., 0., comp_dz)
      temp_mean_dz = np.mean(temp_dz, axis=0)
      pmm_dz = temp_mean_dz * 0.

      pmm_dz[:,:] = prob_match_mean(temp_dz[:,:,:], temp_mean_dz[:,:], neighborhood)

   else: 

###

      wz_0to2                            = np.where((fin.variables["wz_0to2"][:,:,:] > wz_0to2), fin.variables["wz_0to2"][:,:,:], wz_0to2)
      if ((tt > 1) and ((tt % 12) == 1)):
         wz_0to2_hourly                  = fin.variables["wz_0to2"][:,:,:]
      else:
         wz_0to2_hourly                  = np.where((fin.variables["wz_0to2"][:,:,:] > wz_0to2_hourly), fin.variables["wz_0to2"][:,:,:], wz_0to2_hourly)

###

      comp_dz_temp                  = fin.variables["comp_dz"][:,:,:]

      ### Calc probability matched mean for reflectivity only ###
      temp_dz = np.where(comp_dz_temp > 100000., 0., comp_dz_temp)
      temp_mean_dz = np.mean(temp_dz, axis=0)
      pmm_dz = temp_mean_dz * 0.

      pmm_dz[:,:] = prob_match_mean(temp_dz[:,:,:], temp_mean_dz[:,:], neighborhood)

      comp_dz                            = np.where((comp_dz_temp > comp_dz), comp_dz_temp, comp_dz)
      if ((tt > 1) and ((tt % 12) == 1)):
         comp_dz_hourly                  = comp_dz_temp
      else:
         comp_dz_hourly                  = np.where((comp_dz_temp > comp_dz_hourly), comp_dz_temp, comp_dz_hourly)

###

      rain_temp                          = fin.variables["rain"][:,:,:]
      rain                               = rain + rain_temp
      if ((tt > 1) and ((tt % 12) == 1)):
         rain_hourly                     = fin.variables["rain"][:,:,:]
      else:
         rain_hourly                     = rain_hourly + rain_temp

###

      ws_80                              = np.where((fin.variables["ws_80"][:,:,:] > ws_80), fin.variables["ws_80"][:,:,:], ws_80)
      if ((tt > 1) and ((tt % 12) == 1)):
         ws_80_hourly                    = fin.variables["ws_80"][:,:,:]
      else:
         ws_80_hourly                    = np.where((fin.variables["ws_80"][:,:,:] > ws_80_hourly), fin.variables["ws_80"][:,:,:], ws_80_hourly)

   fin.close()
   del fin

##################### Calculate ensemble output variables: ########################

wz_0to2_90, wz_0to2_max, wz_0to2_9km, wz_0to2_15km, wz_0to2_27km = calc_ens_products(wz_0to2, radius_max_9km, radius_max_15km, radius_max_27km, kernel, wz_0to2_thresh_low)
wz_0to2_90_hourly, wz_0to2_max_hourly, wz_0to2_9km_hourly, wz_0to2_15km_hourly, wz_0to2_27km_hourly = calc_ens_products(wz_0to2_hourly, radius_max_9km, radius_max_15km, radius_max_27km, kernel, wz_0to2_thresh_low)

#wz_0to2_90, wz_0to2_max, wz_0to2_9km, wz_0to2_15km, wz_0to2_27km = calc_ens_products(wz_0to2, radius_max_9km, radius_max_15km, radius_max_27km, kernel, wz_0to2_thresh_high)
#wz_0to2_90_hourly, wz_0to2_max_hourly, wz_0to2_9km_hourly, wz_0to2_15km_hourly, wz_0to2_27km_hourly = calc_ens_products(wz_0to2_hourly, radius_max_9km, radius_max_15km, radius_max_27km, kernel, wz_0to2_thresh_high)

comp_dz_90, comp_dz_max, comp_dz_3km, comp_dz_15km, comp_dz_27km = calc_ens_products(comp_dz, 1, radius_max_15km, radius_max_27km, kernel, comp_dz_thresh)
comp_dz_90_hourly, comp_dz_max_hourly, comp_dz_3km_hourly, comp_dz_15km_hourly, comp_dz_27km_hourly = calc_ens_products(comp_dz_hourly, 1, radius_max_15km, radius_max_27km, kernel, comp_dz_thresh)

rain_90, rain_max, rain_two_3km, rain_two_15km, rain_two_27km = calc_ens_products(rain, 1, radius_max_15km, radius_max_27km, kernel, rain_thresh_low)
rain_90_hourly, rain_max_hourly, rain_two_3km_hourly, rain_two_15km_hourly, rain_two_27km_hourly = calc_ens_products(rain_hourly, 1, radius_max_15km, radius_max_27km, kernel, rain_thresh_low)

rain_three_90, rain_three_max, rain_three_3km, rain_three_15km, rain_three_27km = calc_ens_products(rain, 1, radius_max_15km, radius_max_27km, kernel, rain_thresh_high)
rain_three_90_hourly, rain_three_max_hourly, rain_three_3km_hourly, rain_three_15km_hourly, rain_three_27km_hourly = calc_ens_products(rain_hourly, 1, radius_max_15km, radius_max_27km, kernel, rain_thresh_high)

ws_80_90, ws_80_max, ws_80_ts_9km, ws_80_ts_15km, ws_80_ts_27km = calc_ens_products(ws_80, radius_max_9km, radius_max_15km, radius_max_27km, kernel, ws_80_thresh_ts)
ws_80_90_hourly, ws_80_max_hourly, ws_80_ts_9km_hourly, ws_80_ts_15km_hourly, ws_80_ts_27km_hourly = calc_ens_products(ws_80_hourly, radius_max_9km, radius_max_15km, radius_max_27km, kernel, ws_80_thresh_ts)

ws_80_h1_90, ws_80_h1_max, ws_80_h1_9km, ws_80_h1_15km, ws_80_h1_27km = calc_ens_products(ws_80, radius_max_9km, radius_max_15km, radius_max_27km, kernel, ws_80_thresh_h1)
ws_80_h1_90_hourly, ws_80_h1_max_hourly, ws_80_h1_9km_hourly, ws_80_h1_15km_hourly, ws_80_h1_27km_hourly = calc_ens_products(ws_80_hourly, radius_max_9km, radius_max_15km, radius_max_27km, kernel, ws_80_thresh_h1)

ws_80_h3_90, ws_80_h3_max, ws_80_h3_9km, ws_80_h3_15km, ws_80_h3_27km = calc_ens_products(ws_80, radius_max_9km, radius_max_15km, radius_max_27km, kernel, ws_80_thresh_h3)
ws_80_h3_90_hourly, ws_80_h3_max_hourly, ws_80_h3_9km_hourly, ws_80_h3_15km_hourly, ws_80_h3_27km_hourly = calc_ens_products(ws_80_hourly, radius_max_9km, radius_max_15km, radius_max_27km, kernel, ws_80_thresh_h3)

ws_80_50 = np.percentile(ws_80, 50., axis=0) 
ws_80_50_hourly = np.percentile(ws_80_hourly, 50., axis=0) 

##################### Write Summary File: ########################

##################### Get dimensions and attributes of WRFOUT file ########################

### Create file and dimensions: ###

try:
   fout = netCDF4.Dataset(output_path, "w")
except:
   print "Could not create %s!\n" % output_path

fout.createDimension('NX', nx)
fout.createDimension('NY', ny)

### Set Attributes: ###

setattr(fout,'DX',dx)
setattr(fout,'DY',dy)
setattr(fout,'CEN_LAT',cen_lat)
setattr(fout,'CEN_LON',cen_lon)
setattr(fout,'STAND_LON',stand_lon)
setattr(fout,'TRUE_LAT1',true_lat_1)
setattr(fout,'TRUE_LAT2',true_lat_2)
setattr(fout,'PROJECTION','Lambert Conformal')
setattr(fout,'INIT_TIME_SECONDS',init_time_seconds)
setattr(fout,'VALID_TIME_SECONDS',valid_time_seconds)
setattr(fout,'FORECAST_TIME_STEP',t) 

### Create variables ###

xlat1 = fout.createVariable('xlat', 'f4', ('NY','NX',))
xlat1.long_name = "Latitude"
xlat1.units = "degrees North"

xlon1 = fout.createVariable('xlon', 'f4', ('NY','NX',))
xlon1.long_name = "Longitude"
xlon1.units = "degrees West"

hgt1 = fout.createVariable('hgt', 'f4', ('NY','NX',))
hgt1.long_name = "Height AGL"
hgt1.units = "m"

### 90th percentile and max variables ("p" is placeholder to differentiate variable name from one containing data ) ###

wz_0to2_90p = fout.createVariable('wz_0to2_90', 'f4', ('NY','NX',))
wz_0to2_90p.long_name = "Accumulated ensemble 90th percentile value of average 0-2 km vertical vorticity"
wz_0to2_90p.units = "s^-1"

wz_0to2_90p_hourly = fout.createVariable('wz_0to2_90_hourly', 'f4', ('NY','NX',))
wz_0to2_90p_hourly.long_name = "1-hr Accumulated ensemble 90th percentile value of average 0-2 km vertical vorticity"
wz_0to2_90p_hourly.units = "s^-1"

wz_0to2_maxp = fout.createVariable('wz_0to2_max', 'f4', ('NY','NX',))
wz_0to2_maxp.long_name = "Accumulated ensemble max value of average 0-2 km vertical vorticity"
wz_0to2_maxp.units = "s^-1"

wz_0to2_maxp_hourly = fout.createVariable('wz_0to2_max_hourly', 'f4', ('NY','NX',))
wz_0to2_maxp_hourly.long_name = "1-hr Accumulated ensemble max value of average 0-2 km vertical vorticity"
wz_0to2_maxp_hourly.units = "s^-1"

rain_90p = fout.createVariable('rain_90', 'f4', ('NY','NX',))
rain_90p.long_name = "Accumulated ensemble 90th percentile value of accumulated rainfall"
rain_90p.units = "in"

rain_90p_hourly = fout.createVariable('rain_90_hourly', 'f4', ('NY','NX',))
rain_90p_hourly.long_name = "1-hr Accumulated ensemble 90th percentile value of accumulated rainfall"
rain_90p_hourly.units = "in"

rain_maxp = fout.createVariable('rain_max', 'f4', ('NY','NX',))
rain_maxp.long_name = "Accumulated ensemble max value of accumulated rainfall"
rain_maxp.units = "in"

rain_maxp_hourly = fout.createVariable('rain_max_hourly', 'f4', ('NY','NX',))
rain_maxp_hourly.long_name = "1-hr Accumulated ensemble max value of accumulated rainfall"
rain_maxp_hourly.units = "in"

comp_dz_90p = fout.createVariable('comp_dz_90', 'f4', ('NY','NX',))
comp_dz_90p.long_name = "Accumulated ensemble 90th percentile value of composite reflectivity"
comp_dz_90p.units = "dBZ"

comp_dz_90p_hourly = fout.createVariable('comp_dz_90_hourly', 'f4', ('NY','NX',))
comp_dz_90p_hourly.long_name = "1-hr Accumulated ensemble 90th percentile value of composite reflectivity"
comp_dz_90p_hourly.units = "dBZ"

comp_dz_maxp = fout.createVariable('comp_dz_max', 'f4', ('NY','NX',))
comp_dz_maxp.long_name = "Accumulated ensemble max value of composite reflectivity"
comp_dz_maxp.units = "dBZ"

comp_dz_maxp_hourly = fout.createVariable('comp_dz_max_hourly', 'f4', ('NY','NX',))
comp_dz_maxp_hourly.long_name = "1-hr Accumulated ensemble max value of composite reflectivity"
comp_dz_maxp_hourly.units = "dBZ"

ccomp_dz_pmm = fout.createVariable('comp_dz_pmm', 'f4', ('NY','NX',))
ccomp_dz_pmm.long_name = "Probability matched mean composite reflectivity (93 km neighborhood)"
ccomp_dz_pmm.units = "dBZ"

ws_80_50p = fout.createVariable('ws_80_50', 'f4', ('NY','NX',))
ws_80_50p.long_name = "Accumulated ensemble 50th percentile value of 80-m wind speed"
ws_80_50p.units = "kts"

ws_80_50p_hourly = fout.createVariable('ws_80_50_hourly', 'f4', ('NY','NX',))
ws_80_50p_hourly.long_name = "1-hr Accumulated ensemble 50th percentile value of 80-m wind speed"
ws_80_50p_hourly.units = "kts"

ws_80_90p = fout.createVariable('ws_80_90', 'f4', ('NY','NX',))
ws_80_90p.long_name = "Accumulated ensemble 90th percentile value of 80-m wind speed"
ws_80_90p.units = "kts"

ws_80_90p_hourly = fout.createVariable('ws_80_90_hourly', 'f4', ('NY','NX',))
ws_80_90p_hourly.long_name = "1-hr Accumulated ensemble 90th percentile value of 80-m wind speed"
ws_80_90p_hourly.units = "kts"

ws_80_maxp = fout.createVariable('ws_80_max', 'f4', ('NY','NX',))
ws_80_maxp.long_name = "Accumulated ensemble max value of 80-m wind speed"
ws_80_maxp.units = "kts"

ws_80_maxp_hourly = fout.createVariable('ws_80_max_hourly', 'f4', ('NY','NX',))
ws_80_maxp_hourly.long_name = "1-hr Accumulated ensemble max value of 80-m wind speed"
ws_80_maxp_hourly.units = "kts"

### Probability of exceedence ###

wz_0to2_prob_9 = fout.createVariable('wz_0to2_prob_9km', 'f4', ('NY','NX',))
wz_0to2_prob_9.long_name = "Accumulated ensemble probability of 0-2 km average vertical vorticity greater than 0.002 s^-1 (9 km neighborhod)"
wz_0to2_prob_9.units = "%"

wz_0to2_prob_9_hourly = fout.createVariable('wz_0to2_prob_9km_hourly', 'f4', ('NY','NX',))
wz_0to2_prob_9_hourly.long_name = "1-hr Accumulated ensemble probability of 0-2 km average vertical vorticity greater than 0.002 s^-1 (9 km neighborhod)"
wz_0to2_prob_9_hourly.units = "%"

wz_0to2_prob_15 = fout.createVariable('wz_0to2_prob_15km', 'f4', ('NY','NX',))
wz_0to2_prob_15.long_name = "Accumulated ensemble probability of 0-2 km average vertical vorticity greater than 0.002 s^-1 (15 km neighborhod)"
wz_0to2_prob_15.units = "%"

wz_0to2_prob_15_hourly = fout.createVariable('wz_0to2_prob_15km_hourly', 'f4', ('NY','NX',))
wz_0to2_prob_15_hourly.long_name = "1-hr Accumulated ensemble probability of 0-2 km average vertical vorticity greater than 0.002 s^-1 (15 km neighborhod)"
wz_0to2_prob_15_hourly.units = "%"

wz_0to2_prob_27 = fout.createVariable('wz_0to2_prob_27km', 'f4', ('NY','NX',))
wz_0to2_prob_27.long_name = "Accumulated ensemble probability of 0-2 km average vertical vorticity greater than 0.002 s^-1 (27 km neighborhod)"
wz_0to2_prob_27.units = "%"

wz_0to2_prob_27_hourly = fout.createVariable('wz_0to2_prob_27km_hourly', 'f4', ('NY','NX',))
wz_0to2_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of 0-2 km average vertical vorticity greater than 0.002 s^-1 (27 km neighborhod)"
wz_0to2_prob_27_hourly.units = "%"

rain_two_prob_3 = fout.createVariable('rain_two_prob_3km', 'f4', ('NY','NX',))
rain_two_prob_3.long_name = "Accumulated ensemble gridpoint probability of rainfall greater than 2 inches"
rain_two_prob_3.units = "%"

rain_two_prob_3_hourly = fout.createVariable('rain_two_prob_3km_hourly', 'f4', ('NY','NX',))
rain_two_prob_3_hourly.long_name = "1-hr Accumulated ensemble gridpoint probability of rainfall greater than 2 inches"
rain_two_prob_3_hourly.units = "%"

rain_two_prob_15 = fout.createVariable('rain_two_prob_15km', 'f4', ('NY','NX',))
rain_two_prob_15.long_name = "Accumulated ensemble probability of rainfall greater than 2 inches (15 km neighborhood)"
rain_two_prob_15.units = "%"

rain_two_prob_15_hourly = fout.createVariable('rain_two_prob_15km_hourly', 'f4', ('NY','NX',))
rain_two_prob_15_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than 2 inches (15 km neighborhood)"
rain_two_prob_15_hourly.units = "%"

rain_two_prob_27 = fout.createVariable('rain_two_prob_27km', 'f4', ('NY','NX',))
rain_two_prob_27.long_name = "Accumulated ensemble probability of rainfall greater than 2 inches (27 km neighborhood)"
rain_two_prob_27.units = "%"

rain_two_prob_27_hourly = fout.createVariable('rain_two_prob_27km_hourly', 'f4', ('NY','NX',))
rain_two_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than 2 inches (27 km neighborhood)"
rain_two_prob_27_hourly.units = "%"

rain_three_prob_3 = fout.createVariable('rain_three_prob_3km', 'f4', ('NY','NX',))
rain_three_prob_3.long_name = "Accumulated ensemble gridpoint probability of rainfall greater than 3 inches"
rain_three_prob_3.units = "%"

rain_three_prob_3_hourly = fout.createVariable('rain_three_prob_3km_hourly', 'f4', ('NY','NX',))
rain_three_prob_3_hourly.long_name = "1-hr Accumulated ensemble gridpoint probability of rainfall greater than 3 inches"
rain_three_prob_3_hourly.units = "%"

rain_three_prob_15 = fout.createVariable('rain_three_prob_15km', 'f4', ('NY','NX',))
rain_three_prob_15.long_name = "Accumulated ensemble probability of rainfall greater than 3 inches (15 km neighborhood)"
rain_three_prob_15.units = "%"

rain_three_prob_15_hourly = fout.createVariable('rain_three_prob_15km_hourly', 'f4', ('NY','NX',))
rain_three_prob_15_hourly.long_name = "1-hr Accumulated ensemble gridpoint probability of rainfall greater than 3 inches (15 km neighborhood)"
rain_three_prob_15_hourly.units = "%"

rain_three_prob_27 = fout.createVariable('rain_three_prob_27km', 'f4', ('NY','NX',))
rain_three_prob_27.long_name = "Accumulated ensemble probability of rainfall greater than 3 inches (27 km neighborhood)"
rain_three_prob_27.units = "%"

rain_three_prob_27_hourly = fout.createVariable('rain_three_prob_27km_hourly', 'f4', ('NY','NX',))
rain_three_prob_27_hourly.long_name = "1-hr Accumulated ensemble gridpoint probability of rainfall greater than 3 inches (27 km neighborhood)"
rain_three_prob_27_hourly.units = "%"

comp_dz_prob_3 = fout.createVariable('comp_dz_prob_3km', 'f4', ('NY','NX',))
comp_dz_prob_3.long_name = "Accumulated ensemble gridpoint probability of composite reflectivity greater than 40 dBZ"
comp_dz_prob_3.units = "%"

comp_dz_prob_3_hourly = fout.createVariable('comp_dz_prob_3km_hourly', 'f4', ('NY','NX',))
comp_dz_prob_3_hourly.long_name = "1-hr Accumulated ensemble gridpoint probability of composite reflectivity greater than 40 dBZ"
comp_dz_prob_3_hourly.units = "%"

comp_dz_prob_15 = fout.createVariable('comp_dz_prob_15km', 'f4', ('NY','NX',))
comp_dz_prob_15.long_name = "Accumulated ensemble probability of composite reflectivity greater than 40 dBZ (15 km neighborhood)"
comp_dz_prob_15.units = "%"

comp_dz_prob_15_hourly = fout.createVariable('comp_dz_prob_15km_hourly', 'f4', ('NY','NX',))
comp_dz_prob_15_hourly.long_name = "1-hr Accumulated ensemble probability of composite reflectivity greater than 40 dBZ (15 km neighborhood)"
comp_dz_prob_15_hourly.units = "%"

comp_dz_prob_27 = fout.createVariable('comp_dz_prob_27km', 'f4', ('NY','NX',))
comp_dz_prob_27.long_name = "Accumulated ensemble probability of composite reflectivity greater than 40 dBZ (27 km neighborhood)"
comp_dz_prob_27.units = "%"

comp_dz_prob_27_hourly = fout.createVariable('comp_dz_prob_27km_hourly', 'f4', ('NY','NX',))
comp_dz_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of composite reflectivity greater than 40 dBZ (27 km neighborhood)"
comp_dz_prob_27_hourly.units = "%"

ws_80_prob_9_ts = fout.createVariable('ws_80_prob_9km_ts', 'f4', ('NY','NX',))
ws_80_prob_9_ts.long_name = "Accumulated ensemble probability of 80-m wind speed greater than 34 kts (9 km neighborhod)"
ws_80_prob_9_ts.units = "%"

ws_80_prob_9_ts_hourly = fout.createVariable('ws_80_prob_9km_ts_hourly', 'f4', ('NY','NX',))
ws_80_prob_9_ts_hourly.long_name = "1-hr Accumulated ensemble probability of 80-m wind speed greater than 34 kts (9 km neighborhod)"
ws_80_prob_9_ts_hourly.units = "%"

ws_80_prob_15_ts = fout.createVariable('ws_80_prob_15km_ts', 'f4', ('NY','NX',))
ws_80_prob_15_ts.long_name = "Accumulated ensemble probability of 80-m wind speed greater than 34 kts (15 km neighborhod)"
ws_80_prob_15_ts.units = "%"

ws_80_prob_15_ts_hourly = fout.createVariable('ws_80_prob_15km_ts_hourly', 'f4', ('NY','NX',))
ws_80_prob_15_ts_hourly.long_name = "1-hr Accumulated ensemble probability of 80-m wind speed greater than 34 kts (15 km neighborhod)"
ws_80_prob_15_ts_hourly.units = "%"

ws_80_prob_27_ts = fout.createVariable('ws_80_prob_27km_ts', 'f4', ('NY','NX',))
ws_80_prob_27_ts.long_name = "Accumulated ensemble probability of 80-m wind speed greater than 34 kts (27 km neighborhod)"
ws_80_prob_27_ts.units = "%"

ws_80_prob_27_ts_hourly = fout.createVariable('ws_80_prob_27km_ts_hourly', 'f4', ('NY','NX',))
ws_80_prob_27_ts_hourly.long_name = "1-hr Accumulated ensemble probability of 80-m wind speed greater than 34 kts (27 km neighborhod)"
ws_80_prob_27_ts_hourly.units = "%"

ws_80_prob_9_h1 = fout.createVariable('ws_80_prob_9km_h1', 'f4', ('NY','NX',))
ws_80_prob_9_h1.long_name = "Accumulated ensemble probability of 80-m wind speed greater than 64 kts (9 km neighborhod)"
ws_80_prob_9_h1.units = "%"

ws_80_prob_9_h1_hourly = fout.createVariable('ws_80_prob_9km_h1_hourly', 'f4', ('NY','NX',))
ws_80_prob_9_h1_hourly.long_name = "1-hr Accumulated ensemble probability of 80-m wind speed greater than 64 kts (9 km neighborhod)"
ws_80_prob_9_h1_hourly.units = "%"

ws_80_prob_15_h1 = fout.createVariable('ws_80_prob_15km_h1', 'f4', ('NY','NX',))
ws_80_prob_15_h1.long_name = "Accumulated ensemble probability of 80-m wind speed greater than 64 kts (15 km neighborhod)"
ws_80_prob_15_h1.units = "%"

ws_80_prob_15_h1_hourly = fout.createVariable('ws_80_prob_15km_h1_hourly', 'f4', ('NY','NX',))
ws_80_prob_15_h1_hourly.long_name = "1-hr Accumulated ensemble probability of 80-m wind speed greater than 64 kts (15 km neighborhod)"
ws_80_prob_15_h1_hourly.units = "%"

ws_80_prob_27_h1 = fout.createVariable('ws_80_prob_27km_h1', 'f4', ('NY','NX',))
ws_80_prob_27_h1.long_name = "Accumulated ensemble probability of 80-m wind speed greater than 64 kts (27 km neighborhod)"
ws_80_prob_27_h1.units = "%"

ws_80_prob_27_h1_hourly = fout.createVariable('ws_80_prob_27km_h1_hourly', 'f4', ('NY','NX',))
ws_80_prob_27_h1_hourly.long_name = "1-hr Accumulated ensemble probability of 80-m wind speed greater than 64 kts (27 km neighborhod)"
ws_80_prob_27_h1_hourly.units = "%"

ws_80_prob_9_h3 = fout.createVariable('ws_80_prob_9km_h3', 'f4', ('NY','NX',))
ws_80_prob_9_h3.long_name = "Accumulated ensemble probability of 80-m wind speed greater than 94 kts (9 km neighborhod)"
ws_80_prob_9_h3.units = "%"

ws_80_prob_9_h3_hourly = fout.createVariable('ws_80_prob_9km_h3_hourly', 'f4', ('NY','NX',))
ws_80_prob_9_h3_hourly.long_name = "1-hr Accumulated ensemble probability of 80-m wind speed greater than 94 kts (9 km neighborhod)"
ws_80_prob_9_h3_hourly.units = "%"

ws_80_prob_15_h3 = fout.createVariable('ws_80_prob_15km_h3', 'f4', ('NY','NX',))
ws_80_prob_15_h3.long_name = "Accumulated ensemble probability of 80-m wind speed greater than 94 kts (15 km neighborhod)"
ws_80_prob_15_h3.units = "%"

ws_80_prob_15_h3_hourly = fout.createVariable('ws_80_prob_15km_h3_hourly', 'f4', ('NY','NX',))
ws_80_prob_15_h3_hourly.long_name = "1-hr Accumulated ensemble probability of 80-m wind speed greater than 94 kts (15 km neighborhod)"
ws_80_prob_15_h3_hourly.units = "%"

ws_80_prob_27_h3 = fout.createVariable('ws_80_prob_27km_h3', 'f4', ('NY','NX',))
ws_80_prob_27_h3.long_name = "Accumulated ensemble probability of 80-m wind speed greater than 94 kts (27 km neighborhod)"
ws_80_prob_27_h3.units = "%"

ws_80_prob_27_h3_hourly = fout.createVariable('ws_80_prob_27km_h3_hourly', 'f4', ('NY','NX',))
ws_80_prob_27_h3_hourly.long_name = "1-hr Accumulated ensemble probability of 80-m wind speed greater than 94 kts (27 km neighborhod)"
ws_80_prob_27_h3_hourly.units = "%"

### Write variables ###

fout.variables['xlat'][:] = xlat
fout.variables['xlon'][:] = xlon
fout.variables['hgt'][:] = hgt

fout.variables['wz_0to2_90'][:] = wz_0to2_90
fout.variables['wz_0to2_90_hourly'][:] = wz_0to2_90_hourly
fout.variables['wz_0to2_max'][:] = wz_0to2_max
fout.variables['wz_0to2_max_hourly'][:] = wz_0to2_max_hourly

fout.variables['comp_dz_90'][:] = comp_dz_90
fout.variables['comp_dz_90_hourly'][:] = comp_dz_90_hourly
fout.variables['comp_dz_max'][:] = comp_dz_max
fout.variables['comp_dz_max_hourly'][:] = comp_dz_max_hourly
fout.variables['comp_dz_pmm'][:] = pmm_dz

fout.variables['rain_90'][:] = rain_90
fout.variables['rain_90_hourly'][:] = rain_90_hourly
fout.variables['rain_max'][:] = rain_max
fout.variables['rain_max_hourly'][:] = rain_max_hourly

fout.variables['ws_80_50'][:] = ws_80_50
fout.variables['ws_80_50_hourly'][:] = ws_80_50_hourly
fout.variables['ws_80_90'][:] = ws_80_90
fout.variables['ws_80_90_hourly'][:] = ws_80_90_hourly
fout.variables['ws_80_max'][:] = ws_80_max
fout.variables['ws_80_max_hourly'][:] = ws_80_max_hourly

### 

fout.variables['wz_0to2_prob_9km'][:] = wz_0to2_9km
fout.variables['wz_0to2_prob_9km_hourly'][:] = wz_0to2_9km_hourly
fout.variables['wz_0to2_prob_15km'][:] = wz_0to2_15km
fout.variables['wz_0to2_prob_15km_hourly'][:] = wz_0to2_15km_hourly
fout.variables['wz_0to2_prob_27km'][:] = wz_0to2_27km
fout.variables['wz_0to2_prob_27km_hourly'][:] = wz_0to2_27km_hourly
fout.variables['comp_dz_prob_3km'][:] = comp_dz_3km
fout.variables['comp_dz_prob_3km_hourly'][:] = comp_dz_3km_hourly
fout.variables['comp_dz_prob_15km'][:] = comp_dz_15km
fout.variables['comp_dz_prob_15km_hourly'][:] = comp_dz_15km_hourly
fout.variables['comp_dz_prob_27km'][:] = comp_dz_27km
fout.variables['comp_dz_prob_27km_hourly'][:] = comp_dz_27km_hourly
fout.variables['rain_two_prob_3km'][:] = rain_two_3km
fout.variables['rain_two_prob_3km_hourly'][:] = rain_two_3km_hourly
fout.variables['rain_two_prob_15km'][:] = rain_two_15km
fout.variables['rain_two_prob_15km_hourly'][:] = rain_two_15km_hourly
fout.variables['rain_two_prob_27km'][:] = rain_two_27km
fout.variables['rain_two_prob_27km_hourly'][:] = rain_two_27km_hourly

fout.variables['rain_three_prob_3km'][:] = rain_three_3km
fout.variables['rain_three_prob_3km_hourly'][:] = rain_three_3km_hourly
fout.variables['rain_three_prob_15km'][:] = rain_three_15km
fout.variables['rain_three_prob_15km_hourly'][:] = rain_three_15km_hourly
fout.variables['rain_three_prob_27km'][:] = rain_three_27km
fout.variables['rain_three_prob_27km_hourly'][:] = rain_three_27km_hourly

fout.variables['ws_80_prob_9km_ts'][:] = ws_80_ts_9km
fout.variables['ws_80_prob_15km_ts'][:] = ws_80_ts_15km
fout.variables['ws_80_prob_27km_ts'][:] = ws_80_ts_27km

fout.variables['ws_80_prob_9km_ts_hourly'][:] = ws_80_ts_9km_hourly
fout.variables['ws_80_prob_15km_ts_hourly'][:] = ws_80_ts_15km_hourly
fout.variables['ws_80_prob_27km_ts_hourly'][:] = ws_80_ts_27km_hourly

fout.variables['ws_80_prob_9km_h1'][:] = ws_80_h1_9km
fout.variables['ws_80_prob_15km_h1'][:] = ws_80_h1_15km
fout.variables['ws_80_prob_27km_h1'][:] = ws_80_h1_27km

fout.variables['ws_80_prob_9km_h1_hourly'][:] = ws_80_h1_9km_hourly
fout.variables['ws_80_prob_15km_h1_hourly'][:] = ws_80_h1_15km_hourly
fout.variables['ws_80_prob_27km_h1_hourly'][:] = ws_80_h1_27km_hourly

fout.variables['ws_80_prob_9km_h3'][:] = ws_80_h3_9km
fout.variables['ws_80_prob_15km_h3'][:] = ws_80_h3_15km
fout.variables['ws_80_prob_27km_h3'][:] = ws_80_h3_27km

fout.variables['ws_80_prob_9km_h3_hourly'][:] = ws_80_h3_9km_hourly
fout.variables['ws_80_prob_15km_h3_hourly'][:] = ws_80_h3_15km_hourly
fout.variables['ws_80_prob_27km_h3_hourly'][:] = ws_80_h3_27km_hourly

### Close output file ### 

fout.close()
del fout




