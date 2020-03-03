#!/usr/bin/env python
#
# "mrms_aggregator.py"
#
# Script to aggregate interpolated MRMS QPE files at various
#	intervals which will be blended with WoFS model QPF
#	to generate QPE->QPF products. 
#
# Script written by Brian Matilla (CIMMS/NSSL)
# Script latest date: 2019 MAY 31
#
# Load modules

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
#from news_e_post_cbook import *
import shutil

# Add parsing options
parser = OptionParser()
parser.add_option("-m", dest="mrms_dir", type="string", default=None, help="Directory of MRMS Remapped QPE Files")
parser.add_option("-o", dest="out_dir", type="string", default=None, help="Output directory of MRMS aggregated files")
parser.add_option("-t", dest="t", type="int", help="Forecast timestep being processed")
parser.add_option("-e", dest="tot_nt", type="int", help="Total number of forecast timesteps")
parser.add_option("-g", dest="date", type="string", help="Date processed")
parser.add_option("-i", dest="timerun", type="string", help="Time HHMM processed")
parser.add_option("-b", dest="blank_dir", type="string", default=None, help="Directory of blank QPE (in case of missing files)")

(options, args) = parser.parse_args()

if ((options.mrms_dir == None) or (options.out_dir == None) or (options.t == None) or (options.tot_nt == None) or (options.date == None) or (options.timerun == None) or (options.blank_dir == None)):
   print
   parser.print_help()
   print
   sys.exit(1)
else:
   mrms_dir = options.mrms_dir
   out_dir = options.out_dir
   t = options.t
   tot_nt = options.tot_nt
   date = options.date
   timerun = options.timerun
   blank_dir = options.blank_dir

#######

####### Find All Remapped MRMS QPE files #######
start_t = 0

qpe_files = []
qpe_files_temp = os.listdir(mrms_dir)

### Tortured way to extract time in seconds from the runtime. but it works...###
runtime = str(timerun)
runtime_hh = runtime[0:2]
runtime_hh = int(runtime_hh)
runtime_hh = runtime_hh * 3600.
runtime_mm = runtime[2:4]
runtime_mm = int(runtime_mm)
runtime_mm = runtime_mm * 60.

runtime_sec = runtime_hh + runtime_mm
if (runtime_sec < 44099.):
   runtime_sec += 86400.

runtime_ts = (900. * t)

runtime_end = int(runtime_sec + runtime_ts)

if (runtime_end < 86400.):
    runtime_end_hh = '{:02}'.format(runtime_end //3600)
    runtime_end_mm = '{:02}'.format(runtime_end %3600 //60)
    runtime_end_hhmm = runtime_end_hh + runtime_end_mm
    
elif runtime_end >= 86400.:
    runtime_end = int(runtime_end - 86400.)
    runtime_end_hh = '{:02}'.format(runtime_end //3600)
    runtime_end_mm = '{:02}'.format(runtime_end %3600 // 60)
    runtime_end_hhmm = runtime_end_hh + runtime_end_mm
###################################################

for q, qpefile in enumerate(qpe_files_temp):
   qpe_date = qpefile[-18:-10]
   qpe_hh = qpefile[-9:-7]
   qpe_mm = qpefile[-7:-5]

   qpe_time_seconds = int(qpe_hh) * 3600. + int(qpe_mm) * 60.

   if ((qpe_time_seconds < 44099.) and (qpe_date == date)):
      qpe_time_seconds += 86400.
   elif ((qpe_time_seconds < 44099.) and (qpe_date != date)):
      qpe_time_seconds += 86400.
#
   if ((qpe_time_seconds >= runtime_sec) and ((qpe_time_seconds - runtime_sec) <= 21600.)):
      qpe_files.append(qpefile)

   orig_date = str(int(date) - 1) # Hack for now.
#
qpe_files.sort()

print(qpe_time_seconds, qpe_files)
#print 'MRMS AGG STARTED. CURRENT TIMESTEP: ', t     

for tt in range(0, t+1):
   qpe_file = qpe_files[tt]
   infile = os.path.join(mrms_dir, qpe_file)

   if (tt == t):
### Set output path ###
      timestep = str(t)
      if (len(timestep) == 1):
         timestep = '0' + timestep
      outname = date + '_' + timerun + '_agg_t' + timestep + '_' + runtime_end_hhmm + '.nc'
#      outname = 'MRMS_AGG_' + timestep + '_' + date + '_' + timerun + '_' + runtime_end_hhmm + '.nc'
      output_path = out_dir + outname

#   time.sleep(2)
   try:
      fin = netCDF4.Dataset(infile, "r")
      print("Opening %s \n" % infile)
   except:
      print("%s does not exist! \n" % infile)
      sys.exit(1)

   if (tt == 0):
      xlat = fin.variables["XLAT"][:,:]
      xlon = fin.variables["XLON"][:,:]
      qpe = fin.variables["MRMS_QPE"][:,:]
      qpe[np.isnan(qpe)] = 0.
      qpe = 0.
      qpe_15m = fin.variables["MRMS_QPE"][:,:]
      qpe_15m[np.isnan(qpe_15m)] = 0.
      qpe_15m = 0.
      qpe_1hr = fin.variables["MRMS_QPE"][:,:]
      qpe_1hr[np.isnan(qpe_1hr)] = 0.
      qpe_1hr = 0.
      qpe_3hr = fin.variables["MRMS_QPE"][:,:]
      qpe_3hr[np.isnan(qpe_3hr)] = 0.
      qpe_3hr = 0.
      qpe_6hr = fin.variables["MRMS_QPE"][:,:]
      qpe_6hr[np.isnan(qpe_6hr)] = 0.
      qpe_6hr = 0.

      ny = xlat.shape[0]
      nx = xlat.shape[1]

   else:
      qpe_15m = fin.variables["MRMS_QPE"][:,:]
      qpe_15m[np.isnan(qpe_15m)] = 0.
      qpe_temp = fin.variables["MRMS_QPE"][:,:]
      qpe_temp[np.isnan(qpe_temp)] = 0.
      qpe = qpe + qpe_temp

#      if ((tt > 0) and ((tt % 1) == 1)):
#         qpe_15min = fin.variables["MRMS_QPE"][:,:]

      if ((tt > 1) and ((tt % 4) == 1)):
         qpe_1hr = fin.variables["MRMS_QPE"][:,:]
      else:
         qpe_1hr = qpe_1hr + qpe_temp

      if ((tt > 1) and ((tt % 12) == 1)):
         qpe_3hr = fin.variables["MRMS_QPE"][:,:]
      else:
         qpe_3hr = qpe_3hr + qpe_temp

      if ((tt > 1) and ((tt % 24) == 1)):
         qpe_6hr = fin.variables["MRMS_QPE"][:,:]
      else:
         qpe_6hr = qpe_6hr + qpe_temp

   fin.close()
   del fin
##################################################################
#time.sleep(2)
try:
    fout = netCDF4.Dataset(output_path, "w")
except:
    print("Could not create %s!\n" % output_path)

fout.createDimension('NX', nx)
fout.createDimension('NY', ny)

xlat1 = fout.createVariable('XLAT', 'f4', ('NY','NX',))
xlat1.long_name = "Latitude"
xlat1.units = "degrees North"

xlon1 = fout.createVariable('XLON', 'f4', ('NY','NX',))
xlon1.long_name = "Longitude"
xlon1.units = "degrees West"

qpe_1 = fout.createVariable('MRMS_QPE_SWT', 'f4', ('NY', 'NX',))
qpe_1.long_name = "MRMS QPE Swath"
qpe_1.units = "inches"

qpe_15m1 = fout.createVariable('MRMS_QPE_15M', 'f4', ('NY', 'NX',))
qpe_15m1.long_name = "15-minute accumulated MRMS QPE"
qpe_15m1.units = "inches"

qpe_1hr1 = fout.createVariable('MRMS_QPE_1HR', 'f4', ('NY', 'NX',))
qpe_1hr1.long_name = "1-hour accumulated MRMS QPE"
qpe_1hr1.units = "inches"

qpe_3hr1 = fout.createVariable('MRMS_QPE_3HR', 'f4', ('NY', 'NX',))
qpe_3hr1.long_name = "3-hour accumulated MRMS QPE"
qpe_3hr1.units = "inches"

qpe_6hr1 = fout.createVariable('MRMS_QPE_6HR', 'f4', ('NY', 'NX',))
qpe_6hr1.long_name = "6-hour accumulated MRMS QPE"
qpe_6hr1.units = "inches"

fout.variables['XLAT'][:] = xlat
fout.variables['XLON'][:] = xlon
fout.variables['MRMS_QPE_SWT'][:] = qpe
fout.variables['MRMS_QPE_15M'][:] = qpe_15m
fout.variables['MRMS_QPE_1HR'][:] = qpe_1hr
fout.variables['MRMS_QPE_3HR'][:] = qpe_3hr
fout.variables['MRMS_QPE_6HR'][:] = qpe_6hr

fout.close()
del fout
