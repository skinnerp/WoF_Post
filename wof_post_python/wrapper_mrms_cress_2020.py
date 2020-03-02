###################################################################################################

from mpl_toolkits.basemap import Basemap
import matplotlib
import math
from scipy import *
import pylab as P
import numpy as np
import sys, glob
import os
import time
from optparse import OptionParser
import netCDF4
from news_e_post_cbook import *

from multiprocessing import Pool

#sys.path.append("/scratch/software/anaconda/bin")

#=======================================================================================================================
# run_script is a function that runs a system command

def run_script(cmd):
    print "Executing command:  " + cmd
    os.system(cmd)
    print cmd + "  is finished...."
    return

###################################################################################################

parser = OptionParser()
parser.add_option("-m", dest="mrms_dir", type="string", default=None, help="Input Directory (raw MRMS files)")
parser.add_option("-d", dest="exp_dir", type="string", default=None, help="Input Directory (summary files)")
parser.add_option("-o", dest="out_dir", type="string", default=None, help="Output Directory")

(options, args) = parser.parse_args()

if ((options.mrms_dir == None) or (options.exp_dir == None) or (options.out_dir == None)):
    print
    parser.print_help()
    print
    sys.exit(1)
else:
    mrms_dir = options.mrms_dir
    exp_dir = options.exp_dir
    out_dir = options.out_dir

pool = Pool(processes=(30))              # set up a queue to run

################## Make sure summary file exists to get grid info: ########################

files = []
files_temp = os.listdir(exp_dir)
for f, file in enumerate(files_temp):
   if (file[-28:-25] == 'ENS'):
      files.append(file)

files.sort()

if (len(files) > 0): 
   exp_file = os.path.join(exp_dir, files[0])
else: 
   print 'NO SUMMARY FILES YET'
   sys.exit(1)

################## Find files that have already been processed: ########################

cress_files = os.listdir(out_dir)
cress_files.sort()

################## Find available raw MRMS files: ########################

mrms_dz_dir = mrms_dir + '/MergedReflectivityQCComposite/00.25'   ###HARDCODED!!!!
mrms_low_dir = mrms_dir + '/MergedAzShear_0-2kmAGL/00.00'
mrms_mid_dir = mrms_dir + '/MergedAzShear_2-5kmAGL/00.00'
mrms_tot_low_dir = mrms_dir + '/MergedVelocity_Gradient_0-2kmAGL/00.00'
mrms_tot_mid_dir = mrms_dir + '/MergedVelocity_Gradient_2-5kmAGL/00.00'
mrms_mesh_dir = mrms_dir + '/MESH/00.25'

mrms_dz_files = os.listdir(mrms_dz_dir)
mrms_low_files = os.listdir(mrms_low_dir)
mrms_mid_files = os.listdir(mrms_mid_dir)
mrms_tot_low_files = os.listdir(mrms_tot_low_dir)
mrms_tot_mid_files = os.listdir(mrms_tot_mid_dir)
mrms_mesh_files = os.listdir(mrms_mesh_dir)

mrms_dz_files.sort()
mrms_low_files.sort()
mrms_mid_files.sort()
mrms_tot_low_files.sort()
mrms_tot_mid_files.sort()
mrms_mesh_files.sort()

init_dz_date = mrms_dz_files[0][-22:-14]

################## Loop through MRMS DZ files (ties processing to existence of DZ file): ########################

for d, dzfile in enumerate(mrms_dz_files): 
   low_file_found = 0 
   mid_file_found = 0 
   tot_low_file_found = 0 
   tot_mid_file_found = 0 
   mesh_file_found = 0 
   cress_file_found = 0 

   dz_date = dzfile[-22:-14]
   dz_hh = dzfile[-13:-11]
   dz_mm = dzfile[-11:-9]
   dz_ss = dzfile[-9:-7]

   dz_time_seconds = int(dz_hh) * 3600. + int(dz_mm) * 60. + int(dz_ss) 
   if ((dz_time_seconds < 43200.) and (dz_date != init_dz_date)):  #handle times past 00 UTC (assumes no important times before 12 UTC) 
      dz_time_seconds = dz_time_seconds + 86400.  

### Get dz_time_seconds from nearest 5 minute interval for checking other files ### 
   dz_out_time_seconds = int(300. * round(float(dz_time_seconds) / 300.))

### Hacky way of handling date change: ###
   if (dz_out_time_seconds == 0.): 
      dz_out_time_seconds_day_change = 86400. 
   else: 
      dz_out_time_seconds_day_change = -999. 
 
   dz_process_file = os.path.join(mrms_dz_dir, dzfile)

################## Check for existence of matching AWS low/mid files: ########################

   for l, lowfile in enumerate(mrms_low_files): #low 
      low_date = lowfile[-22:-14]
      low_hh = lowfile[-13:-11]
      low_mm = lowfile[-11:-9]
      low_ss = lowfile[-9:-7]

      low_time_seconds = int(low_hh) * 3600. + int(low_mm) * 60. + int(low_ss)
      if (low_time_seconds < 43200.):  #handle times past 00 UTC (assumes no important times before 12 UTC) 
         low_time_seconds = low_time_seconds + 86400.

############### If file within 150 seconds of DZ out time, set low_file_found to 1 (true) #############

      if ((low_date == dz_date) and (low_time_seconds > (dz_out_time_seconds - 180.)) and (low_time_seconds < (dz_out_time_seconds + 180.))): 
         low_file_found = 1
         low_process_file = os.path.join(mrms_low_dir, lowfile)
      elif ((low_time_seconds > (dz_out_time_seconds_day_change - 180.)) and (low_time_seconds < (dz_out_time_seconds_day_change + 180.))): 
         low_file_found = 1
         low_process_file = os.path.join(mrms_low_dir, lowfile)

#############

   for m, midfile in enumerate(mrms_mid_files): #mid
      mid_date = midfile[-22:-14]
      mid_hh = midfile[-13:-11]
      mid_mm = midfile[-11:-9]
      mid_ss = midfile[-9:-7]

      mid_time_seconds = int(mid_hh) * 3600. + int(mid_mm) * 60. + int(mid_ss)
      if (mid_time_seconds < 43200.):  #handle times past 00 UTC (assumes no important times before 12 UTC) 
         mid_time_seconds = mid_time_seconds + 86400.

      if ((mid_date == dz_date) and (mid_time_seconds > (dz_out_time_seconds - 180.)) and (mid_time_seconds < (dz_out_time_seconds + 180.))):
         mid_file_found = 1
         mid_process_file = os.path.join(mrms_mid_dir, midfile)
      elif ((mid_time_seconds > (dz_out_time_seconds_day_change - 180.)) and (mid_time_seconds < (dz_out_time_seconds_day_change + 180.))): 
         mid_file_found = 1
         mid_process_file = os.path.join(mrms_mid_dir, midfile)

#############

   for t, totmidfile in enumerate(mrms_tot_mid_files): #mid
      tot_mid_date = totmidfile[-22:-14]
      tot_mid_hh = totmidfile[-13:-11]
      tot_mid_mm = totmidfile[-11:-9]
      tot_mid_ss = totmidfile[-9:-7]

      tot_mid_time_seconds = int(tot_mid_hh) * 3600. + int(tot_mid_mm) * 60. + int(tot_mid_ss)
      if (tot_mid_time_seconds < 43200.):  #handle times past 00 UTC (assumes no important times before 12 UTC) 
         tot_mid_time_seconds = tot_mid_time_seconds + 86400.

      if ((tot_mid_date == dz_date) and (tot_mid_time_seconds > (dz_out_time_seconds - 180.)) and (tot_mid_time_seconds < (dz_out_time_seconds + 180.))):
         tot_mid_file_found = 1
         tot_mid_process_file = os.path.join(mrms_tot_mid_dir, totmidfile)
      elif ((tot_mid_time_seconds > (dz_out_time_seconds_day_change - 180.)) and (tot_mid_time_seconds < (dz_out_time_seconds_day_change + 180.))):
         tot_mid_file_found = 1
         tot_mid_process_file = os.path.join(mrms_tot_mid_dir, totmidfile)

#############

   for t, totlowfile in enumerate(mrms_tot_low_files): #low
      tot_low_date = totlowfile[-22:-14]
      tot_low_hh = totlowfile[-13:-11]
      tot_low_mm = totlowfile[-11:-9]
      tot_low_ss = totlowfile[-9:-7]

      tot_low_time_seconds = int(tot_low_hh) * 3600. + int(tot_low_mm) * 60. + int(tot_low_ss)
      if (tot_low_time_seconds < 43200.):  #handle times past 00 UTC (assumes no important times before 12 UTC) 
         tot_low_time_seconds = tot_low_time_seconds + 86400.

      if ((tot_low_date == dz_date) and (tot_low_time_seconds > (dz_out_time_seconds - 180.)) and (tot_low_time_seconds < (dz_out_time_seconds + 180.))):
         tot_low_file_found = 1
         tot_low_process_file = os.path.join(mrms_tot_low_dir, totlowfile)
      elif ((tot_low_time_seconds > (dz_out_time_seconds_day_change - 180.)) and (tot_low_time_seconds < (dz_out_time_seconds_day_change + 180.))):
         tot_low_file_found = 1
         tot_low_process_file = os.path.join(mrms_tot_low_dir, totlowfile)

#############

   for t, meshfile in enumerate(mrms_mesh_files): #mesh
      mesh_date = meshfile[-22:-14]
      mesh_hh = meshfile[-13:-11]
      mesh_mm = meshfile[-11:-9]
      mesh_ss = meshfile[-9:-7]

      mesh_time_seconds = int(mesh_hh) * 3600. + int(mesh_mm) * 60. + int(mesh_ss)
      if (mesh_time_seconds < 43200.):  #handle times past 00 UTC (assumes no important times before 12 UTC) 
         mesh_time_seconds = mesh_time_seconds + 86400.

      if ((mesh_date == dz_date) and (mesh_time_seconds > (dz_out_time_seconds - 180.)) and (mesh_time_seconds < (dz_out_time_seconds + 180.))):
         mesh_file_found = 1
         mesh_process_file = os.path.join(mrms_mesh_dir, meshfile)
      elif ((mesh_time_seconds > (dz_out_time_seconds_day_change - 180.)) and (mesh_time_seconds < (dz_out_time_seconds_day_change + 180.))):
         mesh_file_found = 1
         mesh_process_file = os.path.join(mrms_mesh_dir, meshfile)

################## Check for existence of interpolated Cressman file: ########################

   for c, cressfile in enumerate(cress_files): #low 
      cress_date = cressfile[-18:-10]
      cress_hh = cressfile[-9:-7]
      cress_mm = cressfile[-7:-5]

      cress_time_seconds = int(cress_hh) * 3600. + int(cress_mm) * 60.
      if (cress_time_seconds < 43200.):  #handle times past 00 UTC (assumes no important times before 12 UTC) 
         cress_time_seconds = cress_time_seconds + 86400.

      if ((cress_date == dz_date) and (cress_time_seconds > (dz_out_time_seconds - 150.)) and (cress_time_seconds < (dz_out_time_seconds + 150.))):
         cress_file_found = 1

################## If all 6 raw MRMS files exist and output Cress file does not - Process!: ########################

#   if ((cress_file_found == 0)): 
   if ((low_file_found == 1) and (mid_file_found == 1) and (tot_low_file_found == 1) and (tot_mid_file_found == 1) and (mesh_file_found == 1) and (cress_file_found == 0) and (dz_time_seconds >= 64800)): 
      out_hh = str(int(np.floor(dz_out_time_seconds / 3600.))) 
      out_mm = str(int((dz_out_time_seconds - (int(out_hh) * 3600.)) / 60.))

      if (int(out_hh) > 23): 
         out_hh = str(int(out_hh) - 24)
      if (len(out_hh) == 1): 
         out_hh = '0' + out_hh
      if (len(out_mm) == 1): 
         out_mm = '0' + out_mm

      out_file = dz_date + '_' + out_hh + out_mm + '00.nc'

      out_path = os.path.join(out_dir, out_file)
      cmd = "/home/louis.wicker/anaconda2/envs/wof-test/bin/python mrms_post_cress_2020.py -a %s -m %s -z %s -t %s -u %s -v %s -o %s -f %s" % (low_process_file, mid_process_file, dz_process_file, tot_low_process_file, tot_mid_process_file, mesh_process_file, out_path, exp_file)
      pool.apply_async(run_script, (cmd,))

time.sleep(10)

pool.close()
pool.join()

