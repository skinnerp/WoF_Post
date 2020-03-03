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
    print("Executing command:  " + cmd)
    os.system(cmd)
    print(cmd + "  is finished....")
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

pool = Pool(processes=(12))              # set up a queue to run

time.sleep(5)          #Give time for the QPE file to be processed completely
################## Make sure summary file exists to get grid info: ########################

files = []
files_temp = os.listdir(exp_dir)

if os.path.isdir(exp_dir):
   files_temp = os.listdir(exp_dir)
else:
   print("Today's summary file directory not created yet. Trying again in 5 minutes.")
   sys.exit(1)

for f, file in enumerate(files_temp):
   if (file[-28:-25] == 'ENS'):
      files.append(file)

files.sort()

if (len(files) > 0): 
   exp_file = os.path.join(exp_dir, files[0])
else: 
   print('NO SUMMARY FILES YET')
   sys.exit(1)

################## Find files that have already been processed: ########################

qpe_rmp_files = os.listdir(out_dir)
qpe_rmp_files.sort()

################## Find available raw MRMS files: ########################

mrms_qpe_dir = mrms_dir
mrms_qpe_files = os.listdir(mrms_qpe_dir)
mrms_qpe_files.sort()

init_mrms_date = mrms_qpe_files[0][-22:-14]
print('the initial start date for raw QPE data is %s' % (init_mrms_date))
 
################## Loop through MRMS QPE files (ties processing to existence of QPE file): ########################

for q, qpefile in enumerate(mrms_qpe_files): 
   rmp_file_found = 0 

   qpe_date = qpefile[-22:-14]
   qpe_hh = qpefile[-13:-11]
   qpe_mm = qpefile[-11:-9]
   qpe_ss = qpefile[-9:-7]

   qpe_time_seconds = int(qpe_hh) * 3600. + int(qpe_mm) * 60. + int(qpe_ss) 
   if ((qpe_time_seconds < 44099.) and (qpe_date != init_mrms_date)):  #handle times past 00 UTC (assumes no important times before 12 UTC) 
      qpe_time_seconds = qpe_time_seconds + 86400.  

   qpe_out_time_seconds = int(900. * round(float(qpe_time_seconds) / 900.))

   if (qpe_out_time_seconds == 0.):
      qpe_out_time_seconds_day_change = 86400.
   else:
      qpe_out_time_seconds_day_change = -999.

   qpe_process_file = os.path.join(mrms_dir, qpefile)

################## Now check for already remapped files ##########################

   for r, rmpfile in enumerate(qpe_rmp_files):
      
      rmp_date = rmpfile[-22:-14]
      rmp_hh = rmpfile[-13:-11]
      rmp_mm = rmpfile[-11:-9]
      rmp_ss = rmpfile[-9:-7]

      rmp_time_seconds = int(rmp_hh) * 3600. + int(rmp_mm) * 60. + int(rmp_ss)
      if ((rmp_time_seconds < 44099.) and (rmp_date != init_mrms_date)):
         rmp_time_seconds = rmp_time_seconds + 86400.

      if ((rmp_date == qpe_date) and (rmp_time_seconds > (qpe_out_time_seconds - 900.)) and (rmp_time_seconds < (qpe_out_time_seconds + 900.))):
#      if ((rmp_date == qpe_date) and (rmp_time_seconds == qpe_out_time_seconds)):
         print('Remapped QPE file %s found' % (rmpfile))
         rmp_file_found = 1
         break

################## If raw MRMS files exist and output interpolated file does not - Process!: ########################

   if ((rmp_file_found == 0)):
      out_file = qpe_date + '-' + qpe_hh + qpe_mm + '00.nc'
      out_path = os.path.join(out_dir, out_file)
      cmd = "python mrms_remap.py -z %s -o %s -f %s" % (qpe_process_file, out_path, exp_file)
      pool.apply_async(run_script, (cmd,))
   
   elif ((rmp_file_found == 1)):
      pass

time.sleep(10)

pool.close()
pool.join()

