#!/usr/bin/python
###################################################################

import pylab as P
import numpy as np
import sys, glob
import os
import time
from optparse import OptionParser
import netCDF4
from news_e_post_cbook import *

from multiprocessing import Pool

#=======================================================================================================================
# run_script is a function that runs a system command

def run_script(cmd):
    print "Executing command:  " + cmd
    os.system(cmd)
    print cmd + "  is finished...."
    return

########################################################################################################################

parser = OptionParser()
parser.add_option("-m", dest="mrms_dir", type="string", default=None, help="Input Directory of MRMS Remapped QPE files")
parser.add_option("-o", dest="out_dir", type="string", default=None, help="Output Directory With Aggregated Files")
parser.add_option("-a", dest="date", type="string", default=None, help="Day of forecast case (YYYYMMDD)")
(options, args) = parser.parse_args()

if ((options.mrms_dir == None) or (options.out_dir == None) or (options.date == None)):
   print
   parser.print_help()
   print
   sys.exit(1)

else:
   mrms_dir = options.mrms_dir
   out_dir = options.out_dir
   date = options.date

pool = Pool(processes=(1))

####
#### Find files that have already been aggregated ####

agg_files = os.listdir(out_dir)
agg_files.sort()

#### Find MRMS files that have already been remapped ####

qpe_files = os.listdir(mrms_dir)
qpe_files.sort()

print 'Checking for available MRMS remapped files.'

if (len(qpe_files) == 0):
   print 'No MRMS files remapped yet.'
   sys.exit(1)
else:
   print 'We have files! Currently, there are %s files in the directory' % (len(qpe_files))

init_qpe_date = qpe_files[0][-18:-10]

print 'the initial start date for QPE data is %s' % (init_qpe_date)
#### 

for q, qpefile in enumerate(qpe_files):
   agg_file_found = 0

   qpe_date = qpefile[-18:-10]
   qpe_hh = qpefile[-9:-7]
   qpe_mm = qpefile[-7:-5]
   qpe_ss = qpefile[-5:-3]

   qpe_time_seconds = int(qpe_hh) * 3600. + int(qpe_mm) * 60. + int(qpe_ss)
   if ((qpe_time_seconds < 43200.) and (qpe_date != init_qpe_date)):
      qpe_time_seconds = qpe_time_seconds + 86400.

   qpe_out_time_seconds = int(900. * round(float(qpe_time_seconds) / 900.))

   if (qpe_out_time_seconds == 0.):
      qpe_out_time_seconds_day_change = 86400.
   else:
      qpe_out_time_seconds_day_change = -999.

   qpe_process_file = os.path.join(mrms_dir, qpefile)

### Now check for aggregated files ###

   for p, aggfile in enumerate(agg_files):
      
      agg_date = aggfile[-22:-14]
      agg_hh = aggfile[-13:-11]
      agg_mm = aggfile[-11:-9]
      agg_ss = aggfile[-9:-7]

      agg_time_seconds = int(agg_hh) * 3600. + int(agg_mm) * 60. + int(agg_ss)
      if ((agg_time_seconds < 43200.) and (agg_date != init_qpe_date)):
         agg_time_seconds = agg_time_seconds + 86400.

      if ((agg_date == qpe_date) and (agg_time_seconds > (qpe_out_time_seconds - 900.)) and (agg_time_seconds < (qpe_out_time_seconds + 900.))):
#      if ((agg_date == qpe_date) and (agg_time_seconds == qpe_out_time_seconds)):
         print 'Aggregated file %s found' % (aggfile)
         agg_file_found = 1
         break

### If the MRMS QPE Remap file is present but has not been aggregated, what are we waiting for?!? ### 

   if (agg_file_found == 0):
      out_file = qpe_date + '-' + qpe_hh + qpe_mm + qpe_ss + '_agg.nc'
      out_path = os.path.join(out_dir, out_file)
      cmd = "python mrms_aggregator.py -m %s -z %s -o %s " % (mrms_dir, qpe_process_file, out_path)
      pool.apply_async(run_script, (cmd,))
   elif (agg_file_found == 1):
      pass

time.sleep(2)

pool.close()
pool.join()

#### 




      

