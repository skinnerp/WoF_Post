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

sys.path.append("/scratch/software/Anaconda2/bin")

###################################################################################################
# run_script is a function that runs a system command

def run_script(cmd):
    print("Executing command:  " + cmd)
    os.system(cmd)
    print(cmd + "  is finished....")
    return

###################################################################################################

parser = OptionParser()
parser.add_option("-d", dest="summary_dir", type="string", default= None, help="Input Directory (of summary files)")
parser.add_option("-o", dest="outdir", type="string", help = "Output Directory (for summary file)")
parser.add_option("-m", dest="mrms_dir", type="string", help = "Input Directory for MRMS QPE files")
parser.add_option("-q", dest="qpeqc_dir", type="string", default=None, help = "Output directory for qpe_qc files")
#parser.add_option("-m", dest="mapname", type="string", help = "Path to pickled basemap for plotting")
parser.add_option("-e", dest="fcst_nt", type="int", help = "Total number of timesteps in forecast")
parser.add_option("-p", dest="date", type="string", help = "Date of case")
(options, args) = parser.parse_args()

if ((options.summary_dir == None) or (options.outdir == None) or (options.mrms_dir == None) or (options.qpeqc_dir == None) or (options.fcst_nt == None) or (options.date == None)):
   print
   parser.print_help()
   print
   sys.exit(1)
else:
   summary_dir = options.summary_dir
   outdir = options.outdir
   mrms_dir = options.mrms_dir
   qpeqc_dir = options.qpeqc_dir
#   mapname = options.mapname
   fcst_nt = options.fcst_nt
   date = options.date

pool = Pool(processes=(12))              # set up a queue to run

######### Plot swath products #########
fcst_times = np.arange(fcst_nt+1)
process = fcst_times * 0
counter = 0

while (np.min(process) == 0):
   counter = counter + 1
   print('iteration: ', counter)
   for t in fcst_times:
      str_t = str(t)
      if (len(str_t) == 1):
         str_t = '0' + str_t

      summary_files_temp = os.listdir(summary_dir)
      summary_files_temp.sort()

      for f, file in enumerate(summary_files_temp):
         if ((file[-28:-25] == 'PCP') and (file[-24:-22] == '00')):
            init_hour = file[-12:-10]
            init_min = file[-10:-8]
            init_hhmm = init_hour + init_min

         if ((file[-28:-25] == 'PCP') and (file[-24:-22] == str_t) and (process[t] == 0)):
            ens_hhmm = file[-7:-3]
########################################
            mrms_files_temp = os.listdir(mrms_dir)
            mrms_files_temp.sort()

            for m, mfile in enumerate(mrms_files_temp):
               mrms_hhmm = mfile[-7:-3]
#               time.sleep(2)
               if ((mfile[-15:-12] == 'agg') and (ens_hhmm == mrms_hhmm) and (process[t] == 0)):
                  wofs_file = os.path.join(summary_dir, file)
                  process[t] = 1
#                  time.sleep(2)
                  cmd = "python news_e_paintball_qpf_qpe.py -m %s -d %s -o %s -q %s -t %d -n %d -p %s" % (mrms_dir, summary_dir, outdir, qpeqc_dir, t, fcst_nt, date)
                  pool.apply_async(run_script, (cmd,))
                  time.sleep(5)
   if (process[fcst_nt] == 1):
      break
   else:
      time.sleep(5)
      if (counter >= 4350.):
         print("6 hours have passed. Exiting")
         break

time.sleep(10)
print("MRMS QPE-QPF Paintball verification complete for time: %s" % (init_hhmm))

pool.close()
pool.join()


