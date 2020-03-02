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
import datetime
from news_e_post_cbook import *

from multiprocessing import Pool

#sys.path.append("/scratch/software/Anaconda2/bin")

###################################################################################################
# run_script is a function that runs a system command

def run_script(cmd):
    print "Executing command:  " + cmd
    os.system(cmd)
    print cmd + "  is finished...."
    return

###################################################################################################

parser = OptionParser()
parser.add_option("-m", dest="mrms_dir", type="string", default= None, help="Input directory of MRMS Cressman files to plot")
parser.add_option("-d", dest="summary_dir", type="string", default= None, help="Input Directory (of summary files)")
parser.add_option("-o", dest="outdir", type="string", help = "Output Directory (for images)")
parser.add_option("-a", dest="date", type="string", help = "Day of forecast case (YYYYMMDD)")
#parser.add_option("-m", dest="mapname", type="string", default=None, help="Path to Basemap Instance for plotting")
parser.add_option("-e", dest="fcst_nt", type="int", help = "Total number of timesteps in forecast")

(options, args) = parser.parse_args()

if ((options.mrms_dir == None) or (options.summary_dir == None) or (options.outdir == None) or (options.date == None) or (options.fcst_nt == None)):
    print
    parser.print_help()
    print
    sys.exit(1)
else:
   mrms_dir = options.mrms_dir
   summary_dir = options.summary_dir
   outdir = options.outdir
   date = options.date
#    mapname = options.mapname
   fcst_nt = options.fcst_nt

#mapname = '/scratch2/patrick.skinner/images/map.pickle'

pool = Pool(processes=(18))              # set up a queue to run

##### process swath files from completed ens files

fcst_times = np.arange(fcst_nt+1)
process = fcst_times * 0
counter = 0

while (np.min(process) == 0):
   counter = counter + 1
   print 'iteration: ', counter
   for t in fcst_times:
      str_t = str(t)
      if (len(str_t) == 1):
         str_t = '0' + str_t

      summary_files_temp = os.listdir(summary_dir)
      summary_files_temp.sort()

      for f, file in enumerate(summary_files_temp):
         if ((file[-28:-25] == 'ENS') and (file[-24:-22] == '00')):
            init_hour = file[-12:-10]
            init_min = file[-10:-8]
            init_hhmm = init_hour + init_min

         if ((file[-28:-25] == 'ENS') and (file[-24:-22] == str_t) and (process[t] == 0)):
            ens_hhmm = file[-7:-3]
##### look for corresponding MRMS file

            mrms_files = os.listdir(mrms_dir)
            mrms_files.sort()

            for m, mfile in enumerate(mrms_files): 
               mrms_hhmm = mfile[-9:-5] 

               if (ens_hhmm == mrms_hhmm): 
                  wofs_file = os.path.join(summary_dir, file)

                  process[t] = 1

                  cmd = "python mrms_plotter.py -m %s -d %s -o %s -a %s -t %d -s %s" % (mrms_dir, wofs_file, outdir, date, t, init_hhmm)
                  pool.apply_async(run_script, (cmd,))
                  time.sleep(2)

   time.sleep(30)
   if (counter > 760): 
      break

time.sleep(2)

pool.close()
pool.join()

