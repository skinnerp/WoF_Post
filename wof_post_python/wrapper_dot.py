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
parser.add_option("-d", dest="summary_dir", type="string", default= None, help="Input Directory (of summary files)")
parser.add_option("-a", dest="asos_dir", type="string", default= None, help="Input Directory (of ASOS files)")
parser.add_option("-i", dest="image_dir", type="string", default=None, help="Image Directory")
#parser.add_option("-m", dest="mapname", type="string", default=None, help="Path to Basemap Instance for plotting")
parser.add_option("-e", dest="fcst_nt", type="int", help = "Total number of timesteps in forecast")

(options, args) = parser.parse_args()

if ((options.summary_dir == None) or (options.asos_dir == None) or (options.image_dir == None) or (options.fcst_nt == None)):
    print
    parser.print_help()
    print
    sys.exit(1)
else:
    summary_dir = options.summary_dir
    asos_dir = options.asos_dir
    image_dir = options.image_dir
#    mapname = options.mapname
    fcst_nt = options.fcst_nt

#mapname = '/scratch2/patrick.skinner/images/map.pickle'

pool = Pool(processes=(6))              # set up a queue to run

######### Get Fcst initialization time from summary file directory: #############

init_hour = int(summary_dir[-5:-3])
init_minute = int(summary_dir[-3:-1])

if (init_hour < 10):
   init_hour = init_hour + 24

init_time_seconds = init_hour * 3600. + init_minute * 60.

######### Make Dot Plots #########

ens_t = 0
prev_ens_t = 0
iteration = 0

while (ens_t < fcst_nt):

#get current time: 
   current = datetime.datetime.now()
   current_hour = current.hour
   current_minute = current.minute
   if (current_hour < 12):      current_hour = current_hour + 24.

#   current_time_seconds = current_hour * 3600. + current_minute * 60. + (3600. * 5.) #convert to UTC from CDT 
   current_time_seconds = 20000000.  #hack to force plotting in retro mode

   summary_files_temp = os.listdir(summary_dir)
   summary_files_temp.sort()
   for f, file in enumerate(summary_files_temp):
      valid_time_seconds = ens_t * 300. + init_time_seconds
      str_ens_t = str(ens_t)
      if (len(str_ens_t) == 1):
         str_ens_t = '0' + str_ens_t

      if ((file[-28:-25] == 'ENV') and (file[-24:-22] == str_ens_t)):
         if (current_time_seconds > (valid_time_seconds + 900.)): #If it is 15 minutes later than the valid time of the forecast 
            print 'STARTING COMMAND (valid, current, ens_t, iter): ', valid_time_seconds, current_time_seconds, ens_t, iteration
            cmd = "/home/louis.wicker/anaconda2/envs/wof-test/bin/python news_e_dot.py -d %s -a %s -o %s -t %d " % (summary_dir, asos_dir, image_dir, ens_t)
            pool.apply_async(run_script, (cmd,))

            ens_t = ens_t + 1
            iteration = 0
   else: 
      time.sleep(30)
  
   if (ens_t == prev_ens_t): 
      time.sleep(10)
      iteration = iteration + 1
      print 'NOTHING HAPPENED, ', ens_t, iteration
 
   if (iteration > 60): 
      print 'NOTHING HAPPENED FOR 30 MINUTES, GIVING UP'
      break

   prev_ens_t = ens_t
   
time.sleep(300)

cmd = "/home/louis.wicker/anaconda2/envs/wof-test/bin/python news_e_dot.py -d %s -a %s -o %s -t %d " % (summary_dir, asos_dir, image_dir, fcst_nt)
pool.apply_async(run_script, (cmd,))

time.sleep(2)

pool.close()
pool.join()

