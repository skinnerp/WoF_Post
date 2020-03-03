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
#from news_e_post_cbook import *

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
parser.add_option("-r", dest="qpe_dir", type="string", default= None, help="Input Directory (of qpe_qc files)")
parser.add_option("-f", dest="flag_dir", type="string", default= None, help="Input Directory (of flag files)")
parser.add_option("-i", dest="image_dir", type="string", help= "Output directory (for images)")
parser.add_option("-a", dest="date", type="string", help= "Day of forecast case (YYYYMMDD)")
parser.add_option("-e", dest="fcst_nt", type="int", help= "Total forecast time steps")
 
(options, args) = parser.parse_args()

if ((options.qpe_dir == None) or (options.flag_dir == None) or (options.image_dir == None) or (options.date == None) or (options.fcst_nt == None)):
   print
   parser.print_help()
   print
   sys.exit(1)
else:
   qpe_dir = options.qpe_dir
   flag_dir = options.flag_dir
   image_dir = options.image_dir
   date = options.date
   fcst_nt = options.fcst_nt

pool = Pool(processes=(1))              # set up a queue to run

#time.sleep(120)                         # give time to finish writing some of the ensemble files

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

      qpeqc_files_temp = os.listdir(qpe_dir)
      qpeqc_files_temp.sort()

      for f, file in enumerate(qpeqc_files_temp):
         if ((file[-28:-25] == 'QPE') and (file[-24:-22] == '00')):
            init_hour = file[-12:-10]
            init_min = file[-10:-8]
            init_hhmm = init_hour + init_min

         if ((file[-28:-25] == 'QPE') and (file[-24:-22] == str_t) and (process[t] == 0)):
            qpe_hhmm = file[-7:-3]
            process[t] = 1
            cmd = "python news_e_verif_obj.py -r %s -i %s -a %s -v %s -t %d -n %d " % (qpe_dir, image_dir, date, 'rain1', t, fcst_nt)
            pool.apply_async(run_script, (cmd,)) 
            cmd = "python news_e_verif_obj.py -r %s -i %s -a %s -v %s -t %d -n %d " % (qpe_dir, image_dir, date, 'rain2', t, fcst_nt)
            pool.apply_async(run_script, (cmd,)) 
            cmd = "python news_e_verif_obj.py -r %s -i %s -a %s -v %s -t %d -n %d " % (qpe_dir, image_dir, date, 'rain3', t, fcst_nt)
            pool.apply_async(run_script, (cmd,)) 
            time.sleep(5)
########################################
#            mrms_files_temp = os.listdir(mrms_dir)
#            mrms_files_temp.sort()

#            for m, mfile in enumerate(mrms_files_temp):
#               mrms_hhmm = mfile[-7:-3]
#               if ((mfile[-15:-12] == 'agg') and (ens_hhmm == mrms_hhmm) and (process[t] == 0)):
#                  wofs_file = os.path.join(summary_dir, file)
#                  process[t] = 1
#                  cmd = "python news_e_paintball_qpf_qpe.py -m %s -d %s -o %s -q %s -t %d -n %d " % (mrms_dir, summary_dir, outdir, qpeqc_dir, t, fcst_nt)
   if (process[fcst_nt] == 1):
      break
   else:
      time.sleep(60)
      if (counter >= 390.):
         print("6 hours have passed. Exiting")
         break
#ens_t = 0
#iteration = 0
#
#while (ens_t < fcst_nt):
#   summary_files_temp = os.listdir(summary_dir)
#
#   for f, file in enumerate(summary_files_temp):
#      str_ens_t = str(ens_t)
#      if (len(str_ens_t) == 1):
#         str_ens_t = '0' + str_ens_t
#
#      if ((file[-28:-25] == 'PCP') and (file[-24:-22] == str_ens_t)):
#         init_hour = file[-12:-10]
#         init_min = file[-10:-8]
#         ens_hhmm = file[-7:-3]
#         init_hhmm = init_hour + init_min
#
#         print init_hhmm
##         time.sleep(1)
#
#         mrms_files_temp = os.listdir(mrms_dir)
#         mrms_files_temp.sort()
#
#         for m, mfile in enumerate(mrms_files_temp):
#            mrms_hhmm = mfile[-7:-3]
#            if ((mfile[-15:-12] == 'agg') and (ens_hhmm == mrms_hhmm)):
#               cmd = "python news_e_paintball_qpf_qpe.py -m %s -d %s -o %s -t %d -n %d " % (mrms_dir, summary_dir, outdir, (ens_t), fcst_nt)
#               pool.apply_async(run_script, (cmd,))
# 
#            ens_t = ens_t + 1
##   time.sleep(2)
#
#time.sleep(10)
#
#cmd = "python news_e_paintball_qpf_qpe.py -m %s -d %s -o %s -t %d -n %d " % (mrms_dir, summary_dir, outdir, (ens_t-3), fcst_nt)
#pool.apply_async(run_script, (cmd,))
#
#cmd = "python news_e_paintball_qpf_qpe.py -m %s -d %s -o %s -t %d -n %d " % (mrms_dir, summary_dir, outdir, (ens_t-2), fcst_nt)
#pool.apply_async(run_script, (cmd,))
#
#cmd = "python news_e_paintball_qpf_qpe.py -m %s -d %s -o %s -t %d -n %d " % (mrms_dir, summary_dir, outdir, (ens_t-1), fcst_nt)
#pool.apply_async(run_script, (cmd,))
#
#cmd = "python news_e_paintball_qpf_qpe.py -m %s -d %s -o %s -t %d -n %d " % (mrms_dir, summary_dir, outdir, (ens_t), fcst_nt)
#pool.apply_async(run_script, (cmd,))

time.sleep(10)
print("MRMS QPE-QPF object-matching verification complete for time: %s" % (init_hhmm))

pool.close()
pool.join()


