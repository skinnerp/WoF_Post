#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
from scipy import signal
from scipy import *
from scipy import ndimage
import math
from math import radians, tan, sin, cos, pi, atan, sqrt, pow, asin, acos
import pylab as P
import numpy as np
from numpy import NAN
import sys
import netCDF4
from optparse import OptionParser
import pickle as pl
from netcdftime import utime
import os
import time
import ctables
import news_e_post_cbook
from news_e_post_cbook import *
import news_e_plotting_cbook_v2
from news_e_plotting_cbook_v2 import *

####################################### File Variables: ######################################################

parser = OptionParser()
parser.add_option("-m", dest="mrms_dir", type="string", default= None, help="Input directory of MRMS QPE files to plot")
parser.add_option("-d", dest="summary_file", type="string", default= None, help="Input summary file (for grid info)")
parser.add_option("-o", dest="outdir", type="string", help = "Output Directory (for images)")
parser.add_option("-a", dest="date", type="string", help = "Day of forecast case (YYYYMMDD)")
parser.add_option("-v", dest="name", type="string", help = "QPE variable (accum) to process")
parser.add_option("-t", dest="t", type="int", help = "Timestep to process")
parser.add_option("-p", dest="init_hhmm", type="string", help= "Initialization Time from summary file.")

(options, args) = parser.parse_args()

if ((options.mrms_dir == None) or (options.summary_file == None) or (options.outdir == None) or (options.name == None) or (options.date == None) or (options.t == None) or (options.init_hhmm == None)):
   print
   parser.print_help()
   print
   sys.exit(1)
else:
   mrms_dir = options.mrms_dir
   summary_file = options.summary_file
   outdir = options.outdir
   date = options.date
   name = options.name
   t = options.t
   init_hhmm = options.init_hhmm

#################################### User-Defined Variables:  #####################################################

domain                     = 'full'        #vestigial variable that's still needed in the plotting subroutines ... 
edge                       = 7 		#number of grid points to remove from near domain boundaries

radius_max                 = 3                  #grid point radius for maximum value filter (3x3 square neighborhood)
radius_gauss               = 2                  #grid point radius of convolution operator

kernel                     = gauss_kern(radius_gauss)   #convolution operator kernel

neighborhood               = 15                 #15 gridpoint radius for probability matched mean 

plot_alpha          	   = 0.55		#transparency value for filled contour plots

#################################### Basemap Variables:  #####################################################

resolution 	= 'h'
area_thresh 	= 1000.

damage_files = '' #['/Volumes/fast_scr/pythonletkf/vortex_se/2013-11-17/shapefiles/extractDamage_11-17/extractDamagePaths']

#################################### Contour Levels:  #####################################################

#dz_levels_nws       	= np.arange(20.0,80.,5.)		#(dBZ)
qpe_inst_levels		= [0.01, 0.05, 0.10, 0.15, 0.25, 0.35, 0.50, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]
qpe_hrs_levels		= [0.01, 0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0]

qpe_colors_xparent	= [cb_colors.xparent,cb_colors.xparent, cb_colors.xparent,cb_colors.xparent,cb_colors.xparent,cb_colors.xparent,cb_colors.xparent,cb_colors.xparent,cb_colors.xparent,cb_colors.xparent,cb_colors.xparent,cb_colors.xparent,cb_colors.xparent]

if (name == 'rain'):
   var_label       = 'Accumulated Rainfall'
   var_units       = 'inches'
#elif (name == 'MRMS_QPE_1HR'):
#   var_label       = '1-hr Accum.'
#   var_units       = 'inches'
#elif (name == 'MRMS_QPE_3HR'):
#   var_label       = '3-hr Accum.'
#   var_units       = 'inches'
#elif (name == 'MRMS_QPE_6HR'):
#   var_label       = '6-hr Accum.'
#   var_units       = 'inches'
 
#################################### Initialize plot attributes using 'web plot' objects:  #####################################################

mrmsqpe_inst_plot = web_mrms_plot('mrms_qpe',                   \
                   'MRMS Quantitative Precipitation Estimate (in)',                  \
                   qpe_inst_levels,      \
                   qpe_inst_levels,            \
          		   qpe_colors_xparent,			\
         		   cb_colors.red9,		\
         		   'none',				\
                   cb_colors.rain_cmap,              \
                   'max',                \
                   plot_alpha,               \
                   neighborhood)

mrmsqpe_agg_plot = web_mrms_plot('mrms_qpe_agg',                   \
                   'MRMS Quantitative Precipitation Estimate (in)',                  \
                   qpe_hrs_levels,      \
                   qpe_hrs_levels,            \
                   qpe_colors_xparent,                  \
                   cb_colors.red9,              \
                   'none',                              \
                   cb_colors.rain_cmap,              \
                   'max',                \
                   plot_alpha,               \
                   neighborhood)
############################ Find WRFOUT files to process: #################################

try:                                                 #open WRFOUT file
   fin = netCDF4.Dataset(summary_file, "r")
   print("Opening %s \n" % summary_file)
except:
   print("%s does not exist! \n" %summary_file)
   sys.exit(1)

######################### Read Attributes: ####################################

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

######################### Set domain: ####################################

xlat = fin.variables["xlat"][edge:-edge,edge:-edge]                     #latitude (dec deg; Lambert conformal)
xlon = fin.variables["xlon"][edge:-edge,edge:-edge]                     #longitude (dec deg; Lambert conformal)

sw_lat_full = xlat[0,0]
sw_lon_full = xlon[0,0]
ne_lat_full = xlat[-1,-1]
ne_lon_full = xlon[-1,-1]

######################### Read 2-5 km UH and composite reflectivity: ####################################

#      uh_2to5                            = fin.variables["uh_2to5"][:,edge:-edge,edge:-edge]
#      comp_dz                            = fin.variables["comp_dz"][:,edge:-edge,edge:-edge]
#   else:
#      uh_2to5                            = np.where((fin.variables["uh_2to5"][:,edge:-edge,edge:-edge] > uh_2to5), fin.variables["uh_2to5"][:,edge:-edge,edge:-edge], uh_2to5)
#      comp_dz                            = fin.variables["comp_dz"][:,edge:-edge,edge:-edge]

######################### Parse Init time from filename: ####################################

year = summary_file[-21:-17]
month = summary_file[-17:-15]
day = summary_file[-15:-13]
valid_hour = summary_file[-7:-5]
valid_min = summary_file[-5:-3]
init_hour = summary_file[-12:-10]
init_min = summary_file[-10:-8]
#init_hour = init_hhmm[0:2]
#init_min = init_hhmm[2:4]

if ((int(valid_hour) < 12) and (int(init_hour) > 18)):
   temp_day = int(day) + 1
   valid_day = str(temp_day)
   if (len(valid_day) == 1):
      valid_day = '0' + valid_day
else:
   valid_day = day

init_label = 'Init: ' + year + '-' + month + '-' + day + ', ' + init_hour + init_min + ' UTC'
valid_label = 'Valid: ' + year + '-' + month + '-' + valid_day + ', ' + valid_hour + valid_min + ' UTC'

init_time_seconds = int(init_hour) * 3600. + int(init_min) * 60. 
if (init_time_seconds <= 44100.): 
   init_time_seconds = init_time_seconds + 86400. 
#
valid_time_seconds = int(valid_hour) * 3600. + int(valid_min) * 60.

#if (valid_time_seconds > 44099
if (valid_time_seconds <= 44100.):
    valid_time_seconds = valid_time_seconds + 86400.
else:
    valid_time_seconds = valid_time_seconds
######### Parse valid time from t #########
#
#valid_time_seconds = init_time_seconds + (t * 300.)  #assumes 15-min output timesteps (unlike MRMS DZ/AZ at 5 min)
#
#valid_hour = np.floor(valid_time_seconds / 3600.)
#
#valid_min = np.floor((valid_time_seconds - (valid_hour * 3600.)) / 60.)
#
#if (valid_hour > 23.): 
#   valid_hour = str(int(valid_hour - 24))
#   if (int(init_hour) <= 23): 
#      valid_day = int(day) + 1
#   else: 
#      valid_day = day
#else: 
#   valid_day = day
#   valid_hour = str(int(valid_hour))
#
#valid_day = str(valid_day)
#valid_min = str(int(valid_min))
#
#if (len(valid_day) == 1): 
#   valid_day = '0' + valid_day
#
#if (len(valid_hour) == 1): 
#   valid_hour = '0' + valid_hour
#
#if (len(valid_min) == 1): 
#   valid_min = '0' + valid_min
#
#valid_label = 'Valid: ' + year + '-' + month + '-' + valid_day + ', ' + valid_hour + valid_min + ' UTC'

fin.close()
del fin

#################################### Find and process MRMS Files:  #####################################################

mrms_files = os.listdir(mrms_dir)
mrms_files.sort()

for f, file in enumerate(mrms_files):
   temp_time = int(file[22:24]) * 3600. + int(file[24:26]) * 60.
   temp_date = file[0:8]
   if (temp_time <= 44100.):
#      if (temp_date != date):
      temp_time = temp_time + 86400.
#      else:
#         temp_time = 0.
#   if (temp_time == init_time_seconds): 
#      temp_file = os.path.join(mrms_dir, file)
#      try:                                                 #open WRFOUT file
#         fin = netCDF4.Dataset(temp_file, "r")
#         print "Opening %s \n" % temp_file
#      except:
#         print "%s does not exist! \n" %temp_file
#         sys.exit(1)
#
##      qpe_15m = fin.variables["MRMS_QPE_15M"][edge:-edge,edge:-edge]
##      qpe_1hr = fin.variables["MRMS_QPE_1HR"][edge:-edge,edge:-edge]
##      qpe_3hr = fin.variables["MRMS_QPE_3HR"][edge:-edge,edge:-edge]
##      qpe_6hr = fin.variables["MRMS_QPE_6HR"][edge:-edge,edge:-edge]
# 
#   if ((temp_time > init_time_seconds) and (temp_time < valid_time_seconds)):
#      temp_file = os.path.join(mrms_dir, file)
#      try:                                                 #open WRFOUT file
#         fin = netCDF4.Dataset(temp_file, "r")
#         print "Opening %s \n" % temp_file
#      except:
#         print "%s does not exist! \n" %temp_file
#         sys.exit(1)
#
##      qpe_15m = fin.variables["MRMS_QPE_15M"][edge:-edge,edge:-edge]
##      qpe_1hr = fin.variables["MRMS_QPE_1HR"][edge:-edge,edge:-edge]
##      qpe_3hr = fin.variables["MRMS_QPE_3HR"][edge:-edge,edge:-edge]
##      qpe_6hr = fin.variables["MRMS_QPE_6HR"][edge:-edge,edge:-edge]
#      
#   print(temp_time, valid_time_seconds)

   if (temp_time == valid_time_seconds):
      temp_file = os.path.join(mrms_dir, file)

      try:                                                 #open WRFOUT file
         fin = netCDF4.Dataset(temp_file, "r")
         print("Opening %s \n" % temp_file)
      except:
         print("%s does not exist! \n" %temp_file)
         sys.exit(1)

#      qpe_15m = fin.variables["MRMS_QPE_15M"][edge:-edge,edge:-edge]
#      qpe_1hr = fin.variables["MRMS_QPE_1HR"][edge:-edge,edge:-edge]
#      qpe_3hr = fin.variables["MRMS_QPE_3HR"][edge:-edge,edge:-edge]
#      qpe_6hr = fin.variables["MRMS_QPE_6HR"][edge:-edge,edge:-edge]

#print(temp_time, valid_time_seconds)

print('basemap part')

map, fig, ax1, ax2, ax3 = create_fig(sw_lat_full, sw_lon_full, ne_lat_full, ne_lon_full, true_lat_1, true_lat_2, cen_lat, stand_lon, damage_files, resolution, area_thresh)

x, y = map(xlon[:], xlat[:])
#xx, yy = map.makegrid(xlat.shape[1], xlat.shape[0], returnxy=True)[2:4]   #equidistant x/y grid for streamline plots

##########################################################################################################
###################################### Make Plots: #######################################################
##########################################################################################################
######### HURRICANE FLORENCE SPECIFIC - Plot warnings/LSRs: ############

########### LSRs ##########

from datetime import datetime, timedelta

daterun = str(date)

daterun_tm = datetime.strptime(daterun, '%Y%m%d') + timedelta(days=1)
daterun_p1 = '%04d%02d%02d' % ((daterun_tm.year), (daterun_tm.month), (daterun_tm.day))

lsr_file = '/scratch/brian.matilla/WoFS_2019/lsr_wwa/'+daterun+'/lsr_'+daterun+'1600_'+daterun_p1+'1200'
warn_file = '/scratch/brian.matilla/WoFS_2019/lsr_wwa/'+daterun+'/wwa_'+daterun+'1600_'+daterun_p1+'1200'

init_seconds = int(init_hour) * 3600. + int(init_min) * 60.
if (init_seconds <= 44100.):
   init_seconds = init_seconds + 86400.

valid_seconds = int(valid_hour) * 3600. + int(valid_min) * 60.
if (valid_seconds <= 44100.):
   valid_seconds = valid_seconds + 86400.

#print b_seconds, e_seconds
hail, wind, tornado = plot_lsr(map, fig, ax1, ax2, ax3, lsr_file, init_seconds, valid_seconds, plot_h='True', plot_w='True', plot_t='True')

########### Warnings #########

svr, tor, ff = plot_warn(map, fig, ax1, ax2, ax3, warn_file, (valid_seconds-300.), valid_seconds, cb_colors.blue6, cb_colors.red6, cb_colors.green6)

print('plot part')

######################## MRMS Plots: #####################

if (name == 'rain'):
   qpe = fin.variables['MRMS_QPE_SWT'][edge:-edge,edge:-edge]
   mrmsqpe_inst_plot.name = 'mrms_qpe_15m'
   mrmsqpe_inst_plot.var1_title = 'MRMS Accumulated Quantitative Precipitation Estimate (in)'
   mrms_qpe_plot(map, fig, ax1, ax2, ax3, x, y, mrmsqpe_agg_plot, qpe[:,:], t, init_label, valid_label, domain, outdir, 5, 0, showmax='True')

   if ((t > 0) and ((t % 12) == 0)):
      qpe_1hr = fin.variables['MRMS_QPE_1HR'][edge:-edge,edge:-edge]
      mrmsqpe_agg_plot.name = 'mrms_qpe_1hr'
      mrmsqpe_agg_plot.var1_title = 'MRMS 1-hr Accum. Quantitative Precipitation Estimate (in)'
      mrms_qpe_plot(map, fig, ax1, ax2, ax3, x, y, mrmsqpe_agg_plot, qpe_1hr[:,:], (t/12), init_label, valid_label, domain, outdir, 1, 0, showmax='True')

   if ((t > 0) and ((t % 36) == 0)):
      qpe_3hr = fin.variables['MRMS_QPE_3HR'][edge:-edge,edge:-edge]
      mrmsqpe_agg_plot.name = 'mrms_qpe_3hr'
      mrmsqpe_agg_plot.var1_title = 'MRMS 3-hr Accum. Quantitative Precipitation Estimate (in)'
      mrms_qpe_plot(map, fig, ax1, ax2, ax3, x, y, mrmsqpe_agg_plot, qpe_3hr[:,:], (t/12), init_label, valid_label, domain, outdir, 1, 0, showmax='True')

   if ((t > 0) and ((t % 72) == 0)):
      qpe_6hr = fin.variables['MRMS_QPE_6HR'][edge:-edge,edge:-edge]
      mrmsqpe_agg_plot.name = 'mrms_qpe_6hr'
      mrmsqpe_agg_plot.var1_title = 'MRMS 6-hr Accum. Quantitative Precipitation Estimate (in)'
      mrms_qpe_plot(map, fig, ax1, ax2, ax3, x, y, mrmsqpe_agg_plot, qpe_6hr[:,:], (t/12), init_label, valid_label, domain, outdir, 1, 0, showmax='True')
fin.close()
del fin
