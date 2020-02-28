#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
from scipy import signal
from scipy import *
from scipy import ndimage
from skimage.morphology import label
from skimage.measure import regionprops
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
parser.add_option("-m", dest="mrms_dir", type="string", default= None, help="Input directory of MRMS Cressman files to plot")
parser.add_option("-d", dest="summary_file", type="string", default= None, help="Input summary file (for grid info)")
parser.add_option("-o", dest="outdir", type="string", help = "Output Directory (for images)")
parser.add_option("-a", dest="date", type="string", help = "Day of forecast case (YYYYMMDD)")
#parser.add_option("-m", dest="mapname", type="string", help = "Path to Pickled Basemap instance")
parser.add_option("-t", dest="t", type="int", help = "Timestep to process")
parser.add_option("-s", dest="init_hhmm", type="string", help = "HHMM of forecast initialization")

(options, args) = parser.parse_args()

if ((options.mrms_dir == None) or (options.summary_file == None) or (options.outdir == None) or (options.date == None) or (options.t == None) or (options.init_hhmm == None)):
   print
   parser.print_help()
   print
   sys.exit(1)
else:
   mrms_dir = options.mrms_dir
   summary_file = options.summary_file
   outdir = options.outdir
   date = options.date
#   mapname = options.mapname
   t = options.t
   init_hhmm = options.init_hhmm

#path to pickled map instance for plotting: 
#mapname = '/scratch2/patrick.skinner/images/map.pickle'

time.sleep(30)  #wait 30 seconds before starting in case MRMS file is still writing

#################################### User-Defined Variables:  #####################################################

domain                     = 'full'        #vestigial variable that's still needed in the plotting subroutines ... 
edge                       = 7 		#number of grid points to remove from near domain boundaries

radius_max                 = 3                  #grid point radius for maximum value filter (3x3 square neighborhood)
radius_gauss               = 2                  #grid point radius of convolution operator

kernel                     = gauss_kern(radius_gauss)   #convolution operator kernel

neighborhood               = 15                 #15 gridpoint radius for probability matched mean 

plot_alpha          	   = 0.6		#transparency value for filled contour plots

#################################### Basemap Variables:  #####################################################

resolution 	= 'h'
area_thresh 	= 1000.

damage_files = '' #['/Volumes/fast_scr/pythonletkf/vortex_se/2013-11-17/shapefiles/extractDamage_11-17/extractDamagePaths']

#################################### Contour Levels:  #####################################################

dz_levels_nws       	= np.arange(20.0,80.,5.)		#(dBZ)

aws_levels              = [0.005, 0.008, 7000.]                  #az shear
uh_2to5_levels 		= [100., 300., 7000.]				#(dBZ) 
pmm_dz_colors_gray	= [cb_colors.gray5, cb_colors.gray8, 'none']	#gray contours

#################################### Initialize plot attributes using 'web plot' objects:  #####################################################

dz_plot = web_plot('',                   \
                   'MRMS Composite Reflectivity (dBZ)',                  \
                   '',                   \
                   cb_colors.gray6,      \
                   dz_levels_nws,            \
                   aws_levels,        \
                   '',                   \
                   '',                   \
                   pmm_dz_colors_gray,        \
                   cb_colors.c75,               \
                   'none',              \
                   cb_colors.nws_dz_cmap,              \
                   'max',                \
                   0.35,               \
                   neighborhood)

############################ Find WRFOUT files to process: #################################

try:                                                 #open WRFOUT file
   fin = netCDF4.Dataset(summary_file, "r")
   print "Opening %s \n" % summary_file
except:
   print "%s does not exist! \n" %summary_file
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
init_hour = init_hhmm[0:2]
init_min = init_hhmm[2:4]

init_label = 'Init: ' + year + '-' + month + '-' + day + ', ' + init_hour + init_min + ' UTC'

init_time_seconds = int(init_hour) * 3600. + int(init_min) * 60. 
if (init_time_seconds < 36000.): 
   init_time_seconds = init_time_seconds + 86400. 

######## Parse valid time from t #########

valid_time_seconds = init_time_seconds + (t * 300.)  #assumes 5-min output timesteps

valid_hour = np.floor(valid_time_seconds / 3600.)

valid_min = np.floor((valid_time_seconds - (valid_hour * 3600.)) / 60.)

if (valid_hour > 23.): 
   valid_hour = str(int(valid_hour - 24))
   if (int(init_hour) <= 23): 
      valid_day = int(day) + 1
   else: 
      valid_day = day
else: 
   valid_day = day
   valid_hour = str(int(valid_hour))

valid_day = str(valid_day)
valid_min = str(int(valid_min))

if (len(valid_day) == 1): 
   valid_day = '0' + valid_day

if (len(valid_hour) == 1): 
   valid_hour = '0' + valid_hour

if (len(valid_min) == 1): 
   valid_min = '0' + valid_min

valid_label = 'Valid: ' + year + '-' + month + '-' + valid_day + ', ' + valid_hour + valid_min + ' UTC'

fin.close()
del fin

#################################### Find and process MRMS Files:  #####################################################

mrms_files = os.listdir(mrms_dir)
mrms_files.sort()

#swath_files = []

for f, file in enumerate(mrms_files):
   temp_time = int(file[9:11]) * 3600. + int(file[11:13]) * 60.
   temp_date = file[0:8]
   if (temp_time < 36000.):
      if (temp_date != date):
         temp_time = temp_time + 86400.
      else:
         temp_time = 0.
   if (temp_time == init_time_seconds): 
      temp_file = os.path.join(mrms_dir, file)

      try:                                                 #open WRFOUT file
         fin = netCDF4.Dataset(temp_file, "r")
         print "Opening %s \n" % temp_file
      except:
         print "%s does not exist! \n" %temp_file
         sys.exit(1)

      aws_mid = fin.variables["MID_CRESSMAN"][edge:-edge,edge:-edge]
      aws_low = fin.variables["LOW_CRESSMAN"][edge:-edge,edge:-edge]
      fin.close()
      del fin
   if ((temp_time > init_time_seconds) and (temp_time < valid_time_seconds)):  
      temp_file = os.path.join(mrms_dir, file)

      try:                                                 #open WRFOUT file
         fin = netCDF4.Dataset(temp_file, "r")
         print "Opening %s \n" % temp_file
      except:
         print "%s does not exist! \n" %temp_file
         sys.exit(1)

      aws_mid = np.where((fin.variables["MID_CRESSMAN"][edge:-edge,edge:-edge] > aws_mid), fin.variables["MID_CRESSMAN"][edge:-edge,edge:-edge], aws_mid)
      aws_low = np.where((fin.variables["LOW_CRESSMAN"][edge:-edge,edge:-edge] > aws_low), fin.variables["LOW_CRESSMAN"][edge:-edge,edge:-edge], aws_low)

      fin.close()
      del fin
   if (temp_time == valid_time_seconds): 
      temp_file = os.path.join(mrms_dir, file)

      try:                                                 #open WRFOUT file
         fin = netCDF4.Dataset(temp_file, "r")
         print "Opening %s \n" % temp_file
      except:
         print "%s does not exist! \n" %temp_file
         sys.exit(1)

      dz = fin.variables["DZ_CRESSMAN"][edge:-edge,edge:-edge]
#      aws_mid = dz * 0. 
#      aws_low = dz * 0. 
      aws_mid = np.where((fin.variables["MID_CRESSMAN"][edge:-edge,edge:-edge] > aws_mid), fin.variables["MID_CRESSMAN"][edge:-edge,edge:-edge], aws_mid)
      aws_low = np.where((fin.variables["LOW_CRESSMAN"][edge:-edge,edge:-edge] > aws_low), fin.variables["LOW_CRESSMAN"][edge:-edge,edge:-edge], aws_low)

      fin.close()
      del fin

#################################### Calc Swath:  #####################################################

print 'swath part'

########################################################################################################
### If updraft helicity, vertical vorticity, wind speed, graupel max or updraft plot, run maximum value and convolution filters
### over the raw data to spread and smooth the data
########################################################################################################

aws_mid_convolve_temp = aws_mid * 0.
aws_mid_convolve = aws_mid * 0.

aws_low_convolve_temp = aws_low * 0.
aws_low_convolve = aws_low * 0.

kernel = gauss_kern(radius_gauss)

aws_mid_convolve_temp = get_local_maxima2d(aws_mid, radius_max)
aws_mid_convolve = signal.convolve2d(aws_mid_convolve_temp, kernel, 'same')

aws_low_convolve_temp = get_local_maxima2d(aws_low, radius_max)
aws_low_convolve = signal.convolve2d(aws_low_convolve_temp, kernel, 'same')

################################# Make Figure Template: ###################################################

print 'basemap part'

#Load pickled basemap instance for faster plotting: 

#fig, ax1, ax2, ax3 = create_fig_nomap()

#map_temp = pl.load(open(mapname, 'rb'))

#P.sca(ax1)
#map = mymap_boundaries(map_temp, damage_files)

map, fig, ax1, ax2, ax3 = create_fig(sw_lat_full, sw_lon_full, ne_lat_full, ne_lon_full, true_lat_1, true_lat_2, cen_lat, stand_lon, damage_files, resolution, area_thresh)

x, y = map(xlon[:], xlat[:])
#xx, yy = map.makegrid(xlat.shape[1], xlat.shape[0], returnxy=True)[2:4]   #equidistant x/y grid for streamline plots

##########################################################################################################
###################################### Make Plots: #######################################################
##########################################################################################################

print 'plot part'

######################## MRMS Plots: #####################

dz_plot.name = 'mrms_mid'
dz_plot.var2_title = 'MRMS 2-5 km Azimuthal Wind Shear (s$^{-1}$)'
mem_plot(map, fig, ax1, ax2, ax3, x, y, dz_plot, dz, aws_mid_convolve, t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

dz_plot.name = 'mrms_low'
dz_plot.var2_title = 'MRMS 0-2 km Azimuthal Wind Shear (s$^{-1}$)'
mem_plot(map, fig, ax1, ax2, ax3, x, y, dz_plot, dz, aws_low_convolve, t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

