#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import scipy
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
from netcdftime import utime
import os
import time
import datetime
from optparse import OptionParser
import pickle as pl
import news_e_post_cbook
from news_e_post_cbook import *
import news_e_plotting_cbook_v2
from news_e_plotting_cbook_v2 import *

####################################### File Variables: ######################################################

parser = OptionParser()
parser.add_option("-d", dest="summary_dir", type="string", default= None, help="Input directory of summary files to plot")
parser.add_option("-a", dest="asos_dir", type="string", default= None, help="Input Directory (of ASOS files)")
parser.add_option("-o", dest="outdir", type="string", help = "Output Directory (for images)")
#parser.add_option("-m", dest="mapname", type="string", help = "Path to Pickled Basemap instance")
parser.add_option("-t", dest="t", type="int", help = "Timestep to process")

(options, args) = parser.parse_args()

if ((options.summary_dir == None) or (options.asos_dir == None) or (options.outdir == None) or (options.t == None)):
   print
   parser.print_help()
   print
   sys.exit(1)
else:
   summary_dir = options.summary_dir
   asos_dir = options.asos_dir
   outdir = options.outdir
   t = options.t

#################################### Get input file:  #####################################################

summary_files_temp = os.listdir(summary_dir)
timestep = str(t)
if (len(timestep) == 1):
   timestep = '0' + timestep

for f, file in enumerate(summary_files_temp):
   if ((file[-28:-25] == 'ENV') and (file[-24:-22] == timestep)):                   #assumes filename format of: news-e_ENS_20170516_2200_0000.nc
      infile = os.path.join(summary_dir, file)

print 'Matched ENV File: ', timestep, '   ', infile

#################################### Basemap Variables:  #####################################################

resolution      = 'h'
area_thresh     = 1000.

damage_files = '' #['/Volumes/fast_scr/pythonletkf/vortex_se/2013-11-17/shapefiles/extractDamage_11-17/extractDamagePaths']

#################################### User-Defined Variables:  #####################################################

domain          = 'full'        #vestigial variable that's still needed in the plotting subroutines ... 
edge            = 7             #number of grid points to remove from near domain boundaries
thin            = 13             #thinning factor for quiver values (e.g. 6 means slice every 6th grid point)

neighborhood    = 15            #grid point radius of prob matched mean neighborhood

plot_alpha      = 0.55           #transparency value for filled contour plots

#################################### Contour Levels:  #####################################################

temp_levels_ugly        = np.arange(50., 110., 5.)
temp_levels             = np.arange(-20., 125., 5.)
td_levels		= np.arange(-16., 88., 4.)		#(deg F)
td_levels_2             = np.arange(12., 80., 4.)             #(deg F) 
td_levels_ugly          = np.arange(32., 80., 4.)             #(deg F) 
#td_levels_ugly          = np.arange(40., 80., 4.)             #(deg F) 
uv_levels              = np.arange(-30.,36.,6.)               #(kts)

tdd_diff_levels          = np.arange(-10.5,12.,1.5)
t_diff_levels          = np.arange(-14.,16.,2.)
ws_diff_levels         = np.arange(-10.5,12.,1.5)

pmm_dz_levels 		= [35., 50.]				#(dBZ) 
pmm_dz_colors_gray	= [cb_colors.gray8, cb_colors.gray8]	#gray contours

#################################### Initialize plot attributes using 'web plot' objects:  #####################################################

#################################### Initialize plot attributes using 'web plot' objects:  #####################################################

temp_plot = web_plot('tdot',                 \
                   'Ens. Mean 2 m Temperature ($^{\circ}$F)',                   \
                   'Probability Matched Mean - Composite Reflectivity (dBZ)',                   \
                   cb_colors.gray6,      \
                   temp_levels_ugly,          \
                   pmm_dz_levels,        \
                   '',                   \
                   '',                   \
                   pmm_dz_colors_gray,        \
                   cb_colors.purple5,           \
                   cb_colors.blue8,           \
                   cb_colors.temp_cmap_ugly,              \
                   'both',               \
                   plot_alpha,               \
                   neighborhood)

td_plot = web_plot('tddot',                   \
                   'Ens. Mean 2 m Dewpoint Temp ($^{\circ}$F)',                 \
                   'Probability Matched Mean - Composite Reflectivity (dBZ)',                   \
                   cb_colors.gray6,      \
                   td_levels_ugly,            \
                   pmm_dz_levels,        \
                   '',                   \
                   '',                   \
                   pmm_dz_colors_gray,        \
                   cb_colors.purple5,           \
                   cb_colors.orange4,           \
                   cb_colors.td_cmap_ugly,              \
                   'both',               \
                   plot_alpha,               \
                   neighborhood)

u_plot = web_plot('udot',                   \
                   'Ens. Mean 10 m U Component of Wind (kts)',                 \
                   'Probability Matched Mean - Composite Reflectivity (dBZ)',                   \
                   cb_colors.gray6,      \
                   uv_levels,            \
                   pmm_dz_levels,        \
                   '',                   \
                   '',                   \
                   pmm_dz_colors_gray,        \
                   cb_colors.orange7,              \
                   cb_colors.purple7,              \
                   cb_colors.uv_cmap,              \
                   'both',               \
                   plot_alpha,               \
                   neighborhood)

v_plot = web_plot('vdot',                   \
                   'Ens. Mean 10 m V Component of Wind (kts)',                 \
                   'Probability Matched Mean - Composite Reflectivity (dBZ)',                   \
                   cb_colors.gray6,      \
                   uv_levels,            \
                   pmm_dz_levels,        \
                   '',                   \
                   '',                   \
                   pmm_dz_colors_gray,        \
                   cb_colors.orange7,              \
                   cb_colors.purple7,              \
                   cb_colors.uv_cmap,              \
                   'both',               \
                   plot_alpha,               \
                   neighborhood)

######################################################################################################
#################################### Read Data:  #####################################################
######################################################################################################

try:                                                 #open WRFOUT file
   fin = netCDF4.Dataset(infile, "r")
   print "Opening %s \n" % infile
except:
   print "%s does not exist! \n" %infile
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

######################### Parse Init/Valid times from filename: ####################################

year = infile[-21:-17]
month = infile[-17:-15]
day = infile[-15:-13]
init_hour = infile[-12:-10]
init_min = infile[-10:-8]
valid_hour = infile[-7:-5]
valid_min = infile[-5:-3]

if ((int(valid_hour) < 12) and (int(init_hour) > 18)):
   temp_day = int(day) + 1
   valid_day = str(temp_day)
   if (len(valid_day) == 1):
      valid_day = '0' + valid_day
else:
   valid_day = day

epoch_seconds = (datetime.datetime(int(year),int(month),int(valid_day),int(valid_hour),int(valid_min)) - datetime.datetime(1970,1,1)).total_seconds()
 
init_label = 'Init: ' + year + '-' + month + '-' + day + ', ' + init_hour + init_min + ' UTC'
valid_label = 'Valid: ' + year + '-' + month + '-' + valid_day + ', ' + valid_hour + valid_min + ' UTC'

######################### Set domain: ####################################

xlat = fin.variables["xlat"][edge:-edge,edge:-edge]                     #latitude (dec deg; Lambert conformal)
xlon = fin.variables["xlon"][edge:-edge,edge:-edge]                     #longitude (dec deg; Lambert conformal)

#xlat = xlat[edge:-edge,edge:-edge]
#xlon = xlon[edge:-edge,edge:-edge]

sw_lat_full = xlat[0,0]
sw_lon_full = xlon[0,0]
ne_lat_full = xlat[-1,-1]
ne_lon_full = xlon[-1,-1]

######################### Read PMM composite reflectivity: ####################################

pmm_dz = fin.variables["comp_dz_pmm"][edge:-edge,edge:-edge]

######################### Initialize variables: ####################################

tf_2 = fin.variables["t_2"][:,edge:-edge,edge:-edge]
tdf_2 = fin.variables["td_2"][:,edge:-edge,edge:-edge]
u_10 = fin.variables["u_10"][:,edge:-edge,edge:-edge]
v_10 = fin.variables["v_10"][:,edge:-edge,edge:-edge]

mean_tf_2 = np.mean(tf_2, axis=0)
mean_tdf_2 = np.mean(tdf_2, axis=0)
mean_u_10 = np.mean(u_10, axis=0)
mean_v_10 = np.mean(v_10, axis=0)

#mean_tf_2 = fin.variables["t_2"][edge:-edge,edge:-edge]
#mean_tdf_2 = fin.variables["td_2"][edge:-edge,edge:-edge]
#mean_u_10 = fin.variables["u_10"][edge:-edge,edge:-edge]
#mean_v_10 = fin.variables["v_10"][edge:-edge,edge:-edge]

fin.close()
del fin

#################################### Read HF METAR Data:  #############################################

##################### Get list of summary files to process: ##############################

files = []
files_temp = os.listdir(asos_dir)
for f, file in enumerate(files_temp):
   if (file[0] == '2'):
      files.append(file)

files.sort()

############### for each ensemble member summary file: #############################

for f, file in enumerate(files):
   file_month = int(file[-12:-10])
   file_day = int(file[-10:-8])
   file_hour = int(file[-7:-5])
   if (file[-2:] != 'gz'): 
      if ((int(file_month) == int(month)) and (int(file_day) == int(valid_day)) and (int(file_hour) == int(valid_hour))):
         asos_file = os.path.join(asos_dir, file)
         try:
            fin = netCDF4.Dataset(asos_file, "r")
            print "Opening %s \n" % asos_file
         except:
            print "%s does not exist! \n" % asos_file
            sys.exit(1)
         break

ob_type = fin.variables["metarType"][:]
ob_lat = fin.variables["latitude"][:]
ob_lon = fin.variables["longitude"][:]
ob_time = fin.variables["observationTime"][:]
ob_temp = fin.variables["temperature"][:]
ob_td = fin.variables["dewpoint"][:]
ob_dir = fin.variables["windDir"][:]
ob_speed = fin.variables["windSpeed"][:]
#print 'OB STUFF: ', len(ob_type), len(ob_temp), np.max(ob_temp)
#print 'TIME STUFF: ', epoch_seconds, ob_time
ob_time = np.abs(epoch_seconds - ob_time) #get total displacement in seconds between observation time and valid forecast
print 'MORE TIME STUFF: ', np.min(ob_time), epoch_seconds, ob_time
fin.close()
del fin

#find observations within +- 2.5 minutes of valid forecast time: 
ob_type = ob_type[ob_time < 150.]
ob_lat = ob_lat[ob_time < 150.]
ob_lon = ob_lon[ob_time < 150.]
ob_temp = ob_temp[ob_time < 150.]
ob_td = ob_td[ob_time < 150.]
ob_dir = ob_dir[ob_time < 150.]
ob_speed = ob_speed[ob_time < 150.]

print 'number of obs found in time period: ', len(ob_lat)

#use only Federal ASOS observations: 
ob_lat = ob_lat[ob_type == 1]
ob_lon = ob_lon[ob_type == 1]
ob_temp = ob_temp[ob_type == 1]
ob_td = ob_td[ob_type == 1]
ob_dir = ob_dir[ob_type == 1]
ob_speed = ob_speed[ob_type == 1]

#find observations within NEWS-e domain: 

ob_temp = ob_temp[(ob_lat > np.min(xlat)) & (ob_lat < np.max(xlat)) & (ob_lon > np.min(xlon)) & (ob_lon < np.max(xlon))]
ob_td = ob_td[(ob_lat > np.min(xlat)) & (ob_lat < np.max(xlat)) & (ob_lon > np.min(xlon)) & (ob_lon < np.max(xlon))]
ob_dir = ob_dir[(ob_lat > np.min(xlat)) & (ob_lat < np.max(xlat)) & (ob_lon > np.min(xlon)) & (ob_lon < np.max(xlon))]
ob_speed = ob_speed[(ob_lat > np.min(xlat)) & (ob_lat < np.max(xlat)) & (ob_lon > np.min(xlon)) & (ob_lon < np.max(xlon))]
ob_lat_plot = ob_lat[(ob_lat > np.min(xlat)) & (ob_lat < np.max(xlat)) & (ob_lon > np.min(xlon)) & (ob_lon < np.max(xlon))]
ob_lon_plot = ob_lon[(ob_lat > np.min(xlat)) & (ob_lat < np.max(xlat)) & (ob_lon > np.min(xlon)) & (ob_lon < np.max(xlon))]

print 'number of obs found total: ', len(ob_temp)

#convert wind dir (meteorology coordinates) to u and v components of wind: 
ob_u = ob_speed * np.cos(np.pi / 180. * (270 - ob_dir))
ob_v = ob_speed * np.sin(np.pi / 180. * (270 - ob_dir))

#conver u and v to kts: 
ob_u = ob_u *  1.943844
ob_v = ob_v *  1.943844

#mask lat/lon where obs are missing (should remove them from plot): 
for i in range(0, len(ob_lat_plot)):
   if ((np.ma.is_masked(ob_temp[i])) or (np.ma.is_masked(ob_td[i])) or (np.ma.is_masked(ob_u[i])) or (np.ma.is_masked(ob_v[i]))):
      print 'something is masked ... ', i
      ob_lat_plot[i] = -59. #hack so the plots don't show

#find corresponding NEWS-e observations for each ob: 
temp_err = np.zeros((len(ob_temp)))
td_err = np.zeros((len(ob_td)))
u_err = np.zeros((len(ob_u)))
v_err = np.zeros((len(ob_v)))

combo_lat_lon = np.dstack([xlat.ravel(),xlon.ravel()])[0]  #stacked NEWS-e lat/lon values for KDTree

for i in range(0, len(ob_temp)):
   ob_lat_lon = np.dstack([ob_lat_plot[i], ob_lon_plot[i]])
   mytree = scipy.spatial.cKDTree(combo_lat_lon)
   newse_dist, newse_indices = mytree.query(ob_lat_lon)

   temp_err[i] = ((ob_temp[i] - 273.15) * 1.8 + 32)
   td_err[i] = ((ob_td[i] - 273.15) * 1.8 + 32)
   u_err[i] = ob_u[i]
   v_err[i] = ob_v[i]

#################### Thin arrays used for quiver plots (removes extra vectors): #################################

quiv_xlon = xlon[0:-1:thin,0:-1:thin]
quiv_xlat = xlat[0:-1:thin,0:-1:thin]
quiv_u_10 = mean_u_10[0:-1:thin,0:-1:thin]
quiv_v_10 = mean_v_10[0:-1:thin,0:-1:thin]

################################# Make Figure Template: ###################################################

print 'basemap part'

#Load pickled basemap instance for faster plotting: 

#fig, ax1, ax2, ax3 = create_fig_nomap()

#map_temp = pl.load(open(mapname, 'rb'))

#P.sca(ax1)
#map = mymap_boundaries(map_temp, damage_files)

#map.drawcounties(linewidth=0.5, color=cb_colors.gray3)
#map.drawstates(linewidth=1., color=cb_colors.gray5)
#map.drawcoastlines(linewidth=1., color=cb_colors.gray5)
#map.drawcountries(linewidth=1., color=cb_colors.gray5)

map, fig, ax1, ax2, ax3 = create_fig(sw_lat_full, sw_lon_full, ne_lat_full, ne_lon_full, true_lat_1, true_lat_2, cen_lat, stand_lon, damage_files, resolution, area_thresh)

x, y = map(xlon[:], xlat[:])

ob_x, ob_y = map(ob_lon_plot[:], ob_lat_plot[:])

##########################################################################################################
###################################### Make Plots: #######################################################
##########################################################################################################

print 'plot part'

################# Rotate quiver values according to map projection: #####################

q_u = quiv_u_10
q_v = quiv_v_10

######################## Environment Plots: #####################

t_label = 'Temperature Error (NEWS-e - Ob; $^{\circ}$F)'
dot_plot(map, fig, ax1, ax2, ax3, x, y, ob_x, ob_y, temp_err, t_diff_levels, cb_colors.diff_cmap, t_label, temp_plot, mean_tf_2[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, q_u, q_v, 500, 5, 0, spec='False', quiv='True')

td_label = 'Dewpoint Error (NEWS-e - Ob; $^{\circ}$F)'
dot_plot(map, fig, ax1, ax2, ax3, x, y, ob_x, ob_y, td_err, t_diff_levels, cb_colors.diff_cmap, td_label, td_plot, mean_tdf_2[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, q_u, q_v, 500, 5, 0, spec='False', quiv='True')

u_label = 'U Component of Wind Error (NEWS-e - Ob; kts)'
dot_plot(map, fig, ax1, ax2, ax3, x, y, ob_x, ob_y, u_err, ws_diff_levels, cb_colors.diff_cmap, u_label, u_plot, mean_u_10[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, q_u, q_v, 500, 5, 0, spec='False', quiv='True')

v_label = 'V Component of Wind Error (NEWS-e - Ob; kts)'
dot_plot(map, fig, ax1, ax2, ax3, x, y, ob_x, ob_y, v_err, ws_diff_levels, cb_colors.diff_cmap, v_label, v_plot, mean_v_10[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, q_u, q_v, 500, 5, 0, spec='False', quiv='True')


