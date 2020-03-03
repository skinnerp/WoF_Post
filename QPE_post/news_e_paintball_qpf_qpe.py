#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
from scipy import signal
from scipy import *
from scipy import ndimage
import skimage
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
import time as timeit
from optparse import OptionParser
import pickle as pl
import news_e_post_cbook
from news_e_post_cbook import *
import news_e_plotting_cbook_v2
from news_e_plotting_cbook_v2 import *

####################################### File Variables: ######################################################

parser = OptionParser()
parser.add_option("-m", dest="mrms_dir", type="string", default=None, help="Input directory of MRMS QPE plot")
parser.add_option("-d", dest="summary_dir", type="string", default= None, help="Input directory of summary files to plot")
parser.add_option("-o", dest="outdir", type="string", help = "Output Directory (for images)")
parser.add_option("-q", dest="qpe_outdir", type="string", default=None, help="Output directory for qpe_qc files")
#parser.add_option("-m", dest="mapname", type="string", help = "Path to Pickled Basemap instance")
parser.add_option("-t", dest="t", type="int", help = "Timestep to process")
parser.add_option("-n", dest="nt", type="int", help = "Number of timesteps in forecast")
parser.add_option("-p", dest="date", type="string", help="Date of case")

(options, args) = parser.parse_args()

if ((options.mrms_dir == None) or (options.summary_dir == None) or (options.outdir == None) or (options.qpe_outdir == None) or (options.t == None) or (options.nt == None) or (options.date == None)):
   print
   parser.print_help()
   print
   sys.exit(1)
else:
   mrms_dir = options.mrms_dir
   summary_dir = options.summary_dir
   outdir = options.outdir
   qpe_outdir = options.qpe_outdir
   t = options.t
   nt = options.nt
   date = options.date

#path to pickled map instance for plotting: 
#mapname = '/scratch2/patrick.skinner/images/map.pickle'

#################################### User-Defined Variables:  #####################################################

qpf_thresh_1               = 1.
qpf_thresh_2               = 2.
qpf_thresh_3               = 3. 

newse_dz_thresh            = 45.002             

########## 99.95th percentile thresholds for NEWS-e 2017 cases: 

uh_0to2_thresh             = 14.217
uh_2to5_thresh             = 65.790

area_thresh_uh             = 10.            #Minimum area of rotation object
area_thresh_dz             = 12.           #Minimum area of reflectivity object
cont_thresh                = 2             #Must be greater than this number of forecast timesteps in a rotation object to be retained

swath_window               = 3             #+- 3 timesteps for swaths
domain                     = 'full'        #vestigial variable that's still needed in the plotting subroutines ... 
edge                       = 7          #number of grid points to remove from near domain boundaries

neighborhood               = 15                 #15 gridpoint radius for probability matched mean 

plot_alpha                 = 0.55                #transparency value for filled contour plots

#################################### Basemap Variables:  #####################################################

resolution      = 'h'
area_thresh     = 1000.

damage_files = '' #['/scratch/skinnerp/2018_newse_post/damage_files/extractDamage_new/extractDamagePolys'] #['/Volumes/fast_scr/pythonletkf/vortex_se/2013-11-17/shapefiles/extractDamage_11-17/extractDamagePaths']

#################################### Initialize plot attributes using 'web plot' objects:  #####################################################

qpf_paint_plot_1 = v_plot('qpf_qpe_paint_1',                   \
                   'WoFS Member 1 In. Rainfall',                        \
                   '',                  \
                   cb_colors.gray6,   \
                   cb_colors.gray6,     \
                   [0.99, 1000.],              \
                   [0.99, 1000.],     \
                   [cb_colors.gray8],                      \
                   cb_colors.paintball_colors_list,                     \
                   '',  \
                   '',  \
                   'max',  \
                   0.6, \
                   neighborhood)

qpf_paint_plot_2 = v_plot('qpf_qpe_paint_2',                   \
                   'WoFS Member 2 In. Rainfall',                        \
                   '',                  \
                   cb_colors.gray6,   \
                   cb_colors.gray6,     \
                   [1.99, 1000.],              \
                   [1.99, 1000.],     \
                   [cb_colors.gray8],                      \
                   cb_colors.paintball_colors_list,                     \
                   '',  \
                   '',  \
                   'max',  \
                   0.6, \
                   neighborhood)

qpf_paint_plot_3 = v_plot('qpf_qpe_paint_3',                   \
                   'WoFS Member 3 In. Rainfall',                        \
                   '',                  \
                   cb_colors.gray6,   \
                   cb_colors.gray6,     \
                   [2.99, 1000.],              \
                   [2.99, 1000.],     \
                   [cb_colors.gray8],                      \
                   cb_colors.paintball_colors_list,                     \
                   '',  \
                   '',  \
                   'max',  \
                   0.6, \
                   neighborhood)


############################ Find WRFOUT files to process: #################################
timestep = str(t)
if (len(timestep) == 1):
   timestep = '0' + timestep
### Find ENS Summary files ### 

ne = 18

ens_files = []
summary_files_temp = os.listdir(summary_dir)

for f, file in enumerate(summary_files_temp):
   if ((file[-28:-25] == 'PCP')):                               #assumes filename format of: news-e_ENS_20170516_2200_0000.nc
      ens_files.append(file)

ens_files.sort()

ens_file = ens_files[t]

infile = os.path.join(summary_dir, ens_file)

if (t == 0):
   outfile = ens_file[0:7] + 'QPE' + ens_file[10:]
else:
   outfile = ens_file[0:7] + 'QPE_' + timestep + ens_file[13:]

out_path = os.path.join(qpe_outdir, outfile)

try:                                                 #open WRFOUT file
   fin = netCDF4.Dataset(infile, "r")
   print("Opening %s \n" % infile)
except:
   print("%s does not exist! \n" %infile)
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

init_label = 'Init: ' + year + '-' + month + '-' + day + ', ' + init_hour + init_min + ' UTC'
valid_label = 'Valid: ' + year + '-' + month + '-' + valid_day + ', ' + valid_hour + valid_min + ' UTC'
######################### Set domain: ####################################

xlat = fin.variables["xlat"][edge:-edge,edge:-edge]                     #latitude (dec deg; Lambert conformal)
xlon = fin.variables["xlon"][edge:-edge,edge:-edge]                     #longitude (dec deg; Lambert conformal)

sw_lat_full = xlat[0,0]
sw_lon_full = xlon[0,0]
ne_lat_full = xlat[-1,-1]
ne_lon_full = xlon[-1,-1]

########## Read rain_mem variable: #############

rain_mem                            = fin.variables["rain_mem"][:,edge:-edge,edge:-edge]

fin.close()
del fin

###

radmask = rain_mem[0,:,:] * 0. #dummy variable (will be used when MRMS obs are included) 

rain_obj_1 = find_objects_timestep(rain_mem, radmask, qpf_thresh_1, area_thresh_dz)
rain_obj_2 = find_objects_timestep(rain_mem, radmask, qpf_thresh_2, area_thresh_dz)
rain_obj_3 = find_objects_timestep(rain_mem, radmask, qpf_thresh_3, area_thresh_dz)

########## Read MRMS variables: #############

#### Hack to get valid_time_seconds to jump ####

valid_time_seconds = int(valid_hour) * 3600. + int(valid_min) * 60.

#if (valid_time_seconds > 44099
if (valid_time_seconds <= 44100.):
    valid_time_seconds = valid_time_seconds + 86400.
else:
    valid_time_seconds = valid_time_seconds

mrms_files = os.listdir(mrms_dir)
mrms_files.sort()

for f, file in enumerate(mrms_files):
   temp_time = int(file[22:24]) * 3600. + int(file[24:26]) * 60.
   temp_date = file[0:8]
   if (temp_time <= 44100.):
      temp_time += 86400.

   if (temp_time == valid_time_seconds):
      temp_file = os.path.join(mrms_dir, file)
 
      try:                                                 #open WRFOUT file
         fin = netCDF4.Dataset(temp_file, "r")
         print("Opening %s \n" % temp_file)
      except:
         print("%s does not exist! \n" % temp_file)
         sys.exit(1)

      qpe_swt                             = fin.variables["MRMS_QPE_SWT"][edge:-edge,edge:-edge]
#               mrms_dz                             = fin.variables["DZ_CRESSMAN"][edge:-edge,edge:-edge]
      fin.close()
      del fin

rainmask = qpe_swt[:,:] * 0. #dummy variable

mrms_qpe_obj_1 = find_objects_timestep(qpe_swt, rainmask, qpf_thresh_1, area_thresh_dz)
mrms_qpe_obj_2 = find_objects_timestep(qpe_swt, rainmask, qpf_thresh_2, area_thresh_dz)
mrms_qpe_obj_3 = find_objects_timestep(qpe_swt, rainmask, qpf_thresh_3, area_thresh_dz)

################################# Write to qpe_qc file: ###################################################
try:
   fout = netCDF4.Dataset(out_path, "w")
except:
   print("Could not create %s!\n" % out_path)

fout.createDimension('NX', xlat.shape[1])
fout.createDimension('NY', xlat.shape[0])
fout.createDimension('NE', ne)

setattr(fout,'CEN_LAT',cen_lat)
setattr(fout,'CEN_LON',cen_lon)
setattr(fout,'STAND_LON',stand_lon)
setattr(fout,'TRUE_LAT1',true_lat_1)
setattr(fout,'TRUE_LAT2',true_lat_2)

fout.createVariable('XLAT', 'f4', ('NY','NX',))
fout.createVariable('XLON', 'f4', ('NY','NX',))
fout.createVariable('RAINMASK', 'f4', ('NY','NX',))

fout.createVariable('MRMS_QPE_1_QC', 'f4', ('NY','NX',))
fout.createVariable('MRMS_QPE_1_RAW', 'f4', ('NY','NX',))
fout.createVariable('MRMS_QPE_2_QC', 'f4', ('NY','NX',))
fout.createVariable('MRMS_QPE_2_RAW', 'f4', ('NY','NX',))
fout.createVariable('MRMS_QPE_3_QC', 'f4', ('NY','NX',))
fout.createVariable('MRMS_QPE_3_RAW', 'f4', ('NY','NX',))

fout.createVariable('RAIN_1_QC', 'f4', ('NE','NY','NX',))
fout.createVariable('RAIN_1_RAW', 'f4', ('NE','NY','NX',))
fout.createVariable('RAIN_2_QC', 'f4', ('NE','NY','NX',))
fout.createVariable('RAIN_2_RAW', 'f4', ('NE','NY','NX',))
fout.createVariable('RAIN_3_QC', 'f4', ('NE','NY','NX',))
fout.createVariable('RAIN_3_RAW', 'f4', ('NE','NY','NX',))

fout.variables['XLAT'][:] = xlat
fout.variables['XLON'][:] = xlon

fout.variables['MRMS_QPE_1_QC'][:,:] = mrms_qpe_obj_1
fout.variables['MRMS_QPE_1_RAW'][:,:] = qpe_swt
fout.variables['MRMS_QPE_2_QC'][:,:] = mrms_qpe_obj_2
fout.variables['MRMS_QPE_2_RAW'][:,:] = qpe_swt
fout.variables['MRMS_QPE_3_QC'][:,:] = mrms_qpe_obj_3
fout.variables['MRMS_QPE_3_RAW'][:,:] = qpe_swt

fout.variables['RAIN_1_QC'][:,:,:] = rain_obj_1
fout.variables['RAIN_1_RAW'][:,:,:] = rain_mem
fout.variables['RAIN_2_QC'][:,:,:] = rain_obj_2
fout.variables['RAIN_2_RAW'][:,:,:] = rain_mem
fout.variables['RAIN_3_QC'][:,:,:] = rain_obj_3
fout.variables['RAIN_3_RAW'][:,:,:] = rain_mem

fout.close()
del fout
################################# Make Figure Template: ###################################################

print('basemap part')

map, fig, ax1, ax2, ax3 = create_fig(sw_lat_full, sw_lon_full, ne_lat_full, ne_lon_full, true_lat_1, true_lat_2, cen_lat, stand_lon, damage_files, resolution, area_thresh)

x, y = map(xlon[:], xlat[:])
xx, yy = map.makegrid(xlat.shape[1], xlat.shape[0], returnxy=True)[2:4]   #equidistant x/y grid for streamline plots

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
hail, wind, tornado = plot_lsr(map, fig, ax1, ax2, ax3, lsr_file, init_seconds, valid_seconds, plot_h='True', plot_w='True')

########### Warnings #########

svr, tor, ff = plot_warn(map, fig, ax1, ax2, ax3, warn_file, (valid_seconds-300.), valid_seconds, cb_colors.blue6, cb_colors.red6, cb_colors.green6)

print('plot part')

aws_qc = xlat * 0.   ### set aws_qc to 0 since not plotting verification 

### Plot QPE 1" objects:
print('t is:  ', t)

if (t == 0):
   paintqc_plot(map, fig, ax1, ax2, ax3, x, y, x, y, qpf_paint_plot_1, mrms_qpe_obj_1, rain_obj_1, rainmask, t, init_label, valid_label, domain, outdir, 5, 0, blank='True')
else:   
   paintqc_plot(map, fig, ax1, ax2, ax3, x, y, x, y, qpf_paint_plot_1, mrms_qpe_obj_1, rain_obj_1, rainmask, t, init_label, valid_label, domain, outdir, 5, 0, blank='False')

###########################################################################################################
### make new figure for each paint plot ... can't figure out how to remove each members plot from basemap
###########################################################################################################

print('basemap part, part 2 - the basemappening')
#

########### Warnings #########

map2, fig2, ax21, ax22, ax23 = create_fig(sw_lat_full, sw_lon_full, ne_lat_full, ne_lon_full, true_lat_1, true_lat_2, cen_lat, stand_lon, damage_files, resolution, area_thresh)
hail, wind, tornado = plot_lsr(map, fig, ax21, ax22, ax23, lsr_file, init_seconds, valid_seconds, plot_h='True', plot_w='True')
svr, tor, ff = plot_warn(map, fig, ax1, ax2, ax3, warn_file, (valid_seconds-300.), valid_seconds, cb_colors.blue6, cb_colors.red6, cb_colors.green6)

#
### Plot QPE 2" objects: 
if (t == 0):
   paintqc_plot(map2, fig2, ax21, ax22, ax23, x, y, x, y, qpf_paint_plot_2, mrms_qpe_obj_2, rain_obj_2, rainmask, t, init_label, valid_label, domain, outdir, 5, 0, blank='True')
else:
   paintqc_plot(map2, fig2, ax21, ax22, ax23, x, y, x, y, qpf_paint_plot_2, mrms_qpe_obj_2, rain_obj_2, rainmask, t, init_label, valid_label, domain, outdir, 5, 0, blank='False')

print('basemap part, part 3 - So very tired ...')

########### Warnings #########

map3, fig3, ax31, ax32, ax33 = create_fig(sw_lat_full, sw_lon_full, ne_lat_full, ne_lon_full, true_lat_1, true_lat_2, cen_lat, stand_lon, damage_files, resolution, area_thresh)
svr, tor, ff = plot_warn(map, fig, ax31, ax32, ax33, warn_file, (valid_seconds-300.), valid_seconds, cb_colors.blue6, cb_colors.red6, cb_colors.green6)
hail, wind, tornado = plot_lsr(map, fig, ax31, ax32, ax33, lsr_file, init_seconds, valid_seconds, plot_h='True', plot_w='True')

#### Plot QPE 3" objects: 

if (t == 0):
   paintqc_plot(map3, fig3, ax31, ax32, ax33, x, y, x, y, qpf_paint_plot_3, mrms_qpe_obj_3, rain_obj_3, rainmask, t, init_label, valid_label, domain, outdir, 5, 0, blank='True')
else:
   paintqc_plot(map3, fig3, ax31, ax32, ax33, x, y, x, y, qpf_paint_plot_3, mrms_qpe_obj_3, rain_obj_3, rainmask, t, init_label, valid_label, domain, outdir, 5, 0, blank='False')
