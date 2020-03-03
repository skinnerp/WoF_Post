#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import scipy
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
import time as timepy
import datetime
from optparse import OptionParser
import pickle as pl
import news_e_post_cbook
from news_e_post_cbook import *
import news_e_plotting_cbook_v2
from news_e_plotting_cbook_v2 import *

####################################### File Variables: ######################################################

parser = OptionParser()
parser.add_option("-r", dest="qpeqc_dir", type="string", default= None, help="Input Directory (of qpe_qc files)")
parser.add_option("-i", dest="image_dir", type="string", help = "Output Directory (for image files)")
parser.add_option("-a", dest="date", type="string", help = "date (YYYYMMDD)")
parser.add_option("-v", dest="var", type="string", help = "Variable to process")
parser.add_option("-t", dest="t", type="int", help = "Timestep to process")
parser.add_option("-n", dest="nt", type="int", help = "Number of timesteps")

(options, args) = parser.parse_args()

if ((options.qpeqc_dir == None) or (options.image_dir == None) or (options.date == None) or (options.var == None) or (options.t == None) or (options.nt == None)):
   print
   parser.print_help()
   print
   sys.exit(1)
else:
   qpeqc_dir = options.qpeqc_dir
   image_dir = options.image_dir
   date = options.date
   var = options.var
   t = options.t
   nt = options.nt

#   t_label = t #- 6     ###correct for MRMS matching window
#   nt_label = nt #- 12 

#time.sleep(10)

#################################### User-Defined Variables:  #####################################################

ne = 18

dx = 3.                         #horizontal grid spacing in km
merge_thresh = 10.              #maximum value of the shortest distance in km between two objects to be merged
cutoff_min_dis = 40.       #maximum value of the shortest distance between two objects to be matched in km
cutoff_rad = 40.                #maximum distance between two objects to be matched in km
cutoff_time = 1500.             #maximum time between two objects to be matched in s

score_thresh    = 0.2

#################################### Basemap Variables:  #####################################################

domain          = 'full'
resolution      = 'h'
area_thresh     = 1000.

damage_files = '' #['/Volumes/fast_scr/pythonletkf/vortex_se/2013-11-17/shapefiles/extractDamage_11-17/extractDamagePaths']

#################################### Contour Levels:  #####################################################

#obj_prob_levels          = np.arange(0.05,0.65,0.1)               #(%)
obj_prob_levels          = np.arange(0.0,1.1,0.1)               #(%)

#################################### Colormap Names:  #####################################################

match_cmap        = matplotlib.colors.ListedColormap([cb_colors.blue3, cb_colors.blue4, cb_colors.blue5, cb_colors.blue6, cb_colors.blue7, cb_colors.purple7, cb_colors.purple6, cb_colors.purple5, cb_colors.purple4, cb_colors.purple3])
fa_cmap           = matplotlib.colors.ListedColormap([cb_colors.orange3, cb_colors.orange4, cb_colors.orange5, cb_colors.orange6, cb_colors.orange7, cb_colors.red7, cb_colors.red6, cb_colors.red5, cb_colors.red4, cb_colors.red3])

plot_alpha        = 0.5
neighborhood      = 15 

#################################### Initialize plot attributes:  #####################################################

object_plot = ob_plot('',       \
                   '',                  \
                   '',                        \
                   cb_colors.blue6,      \
                   cb_colors.orange6,     \
                   obj_prob_levels,              \
                   obj_prob_levels,     \
                   match_cmap,                      \
                   fa_cmap,                     \
                   [0.001, 100.],  \
                   [cb_colors.gray7, cb_colors.gray7], \
                   'neither',       \
                   plot_alpha,  \
                   neighborhood)

#################################### Read Grid info/QPE QC Data:  #####################################################

qpe_files = os.listdir(qpeqc_dir)
qpe_files.sort()

#################################### Read Grid info from first QPE QC file and initialize MRMS var:  #####################################################
qpe_path = os.path.join(qpeqc_dir, qpe_files[0])

try:
   fin = netCDF4.Dataset(qpe_path, "r")
   print("Opening %s \n" % qpe_path)
except:
   print("%s does not exist! \n" % qpe_path)
   sys.exit(1)

xlat = fin.variables['XLAT'][:]
xlon = fin.variables['XLON'][:]

rainmask = fin.variables['RAINMASK'][:]
rainmask = rainmask * 0.                     ### hack for sls plotting
sw_lat = xlat[0,0]
sw_lon = xlon[0,0]
ne_lat = xlat[-1,-1]
ne_lon = xlon[-1,-1]

cen_lat = fin.CEN_LAT
cen_lon = fin.CEN_LON
stand_lon = fin.STAND_LON
true_lat1 = fin.TRUE_LAT1
true_lat2 = fin.TRUE_LAT2

#print cen_lat, cen_lon, stand_lon, true_lat1, true_lat2

map = Basemap(llcrnrlon=sw_lon, llcrnrlat=sw_lat, urcrnrlon=ne_lon, urcrnrlat=ne_lat, projection='lcc', lat_1=true_lat1, lat_2=true_lat2, lat_0=cen_lat, lon_0=cen_lon, resolution = 'c', area_thresh = area_thresh)

x_offset, y_offset = map(cen_lon, cen_lat)
x, y = map(xlon[:], xlat[:])

x = x - x_offset
y = y - y_offset

fin.close()
del fin

mrms_indices = np.arange(t, t+1)
mrms = np.zeros((len(mrms_indices),xlat.shape[0],xlat.shape[1]))

for q, qfile in enumerate(qpe_files):
   current_qpe_index = int(qfile[-24:-22])
   if (current_qpe_index == t): 
      qpe_path = os.path.join(qpeqc_dir, qfile)

      try:
         fin = netCDF4.Dataset(qpe_path, "r")
         print("Opening %s \n" % qpe_path)
      except:
         print("%s does not exist! \n" % qpe_path)
         sys.exit(1)

      year = qpe_path[-21:-17]
      month = qpe_path[-17:-15]
      day = qpe_path[-15:-13]
      init_hour = qpe_path[-12:-10]
      init_min = qpe_path[-10:-8]
      valid_hour = qpe_path[-7:-5]
      valid_min = qpe_path[-5:-3]

      if ((int(valid_hour) < 12) and (int(init_hour) > 18)):
         temp_day = int(day) + 1
         valid_day = str(temp_day)
         if (len(valid_day) == 1):
            valid_day = '0' + valid_day
      else:
         valid_day = day

      init_label = 'Init: ' + year + '-' + month + '-' + day + ', ' + init_hour + init_min + ' UTC'
      valid_label = 'Valid: ' + year + '-' + month + '-' + valid_day + ', ' + valid_hour + valid_min + ' UTC'

      init_time = int(init_hour) * 3600. + int(init_min) * 60. 
#      time = init_time + t * 300. 
#      fcst_time = time - init_time

#print time, fcst_time

#      valid_hour = np.floor(time / 3600.)
#      valid_min = np.floor((time - valid_hour * 3600.) / 60.)
#
#      if (valid_hour > 23):
#         valid_hour = valid_hour - 24
#
#         if (init_hour > 18):
#            temp_day = int(date[-2:])+1
#            temp_day = str(temp_day)
#
#            if (len(temp_day) == 1):
#               temp_day = '0' + temp_day
#
#            date = date[:-2] + temp_day
#
#      valid_hour = str(int(valid_hour))
#      valid_min = str(int(valid_min))
#
#      if (len(valid_hour) == 1):
#         valid_hour = '0' + valid_hour
#      if (len(valid_min) == 1):
#         valid_min = '0' + valid_min

      if (var == 'rain1'):
         newse = fin.variables['RAIN_1_QC'][:,:,:]   
         object_plot.name = 'qpe1_obmatch'
         object_plot.var1_title = 'Gridpoint Prob. of Matched Rainfall Objects'
         object_plot.var2_title = 'Gridpoint Prob. of False Alarm Rainfall Objects'
      elif (var == 'rain2'):
         newse = fin.variables['RAIN_2_QC'][:,:,:]
         object_plot.name = 'qpe2_obmatch'
         object_plot.var1_title = 'Gridpoint Prob. of Matched Rainfall Objects'
         object_plot.var2_title = 'Gridpoint Prob. of False Alarm Rainfall Objects'
      elif (var == 'rain3'):
         newse = fin.variables['RAIN_3_QC'][:,:,:]
         object_plot.name = 'qpe3_obmatch'
         object_plot.var1_title = 'Gridpoint Prob. of Matched Rainfall Objects'
         object_plot.var2_title = 'Gridpoint Prob. of False Alarm Rainfall Objects'

      if ((t > 0) and (t <= nt)): 
         blank = 0
      elif (t < 1): 
         blank = 1
#
      fin.close()
      del fin
#
   for temp_index, tt in enumerate(mrms_indices): 
      if ((tt >= 0)): 
         tt_str = str(tt) 

      if (len(tt_str) == 1): 
         tt_str = '0' + tt_str

#      print current_qpe_index, tt_str, rfile
      temp_str_current_qpe_index = str(current_qpe_index)
      if (len(temp_str_current_qpe_index) == 1): 
         temp_str_current_qpe_index = '0' + temp_str_current_qpe_index

      if (temp_str_current_qpe_index == tt_str): 
         qpe_path = os.path.join(qpeqc_dir, qfile)

         try:
            fin = netCDF4.Dataset(qpe_path, "r")
            print("Opening %s \n" % qpe_path)
         except:
            print("%s does not exist! \n" % qpe_path)
            sys.exit(1)

         if (var == 'rain1'):
            mrms[temp_index,:,:] = fin.variables['MRMS_QPE_1_QC'][:,:]
         if (var == 'rain2'):
            mrms[temp_index,:,:] = fin.variables['MRMS_QPE_2_QC'][:,:]
         if (var == 'rain3'):
            mrms[temp_index,:,:] = fin.variables['MRMS_QPE_3_QC'][:,:]
         print('in MRMS update loop: ',(mrms))
         fin.close()
         del fin

#print 'BLANK BLANK BLANK BLANK BLANK BLANK: ', blank
#print mrms.shape, 'MRMS SHAPE'
print("able to open the MRMS and SMR files")
if ((blank == 0) or (var == 'rain1') or (var == 'rain2') or (var == 'rain3')): #if rain object or not in blanking region 
#if ((blank == 0)):

#################################### Initialize Plotting Arrays:  #####################################################

   match_plot = newse[0,:,:] * 0. 
   fa_plot = newse[0,:,:] * 0. 
 
   mrms_plot = mrms[0,:,:]  #mrms obs at t

   mrms_plot_binary = np.where(mrms_plot > 0., 1, 0)
   mrms_plot_labels = skimage.measure.label(mrms_plot_binary)
   mrms_plot_props = regionprops(mrms_plot_labels, mrms_plot)

   mrms_count = len(mrms_plot_props)

   print(mrms_plot.shape)
   print(mrms_plot.max(), mrms_plot.min())
#   timepy.sleep(5)

#   sys.exit(0)
#################################### Find/merge MRMS objects:  #####################################################

   ob_time = []
   ob_centroid_x = []
   ob_centroid_y = []
   ob_coords = []

   for tt in range(0, (len(mrms[:,0,0]))):  #skip 15 min blanking periods at beginning and end
      temp_time = init_time + (t + (tt-3)) * 300.   #time of current obs being processed
      aws_temp = mrms[tt,:,:]
      aws_binary = np.where(aws_temp > 0., 1, 0)
      aws_binary = aws_binary.astype(int)            #not sure if needed, but worried 0's aren't being treated as ints in measure.label
      aws_labels = skimage.measure.label(aws_binary)
      aws_labels = aws_labels.astype(int)
      aws_props = regionprops(aws_labels, aws_temp)

##### Merge objects if < merge_thresh (probably 9 km) apart ####

      merge_labels = aws_labels
      for i in range(0, (len(aws_props)-1)):
       i_tree = scipy.spatial.cKDTree(aws_props[i].coords)
       for j in range((i+1), len(aws_props)):
          j_dist, dummy = i_tree.query(aws_props[j].coords)
          j_dist = j_dist * dx    #convert to km
          if (np.min(j_dist) < merge_thresh):
             merge_labels = np.where(aws_labels == aws_props[i].label, aws_props[j].label, merge_labels)
      merge_props = regionprops(merge_labels, aws_temp)

      for i in range(0, len(merge_props)):
         if ((merge_props[i].max_intensity <= 0.003) or (merge_props[i].area < 2)):
            print('what happened?', tt, i, merge_props[i].max_intensity, merge_props[i].area)
         else:
            ob_time.append(temp_time)
            ob_centroid_x.append(merge_props[i].centroid[1])
            ob_centroid_y.append(merge_props[i].centroid[0])
            ob_coords.append(merge_props[i].coords)

   ob_time = np.asarray(ob_time)
   ob_centroid_x = np.asarray(ob_centroid_x)
   ob_centroid_y = np.asarray(ob_centroid_y)
   ob_coords = np.asarray(ob_coords)

   ob_indices = np.arange(0, len(ob_time))

   print(ob_time, ob_centroid_x, ob_centroid_y, ob_coords)
###################################### Find Rainfall Objects: #######################################################

   pod = np.zeros((newse.shape[0]))
   far = np.zeros((newse.shape[0]))
   bias = np.zeros((newse.shape[0]))
   csi = np.zeros((newse.shape[0]))

   var_centroid_x = []
   var_centroid_y = []
   var_coords = []

   for n in range(0, newse.shape[0]):
      newse_match = []
      newse_fa = []
      temp_time = init_time + t * 300.
      var_temp = newse[n,:,:]
      var_binary = np.where(var_temp > 0., 1, 0)
      var_binary = var_binary.astype(int)
      var_labels = skimage.measure.label(var_binary)
      var_labels = var_labels.astype(int)
      var_props = regionprops(var_labels, var_temp)
      merge_labels = var_labels
#      for i in range(0, (len(var_props)-1)):
#       i_tree = scipy.spatial.cKDTree(var_props[i].coords)
#       for j in range((i+1), len(var_props)):
#           j_dist, dummy = i_tree.query(var_props[j].coords)
#           j_dist = j_dist * dx    #convert to km
#           if (np.min(j_dist) < merge_thresh):
#              merge_labels = np.where(var_labels == var_props[i].label, var_props[j].label, merge_labels)
#              print 'merged! ', i, j, var_props[i].coords, var_props[j].coords
      merge_props = regionprops(merge_labels, var_temp)

      for i in range(0, len(merge_props)):
         if ((merge_props[i].max_intensity < 1.0)):# and (var_props[i].area > var_area_thresh)):
            print('what happened (newse)? ', i, merge_props[i].area, merge_props[i].max_intensity)
            for aa in range(0, len(merge_props)):
               print(merge_props[aa].max_intensity)
               print(merge_props[aa].area)
               print(merge_props[aa].coords)
         else:
            var_centroid_x.append(merge_props[i].centroid[1])
            var_centroid_y.append(merge_props[i].centroid[0])
            var_coords.append(merge_props[i].coords)

            match_ti = 0.

            for j in range(0, len(ob_centroid_x)):
               temp_ob_x = gridpoint_interp(x, ob_centroid_x[j], ob_centroid_y[j])
               temp_ob_y = gridpoint_interp(y, ob_centroid_x[j], ob_centroid_y[j])

               temp_var_x = gridpoint_interp(x, var_centroid_x[-1], var_centroid_y[-1])
               temp_var_y = gridpoint_interp(y, var_centroid_x[-1], var_centroid_y[-1])

               temp_dis = np.sqrt(((temp_var_x - temp_ob_x) / 1000.)**2 + ((temp_var_y - temp_ob_y) / 1000.)**2)
               temp_ob_tree = scipy.spatial.cKDTree(ob_coords[j])
               temp_min_dis, dummy = temp_ob_tree.query(var_coords[-1])
               temp_min_dis = np.min(temp_min_dis) * dx

               temp_dis_score = (cutoff_rad - temp_dis) / cutoff_rad
               temp_min_dis_score = (cutoff_min_dis - temp_min_dis) / cutoff_min_dis
               temp_time_score = (cutoff_time - np.abs(temp_time - ob_time[j])) / cutoff_time

               temp_ti = ((temp_dis_score + temp_min_dis_score) / 2.) * temp_time_score

               if (temp_ti > match_ti):
                  match_ti = temp_ti

            if (match_ti > score_thresh):
               newse_match.append(i+1)
               match_plot = np.where(merge_labels == (i+1), match_plot + 1, match_plot)
            else:
               newse_fa.append(i+1)
               fa_plot = np.where(merge_labels == (i+1), fa_plot + 1, fa_plot)


      a = len(newse_match) * 1.
      b = len(newse_fa) * 1.
      c = (mrms_count - a) * 1.
      if (c < 0.):  #temporal matching can create situations where there are more matched forecast objects than observed objects
         c = 0. 

      if ((a+c) == 0.): 
         pod[n] = -999. 
      else:
         pod[n] = a / (a + c)
      if ((a+b) == 0.):
         far[n] = -999. 
      else:  
         far[n] = b / (a + b)
      if ((a+c) == 0.): 
         bias[n] = -999. 
      else: 
         bias[n] = (a + b) / (a + c)
      if ((a+b+c) == 0.): 
         csi[n] = -999. 
      else: 
         csi[n] = a / (a + b + c)

   match_prob = match_plot / ne
   fa_prob = fa_plot / ne

   masked_match_prob = np.ma.masked_array(match_prob, match_prob<0.01)
   masked_fa_prob = np.ma.masked_array(fa_prob, fa_prob<0.01)

   masked_pod = np.ma.masked_array(pod, pod<0.)
   masked_far = np.ma.masked_array(far, far<0.)
   masked_bias = np.ma.masked_array(bias, bias<0.)
   masked_csi = np.ma.masked_array(csi, csi<0.)

   print(pod)
################################# Make Figure Template: ###################################################
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
#   hail, wind, tornado = plot_lsr(map, fig, ax1, ax2, ax3, lsr_file, init_seconds, valid_seconds, plot_h='True', plot_w='True')
#   
#   ########### Warnings #########
#   
#   svr, tor, ff = plot_warn(map, fig, ax1, ax2, ax3, warn_file, (valid_seconds-300.), valid_seconds, cb_colors.blue6, cb_colors.red6, cb_colors.green6)

   if ((t >= 0) and (t <= nt)): 
      print('basemap part')

      map, fig, ax1, ax2, ax3 = create_fig(sw_lat, sw_lon, ne_lat, ne_lon, true_lat1, true_lat2, cen_lat, cen_lon, damage_files, resolution, area_thresh, object='True')

      hail, wind, tornado = plot_lsr(map, fig, ax1, ax2, ax3, lsr_file, init_seconds, valid_seconds, plot_h='True', plot_w='True')
   
      ########### Warnings #########
   
      svr, tor, ff = plot_warn(map, fig, ax1, ax2, ax3, warn_file, (valid_seconds-300.), valid_seconds, cb_colors.blue6, cb_colors.red6, cb_colors.green6)

      x, y = map(xlon[:], xlat[:])
      xx, yy = map.makegrid(xlat.shape[1], xlat.shape[0], returnxy=True)[2:4]   #equidistant x/y grid for streamline plots

      print('MAX MRMS: ', np.max(mrms_plot), t, np.max(mrms), mrms_indices) 
      ob_ver_plot(map, fig, ax1, ax2, ax3, x, y, x, y, object_plot, masked_match_prob, masked_fa_prob, mrms_plot, rainmask, t, np.mean(masked_pod), np.mean(masked_far), np.mean(masked_bias), np.mean(masked_csi), init_label, valid_label, domain, image_dir, 5, 0)

else: ##Plot "No rainfall objects!" plot for first timestep:

   if ((t == 0)): 
      print('basemap part')

      map, fig, ax1, ax2, ax3 = create_fig(sw_lat, sw_lon, ne_lat, ne_lon, true_lat1, true_lat2, cen_lat, cen_lon, damage_files, resolution, area_thresh, object='False')

      hail, wind, tornado = plot_lsr(map, fig, ax1, ax2, ax3, lsr_file, init_seconds, valid_seconds, plot_h='True', plot_w='True')
   
      ########### Warnings #########
   
      svr, tor, ff = plot_warn(map, fig, ax1, ax2, ax3, warn_file, (valid_seconds-300.), valid_seconds, cb_colors.blue6, cb_colors.red6, cb_colors.green6)

      x, y = map(xlon[:], xlat[:])
      xx, yy = map.makegrid(xlat.shape[1], xlat.shape[0], returnxy=True)[2:4]   #equidistant x/y grid for streamline plots

      ob_ver_plot(map, fig, ax1, ax2, ax3, x, y, x, y, object_plot, mrms_plot, mrms_plot, mrms_plot, rainmask, t, 0., 0., 0., 0., init_label, valid_label, domain, image_dir, 5, 0, blank='True')
