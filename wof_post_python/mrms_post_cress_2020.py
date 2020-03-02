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
import time
from optparse import OptionParser
import news_e_post_cbook
from news_e_post_cbook import *
import news_e_plotting_cbook_v2
from news_e_plotting_cbook_v2 import *
import radar_info
from radar_info import *

####################################### Parse Options: ######################################################

parser = OptionParser()
parser.add_option("-a", dest="aws_low_file", type="string", default= None, help="Input Path (of MRMS 0-2 AWS file)")
parser.add_option("-m", dest="aws_mid_file", type="string", default= None, help="Input Path (of MRMS 2-5 AWS file)")
parser.add_option("-z", dest="dz_file", type="string", default= None, help="Input Path (of MRMS dBZ file)")
parser.add_option("-t", dest="tot_low_file", type="string", default= None, help="Input Path (of MRMS TOT Shear 0-2 file)")
parser.add_option("-u", dest="tot_mid_file", type="string", default= None, help="Input Path (of MRMS TOT Shear 2-5 file)")
parser.add_option("-v", dest="mesh_file", type="string", default= None, help="Input Path (of MRMS MESH file)")
parser.add_option("-o", dest="out_file", type="string", help = "Output File Path")
parser.add_option("-f", dest="newse_path", type="string", help = "Path to NEWS-e Summary File")

(options, args) = parser.parse_args()

if ((options.aws_low_file == None) or (options.aws_mid_file == None) or (options.dz_file == None) or (options.tot_low_file == None) or (options.tot_mid_file == None) or (options.mesh_file == None) or (options.out_file == None) or (options.newse_path == None)):
    print
    parser.print_help()
    print
    sys.exit(1)
else:
    aws_low_file = options.aws_low_file
    aws_mid_file = options.aws_mid_file
    dz_file = options.dz_file
    tot_low_file = options.tot_low_file
    tot_mid_file = options.tot_mid_file
    mesh_file = options.mesh_file
    out_file = options.out_file
    newse_path = options.newse_path

time.sleep(60) #wait a minute for MRMS file to finish writing

#################################### Specify OBAN/qc thresholds for DZ/Az shear:  #####################################################

dz_thresh_1 = 5. 
aws_thresh_1 = 0. 
min_obs = 4
min_dz_obs = 8

#################################### Variable names to OBAN/qc:  #####################################################

mrms_low_var = 'MergedAzShear_0-2kmAGL'
mrms_mid_var = 'MergedAzShear_2-5kmAGL'
#mrms_mid_var = 'MergedAzShear_3-6kmAGL'      #just for hurricane harvey!  Change back later
mrms_dz_var = 'MergedReflectivityQCComposite'  #pre-2019 variable name
#mrms_dz_var = 'MergedReflectivityDPQCLightComposite'
mrms_tot_low_var = 'MergedVelocity_Gradient_0-2kmAGL'
mrms_tot_mid_var = 'MergedVelocity_Gradient_2-5kmAGL'
mrms_mesh_var = 'MESH'

#################################### Generic Basemap Variables (to call as quickly as possible):  #####################################################

damage_files     = ''
area_thresh      = 1000.
resolution       = 'c'

#################################### User-Defined Variables:  #####################################################

dz_thresh        = 45.
dz_rad		 = 20000.	 
roi              = 3000.  		#radius of influence for Cressman scheme (m)
dz_roi           = 1000.  		#radius of influence for Cressman scheme (m) -> Smaller since DZ grid spacing 1/5 of AWS spacing
#range_min       = 15000.               #radius near radar to remove
range_min       = 5000.                 #radius near radar to remove
range_max       = 150000.               #radius away radar to remove
radius_max      = 3             	#grid point radius for maximum value filter
radius_gauss    = 2       	        #grid point radius of convolution operator
time_window     = 3
area_thresh2    = 5
eccent_thresh   = 0.4
kernel = gauss_kern(radius_gauss)

######################################################################################################
#################################### Read Data:  #####################################################
######################################################################################################

################### Get Grid info from NEWS-e summary file: ############################

try:
   newse_in = netCDF4.Dataset(newse_path, "r")
   print "Opening %s \n" % newse_path
except:
   print "%s does not exist! \n" %newse_path
   sys.exit(1)

cen_lat = newse_in.CEN_LAT                                   #center of domain latitude (dec deg)
cen_lon = newse_in.CEN_LON                                   #center of domain longitude (dec deg)
stand_lon = newse_in.STAND_LON                               #center lon of Lambert conformal projection
true_lat1 = newse_in.TRUE_LAT1                               #true lat value 1 for Lambert conformal conversion (dec deg)
true_lat2 = newse_in.TRUE_LAT2                               #true lat value 2 for Lambert conformal conversion (dec deg)

newse_lat = newse_in.variables["xlat"][:,:]                     #latitude (dec deg; Lambert conformal)
newse_lon = newse_in.variables["xlon"][:,:]                     #longitude (dec deg; Lambert conformal)

sw_lat_newse = newse_lat[0,0]
sw_lon_newse = newse_lon[0,0]
ne_lat_newse = newse_lat[-1,-1]
ne_lon_newse = newse_lon[-1,-1]

newse_in.close()
del newse_in

################### Create Basemap instance for grid conversion: ############################

newse_map = Basemap(llcrnrlon=sw_lon_newse, llcrnrlat=sw_lat_newse, urcrnrlon=ne_lon_newse, urcrnrlat=ne_lat_newse, projection='lcc', lat_1=true_lat1, lat_2=true_lat2, lat_0=cen_lat, lon_0=cen_lon, resolution = resolution, area_thresh = area_thresh)

newse_x_offset, newse_y_offset = newse_map(cen_lon, cen_lat)
newse_x, newse_y = newse_map(newse_lon[:], newse_lat[:])

newse_x = newse_x - newse_x_offset
newse_y = newse_y - newse_y_offset

newse_x_ravel = newse_x.ravel()
newse_y_ravel = newse_y.ravel()

print 'NEWSE GRID INFO: ', np.min(newse_x), np.max(newse_x), np.min(newse_y), np.max(newse_y)

################### Create KD Tree of NEWS-e x/y gridpoints: ############################

newse_x_y = np.dstack([newse_y_ravel, newse_x_ravel])[0]
newse_tree = scipy.spatial.cKDTree(newse_x_y)

################### Load locations of 88d sites from radar_sites object: ############################

rad_x, rad_y = newse_map(np.asarray(radar_sites.lon), np.asarray(radar_sites.lat))
rad_x = rad_x - newse_x_offset
rad_y = rad_y - newse_y_offset

################### Create KD Tree of radar site x/y locations: ############################

rad_x_y = np.dstack([rad_y, rad_x])[0]
rad_tree = scipy.spatial.cKDTree(rad_x_y)

################### Find points in radar blanking region: ############################

near_rad_points = newse_tree.query_ball_tree(rad_tree, range_min)
far_rad_points = newse_tree.query_ball_tree(rad_tree, range_max)

rad_mask = newse_x_ravel * 0. #initialize radmask

################### Set points in radar blanking region to 1 in radmask: ############################

for i in range(0,len(near_rad_points)):
   if (len(near_rad_points[i]) > 0.):
      rad_mask[i] = 1.

for i in range(0,len(far_rad_points)):
   if (len(far_rad_points[i]) == 0.):
      rad_mask[i] = 1.

rad_mask = rad_mask.reshape(newse_x.shape[0], newse_x.shape[1])

################### Read MRMS DZ/Az shear data: ############################

try:
   low_fin = netCDF4.Dataset(aws_low_file, "r")
   print "Opening %s \n" % aws_low_file
except:
   print "%s does not exist! \n" %aws_low_file
   sys.exit(1)

################### For first AWS file, build MRMS grid from sparse .netcdf file: ############################

nw_lat = low_fin.Latitude
nw_lon = low_fin.Longitude
lat_dy = low_fin.LatGridSpacing
lon_dx = low_fin.LonGridSpacing
lat_length = len(low_fin.dimensions["Lat"])
lon_length = len(low_fin.dimensions["Lon"])

lat_range = np.arange((nw_lat-((lat_length-1)*lat_dy)),(nw_lat+0.00001),lat_dy)
lon_range = np.arange(nw_lon,(nw_lon+((lon_length-0.99999)*lon_dx)),lon_dx)
xlon_full, xlat_full = np.meshgrid(lon_range, lat_range)

lon_range_indices = np.arange(0,lon_length)
lat_range_indices = np.arange(0,lat_length)

xlon_indices, xlat_indices = np.meshgrid(lon_range_indices, lat_range_indices)

min_lat = (np.abs(xlat_full[:,0]-np.min(newse_lat))).argmin()
max_lat = (np.abs(xlat_full[:,0]-np.max(newse_lat))).argmin()
min_lon = (np.abs(xlon_full[0,:]-np.min(newse_lon))).argmin()
max_lon = (np.abs(xlon_full[0,:]-np.max(newse_lon))).argmin()

xlat = xlat_full[min_lat:max_lat,min_lon:max_lon]
xlon = xlon_full[min_lat:max_lat,min_lon:max_lon]

sw_xlat = xlat[0,0]
sw_xlon = xlon[0,0]
ne_xlat = xlat[-1,-1]
ne_xlon = xlon[-1,-1]

print 'WDSSII LAT/LON: ', sw_xlat, sw_xlon, ne_xlat, ne_xlon

map = Basemap(llcrnrlon=sw_xlon, llcrnrlat=sw_xlat, urcrnrlon=ne_xlon, urcrnrlat=ne_xlat, projection='lcc', lat_1=true_lat1, lat_2=true_lat2, lat_0=cen_lat, lon_0=cen_lon, resolution = resolution, area_thresh = area_thresh)

x_offset, y_offset = map(cen_lon, cen_lat)
x, y = map(xlon[:], xlat[:])

x_full, y_full = map(xlon_full[:], xlat_full[:])

x = x - x_offset
y = y - y_offset

print 'WDSII X/Y: ', np.min(x), np.max(x), np.min(y), np.max(y)

x_full = x_full - x_offset
y_full = y_full - y_offset

min_x, min_y = map(lon_range[min_lon], lat_range[min_lat])
max_x, max_y = map(lon_range[max_lon], lat_range[max_lat])

min_x = min_x - x_offset
max_x = max_x - x_offset
min_y = min_y - y_offset
max_y = max_y - y_offset

print 'MAX X/Y: ', min_x, max_x, min_y, max_y

################### Initialize interpolated variables and output .nc file: ############################

low_cress = np.zeros((newse_lat.shape[0],newse_lat.shape[1]))

try:
   fout = netCDF4.Dataset(out_file, "w")
except:
   print "Could not create %s!\n" % out_file

fout.createDimension('NX', newse_lat.shape[1])
fout.createDimension('NY', newse_lat.shape[0])

fout.createVariable('XLAT', 'f4', ('NY','NX',))
fout.createVariable('XLON', 'f4', ('NY','NX',))
fout.createVariable('RADMASK', 'f4', ('NY','NX',))
fout.createVariable('LOW_CRESSMAN', 'f4', ('NY','NX',))
fout.createVariable('MID_CRESSMAN', 'f4', ('NY','NX',))
fout.createVariable('DZ_CRESSMAN', 'f4', ('NY','NX',))
fout.createVariable('TOT_LOW_CRESSMAN', 'f4', ('NY','NX',))
fout.createVariable('TOT_MID_CRESSMAN', 'f4', ('NY','NX',))
fout.createVariable('MESH_CRESSMAN', 'f4', ('NY','NX',))

fout.variables['XLAT'][:] = newse_lat
fout.variables['XLON'][:] = newse_lon
fout.variables['RADMASK'][:] = rad_mask

########## Read 0-2 km MRMS Az. shear data, remove data outside NEWS-e domain for speed, and flip lat coordinates to be compatable: #######
############ Control added for WDSSii files (Az. Shear only) that use 'LatLonGrid' instead of 'SparseLatLonGrid' (11/2016) ########

grid_type = low_fin.DataType

if (grid_type[0:2] == 'Sp'): #If SparseLatLonGrid 
   pixel = len(low_fin.dimensions["pixel"])
   if (pixel > 0):
      pixel_x_full = low_fin.variables["pixel_x"][:]
      pixel_y_full = low_fin.variables["pixel_y"][:]
      pixel_value_full = low_fin.variables[mrms_low_var][:]
      pixel_count_full = low_fin.variables["pixel_count"][:]

      pixel_x_full = pixel_x_full[pixel_value_full > 0.0001]  
      pixel_y_full = pixel_y_full[pixel_value_full > 0.0001]
      pixel_count_full = pixel_count_full[pixel_value_full > 0.0001]
      pixel_value_full = pixel_value_full[pixel_value_full > 0.0001]

      pixel_x_transpose = xlat_indices.shape[0] - pixel_x_full
      pixel_x = pixel_x_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
      pixel_y = pixel_y_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
      pixel_value = pixel_value_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
      pixel_count = pixel_count_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
#      print 'aws low pixel count orig: ', len(pixel_value)

      pixel_value_temp = pixel_value[pixel_count > 1]
      pixel_x_temp = pixel_x[pixel_count > 1]
      pixel_y_temp = pixel_y[pixel_count > 1]
      pixel_count_temp = pixel_count[pixel_count > 1]

      for i in range(0, len(pixel_count_temp)):
         for j in range(1, pixel_count_temp[i]):
            temp_y_index = pixel_y_temp[i] + j
            temp_x_index = pixel_x_temp[i]            
            if (temp_y_index < max_lon):
               pixel_x = np.append(pixel_x, temp_x_index)
               pixel_y = np.append(pixel_y, temp_y_index)
               pixel_value = np.append(pixel_value, pixel_value_temp[i])
#      print 'aws low pixel count new: ', len(pixel_value), len(pixel_count_temp)

      pixel_x = abs(pixel_x - (xlat_full.shape[0] - min_lat)) #... tortured way of flipping lat values, but it works
      pixel_y = pixel_y - min_lon

      mrms_low_pixel_x_val = x[pixel_x, pixel_y]  
      mrms_low_pixel_y_val = y[pixel_x, pixel_y]
#      mrms_low_pixel_x_val = x_full[pixel_x, pixel_y]  
#      mrms_low_pixel_y_val = y_full[pixel_x, pixel_y]
      mrms_low_pixel_value = pixel_value
      mrms_low_x_y = np.dstack([mrms_low_pixel_y_val, mrms_low_pixel_x_val])[0] #KD Tree searchable index of mrms observations
elif (grid_type[0:2] == 'La'): #if LatLonGrid
   pixel_value_full = low_fin.variables[mrms_low_var][:]
#   print 'shapes ... ', pixel_value_full.shape, xlon_indices.shape, xlat_indices.shape 
   pixel_value_full = pixel_value_full.ravel()
   pixel_x_full = xlat_indices.ravel() 
   pixel_y_full = xlon_indices.ravel()

   pixel_x_full = pixel_x_full[pixel_value_full > 0.0001]  
   pixel_y_full = pixel_y_full[pixel_value_full > 0.0001]
   pixel_value_full = pixel_value_full[pixel_value_full > 0.0001]

   pixel_x_transpose = xlat_full.shape[0] - pixel_x_full
#   pixel_x_transpose = pixel_x_full #xlat_full.shape[0] - pixel_x_full
   pixel_x = pixel_x_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
   pixel_y = pixel_y_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
   mrms_low_pixel_value = pixel_value_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
   pixel_x = abs(pixel_x - (xlat_full.shape[0] - min_lat)) #... tortured way of flipping lat values, but it works
   pixel_y = pixel_y - min_lon
   mrms_low_pixel_x_val = x[pixel_x, pixel_y]
   mrms_low_pixel_y_val = y[pixel_x, pixel_y]
   mrms_low_x_y = np.dstack([mrms_low_pixel_y_val, mrms_low_pixel_x_val])[0] #KD Tree searchable index of mrms observations
else: 
   print 'aws low unknown grid type!!!!'

low_fin.close()
del low_fin

########## Read 2-5 km MRMS Az. shear data, remove data outside NEWS-e domain for speed, and flip lat coordinates to be compatable: #######
############ Control added for WDSSii files (Az. Shear only) that use 'LatLonGrid' instead of 'SparseLatLonGrid' (11/2016) ########

try:
   mid_fin = netCDF4.Dataset(aws_mid_file, "r")
   print "Opening %s \n" % aws_mid_file
except:
   print "%s does not exist! \n" %aws_mid_file
   sys.exit(1)

grid_type = mid_fin.DataType

if (grid_type[0:2] == 'Sp'): #If SparseLatLonGrid 
   pixel = len(mid_fin.dimensions["pixel"])
   if (pixel > 0):
      pixel_x_full = mid_fin.variables["pixel_x"][:]
      pixel_y_full = mid_fin.variables["pixel_y"][:]
      pixel_value_full = mid_fin.variables[mrms_mid_var][:]
      pixel_count_full = mid_fin.variables["pixel_count"][:]

      pixel_x_full = pixel_x_full[pixel_value_full > 0.0001]  
      pixel_y_full = pixel_y_full[pixel_value_full > 0.0001]
      pixel_count_full = pixel_count_full[pixel_value_full > 0.0001]
      pixel_value_full = pixel_value_full[pixel_value_full > 0.0001]
    
      pixel_x_transpose = xlat_indices.shape[0] - pixel_x_full

      pixel_x = pixel_x_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
      pixel_y = pixel_y_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
      pixel_value = pixel_value_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
      pixel_count = pixel_count_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]

#      print 'aws mid pixel count orig: ', len(pixel_value)
  
      pixel_value_temp = pixel_value[pixel_count > 1]
      pixel_x_temp = pixel_x[pixel_count > 1]
      pixel_y_temp = pixel_y[pixel_count > 1]
      pixel_count_temp = pixel_count[pixel_count > 1]

      for i in range(0, len(pixel_count_temp)):
         for j in range(1, pixel_count_temp[i]):
            temp_y_index = pixel_y_temp[i] + j
            temp_x_index = pixel_x_temp[i]            
            if (temp_y_index < max_lon):
               pixel_x = np.append(pixel_x, temp_x_index)
               pixel_y = np.append(pixel_y, temp_y_index)
               pixel_value = np.append(pixel_value, pixel_value_temp[i])

#      print 'aws mid pixel count new: ', len(pixel_value), len(pixel_count_temp)

      pixel_x = abs(pixel_x - (xlat_full.shape[0] - min_lat)) #... tortured way of flipping lat values, but it works
      pixel_y = pixel_y - min_lon

      mrms_mid_pixel_x_val = x[pixel_x, pixel_y]  
      mrms_mid_pixel_y_val = y[pixel_x, pixel_y]
#      mrms_mid_pixel_x_val = x_full[pixel_x, pixel_y]  
#      mrms_mid_pixel_y_val = y_full[pixel_x, pixel_y]
      mrms_mid_pixel_value = pixel_value

      mrms_mid_x_y = np.dstack([mrms_mid_pixel_y_val, mrms_mid_pixel_x_val])[0] #KD Tree searchable index of mrms observations

elif (grid_type[0:2] == 'La'): #if LatLonGrid
   pixel_value_full = mid_fin.variables[mrms_mid_var][:]
#   print 'shapes ... ', pixel_value_full.shape, xlon_indices.shape, xlat_indices.shape 

   pixel_value_full = pixel_value_full.ravel()
   pixel_x_full = xlat_indices.ravel() 
   pixel_y_full = xlon_indices.ravel()
   pixel_x_full = pixel_x_full[pixel_value_full > 0.0001]  
   pixel_y_full = pixel_y_full[pixel_value_full > 0.0001]
   pixel_value_full = pixel_value_full[pixel_value_full > 0.0001] 

   pixel_x_transpose = xlat_full.shape[0] - pixel_x_full
#   pixel_x_transpose = pixel_x_full #xlat_full.shape[0] - pixel_x_full
   pixel_x = pixel_x_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
   pixel_y = pixel_y_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
   mrms_mid_pixel_value = pixel_value_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]

   pixel_x = abs(pixel_x - (xlat_full.shape[0] - min_lat)) #... tortured way of flipping lat values, but it works
   pixel_y = pixel_y - min_lon
   mrms_mid_pixel_x_val = x[pixel_x, pixel_y]
   mrms_mid_pixel_y_val = y[pixel_x, pixel_y]

   mrms_mid_x_y = np.dstack([mrms_mid_pixel_y_val, mrms_mid_pixel_x_val])[0] #KD Tree searchable index of mrms observations

else: 
   print 'aws mid unknown grid type!!!!'

mid_fin.close()
del mid_fin

########## Read MRMS Composite Reflectivity data, remove data outside NEWS-e domain for speed, and flip lat coordinates to be compatable: #######
############ Control added for WDSSii files (Az. Shear only) that use 'LatLonGrid' instead of 'SparseLatLonGrid' (11/2016) ########

try:
   dz_fin = netCDF4.Dataset(dz_file, "r")
   print "Opening %s \n" % dz_file
except:
   print "%s does not exist! \n" %dz_file
   sys.exit(1)

################### For first DZ file, build DZ MRMS grid from sparse .netcdf file (DIFFERENT THAN AWS GRID!!!): ############################

dz_nw_lat = dz_fin.Latitude
dz_nw_lon = dz_fin.Longitude
dz_lat_dy = dz_fin.LatGridSpacing
dz_lon_dx = dz_fin.LonGridSpacing
dz_lat_length = len(dz_fin.dimensions["Lat"])
dz_lon_length = len(dz_fin.dimensions["Lon"])

dz_lat_range = np.arange((dz_nw_lat-((dz_lat_length-1)*dz_lat_dy)),(dz_nw_lat+0.00001),dz_lat_dy)
dz_lon_range = np.arange(dz_nw_lon,(dz_nw_lon+((dz_lon_length-0.99999)*dz_lon_dx)),dz_lon_dx)
dz_xlon_full, dz_xlat_full = np.meshgrid(dz_lon_range, dz_lat_range)

dz_lon_range_indices = np.arange(0,dz_lon_length)
dz_lat_range_indices = np.arange(0,dz_lat_length)

dz_xlon_indices, dz_xlat_indices = np.meshgrid(dz_lon_range_indices, dz_lat_range_indices)

dz_min_lat = (np.abs(dz_xlat_full[:,0]-np.min(newse_lat))).argmin()
dz_max_lat = (np.abs(dz_xlat_full[:,0]-np.max(newse_lat))).argmin()
dz_min_lon = (np.abs(dz_xlon_full[0,:]-np.min(newse_lon))).argmin()
dz_max_lon = (np.abs(dz_xlon_full[0,:]-np.max(newse_lon))).argmin()

dz_xlat = dz_xlat_full[dz_min_lat:dz_max_lat,dz_min_lon:dz_max_lon]
dz_xlon = dz_xlon_full[dz_min_lat:dz_max_lat,dz_min_lon:dz_max_lon]

dz_sw_xlat = dz_xlat[0,0]
dz_sw_xlon = dz_xlon[0,0]
dz_ne_xlat = dz_xlat[-1,-1]
dz_ne_xlon = dz_xlon[-1,-1]

dz_map = Basemap(llcrnrlon=dz_sw_xlon, llcrnrlat=dz_sw_xlat, urcrnrlon=dz_ne_xlon, urcrnrlat=dz_ne_xlat, projection='lcc', lat_1=true_lat1, lat_2=true_lat2, lat_0=cen_lat, lon_0=cen_lon, resolution = resolution, area_thresh = area_thresh)

dz_x_offset, dz_y_offset = dz_map(cen_lon, cen_lat)
dz_x, dz_y = dz_map(dz_xlon[:], dz_xlat[:])

dz_x_full, dz_y_full = dz_map(dz_xlon_full[:], dz_xlat_full[:])

dz_x_full = dz_x_full - dz_x_offset
dz_y_full = dz_y_full - dz_y_offset

dz_x = dz_x - dz_x_offset
dz_y = dz_y - dz_y_offset

dz_min_x, dz_min_y = dz_map(dz_lon_range[dz_min_lon], dz_lat_range[dz_min_lat])
dz_max_x, dz_max_y = dz_map(dz_lon_range[dz_max_lon], dz_lat_range[dz_max_lat])

dz_min_x = dz_min_x - dz_x_offset
dz_max_x = dz_max_x - dz_x_offset
dz_min_y = dz_min_y - dz_y_offset
dz_max_y = dz_max_y - dz_y_offset

grid_type = dz_fin.DataType

if (grid_type[0:2] == 'Sp'): #If SparseLatLonGrid 
   pixel = len(dz_fin.dimensions["pixel"])
   if (pixel > 0):
      pixel_x_full = dz_fin.variables["pixel_x"][:]
      pixel_y_full = dz_fin.variables["pixel_y"][:]
      pixel_value_full = dz_fin.variables[mrms_dz_var][:]
      pixel_count_full = dz_fin.variables["pixel_count"][:]

      pixel_x_full = pixel_x_full[pixel_value_full > 0.01]  
      pixel_y_full = pixel_y_full[pixel_value_full > 0.01]
      pixel_count_full = pixel_count_full[pixel_value_full > 0.01]
      pixel_value_full = pixel_value_full[pixel_value_full > 0.01]

      pixel_x_transpose = dz_xlat_indices.shape[0] - pixel_x_full
      pixel_x = pixel_x_full[(pixel_x_transpose > dz_min_lat) & (pixel_x_transpose < dz_max_lat) & (pixel_y_full > dz_min_lon) & (pixel_y_full < dz_max_lon) & (pixel_value_full > dz_thresh_1)]
      pixel_y = pixel_y_full[(pixel_x_transpose > dz_min_lat) & (pixel_x_transpose < dz_max_lat) & (pixel_y_full > dz_min_lon) & (pixel_y_full < dz_max_lon) & (pixel_value_full > dz_thresh_1)]
      pixel_value = pixel_value_full[(pixel_x_transpose > dz_min_lat) & (pixel_x_transpose < dz_max_lat) & (pixel_y_full > dz_min_lon) & (pixel_y_full < dz_max_lon) & (pixel_value_full > dz_thresh_1)]
      pixel_count = pixel_count_full[(pixel_x_transpose > dz_min_lat) & (pixel_x_transpose < dz_max_lat) & (pixel_y_full > dz_min_lon) & (pixel_y_full < dz_max_lon) & (pixel_value_full > dz_thresh_1)]
#      print 'dz pixel count orig: ', len(pixel_value), len(pixel_x_full), len(pixel_value_full), len(pixel_x)
      pixel_value_temp = pixel_value[pixel_count > 1]
      pixel_x_temp = pixel_x[pixel_count > 1]
      pixel_y_temp = pixel_y[pixel_count > 1]
      pixel_count_temp = pixel_count[pixel_count > 1]
      for i in range(0, len(pixel_count_temp)):
         for j in range(1, pixel_count_temp[i]):
            temp_y_index = pixel_y_temp[i] + j
            temp_x_index = pixel_x_temp[i]            
            if (temp_y_index < dz_max_lon):
               pixel_x = np.append(pixel_x, temp_x_index)
               pixel_y = np.append(pixel_y, temp_y_index)
               pixel_value = np.append(pixel_value, pixel_value_temp[i])
#      print 'dz pixel count new: ', len(pixel_value), len(pixel_count_temp)
      pixel_x = abs(pixel_x - (dz_xlat_full.shape[0] - dz_min_lat)) #... tortured way of flipping lat values, but it works
      pixel_y = pixel_y - dz_min_lon

      mrms_dz_pixel_x_val = dz_x[pixel_x, pixel_y]  
      mrms_dz_pixel_y_val = dz_y[pixel_x, pixel_y]
#      mrms_dz_pixel_x_val = dz_x_full[pixel_x, pixel_y]  
#      mrms_dz_pixel_y_val = dz_y_full[pixel_x, pixel_y]
      mrms_dz_pixel_value = pixel_value
      mrms_dz_x_y = np.dstack([mrms_dz_pixel_y_val, mrms_dz_pixel_x_val])[0] #KD Tree searchable index of mrms observations
elif (grid_type[0:2] == 'La'): #if LatLonGrid
   pixel_value_full = dz_fin.variables[mrms_dz_var][:]
   pixel_value_full = pixel_value_full.ravel()
   pixel_x_full = dz_xlat_indices.ravel() 
   pixel_y_full = dz_xlon_indices.ravel()
   pixel_x_full = pixel_x_full[pixel_value_full > 0.01]  
   pixel_y_full = pixel_y_full[pixel_value_full > 0.01]
   pixel_value_full = pixel_value_full[pixel_value_full > 0.01]

   pixel_x_transpose = dz_xlat_full.shape[0] - pixel_x_full
#   pixel_x_transpose = pixel_x_full #xlat_full.shape[0] - pixel_x_full
   pixel_x = pixel_x_full[(pixel_x_transpose > dz_min_lat) & (pixel_x_transpose < dz_max_lat) & (pixel_y_full > dz_min_lon) & (pixel_y_full < dz_max_lon) & (pixel_value_full > dz_thresh_1)]
   pixel_y = pixel_y_full[(pixel_x_transpose > dz_min_lat) & (pixel_x_transpose < dz_max_lat) & (pixel_y_full > dz_min_lon) & (pixel_y_full < dz_max_lon) & (pixel_value_full > dz_thresh_1)]

   mrms_dz_pixel_value = pixel_value_full[(pixel_x_transpose > dz_min_lat) & (pixel_x_transpose < dz_max_lat) & (pixel_y_full > dz_min_lon) & (pixel_y_full < dz_max_lon) & (pixel_value_full > dz_thresh_1)]
   pixel_x = abs(pixel_x - (dz_xlat_full.shape[0] - dz_min_lat)) #... tortured way of flipping lat values, but it works
   pixel_y = pixel_y - dz_min_lon

   mrms_dz_pixel_x_val = dz_x[pixel_x, pixel_y]
   mrms_dz_pixel_y_val = dz_y[pixel_x, pixel_y]
   mrms_dz_x_y = np.dstack([mrms_dz_pixel_y_val, mrms_dz_pixel_x_val])[0] #KD Tree searchable index of mrms observations
else: 
   print 'dz unknown grid type!!!!'

dz_fin.close()
del dz_fin

########## Read 0-2 km Total shear data, remove data outside NEWS-e domain for speed, and flip lat coordinates to be compatable: #######
############ Control added for WDSSii files (Az. Shear only) that use 'LatLonGrid' instead of 'SparseLatLonGrid' (11/2016) ########

try:
   low_tot_fin = netCDF4.Dataset(tot_low_file, "r")
   print "Opening %s \n" % tot_low_file
except:
   print "%s does not exist! \n" %tot_low_file
   sys.exit(1)

grid_type = low_tot_fin.DataType

if (grid_type[0:2] == 'Sp'): #If SparseLatLonGrid 
   pixel = len(low_tot_fin.dimensions["pixel"])
   if (pixel > 0):
      pixel_x_full = low_tot_fin.variables["pixel_x"][:]
      pixel_y_full = low_tot_fin.variables["pixel_y"][:]
      pixel_value_full = low_tot_fin.variables[mrms_tot_low_var][:]
      pixel_count_full = low_tot_fin.variables["pixel_count"][:]

      pixel_x_full = pixel_x_full[pixel_value_full > 0.0001]  
      pixel_y_full = pixel_y_full[pixel_value_full > 0.0001]
      pixel_count_full = pixel_count_full[pixel_value_full > 0.0001]
      pixel_value_full = pixel_value_full[pixel_value_full > 0.0001]

      pixel_x_transpose = xlat_indices.shape[0] - pixel_x_full
      pixel_x = pixel_x_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
      pixel_y = pixel_y_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
      pixel_value = pixel_value_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
      pixel_count = pixel_count_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
#      print 'aws low pixel count orig: ', len(pixel_value)

      pixel_value_temp = pixel_value[pixel_count > 1]
      pixel_x_temp = pixel_x[pixel_count > 1]
      pixel_y_temp = pixel_y[pixel_count > 1]
      pixel_count_temp = pixel_count[pixel_count > 1]

      for i in range(0, len(pixel_count_temp)):
         for j in range(1, pixel_count_temp[i]):
            temp_y_index = pixel_y_temp[i] + j
            temp_x_index = pixel_x_temp[i]            
            if (temp_y_index < max_lon):
               pixel_x = np.append(pixel_x, temp_x_index)
               pixel_y = np.append(pixel_y, temp_y_index)
               pixel_value = np.append(pixel_value, pixel_value_temp[i])
#      print 'aws low pixel count new: ', len(pixel_value), len(pixel_count_temp)

      pixel_x = abs(pixel_x - (xlat_full.shape[0] - min_lat)) #... tortured way of flipping lat values, but it works
      pixel_y = pixel_y - min_lon

      mrms_tot_low_pixel_x_val = x[pixel_x, pixel_y]  
      mrms_tot_low_pixel_y_val = y[pixel_x, pixel_y]
#      mrms_low_pixel_x_val = x_full[pixel_x, pixel_y]  
#      mrms_low_pixel_y_val = y_full[pixel_x, pixel_y]
      mrms_tot_low_pixel_value = pixel_value
      mrms_tot_low_x_y = np.dstack([mrms_tot_low_pixel_y_val, mrms_tot_low_pixel_x_val])[0] #KD Tree searchable index of mrms observations
elif (grid_type[0:2] == 'La'): #if LatLonGrid
   pixel_value_full = low_tot_fin.variables[mrms_tot_low_var][:]
#   print 'shapes ... ', pixel_value_full.shape, xlon_indices.shape, xlat_indices.shape 
   pixel_value_full = pixel_value_full.ravel()
   pixel_x_full = xlat_indices.ravel() 
   pixel_y_full = xlon_indices.ravel()

   pixel_x_full = pixel_x_full[pixel_value_full > 0.0001]  
   pixel_y_full = pixel_y_full[pixel_value_full > 0.0001]
   pixel_value_full = pixel_value_full[pixel_value_full > 0.0001]

   pixel_x_transpose = xlat_full.shape[0] - pixel_x_full
#   pixel_x_transpose = pixel_x_full #xlat_full.shape[0] - pixel_x_full
   pixel_x = pixel_x_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
   pixel_y = pixel_y_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
   mrms_tot_low_pixel_value = pixel_value_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
   pixel_x = abs(pixel_x - (xlat_full.shape[0] - min_lat)) #... tortured way of flipping lat values, but it works
   pixel_y = pixel_y - min_lon
   mrms_tot_low_pixel_x_val = x[pixel_x, pixel_y]
   mrms_tot_low_pixel_y_val = y[pixel_x, pixel_y]
   mrms_tot_low_x_y = np.dstack([mrms_tot_low_pixel_y_val, mrms_tot_low_pixel_x_val])[0] #KD Tree searchable index of mrms observations
else: 
   print 'tot low unknown grid type!!!!'

low_tot_fin.close()
del low_tot_fin

########## Read 2-5 km Total shear data, remove data outside NEWS-e domain for speed, and flip lat coordinates to be compatable: #######
############ Control added for WDSSii files (Az. Shear only) that use 'LatLonGrid' instead of 'SparseLatLonGrid' (11/2016) ########

try:
   mid_tot_fin = netCDF4.Dataset(tot_mid_file, "r")
   print "Opening %s \n" % tot_mid_file
except:
   print "%s does not exist! \n" %tot_mid_file
   sys.exit(1)

grid_type = mid_tot_fin.DataType

if (grid_type[0:2] == 'Sp'): #If SparseLatLonGrid 
   pixel = len(mid_tot_fin.dimensions["pixel"])
   if (pixel > 0):
      pixel_x_full = mid_tot_fin.variables["pixel_x"][:]
      pixel_y_full = mid_tot_fin.variables["pixel_y"][:]
      pixel_value_full = mid_tot_fin.variables[mrms_tot_mid_var][:]
      pixel_count_full = mid_tot_fin.variables["pixel_count"][:]

      pixel_x_full = pixel_x_full[pixel_value_full > 0.0001]  
      pixel_y_full = pixel_y_full[pixel_value_full > 0.0001]
      pixel_count_full = pixel_count_full[pixel_value_full > 0.0001]
      pixel_value_full = pixel_value_full[pixel_value_full > 0.0001]

      pixel_x_transpose = xlat_indices.shape[0] - pixel_x_full
      pixel_x = pixel_x_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
      pixel_y = pixel_y_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
      pixel_value = pixel_value_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
      pixel_count = pixel_count_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
#      print 'aws low pixel count orig: ', len(pixel_value)

      pixel_value_temp = pixel_value[pixel_count > 1]
      pixel_x_temp = pixel_x[pixel_count > 1]
      pixel_y_temp = pixel_y[pixel_count > 1]
      pixel_count_temp = pixel_count[pixel_count > 1]

      for i in range(0, len(pixel_count_temp)):
         for j in range(1, pixel_count_temp[i]):
            temp_y_index = pixel_y_temp[i] + j
            temp_x_index = pixel_x_temp[i]            
            if (temp_y_index < max_lon):
               pixel_x = np.append(pixel_x, temp_x_index)
               pixel_y = np.append(pixel_y, temp_y_index)
               pixel_value = np.append(pixel_value, pixel_value_temp[i])
#      print 'aws low pixel count new: ', len(pixel_value), len(pixel_count_temp)

      pixel_x = abs(pixel_x - (xlat_full.shape[0] - min_lat)) #... tortured way of flipping lat values, but it works
      pixel_y = pixel_y - min_lon

      mrms_tot_mid_pixel_x_val = x[pixel_x, pixel_y]  
      mrms_tot_mid_pixel_y_val = y[pixel_x, pixel_y]
#      mrms_low_pixel_x_val = x_full[pixel_x, pixel_y]  
#      mrms_low_pixel_y_val = y_full[pixel_x, pixel_y]
      mrms_tot_mid_pixel_value = pixel_value
      mrms_tot_mid_x_y = np.dstack([mrms_tot_mid_pixel_y_val, mrms_tot_mid_pixel_x_val])[0] #KD Tree searchable index of mrms observations
elif (grid_type[0:2] == 'La'): #if LatLonGrid
   pixel_value_full = mid_tot_fin.variables[mrms_tot_mid_var][:]
#   print 'shapes ... ', pixel_value_full.shape, xlon_indices.shape, xlat_indices.shape 
   pixel_value_full = pixel_value_full.ravel()
   pixel_x_full = xlat_indices.ravel() 
   pixel_y_full = xlon_indices.ravel()

   pixel_x_full = pixel_x_full[pixel_value_full > 0.0001]  
   pixel_y_full = pixel_y_full[pixel_value_full > 0.0001]
   pixel_value_full = pixel_value_full[pixel_value_full > 0.0001]

   pixel_x_transpose = xlat_full.shape[0] - pixel_x_full
#   pixel_x_transpose = pixel_x_full #xlat_full.shape[0] - pixel_x_full
   pixel_x = pixel_x_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
   pixel_y = pixel_y_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
   mrms_tot_mid_pixel_value = pixel_value_full[(pixel_x_transpose > min_lat) & (pixel_x_transpose < max_lat) & (pixel_y_full > min_lon) & (pixel_y_full < max_lon) & (pixel_value_full > aws_thresh_1)]
   pixel_x = abs(pixel_x - (xlat_full.shape[0] - min_lat)) #... tortured way of flipping lat values, but it works
   pixel_y = pixel_y - min_lon
   mrms_tot_mid_pixel_x_val = x[pixel_x, pixel_y]
   mrms_tot_mid_pixel_y_val = y[pixel_x, pixel_y]
   mrms_tot_mid_x_y = np.dstack([mrms_tot_mid_pixel_y_val, mrms_tot_mid_pixel_x_val])[0] #KD Tree searchable index of mrms observations
else: 
   print 'tot mid unknown grid type!!!!'

mid_tot_fin.close()
del mid_tot_fin

########## Read MESH data, remove data outside NEWS-e domain for speed, and flip lat coordinates to be compatable: #######
############ Control added for WDSSii files (Az. Shear only) that use 'LatLonGrid' instead of 'SparseLatLonGrid' (11/2016) ########

try:
   mesh_fin = netCDF4.Dataset(mesh_file, "r")
   print "Opening %s \n" % mesh_file
except:
   print "%s does not exist! \n" %mesh_file
   sys.exit(1)

grid_type = mesh_fin.DataType

if (grid_type[0:2] == 'Sp'): #If SparseLatLonGrid 
   pixel = len(mesh_fin.dimensions["pixel"])
   if (pixel > 0):
      pixel_x_full = mesh_fin.variables["pixel_x"][:]
      pixel_y_full = mesh_fin.variables["pixel_y"][:]
      pixel_value_full = mesh_fin.variables[mrms_mesh_var][:]
      pixel_count_full = mesh_fin.variables["pixel_count"][:]

      pixel_x_full = pixel_x_full[pixel_value_full > 0.01]
      pixel_y_full = pixel_y_full[pixel_value_full > 0.01]
      pixel_count_full = pixel_count_full[pixel_value_full > 0.01]
      pixel_value_full = pixel_value_full[pixel_value_full > 0.01]

      pixel_x_transpose = dz_xlat_indices.shape[0] - pixel_x_full
      pixel_x = pixel_x_full[(pixel_x_transpose > dz_min_lat) & (pixel_x_transpose < dz_max_lat) & (pixel_y_full > dz_min_lon) & (pixel_y_full < dz_max_lon) & (pixel_value_full > aws_thresh_1)]
      pixel_y = pixel_y_full[(pixel_x_transpose > dz_min_lat) & (pixel_x_transpose < dz_max_lat) & (pixel_y_full > dz_min_lon) & (pixel_y_full < dz_max_lon) & (pixel_value_full > aws_thresh_1)]
      pixel_value = pixel_value_full[(pixel_x_transpose > dz_min_lat) & (pixel_x_transpose < dz_max_lat) & (pixel_y_full > dz_min_lon) & (pixel_y_full < dz_max_lon) & (pixel_value_full > aws_thresh_1)]
      pixel_count = pixel_count_full[(pixel_x_transpose > dz_min_lat) & (pixel_x_transpose < dz_max_lat) & (pixel_y_full > dz_min_lon) & (pixel_y_full < dz_max_lon) & (pixel_value_full > aws_thresh_1)]
#      print 'dz pixel count orig: ', len(pixel_value), len(pixel_x_full), len(pixel_value_full), len(pixel_x)
      pixel_value_temp = pixel_value[pixel_count > 1]
      pixel_x_temp = pixel_x[pixel_count > 1]
      pixel_y_temp = pixel_y[pixel_count > 1]
      pixel_count_temp = pixel_count[pixel_count > 1]
      for i in range(0, len(pixel_count_temp)):
         for j in range(1, pixel_count_temp[i]):
            temp_y_index = pixel_y_temp[i] + j
            temp_x_index = pixel_x_temp[i]
            if (temp_y_index < dz_max_lon):
               pixel_x = np.append(pixel_x, temp_x_index)
               pixel_y = np.append(pixel_y, temp_y_index)
               pixel_value = np.append(pixel_value, pixel_value_temp[i])
#      print 'dz pixel count new: ', len(pixel_value), len(pixel_count_temp)
      pixel_x = abs(pixel_x - (dz_xlat_full.shape[0] - dz_min_lat)) #... tortured way of flipping lat values, but it works
      pixel_y = pixel_y - dz_min_lon

      mrms_mesh_pixel_x_val = dz_x[pixel_x, pixel_y]
      mrms_mesh_pixel_y_val = dz_y[pixel_x, pixel_y]
#      mrms_dz_pixel_x_val = dz_x_full[pixel_x, pixel_y]  
#      mrms_dz_pixel_y_val = dz_y_full[pixel_x, pixel_y]
      mrms_mesh_pixel_value = pixel_value
      mrms_mesh_x_y = np.dstack([mrms_mesh_pixel_y_val, mrms_mesh_pixel_x_val])[0] #KD Tree searchable index of mrms observations
elif (grid_type[0:2] == 'La'): #if LatLonGrid
   pixel_value_full = mesh_fin.variables[mrms_mesh_var][:]
   pixel_value_full = pixel_value_full.ravel()
   pixel_x_full = dz_xlat_indices.ravel()
   pixel_y_full = dz_xlon_indices.ravel()
   pixel_x_full = pixel_x_full[pixel_value_full > 0.01]
   pixel_y_full = pixel_y_full[pixel_value_full > 0.01]
   pixel_value_full = pixel_value_full[pixel_value_full > 0.01]

   pixel_x_transpose = dz_xlat_full.shape[0] - pixel_x_full
#   pixel_x_transpose = pixel_x_full #xlat_full.shape[0] - pixel_x_full
   pixel_x = pixel_x_full[(pixel_x_transpose > dz_min_lat) & (pixel_x_transpose < dz_max_lat) & (pixel_y_full > dz_min_lon) & (pixel_y_full < dz_max_lon) & (pixel_value_full > aws_thresh_1)]
   pixel_y = pixel_y_full[(pixel_x_transpose > dz_min_lat) & (pixel_x_transpose < dz_max_lat) & (pixel_y_full > dz_min_lon) & (pixel_y_full < dz_max_lon) & (pixel_value_full > aws_thresh_1)]

   mrms_mesh_pixel_value = pixel_value_full[(pixel_x_transpose > dz_min_lat) & (pixel_x_transpose < dz_max_lat) & (pixel_y_full > dz_min_lon) & (pixel_y_full < dz_max_lon) & (pixel_value_full > aws_thresh_1)]
   pixel_x = abs(pixel_x - (dz_xlat_full.shape[0] - dz_min_lat)) #... tortured way of flipping lat values, but it works
   pixel_y = pixel_y - dz_min_lon

   mrms_mesh_pixel_x_val = dz_x[pixel_x, pixel_y]
   mrms_mesh_pixel_y_val = dz_y[pixel_x, pixel_y]
   mrms_mesh_x_y = np.dstack([mrms_mesh_pixel_y_val, mrms_mesh_pixel_x_val])[0] #KD Tree searchable index of mrms observations
else:
   print 'mesh unknown grid type!!!!'

mesh_fin.close()
del mesh_fin


##########################################################################################################################

################### If mrms obs are available within 5 min of time bin, Cressman them: ###############

if ((len(mrms_low_x_y) > 0) and (len(mrms_dz_x_y) > 0)):
   mrms_low_tree = scipy.spatial.cKDTree(mrms_low_x_y) 
else: 
   oban_low_var = newse_x_ravel * 0.
   oban_low_var = oban_low_var.reshape(newse_x.shape[0], newse_x.shape[1])

if ((len(mrms_mid_x_y) > 0) and (len(mrms_dz_x_y) > 0)):
   mrms_mid_tree = scipy.spatial.cKDTree(mrms_mid_x_y) 
else: 
   oban_mid_var = newse_x_ravel * 0.
   oban_mid_var = oban_mid_var.reshape(newse_x.shape[0], newse_x.shape[1])

if (len(mrms_dz_x_y) > 0):
   mrms_dz_tree = scipy.spatial.cKDTree(mrms_dz_x_y)
else: 
   oban_dz_var = newse_x_ravel * 0.
   oban_dz_var = oban_dz_var.reshape(newse_x.shape[0], newse_x.shape[1])

if (len(mrms_tot_low_x_y) > 0):
   mrms_tot_low_tree = scipy.spatial.cKDTree(mrms_tot_low_x_y)
else:
   oban_tot_low_var = newse_x_ravel * 0.
   oban_tot_low_var = oban_dz_var.reshape(newse_x.shape[0], newse_x.shape[1])

if (len(mrms_tot_mid_x_y) > 0):
   mrms_tot_mid_tree = scipy.spatial.cKDTree(mrms_tot_low_x_y)
else:
   oban_tot_mid_var = newse_x_ravel * 0.
   oban_tot_mid_var = oban_dz_var.reshape(newse_x.shape[0], newse_x.shape[1])

if (len(mrms_mesh_x_y) > 0):
   mrms_mesh_tree = scipy.spatial.cKDTree(mrms_mesh_x_y)
else:
   oban_mesh_var = newse_x_ravel * 0.
   oban_mesh_var = oban_mesh_var.reshape(newse_x.shape[0], newse_x.shape[1])

################### Find Az. Shear/DZ points to be interpolated: ###############


if ((len(mrms_low_x_y) > 0) and (len(mrms_dz_x_y) > 0)):
   oban_low_points = newse_tree.query_ball_tree(mrms_low_tree, roi)
else: 
   oban_low_points = []

if ((len(mrms_mid_x_y) > 0) and (len(mrms_dz_x_y) > 0)):
   oban_mid_points = newse_tree.query_ball_tree(mrms_mid_tree, roi)
else: 
   oban_mid_points = []

if ((len(mrms_tot_mid_x_y) > 0) and (len(mrms_dz_x_y) > 0)):
   oban_tot_mid_points = newse_tree.query_ball_tree(mrms_tot_mid_tree, roi)
else:
   oban_tot_mid_points = []

if ((len(mrms_tot_low_x_y) > 0) and (len(mrms_dz_x_y) > 0)):
   oban_tot_low_points = newse_tree.query_ball_tree(mrms_tot_low_tree, roi)
else:
   oban_tot_low_points = []

if ((len(mrms_mesh_x_y) > 0) and (len(mrms_dz_x_y) > 0)):
   oban_mesh_points = newse_tree.query_ball_tree(mrms_mesh_tree, roi)
else:
   oban_mesh_points = []

if (len(mrms_dz_x_y) > 0):
   oban_dz_points = newse_tree.query_ball_tree(mrms_dz_tree, roi)
   qc_dz_points = newse_tree.query_ball_tree(mrms_dz_tree, dz_rad)
else: 
   oban_dz_points = []
   qc_dz_points = []

oban_low_var = newse_x_ravel * 0.
oban_mid_var = newse_x_ravel * 0.
oban_tot_mid_var = newse_x_ravel * 0.
oban_tot_low_var = newse_x_ravel * 0.
oban_mesh_var = newse_x_ravel * 0.
oban_dz_var = newse_x_ravel * 0.

################### Interpolate Az. Shear observations applying a proximity threshold to DZ: ###############

############## 0-2 km AWS: ###################
for i in range(0, len(oban_low_points)):
   dz_max_count = 0.
   if (len(qc_dz_points[i]) > min_dz_obs):
      temp = mrms_dz_pixel_value[qc_dz_points[i]]
      temp = temp[temp > dz_thresh]
      dz_max_count = len(temp)

   if (dz_max_count > min_dz_obs):
      if ((len(oban_low_points[i]) > min_obs)):# and (len(near_rad_points[i]) == 0.) and (len(far_rad_points[i]) > 0.)):
         dis = np.sqrt((mrms_low_pixel_x_val[oban_low_points[i]] - newse_x_ravel[i])**2 + (mrms_low_pixel_y_val[oban_low_points[i]] - newse_y_ravel[i])**2)         
         weight = (roi**2 - dis**2) / (roi**2 + dis**2)
         oban_low_var[i] = np.sum(weight * mrms_low_pixel_value[oban_low_points[i]]) / np.sum(weight)

############## 2-5 km AWS: ###################
for i in range(0, len(oban_mid_points)):
   dz_max_count = 0.
   if (len(qc_dz_points[i]) > min_dz_obs):
      temp = mrms_dz_pixel_value[qc_dz_points[i]]
      temp = temp[temp > dz_thresh]
      dz_max_count = len(temp)

   if (dz_max_count > min_dz_obs):
      if ((len(oban_mid_points[i]) > min_obs)):# and (len(near_rad_points[i]) == 0.) and (len(far_rad_points[i]) > 0.)):
         dis = np.sqrt((mrms_mid_pixel_x_val[oban_mid_points[i]] - newse_x_ravel[i])**2 + (mrms_mid_pixel_y_val[oban_mid_points[i]] - newse_y_ravel[i])**2)         
         weight = (roi**2 - dis**2) / (roi**2 + dis**2)
         oban_mid_var[i] = np.sum(weight * mrms_mid_pixel_value[oban_mid_points[i]]) / np.sum(weight)

############## 0-2 km TOT Shear: ###################
for i in range(0, len(oban_tot_low_points)):
   dz_max_count = 0.
   if (len(qc_dz_points[i]) > min_dz_obs):
      temp = mrms_dz_pixel_value[qc_dz_points[i]]
      temp = temp[temp > dz_thresh]
      dz_max_count = len(temp)

   if (dz_max_count > min_dz_obs):
      if ((len(oban_tot_low_points[i]) > min_obs)):# and (len(near_rad_points[i]) == 0.) and (len(far_rad_points[i]) > 0.)):
         dis = np.sqrt((mrms_tot_low_pixel_x_val[oban_tot_low_points[i]] - newse_x_ravel[i])**2 + (mrms_tot_low_pixel_y_val[oban_tot_low_points[i]] - newse_y_ravel[i])**2)
         weight = (roi**2 - dis**2) / (roi**2 + dis**2)
         oban_tot_low_var[i] = np.sum(weight * mrms_tot_low_pixel_value[oban_tot_low_points[i]]) / np.sum(weight)

############## 2-5 km TOT Shear: ###################
for i in range(0, len(oban_tot_mid_points)):
   dz_max_count = 0.
   if (len(qc_dz_points[i]) > min_dz_obs):
      temp = mrms_dz_pixel_value[qc_dz_points[i]]
      temp = temp[temp > dz_thresh]
      dz_max_count = len(temp)

   if (dz_max_count > min_dz_obs):
      if ((len(oban_tot_mid_points[i]) > min_obs)):# and (len(near_rad_points[i]) == 0.) and (len(far_rad_points[i]) > 0.)):
         dis = np.sqrt((mrms_tot_mid_pixel_x_val[oban_tot_mid_points[i]] - newse_x_ravel[i])**2 + (mrms_tot_mid_pixel_y_val[oban_tot_mid_points[i]] - newse_y_ravel[i])**2)
         weight = (roi**2 - dis**2) / (roi**2 + dis**2)
         oban_tot_mid_var[i] = np.sum(weight * mrms_tot_mid_pixel_value[oban_tot_mid_points[i]]) / np.sum(weight)

############## MESH: ###################
for i in range(0, len(oban_mesh_points)):
   dz_max_count = 0.
   if (len(qc_dz_points[i]) > min_dz_obs):
      temp = mrms_dz_pixel_value[qc_dz_points[i]]
      temp = temp[temp > dz_thresh]
      dz_max_count = len(temp)

   if (dz_max_count > min_dz_obs):
      if ((len(oban_mesh_points[i]) > min_obs)):# and (len(near_rad_points[i]) == 0.) and (len(far_rad_points[i]) > 0.)):
         dis = np.sqrt((mrms_mesh_pixel_x_val[oban_mesh_points[i]] - newse_x_ravel[i])**2 + (mrms_mesh_pixel_y_val[oban_mesh_points[i]] - newse_y_ravel[i])**2)
         weight = (roi**2 - dis**2) / (roi**2 + dis**2)
         oban_mesh_var[i] = np.sum(weight * mrms_mesh_pixel_value[oban_mesh_points[i]]) / np.sum(weight)

############## MRMS: ###################
for i in range(0, len(oban_dz_points)):
   if ((len(oban_dz_points[i]) > min_obs)):# and (len(near_rad_points[i]) == 0.) and (len(far_rad_points[i]) > 0.)):
      dis = np.sqrt((mrms_dz_pixel_x_val[oban_dz_points[i]] - newse_x_ravel[i])**2 + (mrms_dz_pixel_y_val[oban_dz_points[i]] - newse_y_ravel[i])**2)         
      weight = (dz_roi**2 - dis**2) / (dz_roi**2 + dis**2)
      oban_dz_var[i] = np.sum(weight * mrms_dz_pixel_value[oban_dz_points[i]]) / np.sum(weight)

oban_low_var = oban_low_var.reshape(newse_x.shape[0], newse_x.shape[1])
oban_mid_var = oban_mid_var.reshape(newse_x.shape[0], newse_x.shape[1])
oban_tot_low_var = oban_tot_low_var.reshape(newse_x.shape[0], newse_x.shape[1])
oban_tot_mid_var = oban_tot_mid_var.reshape(newse_x.shape[0], newse_x.shape[1])
oban_mesh_var = oban_mesh_var.reshape(newse_x.shape[0], newse_x.shape[1])
oban_dz_var = oban_dz_var.reshape(newse_x.shape[0], newse_x.shape[1])


  
################### Write Cressman interpolated data to .nc file: ###############

fout.variables['LOW_CRESSMAN'][:,:] = oban_low_var
fout.variables['MID_CRESSMAN'][:,:] = oban_mid_var
fout.variables['TOT_LOW_CRESSMAN'][:,:] = oban_tot_low_var
fout.variables['TOT_MID_CRESSMAN'][:,:] = oban_tot_mid_var
fout.variables['MESH_CRESSMAN'][:,:] = oban_mesh_var
fout.variables['DZ_CRESSMAN'][:,:] = oban_dz_var

fout.close()
del fout

