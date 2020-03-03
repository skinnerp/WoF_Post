#!/usr/bin/env python
#
# "mrms_remap.py"
#
# Script to perform interpolation of 1km MRMS QPE field
#	onto 3km WoFS grid using ESMPy. This script follows a similar
#	convention to the MRMS Cressman approach for obtaining the
#	lats and lons for both grids and then uses ESMPy to perform
#	the interpolation technique.
#
# Script written by Brian Matilla (CIMMS/NSSL)
# Script date: 2019 APR 5
#
# Load modules
#
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import scipy
from scipy import signal
from scipy import *
from scipy import ndimage
#import skimage
#from skimage.morphology import label
#from skimage.measure import regionprops
import math
from math import radians, tan, sin, cos, pi, atan, sqrt, pow, asin, acos
import pylab as P
import numpy as np
from numpy import NAN
import sys
import netCDF4
from optparse import OptionParser
#from netcdftime import utime
import os
import time as timeit
from optparse import OptionParser
import news_e_post_cbook
from news_e_post_cbook import *
import news_e_plotting_cbook_v2
from news_e_plotting_cbook_v2 import *
import ESMF 		### import name used for ESMPy 

# Add parsing options 
parser = OptionParser()
parser.add_option("-z", dest="qpe_file", type="string", default= None, help="Input Path (of MRMS QPE file)")
parser.add_option("-o", dest="out_file", type="string", help = "Output File Path")
parser.add_option("-f", dest="newse_path", type="string", help = "Path to NEWS-e Summary File")

(options, args) = parser.parse_args()

if ((options.qpe_file == None) or (options.out_file == None) or (options.newse_path == None)):
    print
    parser.print_help()
    print
    sys.exit(1)
else:
    qpe_file = options.qpe_file
    out_file = options.out_file
    newse_path = options.newse_path

#################################### Generic Basemap Variables (to call as quickly as possible):  #####################################################
damage_files     = ''
area_thresh      = 1000.
resolution       = 'c'
######################################################################################################
#################################### Read Data:  #####################################################
######################################################################################################

################### Get Grid info from NEWS-e summary file: ############################

try:
   newse_in = netCDF4.Dataset(newse_path, "r")
   print("Opening %s \n" % newse_path)
except:
   print("%s does not exist! \n" %newse_path)
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

##### Create the WoFS/NEWS-e basemap instance #####

newse_map = Basemap(llcrnrlon=sw_lon_newse, llcrnrlat=sw_lat_newse, urcrnrlon=ne_lon_newse, urcrnrlat=ne_lat_newse, projection='lcc', lat_1=true_lat1, lat_2=true_lat2, lat_0=cen_lat, lon_0=cen_lon, resolution = resolution, area_thresh = area_thresh)

################### Initialize interpolated variables and output .nc file: ############################
try:
   fout = netCDF4.Dataset(out_file, "w")
except:
   print("Could not create %s!\n" % out_file)

fout.createDimension('NX', newse_lat.shape[1])
fout.createDimension('NY', newse_lat.shape[0])

fout.createVariable('XLAT', 'f4', ('NY','NX',))
fout.createVariable('XLON', 'f4', ('NY','NX',))
fout.createVariable('MRMS_QPE', 'f4', ('NY', 'NX',))

fout.variables['XLAT'][:] = newse_lat
fout.variables['XLON'][:] = newse_lon

########## Read MRMS Quantitative Precipitation Estimate (QPE) data, remove data outside NEWS-e domain for speed, and flip lat coordinates to be compatible: #######

try:
   qpe_fin = netCDF4.Dataset(qpe_file, "r")
   print("Opening %s \n" % qpe_file)
except:
   print("%s does not exist! \n" %qpe_file)
   sys.exit(1)

################### For first QPE file, build QPE MRMS grid from sparse .netcdf file: ############################

qpe_nw_lat = qpe_fin.Latitude
qpe_nw_lon = qpe_fin.Longitude
qpe_lat_dy = qpe_fin.LatGridSpacing
qpe_lon_dx = qpe_fin.LonGridSpacing
qpe_lat_length = len(qpe_fin.dimensions["Lat"])
qpe_lon_length = len(qpe_fin.dimensions["Lon"])

qpe_lat_range = np.arange((qpe_nw_lat-((qpe_lat_length-1)*qpe_lat_dy)),(qpe_nw_lat+0.00001),qpe_lat_dy)
qpe_lon_range = np.arange(qpe_nw_lon,(qpe_nw_lon+((qpe_lon_length-0.99999)*qpe_lon_dx)),qpe_lon_dx)
qpe_xlon_full, qpe_xlat_full = np.meshgrid(qpe_lon_range, qpe_lat_range)

qpe_lon_range_indices = np.arange(0,qpe_lon_length)
qpe_lat_range_indices = np.arange(0,qpe_lat_length)

qpe_xlon_indices, qpe_xlat_indices = np.meshgrid(qpe_lon_range_indices, qpe_lat_range_indices)

qpe_min_lat = (np.abs(qpe_xlat_full[:,0]-np.min(newse_lat))).argmin()
qpe_max_lat = (np.abs(qpe_xlat_full[:,0]-np.max(newse_lat))).argmin()
qpe_min_lon = (np.abs(qpe_xlon_full[0,:]-np.min(newse_lon))).argmin()
qpe_max_lon = (np.abs(qpe_xlon_full[0,:]-np.max(newse_lon))).argmin()

qpe_xlat = qpe_xlat_full[qpe_min_lat:qpe_max_lat,qpe_min_lon:qpe_max_lon]
qpe_xlon = qpe_xlon_full[qpe_min_lat:qpe_max_lat,qpe_min_lon:qpe_max_lon]

qpe_sw_xlat = qpe_xlat[0,0]
qpe_sw_xlon = qpe_xlon[0,0]
qpe_ne_xlat = qpe_xlat[-1,-1]
qpe_ne_xlon = qpe_xlon[-1,-1]

#print(qpe_xlat)
#print(qpe_xlon)

##### We need the data from the source, so import it here #####

mrms_qpe_var = 'RadarOnly_QPE_01H'
qpe_vals = qpe_fin.variables[mrms_qpe_var][:]
qpe_vals = np.flipud(qpe_vals)
qpe_vals = qpe_vals[qpe_min_lat:qpe_max_lat,qpe_min_lon:qpe_max_lon]

qpe_vals[qpe_vals < 0] = 0 ### Mask out all of the negative (-99.9) values.
qpe_vals = qpe_vals / 25.4 ### mm to in.

########## Create the QPE basemap instance #############
qpe_map = Basemap(llcrnrlon=qpe_sw_xlon, llcrnrlat=qpe_sw_xlat, urcrnrlon=qpe_ne_xlon, urcrnrlat=qpe_ne_xlat, projection='lcc', lat_1=true_lat1, lat_2=true_lat2, lat_0=cen_lat, lon_0=cen_lon, resolution = resolution, area_thresh = area_thresh)

########## Begin the ESMPy interpolation routine ############

ESMF.Manager(debug=True); #output a log file for troubleshooting if interpolation fails.

##### Obtain the grid shapes as variables

qpe_shape = qpe_xlat.shape
#print(qpe_shape)
newse_shape = newse_lat.shape
#print(newse_shape)

##### Set up the source and destination grid objects ####

sourcegrid = ESMF.Grid(np.array(qpe_shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys= ESMF.CoordSys.SPH_DEG)
destgrid = ESMF.Grid(np.array(newse_shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys= ESMF.CoordSys.SPH_DEG) 

if ESMF.local_pet() == 0:
	print("Source and destination grid objects created.")

##### Generate blank coordinate variables which will act as pointers for the grid object coordinates#####

sourcelon = sourcegrid.get_coords(0)
sourcelat = sourcegrid.get_coords(1)

destlon = destgrid.get_coords(0)
destlat = destgrid.get_coords(1)

##### Now that pointers are ready, fill the lat/lon values into the pointers using [] accessor (to not overwrite pointers!)#####

sourcelon[...] = qpe_xlon
sourcelat[...] = qpe_xlat

destlon[...] = newse_lon
destlat[...] = newse_lat

##### Generate the Field objects

sourcefield = ESMF.Field(sourcegrid, name= 'MRMS QPE Native Grid')
destfield = ESMF.Field(destgrid, name= 'Warn-on-Forecast Grid')

##### Now create the source field data pointer and fill with QPE data #####

sourcefield.data[...] = qpe_vals

##### Create the regridding instance... Needs the source field, destination field, regrid method, and directive on what to do with unmapped points #####

print("Running regridding procedure for %s \n. Please wait..." % qpe_file )

regrid = ESMF.Regrid(sourcefield, destfield, regrid_method=ESMF.RegridMethod.BILINEAR, unmapped_action = ESMF.UnmappedAction.IGNORE)

if ESMF.local_pet() == 0:
	print("Regrid completed successfully.")
elif ESMF.local_pet() != 0:
	print("Regrid unsuccessful.")

##### Compute the regridding weights and final values by running regrid routine on the source and destination fields #####

destfield = regrid(sourcefield, destfield)

qpe_regrid= destfield.data #Output values

##### Now output the regridded QPE data to the created NC file.

fout.variables['MRMS_QPE'][:,:] = qpe_regrid

fout.close()
del fout
