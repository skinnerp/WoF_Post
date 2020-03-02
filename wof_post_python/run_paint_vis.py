#!/usr/bin/env python

import numpy as np
import datetime
import os
import subprocess
import argparse
import time
import glob

###########################################################################
parser = argparse.ArgumentParser()
parser.add_argument('date',type=str,help='Date format YYYMMDD')

###################################################################################################
# run_script is a function that runs a system command

def run_script(cmd):
    print "Executing command:  " + cmd
    os.system(cmd)
    print cmd + "  is finished...."
    return

###########################################################################
args = parser.parse_args()
date = args.date
date = date + '/'
#you can comment out the input argument stuff if you want and hard code date
#I have some tricks we can use to get date if you want to cron this up to use the system datetime info
#date = '20190529'

#cd to directory where job was submitted from
os.chdir('/scratch2/patrick.skinner/python_2019_verif')

base_dir = '/oldscratch/skinnerp/2018_newse_post/'

#hours = ['1600', '1700', '1800', '1900', '2000', '2100', '2200', '2300', '0000', '0100', '0200', '0300', '0400', '0500', '0600']
#nt = [72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72]

hours = ['0030']
nt = [36]

#hours = ['1900', '1930', '2000', '2030', '2100', '2130', '2200', '2230', '2300', '2330', '0000', '0030', '0100', '0130', '0200', '0230', '0300']
#nt = [72, 36, 72, 36, 72, 36, 72, 36, 72, 36, 72, 36, 72, 36, 72, 36, 72]

# Build directories needed: 

summary = "/oldscratch/skinnerp/2018_newse_post/summary_files/{date}".format(date=date)
#summary = "/oldscratch/skinnerp/2018_newse_post/summary_files_vis/{date}".format(date=date)
image = "/oldscratch/skinnerp/2018_newse_post/summary_files_vis/images/"
flagdir = "/oldscratch/skinnerp/2018_newse_post/flags_vis/{date}".format(date=date)

# Run wrapper script for forecast: 
time.sleep(5)

for itx,hr in enumerate(hours):
    sum_folder = '0030/'
    image_folder = '20180502_0030/'
    flag_name = '0030_paint_v6.txt'

    tempsumdir = os.path.join(summary,sum_folder)
    imagedir = os.path.join(image,image_folder)
#    tempsumdir = tempsumdir + '/'
#    imagedir = imagedir + '/'

    print tempsumdir, imagedir

    flagfile = os.path.join(flagdir,flag_name)

    tempnt = nt[itx]
    tempfiles = os.listdir(tempsumdir)
#    tempfiles = glob.glob(tempsumdir)

    print hr, len(tempfiles), flagfile, os.path.isfile(flagfile)

    if len(tempfiles) > 2:
        if os.path.isfile(flagfile) == False:
            #create empty file
            os.mknod(flagfile)
            cmd = '/home/louis.wicker/anaconda2/envs/wof-test/bin/python wrapper_paintball_vis.py -d {tempsumdir} -o {imagedir} -e {tempnt}'.format(tempsumdir=tempsumdir,imagedir=imagedir,tempnt=tempnt)
            print 'Executing command: ' + cmd
            subprocess.call(cmd,shell=True)
            print cmd + ' is finished ...'
            break
