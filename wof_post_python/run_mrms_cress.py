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
parser.add_argument('dates',type=str,nargs='+',help='List of dates in format YYYMMDD')
###########################################################################
args = parser.parse_args()
dates = args.dates
#if you do not want to use the input argument (just list the dates as separate 
#arguments and they will be read in as a list), then just comment out and create a list
#as below
#dates = ['20190529','20190530']
os.chdir('/scratch2/patrick.skinner/python_2019_verif')
for date in dates:
    print(date)
    mrms_dir = '/work1/anthony.reinhart/realtimeVerif/{date}'.format(date=date)
    exp_dir = '/scratch/skinnerp/2019_wofs_post/summary_files/{date}/1900/'.format(date=date)
    out_dir = '/oldscratch/skinnerp/2019_wofs_post/mrms/{date}'.format(date=date)
    mrms_dir = mrms_dir + '/'
    exp_dir = exp_dir + '/' 
    out_dir = out_dir + '/' 

    # Run wrapper script for forecast: 
    time.sleep(5)

    cmd = '/home/louis.wicker/anaconda2/envs/wof-test/bin/python wrapper_mrms_cress.py -m {mrms_dir} -d {exp_dir} -o {out_dir}'.format(mrms_dir=mrms_dir,exp_dir=exp_dir,out_dir=out_dir)
    print 'Executing command: ' + cmd
    subprocess.call(cmd,shell=True)
    print cmd + ' is finished ...'
    time.sleep(2)

time.sleep(5)

