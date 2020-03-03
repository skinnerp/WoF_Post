#!/bin/csh
#==================================================================

source ~/.tcshrc
set echo

# cd to directory where job was submitted from
#cd /scratch2/patrick.skinner/python_2018_verif
#cd /work/brian.matilla/python_2019_verif
cd /home/brian.matilla/work_backup/python_2019_verif

# User defined variables: 
#
#RETRO
#set dates = (20180619)
#set dates2 = (20180620)

#REALTIME 
set dates = (20190502)
set dates2 = (20190503)
set i = 1

foreach date ($dates)

################## 2019 REALTIME DIRECTORIES ############################
#   set exp_dir = '/scratch/skinnerp/2019_wofs_post/summary_files/'$date'/1900/'

   set out_dir = '/scratch/brian.matilla/MRMS_QPE2WOFS/original_files/'$date'/'

###### RETRO ######

#    set exp_dir = '/scratch/brian.matilla/NEWSe_GSI_retro/summary_files_vis/'$date'/1900'
#    set out_dir = '/scratch/brian.matilla/NEWSe_GSI_retro/mrms_qpe/'$date'/'

   set date2 = $dates2[$i]
# Run wrapper script for forecast: 

   sleep 5

   python wrapper_mrms_remap_v1.py -d $exp_dir -o $out_dir -a $date -b $date2
#   python wrapper_mrms_remap2.py -d $exp_dir -o $out_dir -a $date -b $date2
   sleep 2

   @ i++
end 

#sleep 5 

