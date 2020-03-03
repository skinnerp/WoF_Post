#!/bin/csh
#==================================================================

source ~/.tcshrc
set echo

#cd to directory where job was submitted from
#cd /scratch2/patrick.skinner/python_2019_verif
#cd /work/brian.matilla/python_2019_verif
cd /home/brian.matilla/work_backup/python_2019_verif

# User defined variables: 

set date = 20190518
set base_dir = /oldscratch/skinnerp/2019_wofs_post/
set my_dir = /scratch/brian.matilla/

########## WPC Cycling Start ################
#set times = (1600 1700 1800 1900 2000 2100 2200 2300 0000 0100 0200 0300 0400 0500 0600)
#set nt = (24 24 24 24 24 24 24 24 24 24 24 24 24 24 24)

########## SFE Cycling Start ################
set times = (1700 1900 1930 2000 2030 2100 2130 2200 2230 2300 2330 0000 0030 0100 0130 0200 0230 0300)
set nt = (36 72 36 72 36 72 36 72 36 72 36 72 36 72 36 72 36 72)

########## All other runs ###################
#set times = (1900)
#set nt = (3)

# Build directories needed: 

set image = "/www/wof.nssl.noaa.gov/newse_images/"$date"/"
#set image = "/www/wof.nssl.noaa.gov/newse_research/20190622_obj/"
set flagdir = $my_dir"flags/"$date"/"
#set flagdir = $base_dir"flags/"$date"/"
set qpeqc_dir = $my_dir"qpe_qc_2019/"$date"/"

# Run wrapper script for forecast: 

#sleep 5

set i = 1
### Run Wrapper Script ###

sleep 1

foreach dir ($times)
   set tempqpedir = $qpeqc_dir$dir"/"
   echo $tempqpedir

   set temp_flagdir = $flagdir$dir"/"
   echo $temp_flagdir

   set tempimagedir = $image$dir"/"
   echo $tempimagedir

   set tempnt = $nt[$i]
   echo $tempnt

   set flagfile = $flagdir$dir"_obmatch.txt"
   echo $flagfile

   set tempfiles = `ls -a ${tempqpedir} | wc | awk '{print $1}'`
   if ( "${tempfiles}" > 2) then #$ne ) then
      if ( ! -e $flagfile ) then 
         touch $flagfile       
         python wrapper_verif_obj_qpe.py -r $tempqpedir -f $temp_flagdir -i $tempimagedir -a $date -e $tempnt
#         exit(0)
      endif
   endif
   wait
   @ i++
end

sleep 10

