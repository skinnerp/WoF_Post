#!/bin/csh
#==================================================================

source ~/.tcshrc
set echo

# cd to directory where job was submitted from
#cd /work/brian.matilla/python_2019_verif
cd /home/brian.matilla/work_backup/python_2019_verif

# REALTIME RUNS

set date = 20190718
set base_dir = /scratch/brian.matilla/

if ($date >= 20190617) then
########## WPC Cycling Start ################
   set times = (1600 1700 1800 1900 2000 2100 2200 2300 0000 0100 0200 0300 0400 0500 0600)
   set nt = (24 24 24 24 24 24 24 24 24 24 24 24 24 24 24)

else if ($date >= 20190429 && $date <= 20190616) then
########## SFE Cycling Start ################
   set times = (1700 1900 1930 2000 2030 2100 2130 2200 2230 2300 2330 0000 0030 0100 0130 0200 0230 0300)
   set nt = (12 24 12 24 12 24 12 24 12 24 12 24 12 24 12 24 12 24)

else
########## All other runs ###################
#   set times = (1900 2000 2100 2200 2300 0000 0100 0200 0300)
#   #   set nt = (72 72 72 72 72 72 72 72 72)
#   endif

endif
################## 2019 REALTIME DIRECTORIES ############################
set out_base = /scratch/brian.matilla/MRMS_QPE_AGG/
set mrms_base = /scratch/brian.matilla/MRMS_QPE2WOFS/

#set summary = "/scratch/skinnerp/2019_wofs_post/summary_files/"$date"/"
set summary = "/scratch/brian.matilla/WoFS_2019/summary_files/"$date"/"
set flag_base = $base_dir"flags/"
set blank_dir = $mrms_base"test/"

#   sleep 5

set i = 1

set outdir = $out_base$date"/"
set mrmsdir = $mrms_base$date"/"
set flagdir = $flag_base$date"/"

if ( ! -d ${summary} ) then
   echo "Summary forecast directories not created yet"
   exit(0)
endif

### Run Wrapper Script ###

sleep 1

foreach dir ($times)
   set aggdir = $outdir$dir"/"
   set tempsumdir = $summary$dir"/"   
   set tempnt = $nt[$i]

   echo $aggdir
   echo $tempsumdir
   echo $tempnt

   set flagfile = $flagdir$dir"_mrms_agg.txt"
   echo $flagfile

   set tempfiles = `ls -a ${tempsumdir} | wc | awk '{print $1}'`
   if ( "${tempfiles}" > 2) then #$ne ) then
      if ( ! -e $flagfile ) then
         touch $flagfile
         python wrapper_mrms_agg_realtime.py -m $mrmsdir -o $aggdir -d $date -e $tempnt -t $dir -b $blank_dir
#         exit(0)
      endif
   endif
   wait
   @ i++
end

sleep 5
