#!/bin/csh -f
#
set echo

#cd /work/brian.matilla/MRMS_QPE2WOFS/
cd /scratch/brian.matilla/MRMS_QPE2WOFS/

setenv yyyy `date -u +%Y`
setenv mm `date -u +%m`
setenv dd `date -u +%d`

#set currdate = $yyyy$mm$dd
set currdate = $1
set wpc_start_date = 20190617

mkdir $currdate

cd /scratch/brian.matilla/flags/

mkdir $currdate

cd /scratch/brian.matilla/MRMS_QPE_AGG/

mkdir $currdate

#cp -R template_wpc/* $currdate
if ($currdate < $wpc_start_date) then
   cp -R template_hwt/* $currdate
else
   cp -R template_wpc/* $currdate
endif

cd /scratch/brian.matilla/qpe_qc_2019/

mkdir $currdate

if ($currdate < $wpc_start_date) then
   cp -R template_hwt/* $currdate
else
   cp -R template_wpc/* $currdate
endif

cd /scratch/brian.matilla/GLM/OBS_SEQ/

mkdir $currdate
