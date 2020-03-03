#!/bin/csh -f
set echo

#set base_dir = '/work/brian.matilla/python_2019_verif/'
set base_dir = '/home/brian.matilla/work_backup/python_2019_verif/'

set date = `date -u --date='0 minutes ago' +%Y%m%d` 
set prev_date = `date -u --date='1 day ago' +%Y%m%d`

#set prev_date = $1
#set date = $2

#setenv yyyy `date -u +%Y`
#setenv mm `date -u +%m`
#setenv dd `date -u +%d `
#
#setenv ddy `date -d "$date 1 day ago" +%d`
#setenv ddt `date -d "$date 1 day" +%d`
#
##set date = $yyyy$mm$dd
#set prev_date = $yyyy$mm$ddy
#set next_date = $yyyy$mm$ddt

cd $base_dir

#exec sed -i -e "18s/$prev_date/$date/g" run_mrms_remap.csh -i -e "19s/$date/$next_date/g" run_mrms_remap.csh -i -e "s/$prev_date/$date/g" run_mrms_qpe_plot.csh 

exec sed -i -e "s/$prev_date/$date/g" run_mrms_remap_v2.csh -i -e "s/$prev_date/$date/g" run_mrms_qpe_plot_agg.csh -i -e "s/$prev_date/$date/g" run_mrms_agg_realtime.csh -i -e "s/$prev_date/$date/g" run_paint.csh -i -e "s/$prev_date/$date/g" run_verif_obj.csh "s/$prev_date/$date/g" run_mrms_createnans.csh
