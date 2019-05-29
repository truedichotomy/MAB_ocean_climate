% dg_setup_MABclimate_dir
% 2018-04-06: gong@vims.edu, initial version
% 2019-05-29: gong@vims.edu, clean up for public release
%

if ~exist('sodaflag')
	sodaflag = 0
end %if

[status,machinename] = system('uname -n');

rootdir = ['./MAB_ocean_climate/']
bathydir = [rootdir './']
datadir = [rootdir 'NOAA_NEFSC_data/hydro/']
homedir = rootdir
workdir2 = [rootdir 'forsyth/']
workdir3 = [rootdir 'sensitivity/']
workdir = rootdir
workdirlocal = ['.']
SODAdir = '/dir/to/SODA/';
gissdatadir = [rootdir 'GISStemp/'];

if sodaflag == 0
	figoutdir = [workdir 'pngs/']
elseif sodaflag == 1
	figoutdir = [workdir 'pngs/SODA/']
end %if

