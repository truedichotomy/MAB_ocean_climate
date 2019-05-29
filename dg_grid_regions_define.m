function region = dg_grid_regions_define(LON,LAT,Z,depthlim,sodaflag)
% synposis: function region = dg_grid_regions_define(LON,LAT,Z,depthlim)
%
% this function is called by dg_grid_regions.m or dg_grid_regions_analysis.m
% DG 2017-11-01

if ~exist('sodaflag')
	sodaflag = 0
end %if

[status,machinename] = system('uname -n');

if strcmp(lower(machinename(1:12)),'scarletnight') == 1
	rootdirE = ['/Volumes/Passport16/'];
	rootdir = ['/Users/gong/Documents/'];
	bathydir = [rootdir 'Research/bathy/'];
	datadir = [rootdir 'Research/Data/NOAA/NEFSC/hydro/'];
	homedir = [rootdir 'Research/Projects/Wallace_REU/Donglai/'];
	workdir = [rootdir 'Research/MAB_climate/'];
    if sodaflag == 0
	    figoutdir = [workdir 'pngs/'];
    elseif sodaflag == 1
        figoutdir = [workdir 'pngs/SODA/'];
    end %if
elseif strcmp(lower(machinename(1:5)),'orion') == 1
	rootdir = ['/Users/c2po/Research/MAB_climate/'];
	bathydir = [rootdir 'data/'];
	datadir = [rootdir 'data/NEFSC/hydro/'];
	homedir = [rootdir];
	workdir = [rootdir 'data/'];
    if sodaflag == 0
	    figoutdir = [workdir 'pngs/'];
    elseif sodaflag == 1
        figoutdir = [workdir 'pngs/SODA/'];
    end %if
else
	rootdir = ['/Users/gong/Documents/'];
	bathydir = [rootdir 'Research/bathy/'];
	datadir = [rootdir 'Research/Data/NOAA/NEFSC/hydro/'];
	homedir = [rootdir 'Research/Projects/Wallace_REU/Donglai/'];
	workdir = [rootdir 'Research/MAB_climate/'];
    if sodaflag == 0
	    figoutdir = [workdir 'pngs/'];
    elseif sodaflag == 1
        figoutdir = [workdir 'pngs/SODA/'];
    end %if
end %if

if nargin == 0
    load([bathydir 'gebco_MAB_30arcsec.mat']);
end %if

if nargin < 4
    depthlim = [6 1000]
end %if

region = [];

region(1).label = 'SNE';
%region(1).lon = [-70, -70, -70.544, -71.5034, -72.2848, -71.5034];
%region(1).lat = [40, 40.7537, 41.1418, 41.1418, 40.7537, 39.8582];
region(1).lon = [-70, -70, -71.02379, -71.803099, -71.88708, -72.4060, -71.2699, -70];
region(1).lat = [39.55, 41.2697, 41.55261, 41.44175, 41.074622, 40.8968, 39.5765, 39.55];

region(2).label = 'NYB1';
%region(2).lon = [-71.5034, -72.2848, -73.0485, -74.0253, -72.6933];
%region(2).lat = [39.8582, 40.7537, 40.574, 39.7836, 38.8881];
region(2).lon = [-71.2699, -72.4060, -73.9595, -74.1005, -72.4335, -71.2699];
region(2).lat = [39.5765, 40.8968, 40.5902, 39.8399, 38.7224, 39.5765];

region(3).label = 'NYB2';
%region(3).lon = [-72.6933, -74.0253, -75.0020, -73.7589];
%region(3).lat = [38.8881, 39.7836, 38.7090, 38.0075];
region(3).lon = [-72.4335, -74.1005, -75.068616, -73.37602, -72.4335];
region(3).lat = [38.7224, 39.8399, 38.75436, 37.7167, 38.7224];

region(4).label = 'SS1';
%region(4).lon = [-73.7589, -75.0020, -75.0020, -75.8367, -74.4693, -74.3627];
%region(4).lat = [38.0075, 38.7090, 38.3358, 37.2612, 37.0075, 37.3657];
region(4).lon = [-73.3760, -75.068616, -75.8367, -73.98192, -73.3760];
region(4).lat = [37.7167, 38.75436, 37.2612, 36.87062, 37.7167];

region(5).label = 'SS2';
%region(5).lon = [-74.493, -75.8367, -75.9966, -75.2862, -75.1974, -74.5403];
%region(5).lat = [37.0075, 37.2612, 36.9925, 35.3358, 35.3358, 35.9925];
region(5).lon = [-73.98192, -75.8367, -76.02073, -75.36677, -74.3666, -73.98192];
region(5).lat = [36.87062, 37.2612, 36.90654, 35.230329, 35.61346, 36.87062];

region(6).label = 'MAB';
region(6).lon = [-70, -70, -71.02379, -71.803099, -71.88708, -72.4060, -73.9595, -74.1005, -75.068616, -75.8367, -76.02073, -75.36677, -74.3666, -73.98192, -73.3760, -72.4335, -71.2699, -70];
region(6).lat = [39.55, 41.2697, 41.55261, 41.44175, 41.074622, 40.8968, 40.5902, 39.8399, 38.75436, 37.2612, 36.90654, 35.230329, 35.61346, 36.87062, 37.7167, 38.7224, 39.5765, 39.55];

region(7).label = 'GB';
%region(6).lon = [-69, -69, -66.5, -66.5];
%region(6).lat = [40.6, 42.2, 42.2, 40.6];
region(7).lon = [-68.9971, -68.9971, -68.39713, -67.53589,  -66.17145, -65.416655, -66.500464, -67.56492, -68.51325, -68.9971];
region(7).lat = [40.61436, 41.34471, 41.91941, 42.20676, 42.24268,  41.71587, 40.530550, 40.063606, 39.84809, 40.61436];

region(8).label = 'ENE';
region(8).lon = [-68.51325, -68.9971, -68.9971, -70, -70, -68.51325];
region(8).lat = [39.84809, 40.61436, 41.34471, 41.2697, 39.55, 39.84809];

region(9).label = 'GOM';
region(9).lon = [-68.9971, -70, -71.1572, -70.1955, -67.9136, -67.0905, -66.51, -66.17145, -67.53589, -68.39713, -68.9971];
region(9).lat = [41.34471, 41.7, 42.4879, 43.8169, 44.6210, 44.7885, 43.6, 42.24268, 42.20676, 41.91941, 41.34471];

% find all grid locations within a region
for ii = 1:length(region)
    inregion = inpolygon(LON,LAT,region(ii).lon,region(ii).lat);
    regionalind = find(inregion == logical(1));
    region(ii).longrid = LON(regionalind);
    region(ii).latgrid = LAT(regionalind);
    region(ii).zgrid = Z(regionalind);
end %for

for ri = 1:length(region)
    region(ri).label
    region(ri).zind = find(region(ri).zgrid < -depthlim(1) & region(ri).zgrid >= -depthlim(2));
    %region(ri).zind = find(region(ri).zgrid < -6 & region(ri).zgrid >= -1000);
    zind = region(ri).zind;
    region(ri).llon = region(ri).longrid(zind);
    region(ri).llat = region(ri).latgrid(zind);
    region(ri).zz = region(ri).zgrid(zind);
end %for
