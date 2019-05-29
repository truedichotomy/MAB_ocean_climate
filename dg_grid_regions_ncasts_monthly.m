% DG 2017-12-26, 2018-04-06
% this script calculates the number of casts per month per region for all 40 years together
% run this script after dg_grid_regions_analysis2D_nobias.m

loadregflag = 0
loadgridflag = 0
loadcastflag = 0
loadcccflag = 0

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

if loadregflag == 1
    reg = dg_grid_regions_define_basic;
end %if

if loadgridflag == 1
	display('loading hydroMABgrid: LON, LAT, Z, casts, region')
	tic
	load([workdir 'hydroMABgrid.mat']);
	toc
end %if

if loadcastflag == 1
	display('loading hydroMABcasts: casts, sodacast')
	tic
	load([workdir 'hydroMABcasts.mat']);
	toc
end %if

if loadcccflag == 1
    load([workdir 'MABclimate2D.mat'])
end %if

yyyylist = [1977:2016];
yyyy = [casts.yr]';
dyd = [casts.dyd]';
yd = [casts.yd]';
mm = str2num(datestr(datenum(yyyy,1,dyd),5));
londata = [casts.lon]';
latdata = [casts.lat]';
depth = -double([casts.zdeep])';
tempda = [casts.tempda]';
saltda = [casts.saltda]';
tempbott = [casts.tempbott]';
saltbott = [casts.saltbott]';

ccc.ncastsmonthly = repmat(NaN,[length(ccc.month), length(region)]);
ccc.ncasts = repmat(NaN,[length(yyyylist), length(ccc.month), length(region)]);

% find all grid locations within a region
for ri = 1:length(reg)
    inregion = inpolygon(londata,latdata,reg(ri).lon,reg(ri).lat);
    regionalind = find(inregion == logical(1));
    reg(ri).clon = londata(regionalind);
    reg(ri).clat = latdata(regionalind);
    reg(ri).cdepth = depth(regionalind);
    reg(ri).cyyyy = yyyy(regionalind);
    reg(ri).cmm = mm(regionalind);
    reg(ri).ind = [];
    for mi = 1:12
        %reg(ri).ind{mi} = find(cmm == mi & cdepth <= 1500);
        ccc.ncastsmonthly(mi,ri) = length(find(reg(ri).cmm == mi & reg(ri).cdepth <= 1500));
    end %for

    for yi = 1:length(yyyylist)
        for mi = 1:12
            %reg(ri).ind{mi} = find(cmm == mi & cdepth <= 1500);
            ccc.ncasts(yi,mi,ri) = length(find(reg(ri).cmm == mi & reg(ri).cyyyy == yyyylist(yi) & reg(ri).cdepth <= 1500));
        end %for
    end %for
end %for
