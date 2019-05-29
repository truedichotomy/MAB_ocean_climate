% this script uses NOAA NEFSC hydrographic data to conduct regional climate analysis of temperature and salinity for each region in the U.S. Northeast.
% this script should be run first. it calls dg_grid_regions_define.m and it should be followed by dg_grid_regions_hydroavg.m
%
% Note this script takes a LONG time (10-16 hours) to run on a quad core machine with 16 GB of RAM!!!
%
% Donglai Gong 2017-09-22, 2017-10-09

dbstop if error

loadbathyflag = 0
calcregionflag = 0
loadctdflag = 0
savegridflag = 0
loadgridflag = 1 %
savecastflag = 0
loadcastflag = 1 %
calcflag = 1 %
hydromoflag = 'all'
testplot = 0

% parameters for inverse-distance search algorithm
h = 0 %meters
searchradius = 40000 % in meters, default is 30000
p = 4 % power parameter, default is 2

% set depth limits for
depthlim = [6 1000]; % looking at shelf regions in between these two isobaths
%depthlim = []

% setup directories needed for the script
dg_setup_MABclimate_dir

% load GEBCO bathymetry
if loadbathyflag == 1
	load([bathydir 'gebco_MAB_30arcsec.mat']);
    display('Done loading bathy data.')
end %if

%% define regions
if calcregionflag == 1
    region = dg_grid_regions_define(LON,LAT,Z);
end %if calcregionflag

%% loading NEFSC CTD data into one data structure array
if loadctdflag == 1
    casts1yr = struct('cast',NaN,'lat',NaN,'lon',NaN,'yr',NaN,'yd',NaN,'dyd',NaN,'depth',NaN,'pc',NaN,'vn',NaN,'np',NaN,'s',[],'t',[],'p',[],'salt',[],'temp',[],'pres',[],'cru',{},'opsid',{},'gear','');
    casts = [];
    datafiles = dir([datadir '*.mat']);
    for ii = 1:length(datafiles)
        hydrodata1 = load([datadir datafiles(ii).name]);
        fields = fieldnames(hydrodata1);
        casts1yr = hydrodata1.(fields{1});
        casts = [casts, casts1yr]; % concatinate data structure containing the hydro cast data
    end %for

    % setting longitude to based on East
    for ii = 1:length(casts)
        casts(ii).lon = -casts(ii).lon;
        if strcmp(lower(casts(ii).gear),'bottle') == 1
            casts(ii).z = -casts(ii).p;
        else
            casts(ii).z = gsw_z_from_p(casts(ii).p, casts(ii).lat);
        end %if
    end %for

    % determine the deepest depth at any given lon/lat
    zdeep = interp2(LON,LAT,Z,[casts.lon],[casts.lat],'linear');
    for ii = 1:length(casts)
        casts(ii).zdeep = min([zdeep(ii) casts(ii).z]);
        casts(ii).zz = [-1:-1:casts(ii).zdeep];
    end %for

    %zdiff = [casts.zmax] + [casts.depth];
    %scatter([casts.lon], [casts.lat],10, log10(abs(zmax)),'o','filled'); colorbar;
    %scatter([casts.lon], [casts.lat],10, zdiff,'o','filled'); colorbar;

    %bottind = find(strcmp({casts.gear},'bottle') == 1);
    %instind = find(strcmp({casts.gear},'bottle') == 0);

    for ii = 1:length(casts)
        if length(find(~isnan(casts(ii).t))) >= 2
            casts(ii).temp = interp1(casts(ii).z,casts(ii).t,casts(ii).zz,'linear',NaN);
        else
            casts(ii).temp = repmat(NaN,size(casts(ii).zz));
        end %if

        if length(find(~isnan(casts(ii).s))) >= 2
            casts(ii).salt = interp1(casts(ii).z,casts(ii).s,casts(ii).zz,'linear',NaN);
        else
            casts(ii).salt = repmat(NaN,size(casts(ii).zz));
        end %if

        % fill the bottom boundary layer with value from the deepest measurement of that cast
        if abs(casts(ii).zdeep - min(casts(ii).z)) <= 15
            zbblind = find(casts(ii).zz < min(casts(ii).z)); % find indices of the BBL
            ztempind = find(~isnan(casts(ii).temp)); % find temperature indices for the cast
            if ~isempty(zbblind) & ~isempty(ztempind)
                casts(ii).temp(zbblind) = casts(ii).temp(ztempind(end));
            end %if
            zsaltind = find(~isnan(casts(ii).salt)); % find salinity indices for the cast
            if ~isempty(zbblind) & ~isempty(zsaltind)
                casts(ii).salt(zbblind) = casts(ii).salt(zsaltind(end));
            end %if
        end %if

        % set the value of the top mixed layers (tml) to that of the first measurement at the top of the cast
        if abs(max(casts(ii).zz) - max(casts(ii).z)) <= 5
            tmlind = find(casts(ii).zz > max(casts(ii).z));
            casts(ii).temp(tmlind) = casts(ii).t(1);
            casts(ii).salt(tmlind) = casts(ii).s(1);
        end %if

        zind = find(casts(ii).zz >= -200);
        casts(ii).tempda = nanmean(casts(ii).temp(zind));
        casts(ii).saltda = nanmean(casts(ii).salt(zind));

        zindbott = find(casts(ii).zz <= min(casts(ii).zz) + 15);
        casts(ii).tempbott = nanmean(casts(ii).temp(zindbott));
        casts(ii).saltbott = nanmean(casts(ii).salt(zindbott));
    end %for

	if sodaflag == 1
		% run dg_load_SODA.m if SODAcasts.mat does not exist
		if ~exist('sodacast')
			load([workdir 'SODAcasts.mat']);
		end %if
	end %if

    if testplot == 1
        jj = 102
        figure(3)
        hold off
        plot(casts(jj).temp,casts(jj).zz,'x-');
        hold on
        plot(casts(jj).t,casts(jj).z,'ro');
        casts(jj).depth
    end %if

    display('Done loading NEFSC CTD data.')
end %if loadctdflag

%% save grid file to disk
if savegridflag == 1
    display('saving hydroMABgrid to local drive')
    tic
    save([workdir 'hydroMABgrid.mat'],'region','LON','LAT','Z','-v7.3');
    toc
end %if

if loadgridflag == 1
	display('loading hydroMABgrid: LON, LAT, Z, region')
	tic
	load([workdir 'hydroMABgrid.mat']);
	toc
end %if

if savecastflag == 1
	display('saving hydroMABcasts to local drive')
    tic
    save([workdir 'hydroMABcasts.mat'],'casts','sodacast','-v7.3');
    toc
end %if

if loadcastflag == 1
	display('loading hydroMABcasts: casts, sodacast')
	tic
	load([workdir 'hydroMABcasts.mat']);
	toc
end %if

%% calculate temp & salt on grid
if calcflag == 1
    % calculating the nearest CTD station distances (within 30 km) from each grid point location

	if sodaflag == 1
		casts = sodacast;
	end %if

    display('Calculating grid points proximity to sampling points')

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

    % lon, lat locations of grid points from gebco_MAB_30arcsec
    %llon = reshape(LON,[size(LON,1)*size(LON,2), 1]);
    %llat = reshape(LAT,[size(LAT,1)*size(LAT,2), 1]);
    %zz = reshape(Z,[size(Z,1)*size(Z,2), 1]);

    % select only points on the shelf, shelf-break & slope
    %zind = find(zz < -6 & zz >= -1000);
    %lon = llon(zind);
    %lat = llat(zind);
    %z = zz(zind);

    if testplot == 1
        ri = 1
        lon = region(ri).lon;
        lat = region(ri).lat;
        z = region(ri).z;

        figure(2)
        pcolor(LON,LAT,log10(abs(Z))); shading flat; colorbar;
        hold on;
        plot(lon,lat,'.');
        plot(lontrain,lattrain,'kx');
        plot(region(ri).lon, region(ri).lat, 'bo-')
        hold off;
    end %if

    % Finding indices CTD casts for specific year, time of year, and depth range
    %hydro(39).season(3).region(9).temp = [];

    %yyyylist = unique(yyyy);
    yyyylist = [1977:2016];

	hydromo = [];
    %hydromotemp = [];
	%hydromosalt = [];
	%hydromoz = [];

    for mi = 1:12
        hydromo.ind{mi} = find(mm == mi & depth <= 1500);
		%hydromotemp.ind{mi} = find(mm == mi & depth <= 1500);
		%hydromosalt.ind{mi} = find(mm == mi & depth <= 1500);
		%hydromoz.ind{mi} = find(mm == mi & depth <= 1500);
    end %for

    hydromo.yyyy = yyyylist;
	%hydromotemp.yyyy = yyyylist;
	%hydromosalt.yyyy = yyyylist;
	%hydromoz.yyyy = yyyylist;
    for ii = 1:12
        for ri = 1:length(region)
            hydromo.month(ii).region(ri).temp = [];
            hydromo.month(ii).region(ri).salt = [];
            hydromo.month(ii).region(ri).z = [];

			%hydromotemp.month(ii).region(ri).temp = [];
			%hydromotemp.month(ii).region(ri).salt = [];
			%hydromotemp.month(ii).region(ri).z = [];

            %hydromosalt.month(ii).region(ri).temp = [];
			%hydromosalt.month(ii).region(ri).salt = [];
			%hydromosalt.month(ii).region(ri).z = [];

			%hydromoz.month(ii).region(ri).temp = [];
			%hydromoz.month(ii).region(ri).salt = [];
            %hydromoz.month(ii).region(ri).z = [];
        end %for ri
    end %for ii

    if exist('mc') ~= 1
        %mc = [];
        mc = parpool(4)
    end %if

    hind = hydromo.ind;
    mdata = struct();
    for ii = 1:12
        mdata(ii) = struct();
    end %for

	switch hydromoflag
	case 'all'
	    % loop through all the months
	    parfor ii = 1:length(hind)
        %for ii = 7:7
            ii
	        if ~isempty(hind{ii})
	            lontrain = londata(hind{ii});
	            lattrain = latdata(hind{ii});
	            hdoi = casts(hind{ii}); % hydrographic data of interest

	            % constructing KDtree search model using
	            wgs84 = wgs84Ellipsoid('meters');

	            [extrain,eytrain,eztrain] = geodetic2ecef(wgs84,lattrain,lontrain,h);
	            ECEFxyztrain = [extrain,eytrain,eztrain];
	            dModel = createns(ECEFxyztrain,'NSMethod','kdtree','Distance','Euclidean');
	            display('Done calculating KD Tree model.')

                nnindxri = dg_nnindxfind(wgs84,region,h,dModel,searchradius);
                mdata(ii).nnindx = nnindxri;
	            display('Done nearest neighbor calcuation for all the regions.')
                tic
                [ttemp, ssalt, zz, ttempda, ssaltda] =  dg_varfitregion(region,hdoi,nnindxri,searchradius,p);
                toc
	        else
                [ttemp, ssalt, zz, ttempda, ssaltda] =  dg_varempty(region);
	        end % if isempty ind

            mdata(ii).tempfit = ttemp;
            mdata(ii).saltfit = ssalt;
            mdata(ii).zfit = zz;
            mdata(ii).tempdafit = ttempda;
            mdata(ii).saltdafit = ssaltda;
	    end %parfor ind

        for ii = 1:length(hind)
        %for ii = 7:7
            for ri = 1:length(region)
                hydromo.month(ii).region(ri).temp = mdata(ii).tempfit{ri}';
				mdata(ii).tempfit{ri} = [];
                hydromo.month(ii).region(ri).salt = mdata(ii).saltfit{ri}';
				mdata(ii).saltfit{ri} = [];
                hydromo.month(ii).region(ri).z = mdata(ii).zfit{ri}';
				mdata(ii).zfit{ri} = [];
                hydromo.month(ii).region(ri).tempda = mdata(ii).tempdafit{ri}';
				mdata(ii).tempdafit{ri} = [];
                hydromo.month(ii).region(ri).saltda = mdata(ii).saltdafit{ri}';
				mdata(ii).saltdafit{ri} = [];
                hydromo.month(ii).region(ri).nnindx = mdata(ii).nnindx{ri};
				mdata(ii).nnindx{ri} = [];
            end %for
        end %for

        mdata = [];
		hydromo.r = searchradius;
		hydromo.p = p;

		tic
		display(['saving calculated monthly climatology values to disk.'])
		if sodaflag == 0
			save([workdirlocal 'hydroMABmonthlyDA_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydromo','-v7.3');
		elseif sodaflag == 1
			save([workdirlocal 'hydroMABmonthlyDA_SODA_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydromo','-v7.3');
		end %if
		toc
	end % switch

    display('Done interplating temperature and salinity onto bathymetric grid.')
    delete(mc); clear mc;
end %calcflag

%% test plots
if testplot == 1
    % plotting the regions
    figure(2)
    pcolor(LON,LAT,log10(abs(Z))); shading flat; colorbar;
    hold on;
    ri = 2
    plot(region(ri).lon, region(ri).lat, 'bx-')

    % testing & visualizing
    nnindlen = repmat(NaN,size(nnindx));
    for ii = 1:length(nnindx)
        nnindlen(ii) = length(nnindx{ii});
    end %for

    ind1 = find(nnindlen == 5);
    indx = nnindx{ind1(20000)};
    %[lontrain(indx), lattrain(indx)]

    figure(1)
    plot(lontrain,lattrain,'b.');
    hold on
    plot(lontrain(indx), lattrain(indx), 'rx');
end %if

function nnindxri = dg_nnindxfind(wgs84,region,h,dModel,searchradius)
	nnindxri = cell(length(region),1);
    for ri = 1:length(region)
        ri
        [ex,ey,ez] = geodetic2ecef(wgs84,region(ri).llat,region(ri).llon,h);
        %hydromo.month(ii).region(ri).nnindx = rangesearch(dModel, [ex,ey,ez], searchradius);
        nnindxri{ri} = rangesearch(dModel, [ex,ey,ez], searchradius);
        length(nnindxri{ri})
    end %for
    nnindxri
end

function [ttemp, ssalt, zz, ttempda, ssaltda] =  dg_varfitregion(region,hdoi,nnindx,searchradius,p)
	ttemp = cell(1,length(region));
	ssalt = cell(1,length(region));
	zz = cell(1,length(region));
    for ri = 1:length(region)
        tic
       [tempdafit, saltdafit, zfit, tempfit, saltfit] = dg_idw_casts2d(region(ri).llon, region(ri).llat, region(ri).zz, hdoi, nnindx{ri},searchradius,p);
        %hydromo.month(ii).region(ri).temp = tempfit;
        %hydromo.month(ii).region(ri).salt = saltfit;
        %hydromo.month(ii).region(ri).z = zfit;
        ttemp{ri} = tempfit;
        ssalt{ri} = saltfit;
        zz{ri} = zfit;
        ttempda{ri} = tempdafit;
        ssaltda{ri} = saltdafit;
        toc
    end %for
end

function [ttemp, ssalt, zz, ttempda, ssaltda] =  dg_varempty(region)
    for ri = 1:length(region)
        ttemp{ri} = [];
        ssalt{ri} = [];
        zz{ri} = [];
        ttempda{ri} = [];
        ssaltda{ri} = [];
    end %for
end
