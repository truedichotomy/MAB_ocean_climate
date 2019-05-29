% this script uses NOAA NEFSC hydrographic data to conduct regional climate analysis of temperature and salinity for each region in the U.S. Northeast.
% this script should be run first. it calls dg_grid_regions_define.m and it should be followed by dg_grid_regions_hydroavg.m
%
% Note this script takes a LONG time (10-16 hours) to run on a quad core machine with 16 GB of RAM!!!
%
% Donglai Gong 2017-09-22, 2017-10-09, 2019-01-13

dbstop if error

testplot = 0

loadbathyflag = 1
calcregionflag = 1
loadctdflag = 1
sodaflag = 0
savegridflag = 1
loadgridflag = 1
savecastflag = 1
loadcastflag = 1
calcflag = 1
dataflag = 'depth-averaged' % 'bottom' or 'depth-averaged' or 'bottomSODA' or 'depth-averagedSODA', only the first and last letters matter

dg_setup_MABclimate_dir

% parameters for inverse-distance search algorithm
h = 0; %meters
searchradius = 25000; % 25 km
p = 2; % power parameter, default is 2

depthlim = [6 1000]; % looking at shelf regions in between these two isobaths
%depthlim = []

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
        if abs(casts(ii).zdeep - min(casts(ii).z)) <= 30
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

        % calculating depth-averaged temperature and salinity
        zind = find(casts(ii).zz >= -200);
        casts(ii).tempda = nanmean(casts(ii).temp(zind));
        casts(ii).saltda = nanmean(casts(ii).salt(zind));

        % calculating bottom 15 m temperature and salinity
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
end %if loaddataflag

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
	hydro = [];

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

    for yi = 1:length(yyyylist)

        % winter/spring
        %hydro(yi).ind{1} = find(yyyy == yyyylist(yi) & dyd <= 125 & dyd >= 0 & depth <= 1500);
        % spring/summer
        %hydro(yi).ind{2} = find(yyyy == yyyylist(yi) & dyd <= 250 & dyd >= 126 & depth <= 1500);
        % fall/winter
        %hydro(yi).ind{3} = find(yyyy == yyyylist(yi) & dyd <= 365 & dyd >= 251 & depth <= 1500);

		% winter/spring new
		hydro(yi).ind{1} = find(yyyy == yyyylist(yi) & dyd <= 120 & dyd >= 0 & depth <= 1500);
		% spring/summer new
		hydro(yi).ind{2} = find(yyyy == yyyylist(yi) & dyd <= 243 & dyd >= 121 & depth <= 1500);
		% fall/winter new
		hydro(yi).ind{3} = find(yyyy == yyyylist(yi) & dyd <= 365 & dyd >= 244 & depth <= 1500);

        hydro(yi).yyyy = yyyylist(yi);

        % loop through the seasons to intialize hydro
		switch lower(dataflag(1))
		case 'd'
	        for ii = 1:length(hydro(yi).ind)
	            for ri = 1:length(region)
	                hydro(yi).season(ii).region(ri).tempda = [];
	                hydro(yi).season(ii).region(ri).saltda = [];
	                hydro(yi).season(ii).region(ri).zdeep = [];
	            end %for ri
	        end %for ii
		case 'b'
			for ii = 1:length(hydro(yi).ind)
	            for ri = 1:length(region)
					hydro(yi).season(ii).region(ri).tempbott = [];
	                hydro(yi).season(ii).region(ri).saltbott = [];
	                hydro(yi).season(ii).region(ri).zdeep = [];
	            end %for ri
	        end %for ii
		end %if
    end %for yi

    if exist('mc') ~= 1
        mc = parpool(6)
    end %if

    %hydro1 = hydro(1);
    % loop through the years
    parfor yi = 1:length(yyyylist)
        yyyylist(yi)
        % loop through all the seasons
        for ii = 1:length(hydro(yi).ind)
            if ~isempty(hydro(yi).ind{ii})
                lontrain = londata(hydro(yi).ind{ii});
                lattrain = latdata(hydro(yi).ind{ii});
                hdoi = casts(hydro(yi).ind{ii}); % hydrographic data of interest

                % constructing KDtree search model using
                tic
                wgs84 = wgs84Ellipsoid('meters');

                [extrain,eytrain,eztrain] = geodetic2ecef(wgs84,lattrain,lontrain,h);
                ECEFxyztrain = [extrain,eytrain,eztrain];
                dModel = createns(ECEFxyztrain,'NSMethod','kdtree','Distance','Euclidean');
                display('Done calculating KD Tree model.')

                nnindxri = dg_nnindxfind(wgs84,region,h,dModel,searchradius)
                %hydro(yi).season(ii).nnindx = nnindxri;

                display('Done nearest neighbor calcuation for all the regions.')
                toc

                tic
                %[ttempda, ssaltda, zzdeep, ttempda, ssaltda, zz] =  dg_varfitregion(region,hdoi,nnindxri,searchradius);
				switch lower(dataflag(1))
				case 'd'
                	[ttempda, ssaltda, zzdeep] =  dg_varfitregion(region,hdoi,nnindxri,searchradius,p,dataflag);
				case 'b'
					[ttempbott, ssaltbott, zzdeep] =  dg_varfitregion(region,hdoi,nnindxri,searchradius,p,dataflag);
				end %switch
                toc
            else
				switch lower(dataflag(1))
				case 'd'
                	[ttempda, ssaltda, zzdeep] =  dg_varempty(region);
				case 'b'
					[ttempbott, ssaltbott, zzdeep] =  dg_varempty(region);
				end %switch
            end % if isempty ind

			switch lower(dataflag(1))
			case 'd'
	            for ri = 1:length(region)
	                hydro(yi).season(ii).region(ri).tempda = ttempda{ri};
	                hydro(yi).season(ii).region(ri).saltda = ssaltda{ri};
	                hydro(yi).season(ii).region(ri).zdeep = zzdeep{ri};
	            end %for
	            ttempda = []; ssaltda = []; zzdeep = [];
			case 'b'
				for ri = 1:length(region)
	                hydro(yi).season(ii).region(ri).tempbott = ttempbott{ri};
	                hydro(yi).season(ii).region(ri).saltbott = ssaltbott{ri};
	                hydro(yi).season(ii).region(ri).zdeep = zzdeep{ri};
	            end %for
	            ttempbott = []; ssaltbott = []; zzdeep = [];
			end %switch
        end %for ind
    end %parfor yyyylist
    display(['saving calculated temp and salt values to disk.'])
    tic
    dg_savehydro(hydro, workdirlocal, dataflag, searchradius, p)
    toc

    display('Done interplating temperature and salinity onto bathymetric grid.')
    if exist('mc')
        delete(mc)
    end %if
end %calcflag

%% define a save gridded data function
function dg_savehydro(hydro, workdir,dataflag,searchradius,p)
	if ~exist('searchradius')
		searchradius = 30000 %m
	end

	if ~exist('p')
		p = 2
	end %if

	if lower(dataflag(1)) == 'd' & dataflag(end) ~= 'A'
		save([workdir 'hydroMAB2Da_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydro','-v7.3');
	elseif lower(dataflag(1)) == 'b' & dataflag(end) ~= 'A'
		save([workdir 'hydroMAB2Db_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydro','-v7.3');
	elseif lower(dataflag(1)) == 'd' & dataflag(end) == 'A'
		save([workdir 'hydroMAB2DaSODA_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydro','-v7.3');
	elseif lower(dataflag(1)) == 'b' & dataflag(end) == 'A'
		save([workdir 'hydroMAB2DbSODA_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydro','-v7.3');
	end %if
end

function nnindxri = dg_nnindxfind(wgs84,region,h,dModel,searchradius)
	nnindxri = cell(length(region),1);
    for ri = 1:length(region)
        %ri
        [ex,ey,ez] = geodetic2ecef(wgs84,region(ri).llat,region(ri).llon,h);
        %hydromo.month(ii).region(ri).nnindx = rangesearch(dModel, [ex,ey,ez], searchradius);
        nnindxri{ri} = rangesearch(dModel, [ex,ey,ez], searchradius);
        %length(nnindxri{ri})
    end %for
    %nnindxri
end

function [ttempda, ssaltda, zzdeep, ttemp, ssalt, zz] =  dg_varfitregion(region,hdoi,nnindx,searchradius,p,datatype)
	if ~exist('datatype')
		datetype = 'depth-averaged';
	end %if
    ttempda = cell(1,length(region));
	ssaltda = cell(1,length(region));
    zzdeep = cell(1,length(region));

    if nargout > 3
        ttemp = cell(1,length(region));
        ssalt = cell(1,length(region));
        zz = cell(1,length(region));
    end %if

    for ri = 1:length(region)
        ri
        tic
        if nargout < 4
			switch lower(datatype(1))
			case 'd'
            	[tempdafit, saltdafit] = dg_idw_casts2d(region(ri).llon, region(ri).llat, region(ri).zz, hdoi, nnindx{ri},searchradius, p);
			case 'b'
				[tempdafit, saltdafit] = dg_idw_casts2dbott(region(ri).llon, region(ri).llat, region(ri).zz, hdoi, nnindx{ri},searchradius, p);
			end %switch
            ttempda{ri} = tempdafit;
            ssaltda{ri} = saltdafit;
            zzdeep{ri} = region(ri).zz;
        elseif nargout >= 4
			switch lower(datatype(1))
			case 'd'
            	[tempdafit, saltdafit, zfit, tempfit, saltfit] = dg_idw_casts2d(region(ri).llon, region(ri).llat, region(ri).zz, hdoi, nnindx{ri}, searchradius, p);
			case 'b'
				[tempdafit, saltdafit, zfit, tempfit, saltfit] = dg_idw_casts2dbott(region(ri).llon, region(ri).llat, region(ri).zz, hdoi, nnindx{ri}, searchradius, p);
			end %switch
            zz{ri} = zfit;
            ttemp{ri} = tempfit;
            ssalt{ri} = saltfit;
            zzdeep{ri} = region(ri).zz;
            ttempda{ri} = tempdafit;
            ssaltda{ri} = saltdafit;
        end %if nargout
        toc
    end %for
end

function [ttempda, ssaltda, zzdeep, ttemp, ssalt, zz] =  dg_varempty(region)
    if nargout < 4
        for ri = 1:length(region)
            ttempda{ri} = [];
            ssaltda{ri} = [];
            zzdeep{ri} = region(ri).zz;
        end %for
    elseif nargout >= 4
        for ri = 1:length(region)
            ttemp{ri} = [];
            ssalt{ri} = [];
            zz{ri} = [];
            ttempda{ri} = [];
            ssaltda{ri} = [];
            zzdeep{ri} = region(ri).zz;
        end %for
    end
end
