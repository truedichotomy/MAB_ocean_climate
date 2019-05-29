% this script loads SODA 3.4.2 data files
% DG 2018-03-26, 2018-04-01, 2018-04-03

testflag = 0
loadbathyflag = 1
loadctdflag = 1
loadsodaflag = 1

dg_setup_MABclimate_dir

if loadbathyflag == 1
	load([bathydir 'gebco_MAB_30arcsec.mat']);
    display('Done loading bathy data.')
    negLON = find(LON < 0);
    LON(negLON) = LON(negLON) + 360;
end %if

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

        % calculating depth-averaged temperature and salinity
        zind = find(casts(ii).zz >= -200);
        casts(ii).tempda = nanmean(casts(ii).temp(zind));
        casts(ii).saltda = nanmean(casts(ii).salt(zind));

        % calculating bottom 15 m temperature and salinity
        zindbott = find(casts(ii).zz <= min(casts(ii).zz) + 15);
        casts(ii).tempbott = nanmean(casts(ii).temp(zindbott));
        casts(ii).saltbott = nanmean(casts(ii).salt(zindbott));

    end %for

    if testflag == 1
        jj = 102
        figure(3)
        hold off
        plot(casts(jj).temp,casts(jj).zz,'x-');
        hold on
        plot(casts(jj).t,casts(jj).z,'ro');
        casts(jj).depth
    end %if
    save([workdir 'MABctdcasts.mat'],'casts')
    display('Done loading NEFSC CTD data.')
end %if loaddataflag

if loadsodaflag == 1
    % define soda data structure load file path and file names
    soda.readme = 'http://www.atmos.umd.edu/%7Eocean/index_files/soda3_readme.htm';
    soda.version = '3.4.2';
    SODAdir = '/Users/c2po/Research/SODA/';
    soda.files = dir([SODAdir '*.nc']);
    soda.yyyy = [1980:2015];

    % create cast list to match in SODA
    %tlonlatlist = [[datenum(1985,12,1), 285, 39]; [datenum(2013,12,1), 292.5 42]];
    clon = [casts.lon];
    clat = [casts.lat];
    neglon = find(clon < 0);
    clon(neglon) = clon(neglon) + 360;
    ct = datenum([casts.yr],1,[casts.dyd]);
    tlonlatlist = [ct',clon',clat'];
    toutind = find(ct < datenum(soda.yyyy(1),1,1) | ct >= datenum(soda.yyyy(end),13,1));
    tintind = find(ct >= datenum(soda.yyyy(1),1,1) & ct < datenum(soda.yyyy(end),13,1));
    tlonlatlist(toutind,:) = [];

    yyyyin = str2num(datestr(tlonlatlist(:,1),10));
    mmin = str2num(datestr(tlonlatlist(:,1),5));

    % load SODA file objecdts using nctoolbox
    clear sds vart varxt varyt varst vartemp varsalt
    for fi = 1:length(soda.files)
        sodafile = [soda.files(fi).folder '/' soda.files(fi).name]
        sds{fi} = cfdataset(sodafile);
    end

    % extract variable sizes, they don't change
    stemp = sds{1}.size('temp');
    sx = sds{1}.size('xt_ocean');
    sy = sds{1}.size('yt_ocean');
    sz = sds{1}.size('st_ocean');

    % define LON and LAT SODA
    varxt = sds{1}.variable('xt_ocean');
    varyt = sds{1}.variable('yt_ocean');
    varst = sds{1}.variable('st_ocean');

    xt = varxt.data(1,sx);
    yt = varyt.data(1,sy);
    [soda.LON, soda.LAT] = meshgrid(xt,yt);

    % extract SODA temp and salt cast values for tlonlatlist locations
    tic
    sodacast = [];
    nctdcasts = size(tlonlatlist,1);
    for ii = 1:nctdcasts
        [ii, nctdcasts]
        lat0 = tlonlatlist(ii,3);
        lon0 = tlonlatlist(ii,2);
        [lattmp, lati] = min(abs(lat0-yt));
        [lontmp, loni] = min(abs(lon0-xt));
        fi = find(yyyyin(ii) == soda.yyyy);
        mi = mmin(ii);

        %vart = sds{fi}.variable('time');
        %varxt = sds{fi}.variable('xt_ocean');
        %varyt = sds{fi}.variable('yt_ocean');
        varst = sds{fi}.variable('st_ocean');
        vartemp = sds{fi}.variable('temp');
        varsalt = sds{fi}.variable('salt');

        sodacast(ii).time = tlonlatlist(ii,1);
        sodacast(ii).yr = casts(tintind(ii)).yr;
        sodacast(ii).yd = casts(tintind(ii)).yd;
        sodacast(ii).dyd = casts(tintind(ii)).dyd;
        sodacast(ii).tempSODA = vartemp.data([mi,1,lati,loni],[mi,sz,lati,loni])';
        sodacast(ii).saltSODA = varsalt.data([mi,1,lati,loni],[mi,sz,lati,loni])';
        sodacast(ii).zSODA = -varst.data(1,sz);
        sodacast(ii).lonSODA = xt(loni);
        sodacast(ii).latSODA = yt(lati);
        sodacast(ii).depthSODA = max(sodacast(ii).zSODA(find(~isnan(sodacast(ii).tempSODA))));
        sodacast(ii).zdeepSODA = -sodacast(ii).depthSODA;

        sodacast(ii).lon = lon0;
        sodacast(ii).lat = lat0;
        sodacast(ii).zdeep = floor(interp2(LON,LAT,Z,[sodacast(ii).lon],[sodacast(ii).lat],'linear'));
        sodacast(ii).zz = [-1:-1:sodacast(ii).zdeep];
        sodacast(ii).temp = repmat(NaN,size(sodacast(ii).zz));
        sodacast(ii).salt = repmat(NaN,size(sodacast(ii).zz));

        % calculate the interpolated temperature and salinity from SODA grid to a 1 m grid
        nnanindSODA = find(~isnan(sodacast(ii).tempSODA));
        if length(nnanindSODA) >= 2
            sodacast(ii).temp = interp1(sodacast(ii).zSODA,sodacast(ii).tempSODA,sodacast(ii).zz,'linear',NaN);
            sodacast(ii).salt = interp1(sodacast(ii).zSODA,sodacast(ii).saltSODA,sodacast(ii).zz,'linear',NaN);
        elseif length(nnanindSODA) == 1
            sodacast(ii).temp = sodacast(ii).tempSODA(nnanindSODA);
            sodacast(ii).salt = sodacast(ii).saltSODA(nnanindSODA);
        else
            sodacast(ii).temp = repmat(NaN,size(sodacast(ii).zz));
            sodacast(ii).salt = repmat(NaN,size(sodacast(ii).zz));
        end %if

        % fill in the bottom boundary layer temp and salinity with the deepest SODA value if the bathy grid is deeper than SODA grid
        diffdepth = sodacast(ii).zdeep - sodacast(ii).zdeepSODA; % positive if shallower than SODA, negative if deeper than SODA
        if diffdepth >= -30 & diffdepth < 0
            zbblind = find(sodacast(ii).zz < sodacast(ii).zdeepSODA & sodacast(ii).zz > sodacast(ii).zdeep); % find indices of the BBL
            sodacast(ii).temp(zbblind) = sodacast(ii).tempSODA(end);
            sodacast(ii).salt(zbblind) = sodacast(ii).saltSODA(end);
        end %if

        % fill in the surface 5 m temp and salinity with the shallowest SODA value
        sodacast(ii).temp(1:5) = sodacast(ii).tempSODA(1);
        sodacast(ii).salt(1:5) = sodacast(ii).saltSODA(1);

        % calculating depth-averaged temperature and salinity
        zind = find(sodacast(ii).zz >= -200);
        sodacast(ii).tempda = nanmean(sodacast(ii).temp(zind));
        sodacast(ii).saltda = nanmean(sodacast(ii).salt(zind));

        % calculating bottom 15 m temperature and salinity
        zindbott = find(sodacast(ii).zz <= min(sodacast(ii).zz) + 15);
        sodacast(ii).tempbott = nanmean(sodacast(ii).temp(zindbott));
        sodacast(ii).saltbott = nanmean(sodacast(ii).salt(zindbott));

    end %for
    toc
    save([workdir 'SODAcasts.mat'],'soda','sodacast')
end %loadsodaflag

%for ii = 1:nctdcasts
%    sodacast(ii).zz = sodacast(ii).zz';
%    sodacast(ii).temp = sodacast(ii).temp';
%    sodacast(ii).salt = sodacast(ii).salt';
%end %if

% load temperature and salinity
%vart = sds.variable('time');
%varxt = sds.variable('xt_ocean');
%varyt = sds.variable('yt_ocean');
%varst = sds.variable('st_ocean');
%vartemp = sds.variable('temp');
%varsalt = sds.variable('salt');

if testflag == 1
    fi = 1
    mi = 6
    zi = 10
    temp = vartemp.data([mi,zi,1,1],[mi,zi,sy,sx]);
    salt = varsalt.data([mi,zi,1,1],[mi,zi,sy,sx]);
    xt = varxt.data(1,sx);
    yt = varyt.data(1,sy);
    zt = varst.data(1,sz);

    [soda.LON, soda.LAT] = meshgrid(xt,yt);
    soda.z = zt;
    soda.zind = find(soda.z <= 210);
    soda.z(soda.zind)
    % data structure: time (12 months), depth (50 layers), lat (330 half degree), lon (720 half degree)

    %pcolor(soda.LON,soda.LAT,squeeze(temp(1,1,:,:))); shading flat; colorbar;

    pcolor(soda.LON,soda.LAT,squeeze(salt(1,1,:,:))); shading flat; colorbar; caxis([28 38]);
end %if
