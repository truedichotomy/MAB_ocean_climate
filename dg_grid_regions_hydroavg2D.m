% this script calculates the average temperature and salinity for each MAB region during each season for each year
% Note this script also takes a LONG time (10-16 hours) to run on a quad core machine with 16 GB of RAM!!!
% DG 2017-10-10, 2017-10-14, 2017-10-16, 2017-11-01, 2017-12-15

dbstop if error

dataflag = 'depth-averaged' % 'bottom' or 'depth-averaged' or 'bottomSODA' or 'daSODA', only the first and last letters matter

loadgridflag = 1
loaddataflag = 1

initflag = 1

searchradius = 40000
p = 1

% setup directories needed for the script
dg_setup_MABclimate_dir

% Define years
%workdir = '/Volumes/Passport16/MAB_climate/';
%yyyylist = [1980];
%yyyy = [casts.yr];

% load grid data
if loadgridflag == 1
    gridfile = 'hydroMABgrid.mat';
    gridpath = [workdir gridfile];
    load(gridpath);
end %if

% load data for each year if needed
if loaddataflag == 1
    if lower(dataflag(1)) == 'd' & dataflag(end) ~= 'A'
        datafile = ['hydroMAB2Da_r' num2str(searchradius) '_p' num2str(p) '.mat'];
    elseif lower(dataflag(1)) == 'b' & dataflag(end) ~= 'A'
        datafile = ['hydroMAB2Db_r' num2str(searchradius) '_p' num2str(p) '.mat'];
    elseif lower(dataflag(1)) == 'd' & dataflag(end) == 'A'
        datafile = ['hydroMAB2DaSODA_r' num2str(searchradius) '_p' num2str(p) '.mat'];
    elseif lower(dataflag(1)) == 'b' & dataflag(end) == 'A'
        datafile = ['hydroMAB2DbSODA_r' num2str(searchradius) '_p' num2str(p) '.mat'];
    end %case

    datapath = [workdirlocal datafile];
    load(datapath);
end %if

%if exist('mc') ~= 1
%    mc = parpool(4)
%end %if

hydrom = []; hydrob = [];

if initflag == 1
    %yyyylist = unique(yyyy);
    yyyylist = [1977:2016];

    for yi = 1:length(yyyylist)

        hydrom(yi).yyyy = yyyylist(yi);

        for ii = 1:3
            for ri = 1:length(region)
                hydrom(yi).season(ii).region(ri).temp = [];
                hydrom(yi).season(ii).region(ri).salt = [];
                hydrom(yi).season(ii).region(ri).voltemp = [];
                hydrom(yi).season(ii).region(ri).volsalt = [];
                hydrom(yi).season(ii).region(ri).tempshw = [];
                hydrom(yi).season(ii).region(ri).saltshw = [];
                hydrom(yi).season(ii).region(ri).volshw = [];
                hydrom(yi).season(ii).region(ri).tempcp = [];
                hydrom(yi).season(ii).region(ri).saltcp = [];
                hydrom(yi).season(ii).region(ri).volcp = [];
                hydrom(yi).season(ii).region(ri).volume = [];
                hydrom(yi).season(ii).region(ri).z = [];
                hydrom(yi).season(ii).region(ri).gridarea = [];
                hydrom(yi).season(ii).region(ri).area = [];
            end %for ri
        end %for ii
    end %for yi
end %if

% loop through all the years of interest to calculate regional mean
for yi = 1:length(yyyylist) % years
    yyyy = yyyylist(yi);

    for si = 1:3 % seasons
        for ri = 1:length(region)
            [yyyy si ri]
            tic
            tempvol = 0; saltvol = 0; tempshwvol = 0; saltshwvol = 0;
            voltemp = 0; volsalt = 0; volshw = 0; volume = 0;

            % load temperature, salinity, and depth data from interpolated data structure
            switch lower(dataflag(1))
            case 'd'
                ndata = length(hydro(yi).season(si).region(ri).tempda);
                temp = hydro(yi).season(si).region(ri).tempda;
                salt = hydro(yi).season(si).region(ri).saltda;
            case 'b'
                ndata = length(hydro(yi).season(si).region(ri).tempbott);
                temp = hydro(yi).season(si).region(ri).tempbott;
                salt = hydro(yi).season(si).region(ri).saltbott;
            end %switch
            zdeep = hydro(yi).season(si).region(ri).zdeep;

            temp = reshape(temp,[length(temp), 1]);
            salt = reshape(salt,[length(salt), 1]);
            zdeep = reshape(zdeep,[length(zdeep), 1]);

            % find bathy grid points with temperature, salinity, or SHW data
            tempind = find(~isnan(temp)); % grid points with temp data
            saltind = find(~isnan(salt)); % grid points with salt data
            shwind = find(salt <= 34); % find shelf water volume

            % find volumes
            %gridarea = 0.5 * 0.5 * 1852*1852 ./ cosd(region(ri).llat);
            gridarea = 0.5 * 0.5 * 1852*1852 .* cosd(region(ri).llat);
            gridarea = reshape(gridarea,[length(gridarea), 1]);

            % calculate the volume of each grid point for Depth averaged or bottom (15 m) cases
            switch lower(dataflag(1))
            case 'd'
                volume = abs(zdeep) .* gridarea; % volume of each grid point
            case 'b'
                volume = min([abs(zdeep),repmat(15,[length(zdeep), 1])],[],2) .* gridarea; % volume of each grid point
            end %switch

            if ~isempty(tempind)
                voltemp = nansum(volume(tempind));
                tempavg = nansum(temp(tempind) .* volume(tempind))/voltemp;
            else
                voltemp = NaN;
                tempavg = NaN;
            end %if

            if ~isempty(saltind)
                volsalt = nansum(volume(saltind));
                saltavg = nansum(salt(saltind) .* volume(saltind))/volsalt;
            else
                volsalt = NaN;
                saltavg = NaN;
            end %if

            if ~isempty(shwind)
                volshw = nansum(volume(shwind));
                tempavgshw = nansum(temp(shwind) .* volume(shwind))/volshw;
                saltavgshw = nansum(salt(shwind) .* volume(shwind))/volshw;
            else
                volshw = NaN;
                tempavgshw = NaN;
                saltavgshw = NaN;
            end %if
            % calculate average temperature and salinity for each region of each season of each year

            hydrom(yi).season(si).region(ri).voltemp = voltemp;
            hydrom(yi).season(si).region(ri).volsalt = volsalt;
            hydrom(yi).season(si).region(ri).temp = tempavg;
            hydrom(yi).season(si).region(ri).salt = saltavg;
            hydrom(yi).season(si).region(ri).tempshw = tempavgshw;
            hydrom(yi).season(si).region(ri).saltshw = saltavgshw;
            hydrom(yi).season(si).region(ri).volshw = volshw;
            hydrom(yi).season(si).region(ri).volume = nansum(volume);
            hydrom(yi).season(si).region(ri).z = nanmean(zdeep);
            hydrom(yi).season(si).region(ri).gridarea = gridarea;
            hydrom(yi).season(si).region(ri).area = nansum(gridarea);
            toc
        end %for ri
    end %for si
end %for yi

if lower(dataflag(1)) == 'd' & dataflag(end) ~= 'A'
    save([workdirlocal 'hydroMABavg2Da_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydrom','-v7.3');
elseif lower(dataflag(1)) == 'b' & dataflag(end) ~= 'A'
    save([workdirlocal 'hydroMABavg2Db_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydrob','-v7.3');
elseif lower(dataflag(1)) == 'd' & dataflag(end) == 'A'
    save([workdirlocal 'hydroMABavg2DaSODA_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydrom','-v7.3');
elseif lower(dataflag(1)) == 'b' & dataflag(end) == 'A'
    save([workdirlocal 'hydroMABavg2DbSODA_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydrom','-v7.3');
end %if
