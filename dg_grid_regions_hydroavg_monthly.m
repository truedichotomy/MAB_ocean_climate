% this script calculates the average temperature and salinity for each MAB region during each month for each year
% Note this script also takes a LONG time (10-16 hours) to run on a quad core machine with 16 GB of RAM!!!
% run dg_grid_regions_monthly2D.m before this script
% DG 2017-10-10, 2017-10-14, 2017-10-16, 2017-11-01, 2018-04-06

dbstop if error

loadgridflag = 0
initflag = 1 %
loaddataflag = 1 %
sodaflag = 0

searchradius = 20000
p = 4

dg_setup_MABclimate_dir

% load grid data
if loadgridflag == 1
    gridfile = 'hydroMABgrid.mat';
    gridpath = [workdir gridfile];
    load(gridpath);
end %if

%if exist('mc') ~= 1
%    mc = parpool(4)
%end %if

if initflag == 1
    yyyylist = [1977:2016];
    hydromom.yyyy = yyyylist;
    for ii = 1:12
        for ri = 1:length(region)
            hydromom.month(ii).region(ri).mtemp = [];
            hydromom.month(ii).region(ri).msalt = [];
            hydromom.month(ii).region(ri).voltemp = [];
            hydromom.month(ii).region(ri).volsalt = [];
            hydromom.month(ii).region(ri).mtempshw = [];
            hydromom.month(ii).region(ri).msaltshw = [];
            hydromom.month(ii).region(ri).volshw = [];
            hydromom.month(ii).region(ri).tempcp = [];
            hydromom.month(ii).region(ri).saltcp = [];
            hydromom.month(ii).region(ri).volcp = [];
            hydromom.month(ii).region(ri).volume = [];
            hydromom.month(ii).region(ri).z = [];
        end %for ri
    end %for ii
end %if

% load data for each year if needed
if loaddataflag == 1
    if sodaflag == 0
        datafile = ['hydroMABmonthlyDA_r' num2str(searchradius) '_p' num2str(p) '.mat'];
    elseif sodaflag == 1
        datafile = ['hydroMABmonthlyDA_SODA_r' num2str(searchradius) '_p' num2str(p) '.mat'];
    end %if
    datapath = [workdirlocal datafile];
    load(datapath);
elseif loaddataflag == 2
	if sodaflag == 0
        datafile = ['hydroMABmonthlyDA_r' num2str(searchradius) '_p' num2str(p) '.mat'];
    elseif sodaflag == 1
        datafile = ['hydroMABmonthlyDA_SODA_r' num2str(searchradius) '_p' num2str(p) '.mat'];
    end %if
    datapath = [workdir datafile];
    load(datapath);
end %if

for mi = 1:12 % months of a year
    for ri = 1:length(region)
        [mi ri]
        tic
        tempvol = 0; saltvol = 0; tempshwvol = 0; saltshwvol = 0;
        ndata = length(hydromo.month(mi).region(ri).temp);
        voltemp = 0; volsalt = 0; volshw = 0; volume = 0;
        for ii = 1:ndata
            if ~isempty(hydromo.month(mi).region(ri).z{ii})
                gridarea = 0.5 * 0.5 * 1852*1852/cosd(region(ri).llat(ii));
                volume = volume + max(abs(hydromo.month(mi).region(ri).z{ii}))*gridarea;

                % determine min vector length in z
                mlztemp = min(length(hydromo.month(mi).region(ri).temp{ii}),length(hydromo.month(mi).region(ri).z{ii}));
                mlzsalt = min(length(hydromo.month(mi).region(ri).salt{ii}),length(hydromo.month(mi).region(ri).z{ii}));
                tempshwflag = 0;

                % do the temp calculation if there are two or more not NaN elements in the profile
                if length(find(~isnan(hydromo.month(mi).region(ri).temp{ii}(1:mlztemp)) == 1)) > 1
                    % min vector length in z

                    zz = hydromo.month(mi).region(ri).z{ii}(1:mlztemp);
                    z1 = [-diff(zz);1];
                    voli = z1(1:mlztemp) * gridarea;

                    alltempind = find(~isnan(hydromo.month(mi).region(ri).temp{ii}(1:mlztemp))); % find z indices that are not NaN

                    temp = hydromo.month(mi).region(ri).temp{ii}(1:mlztemp);
                    temp = reshape(temp,length(temp),1);
                    tempfit = interp1(zz(alltempind),temp(alltempind),zz,'linear');
                    %hydromo.month(mi).region(ri).temp{ii} = tempfit;

                    voltemp = voltemp + nansum(voli(alltempind));

                    tempvol = tempvol + nansum(temp(alltempind) .* voli(alltempind));
                    tempshwflag = 1;
                end %if

                % do the salt calculation if there are two or more not NaN elements in the profile
                if length(find(~isnan(hydromo.month(mi).region(ri).salt{ii}(1:mlzsalt)) == 1)) > 1

                    zz = hydromo.month(mi).region(ri).z{ii}(1:mlzsalt);
                    z1 = [-diff(zz);1];
                    voli = z1(1:mlzsalt) * gridarea;

                    shwind = find(hydromo.month(mi).region(ri).salt{ii}(1:mlzsalt) <= 34);
                    allsaltind = find(~isnan(hydromo.month(mi).region(ri).salt{ii}(1:mlzsalt))); % find z indices that are not NaN

                    salt = hydromo.month(mi).region(ri).salt{ii}(1:mlzsalt);
                    salt = reshape(salt,length(salt),1);
                    saltfit = interp1(zz(allsaltind),salt(allsaltind),zz,'linear');
                    %hydromo.month(mi).region(ri).salt{ii} = saltfit;

                    volsalt = volsalt + nansum(voli(allsaltind));
                    volshw = volshw + nansum(voli(shwind));

                    saltvol = saltvol + nansum(salt(allsaltind) .* voli(allsaltind));
                    saltshwvol = saltshwvol + nansum(salt(shwind) .* voli(shwind));
                    if tempshwflag == 1
                        tempshwvol = tempshwvol + nansum(temp(shwind) .* voli(shwind));
                    end %if
                end %if
            end %if
        end %for ii
        %dbstop
        hydromom.month(mi).region(ri).voltemp = voltemp;
        hydromom.month(mi).region(ri).volsalt = volsalt;
        hydromom.month(mi).region(ri).mtemp = tempvol / voltemp;
        hydromom.month(mi).region(ri).msalt = saltvol / volsalt;
        hydromom.month(mi).region(ri).mtempshw = tempshwvol / volshw;
        hydromom.month(mi).region(ri).msaltshw = saltshwvol / volshw;
        hydromom.month(mi).region(ri).volshw = volshw;
        hydromom.month(mi).region(ri).volume = volume;

        %hydromo.month(mi).region(ri).voltemp = voltemp;
        %hydromo.month(mi).region(ri).volsalt = volsalt;
        %hydromo.month(mi).region(ri).mtemp = tempvol / voltemp;
        %hydromo.month(mi).region(ri).msalt = saltvol / volsalt;
        %hydromo.month(mi).region(ri).mtempshw = tempshwvol / volshw;
        %hydromo.month(mi).region(ri).msaltshw = saltshwvol / volshw;
        %hydromo.month(mi).region(ri).volshw = volshw;
        %hydromo.month(mi).region(ri).volume = volume;

        toc
    end %for ri
end %for mi

hydromo = [];

if sodaflag == 0
    save([workdirlocal 'hydroMABavgmonthly_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydromom','-v7.3');
elseif sodaflag == 1
    save([workdirlocal 'hydroMABavgmonthlySODA_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydromom','-v7.3');
end %if
