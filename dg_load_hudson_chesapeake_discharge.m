% this script loads river discharge data from tributaries of the Hudson and Chesapeake Bay estuaries and calculates the discharge for each estuary over the time period specified in t1 and t2.
% data is organized into daily, 3 seasons for each year, monthly for each year, and daily climatology
% the script gracefully hand 'Ice' and no data situations for each river
% DG 20180307 gong@vims.edu

workdir = ['./']
loadflag = 1

t1 = datenum(1977,1,1);
t2 = datenum(2016,12,31);

if loadflag == 1
    % Hudson River main branch
    [usgc,staid,yyyy,mm,dd,flow_ftedward_ft3pm_str,stastat] = textread([workdir 'USGS_01327750_FortEdwardNY_197612-201803.txt'],'%s%s%n-%n-%n%s%s','delimiter','\t ','commentstyle','shell','headerlines',31);
    niceind = find(~strcmp(flow_ftedward_ft3pm_str, 'Ice') & ~strcmp(flow_ftedward_ft3pm_str, ''));
    flow_ftedward_ft3pm = repmat(NaN,size(flow_ftedward_ft3pm_str));
    flow_ftedward_ft3pm(niceind) = str2num(char(string(flow_ftedward_ft3pm_str(niceind))));
    t_ftedward = datenum(yyyy,mm,dd);
    flow_ftedward = flow_ftedward_ft3pm * 0.0283168; % convert to cubic meters per sec.
    stastat_ftedward = stastat;
    staid_ftedward = staid(1);

    % Mohawk River that feeds into the Hudson
    [usgc,staid,yyyy,mm,dd,flow_cohoes_ft3pm_str,stastat] = textread([workdir 'USGS_01357500_CohoesNY_191712-201803.txt'],'%s%s%n-%n-%n%s%s','delimiter','\t ','commentstyle','shell','headerlines',31);
    niceind = find(~strcmp(flow_cohoes_ft3pm_str, 'Ice') & ~strcmp(flow_cohoes_ft3pm_str, ''));
    flow_cohoes_ft3pm = repmat(NaN,size(flow_cohoes_ft3pm_str));
    flow_cohoes_ft3pm(niceind) = str2num(char(string(flow_cohoes_ft3pm_str(niceind))));
    t_cohoes = datenum(yyyy,mm,dd);
    flow_cohoes = flow_cohoes_ft3pm * 0.0283168;
    stastat_cohoes = stastat;
    staid_cohoes = staid(1);

    % Susquehanna River that flows in the Chesapeake
    [usgc,staid,yyyy,mm,dd,flow_susquehanna_ft3pm_str,stastat] = textread([workdir 'USGS_01578310_Susquehanna_196710-201803.txt'],'%s%s%n-%n-%n%s%s','delimiter','\t ','commentstyle','shell','headerlines',31);
    niceind = find(~strcmp(flow_susquehanna_ft3pm_str, 'Ice') & ~strcmp(flow_susquehanna_ft3pm_str, ''));
    flow_susquehanna_ft3pm = repmat(NaN,size(flow_susquehanna_ft3pm_str));
    flow_susquehanna_ft3pm(niceind) = str2num(char(string(flow_susquehanna_ft3pm_str(niceind))));
    t_susquehanna = datenum(yyyy,mm,dd);
    flow_susquehanna = flow_susquehanna_ft3pm * 0.0283168;
    stastat_susquehanna = stastat;
    staid_susquehanna = staid(1);

    % Potomac River that flows in the Chesapeake
    [usgc,staid,yyyy,mm,dd,flow_potomac_ft3pm_str,stastat] = textread([workdir 'USGS_01646500_Potomac_193003-201803.txt'],'%s%s%n-%n-%n%s%s','delimiter','\t ','commentstyle','shell','headerlines',31);
    niceind = find(~strcmp(flow_potomac_ft3pm_str, 'Ice') & ~strcmp(flow_potomac_ft3pm_str, ''));
    flow_potomac_ft3pm = repmat(NaN,size(flow_potomac_ft3pm_str));
    flow_potomac_ft3pm(niceind) = str2num(char(string(flow_potomac_ft3pm_str(niceind))));
    t_potomac = datenum(yyyy,mm,dd);
    flow_potomac = flow_potomac_ft3pm * 0.0283168;
    stastat_potomac = stastat;
    staid_potomac = staid(1);

    % James River that flows in the Chesapeake
    [usgc,staid,yyyy,mm,dd,flow_james_ft3pm_str,stastat] = textread([workdir 'USGS_02037500_James_193410-201803.txt'],'%s%s%n-%n-%n%s%s','delimiter','\t ','commentstyle','shell','headerlines',31);
    niceind = find(~strcmp(flow_james_ft3pm_str, 'Ice') & ~strcmp(flow_james_ft3pm_str, ''));
    flow_james_ft3pm = repmat(NaN,size(flow_james_ft3pm_str));
    flow_james_ft3pm(niceind) = str2num(char(string(flow_james_ft3pm_str(niceind))));
    t_james = datenum(yyyy,mm,dd);
    flow_james = flow_james_ft3pm * 0.0283168;
    stastat_james = stastat;
    staid_james = staid(1);
end %if

t_river = [t1:1:t2]';
t_hudson = t_river;
t_chesapeake = t_river;

yyyy = str2num(datestr(t_river,10));
mm = str2num(datestr(t_river,5));
dd = str2num(datestr(t_river,7));
yd = datenum(yyyy,mm,dd) - datenum(yyyy,1,0);
yyyylist = unique(yyyy);

tind_ftedward = find(t_ftedward >= t1 & t_ftedward <= t2);
tind_cohoes = find(t_cohoes >= t1 & t_cohoes <= t2);
tind_susquehanna = find(t_susquehanna >= t1 & t_susquehanna <= t2);
tind_potomac = find(t_potomac >= t1 & t_potomac <= t2);
tind_james = find(t_james >= t1 & t_james <= t2);

% discharge of Hudson River (cubic meters/s): sum of flows at Fort Edwards and Mohawk River at Cohoes, total * 1.3 according to Choi and Wilkin 2007
flow_hudson = (flow_ftedward(tind_ftedward) + flow_cohoes(tind_cohoes))*1.3;
flow_hudson_JFMA = repmat(NaN,[length(yyyylist),1]);
flow_hudson_MJJA = repmat(NaN,[length(yyyylist),1]);
flow_hudson_SOND = repmat(NaN,[length(yyyylist),1]);

% discharge from the Chesapeake (cubic meters/s): sum of flows of Susquehanna, Potomac, and James, total * 1.25 according to REF???
flow_chesapeake = (flow_susquehanna(tind_susquehanna) + flow_potomac(tind_potomac) + flow_james(tind_james))*1.25;
flow_chesapeake_JFMA = repmat(NaN,[length(yyyylist),1]);
flow_chesapeake_MJJA = repmat(NaN,[length(yyyylist),1]);
flow_chesapeake_SOND = repmat(NaN,[length(yyyylist),1]);

flow_hudson_monthly = repmat(NaN,[length(yyyylist),12]);
flow_chesapeake_monthly = repmat(NaN,[length(yyyylist),12]);

for yi = 1:length(yyyylist)
%for yi = 13:13
    t1ind = find(yyyy == yyyylist(yi) & mm >= 1 & mm <= 4);
    t2ind = find(yyyy == yyyylist(yi) & mm >= 5 & mm <= 8);
    t3ind = find(yyyy == yyyylist(yi) & mm >= 9 & mm <= 12);

    flow_hudson_JFMA(yi) = nanmean(flow_hudson(t1ind));
    flow_hudson_MJJA(yi) = nanmean(flow_hudson(t2ind));
    flow_hudson_SOND(yi) = nanmean(flow_hudson(t3ind));

    flow_chesapeake_JFMA(yi) = nanmean(flow_chesapeake(t1ind));
    flow_chesapeake_MJJA(yi) = nanmean(flow_chesapeake(t2ind));
    flow_chesapeake_SOND(yi) = nanmean(flow_chesapeake(t3ind));

    for mi = 1:12
        tind = find(yyyy == yyyylist(yi) & mm == mi);
        flow_hudson_monthly(yi,mi) = nanmean(flow_hudson(tind));
        flow_chesapeake_monthly(yi,mi) = nanmean(flow_chesapeake(tind));
    end %for
end %for

flow_hudson_dailyclim = repmat(NaN,[365,1]);
flow_chesapeake_dailyclim = repmat(NaN,[365,1]);
flowstd_hudson_dailyclim = repmat(NaN,[365,1]);
flowstd_chesapeake_dailyclim = repmat(NaN,[365,1]);
yday = [1:365]';
for di = 1:length(yday)
    tind = find(yd == yday(di));
    flow_hudson_dailyclim(di) = nanmean(flow_hudson(tind));
    flow_chesapeake_dailyclim(di) = nanmean(flow_chesapeake(tind));
    flowstd_hudson_dailyclim(di) = nanstd(flow_hudson(tind));
    flowstd_chesapeake_dailyclim(di) = nanstd(flow_chesapeake(tind));
end %for

save([workdir 'riverdischarge.mat'],'t_river','yyyylist','yday','flow_hudson_JFMA','flow_hudson_MJJA','flow_hudson_SOND','flow_chesapeake_JFMA','flow_chesapeake_MJJA','flow_chesapeake_SOND','flow_hudson_monthly','flow_chesapeake_monthly','flow_hudson_dailyclim','flowstd_hudson_dailyclim','flow_chesapeake_dailyclim','flowstd_chesapeake_dailyclim','-v7.3');

%plot(yday,flow_chesapeake_dailyclim); hold on;
%plot(yday,flow_chesapeake_dailyclim+flowstd_chesapeake_dailyclim)
%plot(yday,flow_chesapeake_dailyclim-flowstd_chesapeake_dailyclim)
