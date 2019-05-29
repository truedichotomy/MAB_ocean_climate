% setup directory needed to run the script
dg_setup_MABclimate_dir

yyyy1 = 1875;
yyyy2 = 2007;

%yyyy1 = 1977;
%yyyy2 = 2016;


gissdatadir = '/Users/gong/Documents/Research/MAB_climate/GISStemp/';
%gissdatafile = 'GISS_temp_NHstations_J-D.txt';
%gissdatafile = 'GISS_temp_NHstations_DJF.txt';
%gissdatafile = 'GISS_temp_NHstations_MAM.txt';
%gissdatafile = 'GISS_temp_NHstations_JJA.txt';
%gissdatafile = 'GISS_temp_NHstations_SON.txt';
gissdatafile = 'GISS_temp_US_J-D.txt';

[yyyy,temp,tempsm] = textread([gissdatadir gissdatafile],'%n%n%n','delimiter',' \t','commentstyle','matlab');
s = polyfit(yyyy,temp,1)
s_smooth = polyfit(yyyy,tempsm,1)
%tind= find(yyyy <= 2016 & yyyy >= 1977);  
tind= find(yyyy <= yyyy2 & yyyy >= yyyy1);  
s_study = polyfit(yyyy(tind),temp(tind),1)         
s_study_smooth = polyfit(yyyy(tind),tempsm(tind),1)
s_fit = s_study;    
temp_calc = s_fit(1) .* yyyy(tind) + s_fit(2);

fgiss = figure(20)
set(fgiss,'unit','inches')
set(fgiss,'paperposition',[0 0 8 6]);
hold off
hpt = plot(yyyy,temp,'k*-'); hold on;
hptsm = plot(yyyy,tempsm,'r-');
hptfit = plot(yyyy(tind),temp_calc,'b-');
hl = legend('NH Stations Avg. Temperature','5 Year Lowess Smoothing',[num2str(yyyy(tind(1))) '-' num2str(yyyy(tind(end))) ' fit: ' num2str(round(s_fit(1)*10,3)) char(176) 'C/decade'],'location','northwest');
%hl = legend('U.S. Stations Avg. Temperature','5 Year Lowess Smoothing',[num2str(yyyy(tind(1))) '-' num2str(yyyy(tind(end))) ' fit: ' num2str(round(s_fit(1)*10,3)) char(176) 'C/decade'],'location','northwest');
hx = xlabel('Year');
hy = ylabel(['Temperature (' char(176) 'C)']);
ht = title(['NASA GISS Surface Temperature Analysis (' gissdatafile(end-6:end-4) ')']);
set(gca,'xgrid','on','ygrid','on');
set(hpt,'linewidth',1);
set(hptsm,'linewidth',2);
set(hptfit,'linewidth',2);
set(gca,'fontsize',15,'fontweight','bold');
set(hx,'fontsize',16,'fontweight','bold');
set(hy,'fontsize',16,'fontweight','bold');
set(ht,'fontsize',16,'fontweight','bold');

gissfigfile = [gissdatafile(1:end-4) '_' num2str(yyyy(tind(1))) '-' num2str(yyyy(tind(end)))];
eval(['print -depsc -r300 ' gissdatadir gissfigfile '.eps' ]);