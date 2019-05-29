% dg_mab_fig2.m
% this script makes figure 2 for the MAB manuscript showing how the interpolation scheme works.
% 
% DG 2018-06-10

sodaflag = 0
dg_setup_MABclimate_dir
if ~exist('region')
    load([workdir 'hydroMABgrid.mat']);
end %if

if ~exist('hydro')
    load([workdir 'hydroMAB2Da.mat']);
end %if

if ~exist('casts')
    load([workdir 'hydroMABcasts.mat']);
end %if

yi = 33 % 38 is 2014
si = 2 % 1 is winter-spring, 2 is summer
ri = 2 % NYB1

season(1).label = 'Jan-Apr';
season(2).label = 'May-Aug';
season(3).label = 'Sep-Dec';

lon = region(ri).longrid;
lat = region(ri).latgrid;
londata = [casts.lon]';
latdata = [casts.lat]';
depth = -double([casts.zdeep])';
ind = hydro(yi).ind{si};

lontrain = londata(hydro(yi).ind{si});
lattrain = latdata(hydro(yi).ind{si});
hdoi = casts(hydro(yi).ind{si}); % hydrographic data of interest
h = 0; %meters
searchradius = 30000;
wgs84 = wgs84Ellipsoid('meters');

[extrain,eytrain,eztrain] = geodetic2ecef(wgs84,lattrain,lontrain,h);
ECEFxyztrain = [extrain,eytrain,eztrain];
dModel = createns(ECEFxyztrain,'NSMethod','kdtree','Distance','Euclidean');
nnindxri = cell(length(region),1);
[ex,ey,ez] = geodetic2ecef(wgs84,region(ri).llat,region(ri).llon,h);
nnindxri{ri} = rangesearch(dModel, [ex,ey,ez], searchradius);

% figure 2a
f2a = figure('unit','inches')
set(gcf,'paperposition',[0 0 10 8]);

hpgrid = plot(lon,lat,'y.'); hold on;
hpcast = plot(londata(ind),latdata(ind),'k.');
set(hpcast,'markersize',15);
xlim([-74.25 -70.25]);
ylim([38.7 41.25]);
set(hpgrid,'markersize',15)

% find the minimum distance to an user specified lon/lat
%inputll0 = [-72.6817396313364 39.3053935860058];
inputll0 = [-72.8076036866359 39.7222303206997];
%inputll0 =  [-72.3375576036866 40.1013848396501];

inputll0 = ginput;

dist0 = distance(inputll0(2),inputll0(1),region(ri).llat,region(ri).llon);
[tmp,md0] = min(dist0);
lon0 = region(ri).llon(md0);
lat0 = region(ri).llat(md0);

% find all the grid points within the search radius;
dist0 = distance(lat0,lon0,region(ri).llat,region(ri).llon)*1852*60;
ind30 = find(dist0 <= searchradius);


%[tempdafit, saltdafit] = dg_idw_casts2d(region(ri).llon, region(ri).llat, region(ri).zz, hdoi, nnindx{ri},searchradius);
maxR = searchradius;
nnindx = nnindxri{ri};

% find CTD casts used to compute average for specific grid location
clon = [hdoi(nnindx{md0}).lon]';
clat = [hdoi(nnindx{md0}).lat]';

% compute the weight for each CTD station used in IDW weighting
distgrid = deg2km(distance([lat0,lon0],[clat,clon]))*1000; % distance in meters
tempind = find(~isnan([hdoi(nnindx{md0}).tempda]) == 1);
tempweight = (max(0,maxR - distgrid(tempind)) ./ (maxR*distgrid(tempind))).^2;
tempweightrel = tempweight ./ (sum(tempweight));

green = [0.8 0.8 0.8];
blue = [0 0 1];

% calculate the colormap for the IDW weight plot
ptempidw = repmat(NaN,[size(tempweightrel,1),3]);
for ii = 1:length(tempweightrel)
    ptempidw(ii,:) = tempweightrel(ii) .* blue + (1-tempweightrel(ii)) .* green;
end %for

% define a reference scale for colorbar plotting;
idwscl = [0:0.01:1];
prefidw = repmat(NaN,[size(idwscl,1),3]);
for ii = 1:length(idwscl)
    prefidw(ii,:) = idwscl(ii) .* blue + (1-idwscl(ii)) .* green;
end %for

%fcb = figure;
%close(fcb)

% plot the IDW CTD casts selected with corresponding weights
hpsearch = plot(region(ri).llon(ind30),region(ri).llat(ind30),'g.');
hs1 = scatter(clon,clat,30,ptempidw,'o','filled'); cmidw = colormap(prefidw); hcb = colorbar;
hs2 = plot(clon,clat,'bo');
hp0 = plot(lon0,lat0,'r*');

set(hs1,'sizedata',40)
set(hs2,'markersize',7)
set(hp0,'markersize',15)
set(hpsearch,'markersize',15)

hold off

hx = xlabel('Longitude');
hy = ylabel('Latitude');
ht = title(['Inverse Distance Weighting Scheme (' region(ri).label ' ' season(si).label ' ' num2str(hydro(yi).yyyy) ')'])
hl = legend('Bathymetric grid','CTD Stations','Search domain (radius=30km)','IDW weights (grey=0,blue=1)','IDW stations','Grid computed');
set(gca,'fontweight','bold','fontsize',16);
set(hx,'fontweight','bold','fontsize',18);
set(hy,'fontweight','bold','fontsize',18);
set(ht,'fontweight','bold','fontsize',18);

timenowstr = datestr(now,30);

eval(['print -depsc -r300 ' figoutdir 'MABclimate_fig2a_' timenowstr '.eps'])
close(f2a)

%figure 2b temp
ind1 = nnindx{md0}(1);

f2b = figure('unit','inches')
set(gcf,'paperposition',[0 0 8 8]);
hp1 = plot(hdoi(ind1).temp,hdoi(ind1).zz,'b.-'); hold on;
mtemp = repmat(nanmean(hdoi(ind1).temp),size(hdoi(ind1).zz));
hp2 = plot(mtemp,hdoi(ind1).zz,'b--');
set(gca,'box','on')
set(hp1,'markersize',20);
set(hp1,'linewidth',3);
set(hp2,'linewidth',3);

hx = xlabel('Temperature (C)');
hy = ylabel('Depth (m)');
%ht = title('Inverse Distance Weighting Scheme (NYB1 Jan-Apr 2015)')
hl = legend('CTD temperature profile','Depth-averaged temperature','location','southwest');
set(gca,'fontweight','bold','fontsize',16);
set(hx,'fontweight','bold','fontsize',18);
set(hy,'fontweight','bold','fontsize',18);
set(hl,'fontweight','bold','fontsize',13);
%xlim([9.5 12])
%set(ht,'fontweight','bold','fontsize',18);

eval(['print -depsc -r300 ' figoutdir 'MABclimate_fig2b_temp_' timenowstr '.eps'])
close(f2b)

f2bsalt = figure('unit','inches')
set(gcf,'paperposition',[0 0 8 8]);
hp1 = plot(hdoi(ind1).salt,hdoi(ind1).zz,'b.-'); hold on;
msalt = repmat(nanmean(hdoi(ind1).salt),size(hdoi(ind1).zz));
hp2 = plot(msalt,hdoi(ind1).zz,'b--');
set(gca,'box','on')
set(hp1,'markersize',20);
set(hp1,'linewidth',3);
set(hp2,'linewidth',3);

hx = xlabel('Salinity');
hy = ylabel('Depth (m)');
%ht = title(['Inverse Distance Weighting Scheme (' region(ri).label ' Jan-Apr 2015)'])
hl = legend('CTD salinity profile','Depth-averaged salinity','location','southwest');
set(gca,'fontweight','bold','fontsize',16);
set(hx,'fontweight','bold','fontsize',18);
set(hy,'fontweight','bold','fontsize',18);
set(hl,'fontweight','bold','fontsize',14);
%xlim([34.7 35.5])
%set(ht,'fontweight','bold','fontsize',18);

eval(['print -depsc -r300 ' figoutdir 'MABclimate_fig2b_salt_' timenowstr '.eps'])
close(f2bsalt)