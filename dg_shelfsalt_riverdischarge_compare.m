dg_setup_MABclimate_dir;
loadflag = 1

close all

workdirsalt = [rootdir 'MAB_river_salt_analysis/'];

if loadflag == 1
    load([workdir 'MABclimate2Da_nobiasDG.mat']);
    load([workdir 'riverdischarge.mat']);
end %if

% season (3): 1-Winter-Spring, 2-Spring-Summer, 3-Fall-Winter
% region (9): 1-SNE, 2-NYB1, 3-NYB2, 4-SS1, 5-SS2, 6-MAB, 7-GB, 8-ENE, 9-GOM
% variable (4): 1-temp, 2-tempshw, 3-salt, 4-saltshw
% time period (3): 1-1977-1999, 2-1999-2017, 3-1977-2016

% season selection: si = 1: winter/spring, si = 2: spring/summer, si = 3: fall/winter
si = 1
sir = 1 % seasonal index for river

switch sir
case 1,
    yriver = flow_chesapeake_JFMA;
case 2,
    yriver = flow_chesapeake_MJJA;
case 3,
    yriver = flow_chesapeake_SOND;
otherwise,
    yriver = flow_chesapeake_JFMA;
    display('assuming winter for river...')
end

% chesapeake discharge & SS2 salt comparison
ri = 5;
x = ccc.yyyy;
ysalt = ccc.salt(:,si,ri);

%yriver = flow_chesapeake_JFMA;
%yriver = flow_chesapeake_MJJA;

ysalta = (ysalt-nanmean(ysalt))/nanstd(ysalt);
yrivera = (yriver-nanmean(yriver))/nanstd(yriver);
ysalta_ss2 = ysalta;
yrivera_ss2 = yrivera;

nnind = find(~isnan(ysalta));
[r,p]= corrcoef(ysalta(nnind),yrivera(nnind));
chesapeake_ss2_corr = r(2,1)
chesapeake_ss2_p = p(2,1)

fp_ches = polyfit(x(nnind),yriver(nnind),1);
fp_ss2 = polyfit(x(nnind),ysalt(nnind),1);

% delta discharge / shelf vol = delta salt / shelf salt * x
D1fresh = fp_ches(1)*3600*24*120; % cubic meters
SS2vol = 20e9*50; % cubic meters
%SS2vol = vol
D1salt = fp_ss2(1);
SS2salt = 34;

x_ss2 = -(D1fresh/SS2vol) / (D1salt/SS2salt)



% hudson discharge & NYB1 salt comparison
ri = 2;
x = ccc.yyyy;
ysalt = ccc.salt(:,si,ri);

switch sir
case 1,
    yriverh = flow_hudson_JFMA;
case 2,
    yriverh = flow_hudson_MJJA;
case 3,
    yriverh = flow_hudson_SOND;
otherwise,
    yriverh = flow_hudson_JFMA;
    display('assuming winter for river...')
end

%yriverh = flow_hudson_JFMA;
%yriverh = flow_hudson_MJJA;

ysalta = (ysalt-nanmean(ysalt))/nanstd(ysalt);
yriverah = (yriverh-nanmean(yriverh))/nanstd(yriverh);
ysalta_nyb1 = ysalta;
yrivera_nyb1 = yriverah;

nnind = find(~isnan(ysalta));
[r,p]= corrcoef(ysalta(nnind),yriverah(nnind));
hudson_nyb1_corr = r(2,1)
hudson_nyb1_p = p(2,1)

fp_huds = polyfit(x(nnind),yriverh(nnind),1);
fp_nyb1 = polyfit(x(nnind),ysalt(nnind),1);

% delta_discharge / shelf_vol = delta_salt / shelf_salt * fraction_attributable_to_river
D1fresh = fp_huds(1)*3600*24*120; % cubic meters
NYB1vol = 30e9*50; % cubic meters
%NYB1vol = vol
D1salt = fp_nyb1(1);
NYB1salt = 34;

x_nyb1 = -(D1fresh/NYB1vol) / (D1salt/NYB1salt)

figure(2)
set(gcf,'unit','inches')
set(gcf,'paperposition',[0 0 10 10])

subplot(2,1,1)
hold on
hpsalt_nyb1 = plot(x,ysalta_nyb1,'k-');
hpriver_nyb1 = plot(x,yrivera_nyb1,'k:');
set(hpsalt_nyb1,'linewidth',2);
set(hpriver_nyb1,'linewidth',2)
ht_nyb1 = title(['NYB1 Salinity (' season(si).name ') vs. Hudson Discharge (' season(sir).name ') (1977-2016)']);
set(gca,'box','on','xgrid','on','ygrid','on');
hy = ylabel('Normalized Anomaly');
ylim([-3 3]);
set(gca,'fontsize',16,'fontweight','bold');
set(hy,'fontsize',16,'fontweight','bold');
hl_nyb1 = legend('NYB1 Salinity','Hudson Discharge');
htxt_nyb1 = annotation(gcf,'textbox',[0.67 0.59 0.2 0.04],'string',{['R=' num2str(round(hudson_nyb1_corr,2)) ', P=' num2str(round(hudson_nyb1_p,3))]});
set(htxt_nyb1,'fontsize',16,'fontweight','bold')
hold off

subplot(2,1,2)
hold on
hpsalt_ss2 = plot(x,ysalta_ss2,'k-');
hpriver_ss2 = plot(x,yrivera_ss2,'k:');
set(hpsalt_ss2,'linewidth',2);
set(hpriver_ss2,'linewidth',2);
ht_ss2 = title(['SS2 Salinity (' season(si).name ') vs. Chesapeake Discharge (' season(sir).name ') (1977-2016)']);
set(gca,'box','on','xgrid','on','ygrid','on');
hx = xlabel('Year');
hy = ylabel('Normlized Anomaly');
ylim([-3 3])
set(gca,'fontsize',16,'fontweight','bold');
hl_ss2 = legend('SS2 Salinity','Chesapeake Discharge');
htxt_nyb1 = annotation(gcf,'textbox',[0.67 0.12 0.2 0.04],'string',{['R=' num2str(round(chesapeake_ss2_corr,2)) ', P=' num2str(round(chesapeake_ss2_p,3))]});
set(htxt_nyb1,'fontsize',16,'fontweight','bold')
hold off

eval(['print -deps -r300 ' figoutdir 'river_salinity_r[' season(sir).name ']_s[' season(si).name '].eps'])


% River correlation

f = figure;
set(gcf,'unit','inches')
set(gcf,'paperposition',[0 0 10 10])

for sir = 1:3

switch sir
case 1,
    yriverh = flow_hudson_JFMA;
    yriver = flow_chesapeake_JFMA;
case 2,
    yriverh = flow_hudson_MJJA;
    yriver = flow_chesapeake_MJJA;
case 3,
    yriverh = flow_hudson_SOND;
    yriver = flow_chesapeake_SOND;
otherwise,
    yriverh = flow_hudson_JFMA;
    yriver = flow_chesapeake_JFMA;
    display('assuming winter for river...')
end

subplot(3,1,sir)

yrivera_nyb1 = (yriverh-nanmean(yriverh))/nanstd(yriverh);
yrivera_ss2 = (yriver-nanmean(yriver))/nanstd(yriver);

hold on;
hp_0 = plot([1960 2030],[0 0],'k.');
hp_nyb1 = plot(x,yrivera_nyb1);
hp_ss2 = plot(x,yrivera_ss2);
hl_river = legend(['R = ' num2str(hudson_chesapeake_corr,3)],'location','southeast');
set(hp_nyb1,'linewidth',2);
set(hp_ss2,'linewidth',2)
[rr,pr]= corrcoef(yrivera_nyb1,yrivera_ss2);
hudson_chesapeake_corr = rr(2,1)
hudson_chesapeake_p = pr(2,1)
ht_hudson_chesapeake = title(['Hudson vs. Chesapeake (' season(sir).name ') (1977-2016)']);
set(gca,'box','on','xgrid','on','ygrid','on');
%hx = xlabel('Year');
hy = ylabel('Normlized Anomaly');
xlim([1975 2020])
ylim([-2.5 2.5])
set(gca,'fontsize',16,'fontweight','bold');
hold off
end %for
eval(['print -depsc -r300 ' figoutdir 'hudson_chesapeake_correlation.eps'])

close(f)


fMAB = figure;
set(gcf,'unit','inches')
set(gcf,'paperposition',[0 0 10 10])

si = 3

subplot(2,1,1)
hold on;
hp1 = plot(x,ccc.temp_nobias(:,si,1)); % SNE
hp2 = plot(x,ccc.temp_nobias(:,si,2)); % NYB1
hp3 = plot(x,ccc.temp_nobias(:,si,3)); % NYB2
hp4 = plot(x,ccc.temp_nobias(:,si,4)); % SS1
hp5 = plot(x,ccc.temp_nobias(:,si,5)); % SS2

nnan12 = find(~isnan(ccc.temp_nobias(:,si,1)) & ~isnan(ccc.temp_nobias(:,si,2)));
[rt12,pt12]= corrcoef(ccc.temp_nobias(nnan12,si,1),ccc.temp_nobias(nnan12,si,2));

nnan23 = find(~isnan(ccc.temp_nobias(:,si,2)) & ~isnan(ccc.temp_nobias(:,si,3)));
[rt23,pt23]= corrcoef(ccc.temp_nobias(nnan23,si,2),ccc.temp_nobias(nnan23,si,3));

nnan34 = find(~isnan(ccc.temp_nobias(:,si,3)) & ~isnan(ccc.temp_nobias(:,si,4)));
[rt34,pt34]= corrcoef(ccc.temp_nobias(nnan34,si,3),ccc.temp_nobias(nnan34,si,4));

nnan45 = find(~isnan(ccc.temp_nobias(:,si,4)) & ~isnan(ccc.temp_nobias(:,si,5)));
[rt45,pt45]= corrcoef(ccc.temp_nobias(nnan45,si,4),ccc.temp_nobias(nnan45,si,5));

set(hp1,'linewidth',2);
set(hp2,'linewidth',2);
set(hp3,'linewidth',2);
set(hp4,'linewidth',2);
set(hp5,'linewidth',2);

xlim([1975 2030]);
if si == 1
    ylim([6 18]);
elseif si == 2
    ylim([8 18]);
elseif si == 3
    ylim([11 20])
end %if


ht = title(['Temperature for ' season(si).name])
hx = xlabel('Year')
hy = ylabel('Temperature')
set(gca,'fontsize',16,'fontweight','bold','box','on','xgrid','on','ygrid','on');
set(hy,'fontsize',16,'fontweight','bold');
hl = legend('SNE',['NYB1 R=' num2str(rt12(2,1),2)],['NYB2 R=' num2str(rt23(2,1),2)],['SS1 R=' num2str(rt34(2,1),2)],['SS2 R=' num2str(rt45(2,1),2)],'location','east');


subplot(2,1,2)
hold on;
hp1 = plot(x,ccc.salt_nobias(:,si,1)); % SNE
hp2 = plot(x,ccc.salt_nobias(:,si,2)); % NYB1
hp3 = plot(x,ccc.salt_nobias(:,si,3)); % NYB2
hp4 = plot(x,ccc.salt_nobias(:,si,4)); % SS1
hp5 = plot(x,ccc.salt_nobias(:,si,5)); % SS2

nnan12 = find(~isnan(ccc.salt_nobias(:,si,1)) & ~isnan(ccc.salt_nobias(:,si,2)));
[rs12,ps12]= corrcoef(ccc.salt_nobias(nnan12,si,1),ccc.salt_nobias(nnan12,si,2));

nnan23 = find(~isnan(ccc.salt_nobias(:,si,2)) & ~isnan(ccc.salt_nobias(:,si,3)));
[rs23,ps23]= corrcoef(ccc.salt_nobias(nnan23,si,2),ccc.salt_nobias(nnan23,si,3));

nnan34 = find(~isnan(ccc.salt_nobias(:,si,3)) & ~isnan(ccc.salt_nobias(:,si,4)));
[rs34,ps34]= corrcoef(ccc.salt_nobias(nnan34,si,3),ccc.salt_nobias(nnan34,si,4));

nnan45 = find(~isnan(ccc.salt_nobias(:,si,4)) & ~isnan(ccc.salt_nobias(:,si,5)));
[rs45,ps45]= corrcoef(ccc.salt_nobias(nnan45,si,4),ccc.salt_nobias(nnan45,si,5));

set(hp1,'linewidth',2);
set(hp2,'linewidth',2);
set(hp3,'linewidth',2);
set(hp4,'linewidth',2);
set(hp5,'linewidth',2);

xlim([1975 2030]);
if si == 1
    ylim([32 36])
elseif si == 2
    ylim([31.5 35.5])
elseif si == 3
    ylim([32 35.5])
end %if
ht = title(['Salinity for ' season(si).name])
hx = xlabel('Year')
hy = ylabel('Salinity')
set(gca,'fontsize',16,'fontweight','bold','box','on','xgrid','on','ygrid','on');
set(hy,'fontsize',16,'fontweight','bold');
%hl = legend('SNE','NYB1','NYB2','SS1','SS2');
hl = legend('SNE',['NYB1 R=' num2str(rs12(2,1),2)],['NYB2 R=' num2str(rs23(2,1),2)],['SS1 R=' num2str(rs34(2,1),2)],['SS2 R=' num2str(rs45(2,1),2)],'location','east');

eval(['print -depsc -r300 ' figoutdir 'MAB_regional_' season(si).name '_correlation.eps'])
close(fMAB);