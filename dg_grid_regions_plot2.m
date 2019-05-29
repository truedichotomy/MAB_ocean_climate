% dg_grid_region_plot: this script plots the calculated temperature and salinity for each region of each season in the MAB, GOM and GB domain
%
% NOTE: should run dg_grid_regions_analysis.m first if the ccc regression data are not yet saved in MABclimate.mat.
%
% DG 2017-11-01

loadflag = 0
displayflag = 0
printflag = 1
plotmonthlyflag = 0
seasonbiasflag = 1
sodaflag = 0

psig = 3;

% setup directories needed for the script
dg_setup_MABclimate_dir

if loadflag == 1
    if sodaflag == 0
        load([workdir 'MABclimate2Da_nobiasDG.mat']);
    elseif sodaflag == 1
        load([workdir 'MABclimate2DaSODA_nobiasDG.mat']);
    end %if
end %if

ns = 3;
nr = 9;
nv = 4;
nt = 3;

yyyy1 = yyyy(1:23);
yyyy2 = yyyy(23:end);

tempi = 1;
tempshwi = 2;
salti = 3;
saltshwi = 4;
p1 = 1;
p2 = 2;
p3 = 3;

if displayflag == 1

    si = 3
    ri = 6

    if seasonalbiasflag == 0
        tempfilter = ccc.tempfilt(:,si,ri);
        saltfilter = ccc.saltfilt(:,si,ri);
        tempshwfilter = ccc.tempshwfilt(:,si,ri);
        saltshwfilter = ccc.saltshwfilt(:,si,ri);

        tempavg = nanmean(tempfilter);
        tempstd = nanstd(tempfilter);
        saltavg = nanmean(saltfilter);
        saltstd = nanstd(saltfilter);
        tempshwavg = nanmean(tempshwfilter);
        tempshwstd = nanstd(tempshwfilter);
        saltshwavg = nanmean(saltshwfilter);
        saltshwstd = nanstd(saltshwfilter);
    elseif seasonalbiasflag == 1
        tempfilter = ccc.tempfilt_nobias(:,si,ri);
        saltfilter = ccc.saltfilt_nobias(:,si,ri);
        tempshwfilter = ccc.tempshwfilt_nobias(:,si,ri);
        saltshwfilter = ccc.saltshwfilt_nobias(:,si,ri);

        tempavg = nanmean(tempfilter);
        tempstd = nanstd(tempfilter);
        saltavg = nanmean(saltfilter);
        saltstd = nanstd(saltfilter);
        tempshwavg = nanmean(tempshwfilter);
        tempshwstd = nanstd(tempshwfilter);
        saltshwavg = nanmean(saltshwfilter);
        saltshwstd = nanstd(saltshwfilter);
    end %if
    [region(ri).label ' ' season(si).label]

    figure(1)
    subplot(2,1,1)
    plot(yyyy,ccc.temp(:,si,ri),'ko-');
    %plot(yyyy,temp(:,si,ri),'ko-',yyyy,tempshw(:,si,ri),'bo-');
    ht1 = title([region(ri).label ' Temperature for ' season(si).label]);
    ylim([round(tempavg-3.5*tempstd) round(tempavg+3.5*tempstd)]);
    set(gca,'xgrid','on','ygrid','on')
    subplot(2,1,2)
    plot(yyyy,ccc.salt(:,si,ri),'ko-');
    %plot(yyyy,salt(:,si,ri),'ko-',yyyy,saltshw(:,si,ri),'bo-');
    ht2 = title([region(ri).label ' Salinity for ' season(si).label]);
    ylim([round(saltavg-psig*saltstd) round(saltavg+psig*saltstd)]);
    set(gca,'xgrid','on','ygrid','on')

    figure(2)
    subplot(2,1,1)
    %plot(yyyy,temp(:,si,ri),'ko-');
    plot(yyyy,ccc.temp(:,si,ri),'ko-',yyyy,ccc.tempshw(:,si,ri),'bo-');
    ht1 = title([region(ri).label ' Temperature for ' season(si).label]);
    ylim([round(tempshwavg-psig*tempshwstd) round(tempavg+psig*tempstd)]);
    set(gca,'xgrid','on','ygrid','on')
    subplot(2,1,2)
    %plot(yyyy,salt(:,si,ri),'ko-');
    plot(yyyy,ccc.salt(:,si,ri),'ko-',yyyy,ccc.saltshw(:,si,ri),'bo-');
    ht2 = title([region(ri).label ' Salinity for ' season(si).label]);
    ylim([round(saltshwavg-psig*saltshwstd) round(saltavg+psig*saltstd)]);
    set(gca,'xgrid','on','ygrid','on')

    figure(3)
    subplot(2,1,1)
    %plot(yyyy,voltemp(:,si,ri),'ko-');
    plot(yyyy,ccc.voltemp(:,si,ri),'ko-',yyyy,ccc.volshw(:,si,ri),'bo-');
    ht1 = title([region(ri).label ' volume temperature for ' season(si).label]);
    subplot(2,1,2)
    %plot(yyyy,volsalt(:,si,ri),'ko-');
    plot(yyyy,ccc.volsalt(:,si,ri),'ko-',yyyy,ccc.volshw(:,si,ri),'bo-');
    ht2 = title([region(ri).label ' volume salinity for ' season(si).label]);
end %if

if printflag == 1
    for si = 1:ns
        for ri = 1:nr
            [region(ri).label ' ' season(si).label]

            tempfilter = ccc.tempfilt(:,si,ri);
            saltfilter = ccc.saltfilt(:,si,ri);
            tempshwfilter = ccc.tempshwfilt(:,si,ri);
            saltshwfilter = ccc.saltshwfilt(:,si,ri);

            tempavg = nanmean(tempfilter);
            tempstd = nanstd(tempfilter);
            saltavg = nanmean(saltfilter);
            saltstd = nanstd(saltfilter);
            tempshwavg = nanmean(tempshwfilter);
            tempshwstd = nanstd(tempshwfilter);
            saltshwavg = nanmean(saltshwfilter);
            saltshwstd = nanstd(saltshwfilter);

            mtemp = [ccc.b(si,ri,tempi,p3) ccc.m(si,ri,tempi,p3)];
            mtemp1 = [ccc.b(si,ri,tempi,p1) ccc.m(si,ri,tempi,p1)];
            mtemp2 = [ccc.b(si,ri,tempi,p2) ccc.m(si,ri,tempi,p2)];

            mtempshw = [ccc.b(si,ri,tempshwi,p3) ccc.m(si,ri,tempshwi,p3)];
            mtempshw1 = [ccc.b(si,ri,tempshwi,p1) ccc.m(si,ri,tempshwi,p1)];
            mtempshw2 = [ccc.b(si,ri,tempshwi,p2) ccc.m(si,ri,tempshwi,p2)];

            msalt = [ccc.b(si,ri,salti,p3) ccc.m(si,ri,salti,p3)];
            msalt1 = [ccc.b(si,ri,salti,p1) ccc.m(si,ri,salti,p1)];
            msalt2 = [ccc.b(si,ri,salti,p2) ccc.m(si,ri,salti,p2)];

            msaltshw = [ccc.b(si,ri,saltshwi,p3) ccc.m(si,ri,saltshwi,p3)];
            msaltshw1 = [ccc.b(si,ri,saltshwi,p1) ccc.m(si,ri,saltshwi,p1)];
            msaltshw2 = [ccc.b(si,ri,saltshwi,p2) ccc.m(si,ri,saltshwi,p2)];

            temphat = mtemp(1) + yyyy.*mtemp(2);
            temp1hat = mtemp1(1) + yyyy1.*mtemp1(2);
            temp2hat = mtemp2(1) + yyyy2.*mtemp2(2);

            tempshwhat = mtempshw(1) + yyyy.*mtempshw(2);
            tempshw1hat = mtempshw1(1) + yyyy1.*mtempshw1(2);
            tempshw2hat = mtempshw2(1) + yyyy2.*mtempshw2(2);

            salthat = msalt(1) + yyyy.*msalt(2);
            salt1hat = msalt1(1) + yyyy1.*msalt1(2);
            salt2hat = msalt2(1) + yyyy2.*msalt2(2);

            saltshwhat = msaltshw(1) + yyyy.*msaltshw(2);
            saltshw1hat = msaltshw1(1) + yyyy1.*msaltshw1(2);
            saltshw2hat = msaltshw2(1) + yyyy2.*msaltshw2(2);


            f1 = figure(1);
            set(gcf,'visible','off','unit','inches')
            set(gcf,'paperposition',[0 0 6 7.5])

            subplot(2,1,1)
            symcolor = 'kkk';
            sigind = find(ccc.hypothesis(si,ri,tempi,:)==1);
            if ~isempty(sigind)
                symcolor(sigind) = 'r';
            end %if

            switch seasonbiasflag
            case 0
                plot(yyyy,ccc.temp(:,si,ri),'ko-',yyyy,temphat,[symcolor(3) '--'], yyyy1,temp1hat,[symcolor(1) '-'],yyyy2,temp2hat,[symcolor(2) '-']);
            case 1
                plot(yyyy,ccc.temp_nobias(:,si,ri),'ko-',yyyy,temphat,[symcolor(3) '--'], yyyy1,temp1hat,[symcolor(1) '-'],yyyy2,temp2hat,[symcolor(2) '-']);
            end %switch
            %plot(yyyy,temp(:,si,ri),'ko-',yyyy,tempshw(:,si,ri),'bo-');
            ht1 = title([region(ri).label ' Temperature for ' season(si).label]);
            ylim([round(tempavg-psig*tempstd) round(tempavg+psig*tempstd)]);
            set(gca,'xgrid','on','ygrid','on')
            legend('temp', [num2str(round(ccc.m(si,ri,tempi,p3),2,'significant')) '/yr'], [num2str(round(ccc.m(si,ri,tempi,p1),2,'significant')) '/yr'], [num2str(round(ccc.m(si,ri,tempi,p2),2,'significant')) '/yr'], 'location','southeast')

            subplot(2,1,2)
            symcolor = 'kkk';
            sigind = find(ccc.hypothesis(si,ri,salti,:)==1);
            if ~isempty(sigind)
                symcolor(sigind) = 'r';
            end %if

            switch seasonbiasflag
            case 0
                plot(yyyy,ccc.salt(:,si,ri),'ko-',yyyy,salthat,[symcolor(3) '--'], yyyy1,salt1hat,[symcolor(1) '-'], yyyy2,salt2hat,[symcolor(2) '-']);
            case 1
                plot(yyyy,ccc.salt_nobias(:,si,ri),'ko-',yyyy,salthat,[symcolor(3) '--'], yyyy1,salt1hat,[symcolor(1) '-'], yyyy2,salt2hat,[symcolor(2) '-']);
            end %if

            %plot(yyyy,salt(:,si,ri),'ko-',yyyy,saltshw(:,si,ri),'bo-');
            ht2 = title([region(ri).label ' Salinity for ' season(si).label]);
            ylim([floor(saltavg-psig*saltstd) ceil(saltavg+psig*saltstd)]);
            set(gca,'xgrid','on','ygrid','on')
            legend('salt', [num2str(round(ccc.m(si,ri,salti,p3),2,'significant')) '/yr'], [num2str(round(ccc.m(si,ri,salti,p1),2,'significant')) '/yr'], [num2str(round(ccc.m(si,ri,salti,p2),2,'significant')) '/yr'], 'location','southeast')
            eval(['print -depsc -r300 ' figoutdir 'tempsalt_' region(ri).label '_' season(si).name '.eps']);
            close(f1)


            f2 = figure(2);
            set(gcf,'visible','off','unit','inches')
            set(gcf,'paperposition',[0 0 6 7.5])

            subplot(2,1,1)
            symcolor = 'kkk';
            sigind = find(ccc.hypothesis(si,ri,tempshwi,:)==1);
            if ~isempty(sigind)
                symcolor(sigind) = 'r';
            end %if

            switch seasonbiasflag
            case 0
                plot(yyyy,ccc.tempshw(:,si,ri),'ko-',yyyy,tempshwhat,[symcolor(3) '--'], yyyy1,tempshw1hat,[symcolor(1) '-'],yyyy2,tempshw2hat,[symcolor(2) '-']);
            case 1
                plot(yyyy,ccc.tempshw_nobias(:,si,ri),'ko-',yyyy,tempshwhat,[symcolor(3) '--'], yyyy1,tempshw1hat,[symcolor(1) '-'],yyyy2,tempshw2hat,[symcolor(2) '-']);
            end %switch
            %plot(yyyy,temp(:,si,ri),'ko-',yyyy,tempshw(:,si,ri),'bo-');
            ht1 = title([region(ri).label ' SHW Temperature for ' season(si).label]);
            ylim([round(tempshwavg-psig*tempshwstd) round(tempshwavg+psig*tempshwstd)]);
            set(gca,'xgrid','on','ygrid','on')
            legend('SHW temp', [num2str(round(ccc.m(si,ri,tempshwi,p3),2,'significant')) '/yr'], [num2str(round(ccc.m(si,ri,tempshwi,p1),2,'significant')) '/yr'], [num2str(round(ccc.m(si,ri,tempshwi,p2),2,'significant')) '/yr'], 'location','southeast')

            subplot(2,1,2)
            symcolor = 'kkk';
            sigind = find(ccc.hypothesis(si,ri,saltshwi,:)==1);
            if ~isempty(sigind)
                symcolor(sigind) = 'r';
            end %if

            switch seasonbiasflag
            case 0
                plot(yyyy,ccc.saltshw(:,si,ri),'ko-',yyyy,saltshwhat,[symcolor(3) '--'], yyyy1,saltshw1hat,[symcolor(1) '-'], yyyy2,saltshw2hat,[symcolor(2) '-']);
            case 1
                plot(yyyy,ccc.saltshw_nobias(:,si,ri),'ko-',yyyy,saltshwhat,[symcolor(3) '--'], yyyy1,saltshw1hat,[symcolor(1) '-'], yyyy2,saltshw2hat,[symcolor(2) '-']);
            end %switch

            %plot(yyyy,salt(:,si,ri),'ko-',yyyy,saltshw(:,si,ri),'bo-');
            ht2 = title([region(ri).label ' SHW Salinity for ' season(si).label]);
            ylim([floor(saltshwavg-psig*saltshwstd) ceil(saltshwavg+psig*saltshwstd)]);
            set(gca,'xgrid','on','ygrid','on')
            legend('SHW salt', [num2str(round(ccc.m(si,ri,saltshwi,p3),2,'significant')) '/yr'], [num2str(round(ccc.m(si,ri,saltshwi,p1),2,'significant')) '/yr'], [num2str(round(ccc.m(si,ri,saltshwi,p2),2,'significant')) '/yr'], 'location','southeast')

            eval(['print -depsc -r300 ' figoutdir 'tempsaltshw_' region(ri).label '_' season(si).name '.eps']);
            close(f2);


            f3 = figure(3);
            set(gcf,'visible','off')

            plot(yyyy,ccc.voltemp(:,si,ri),'bo-',yyyy,ccc.volsalt(:,si,ri),'ro-',yyyy,ccc.volshw(:,si,ri),'ko-');
            ht1 = title([region(ri).label ' volume temperature for ' season(si).label]);
            set(gca,'xgrid','on','ygrid','on')
            legend('vol temp','vol salt','vol shw', 'location','southeast')

            eval(['print -depsc -r300 ' figoutdir 'voltempsalt_' region(ri).label '_' season(si).name '.eps']);
            close(f3)
        end %for
    end %for
elseif printflag == 2
    gg(si,ri).title = [];

    gg(si,ri).tempfilter = [];
    gg(si,ri).saltfilter = [];
    gg(si,ri).tempshwfilter = [];
    gg(si,ri).saltshwfilter = [];

    gg(si,ri).tempavg = [];
    gg(si,ri).tempstd = [];
    gg(si,ri).saltavg = [];
    gg(si,ri).saltstd = [];
    gg(si,ri).tempshwavg = [];
    gg(si,ri).tempshwstd = [];
    gg(si,ri).saltshwavg = [];
    gg(si,ri).saltshwstd = [];

    gg(si,ri).mtemp = [];
    gg(si,ri).mtemp1 = [];
    gg(si,ri).mtemp2 = [];

    gg(si,ri).mtempshw = [];
    gg(si,ri).mtempshw1 = [];
    gg(si,ri).mtempshw2 = [];

    gg(si,ri).msalt = [];
    gg(si,ri).msalt1 = [];
    gg(si,ri).msalt2 = [];

    gg(si,ri).msaltshw = [];
    gg(si,ri).msaltshw1 = [];
    gg(si,ri).msaltshw2 = [];

    gg(si,ri).temphat = [];
    gg(si,ri).temp1hat = [];
    gg(si,ri).temp2hat = [];

    gg(si,ri).tempshwhat = [];
    gg(si,ri).tempshw1hat = [];
    gg(si,ri).tempshw2hat = [];

    gg(si,ri).salthat = [];
    gg(si,ri).salt1hat = [];
    gg(si,ri).salt2hat = [];

    gg(si,ri).saltshwhat = [];
    gg(si,ri).saltshw1hat = [];
    gg(si,ri).saltshw2hat = [];
 
    for si = 1:ns
        for ri = 1:nr
            gg(si,ri).title = [region(ri).label ' ' season(si).label];

            gg(si,ri).tempfilter = ccc.tempfilt(:,si,ri);
            gg(si,ri).saltfilter = ccc.saltfilt(:,si,ri);
            gg(si,ri).tempshwfilter = ccc.tempshwfilt(:,si,ri);
            gg(si,ri).saltshwfilter = ccc.saltshwfilt(:,si,ri);

            gg(si,ri).tempavg = nanmean(tempfilter);
            gg(si,ri).tempstd = nanstd(tempfilter);
            gg(si,ri).saltavg = nanmean(saltfilter);
            gg(si,ri).saltstd = nanstd(saltfilter);
            gg(si,ri).tempshwavg = nanmean(tempshwfilter);
            gg(si,ri).tempshwstd = nanstd(tempshwfilter);
            gg(si,ri).saltshwavg = nanmean(saltshwfilter);
            gg(si,ri).saltshwstd = nanstd(saltshwfilter);

            gg(si,ri).mtemp = [ccc.b(si,ri,tempi,p3) ccc.m(si,ri,tempi,p3)];
            gg(si,ri).mtemp1 = [ccc.b(si,ri,tempi,p1) ccc.m(si,ri,tempi,p1)];
            gg(si,ri).mtemp2 = [ccc.b(si,ri,tempi,p2) ccc.m(si,ri,tempi,p2)];

            gg(si,ri).mtempshw = [ccc.b(si,ri,tempshwi,p3) ccc.m(si,ri,tempshwi,p3)];
            gg(si,ri).mtempshw1 = [ccc.b(si,ri,tempshwi,p1) ccc.m(si,ri,tempshwi,p1)];
            gg(si,ri).mtempshw2 = [ccc.b(si,ri,tempshwi,p2) ccc.m(si,ri,tempshwi,p2)];

            gg(si,ri).msalt = [ccc.b(si,ri,salti,p3) ccc.m(si,ri,salti,p3)];
            gg(si,ri).msalt1 = [ccc.b(si,ri,salti,p1) ccc.m(si,ri,salti,p1)];
            gg(si,ri).msalt2 = [ccc.b(si,ri,salti,p2) ccc.m(si,ri,salti,p2)];

            gg(si,ri).msaltshw = [ccc.b(si,ri,saltshwi,p3) ccc.m(si,ri,saltshwi,p3)];
            gg(si,ri).msaltshw1 = [ccc.b(si,ri,saltshwi,p1) ccc.m(si,ri,saltshwi,p1)];
            gg(si,ri).msaltshw2 = [ccc.b(si,ri,saltshwi,p2) ccc.m(si,ri,saltshwi,p2)];

            gg(si,ri).temphat = mtemp(1) + yyyy.*mtemp(2);
            gg(si,ri).temp1hat = mtemp1(1) + yyyy1.*mtemp1(2);
            gg(si,ri).temp2hat = mtemp2(1) + yyyy2.*mtemp2(2);

            gg(si,ri).tempshwhat = mtempshw(1) + yyyy.*mtempshw(2);
            gg(si,ri).tempshw1hat = mtempshw1(1) + yyyy1.*mtempshw1(2);
            gg(si,ri).tempshw2hat = mtempshw2(1) + yyyy2.*mtempshw2(2);

            gg(si,ri).salthat = msalt(1) + yyyy.*msalt(2);
            gg(si,ri).salt1hat = msalt1(1) + yyyy1.*msalt1(2);
            gg(si,ri).salt2hat = msalt2(1) + yyyy2.*msalt2(2);

            gg(si,ri).saltshwhat = msaltshw(1) + yyyy.*msaltshw(2);
            gg(si,ri).saltshw1hat = msaltshw1(1) + yyyy1.*msaltshw1(2);
            gg(si,ri).saltshw2hat = msaltshw2(1) + yyyy2.*msaltshw2(2);
        end %for
    end %for

    f1 = figure(1);
    set(gcf,'visible','off')

    subplot(3,2,1)
    ri = 2 %NYB1
    si = 1 %winter-spring
    symcolor = 'kkk';
    sigind = find(ccc.hypothesis(si,ri,tempi,:)==1);
    if ~isempty(sigind)
        symcolor(sigind) = 'r';
    end %if

    switch seasonbiasflag
    case 0
        plot(yyyy,ccc.temp(:,si,ri),'ko-',yyyy,gg(si,ri).temphat,[symcolor(3) '--'], yyyy1,gg(si,ri).temp1hat,[symcolor(1) '-'],yyyy2,gg(si,ri).temp2hat,[symcolor(2) '-']);
    case 1
        plot(yyyy,ccc.temp_nobias(:,si,ri),'ko-',yyyy,gg(si,ri).temphat,[symcolor(3) '--'], yyyy1,gg(si,ri).temp1hat,[symcolor(1) '-'],yyyy2,gg(si,ri).temp2hat,[symcolor(2) '-']);
    end %switch
    %plot(yyyy,temp(:,si,ri),'ko-',yyyy,tempshw(:,si,ri),'bo-');
    ht1 = title([region(ri).label ' Temperature for ' season(si).label]);
    ylim([round(gg(si,ri).tempavg-psig*gg(si,ri).tempstd) round(gg(si,ri).tempavg+psig*gg(si,ri).tempstd)]);
    set(gca,'xgrid','on','ygrid','on')
    %legend('temp', [num2str(round(ccc.m(si,ri,tempi,p3),2,'significant')) '/yr'], [num2str(round(ccc.m(si,ri,tempi,p1),2,'significant')) '/yr'], [num2str(round(ccc.m(si,ri,tempi,p2),2,'significant')) '/yr'], 'location','southeastoutside')


    subplot(3,2,4)
    symcolor = 'kkk';
    sigind = find(ccc.hypothesis(si,ri,salti,:)==1);
    if ~isempty(sigind)
        symcolor(sigind) = 'r';
    end %if

    switch seasonbiasflag
    case 0
        plot(yyyy,ccc.salt(:,si,ri),'ko-',yyyy,gg(si,ri).salthat,[symcolor(3) '--'], yyyy1,gg(si,ri).salt1hat,[symcolor(1) '-'], yyyy2,gg(si,ri).salt2hat,[symcolor(2) '-']);
    case 1
        plot(yyyy,ccc.salt_nobias(:,si,ri),'ko-',yyyy,gg(si,ri).salthat,[symcolor(3) '--'], yyyy1,gg(si,ri).salt1hat,[symcolor(1) '-'], yyyy2,gg(si,ri).salt2hat,[symcolor(2) '-']);
    end %if

    %plot(yyyy,salt(:,si,ri),'ko-',yyyy,saltshw(:,si,ri),'bo-');
    ht2 = title([region(ri).label ' Salinity for ' season(si).label]);
    ylim([floor(gg(si,ri).saltavg-psig*gg(si,ri).saltstd) ceil(gg(si,ri).saltavg+psig*gg(si,ri).saltstd)]);
    set(gca,'xgrid','on','ygrid','on')
    %legend('salt', [num2str(round(ccc.m(si,ri,salti,p3),2,'significant')) '/yr'], [num2str(round(ccc.m(si,ri,salti,p1),2,'significant')) '/yr'], [num2str(round(ccc.m(si,ri,salti,p2),2,'significant')) '/yr'], 'location','southeastoutside')
    eval(['print -depsc -r300 ' figoutdir 'tempsalt_6panel.eps']);
    close(f1)

end %if

if plotmonthlyflag == 1
    fmonthly = figure(1);
    hp = [];
    set(gcf,'visible','on','unit','inches')
    set(gcf,'paperposition',[0 0 10 7]);
    hold on

    pmonth = datenum(0, ccc.month, 15);
    pregions = [9 7 8 1 2 3 4 5 6]; % north to south, then MAB
    for ri = 1:length(pregions)
        hp{ri} = plot(pmonth, ccc.tempshwmonthly(:,pregions(ri)));
        set(hp{ri},'color',[ri*0.1 0.3 1-ri*0.1],'linewidth',1);
    end %for
    set(gca,'box','on','xgrid','on','ygrid','on')
    datetick('x',5);
    set(hp{9},'color',[0 0 0],'linewidth',2)
    ht = title('Monthly Temperature Climatology')
    hx = xlabel('Month');
    hy = ylabel('Temperature');
    set(gca,'fontsize',18,'fontweight','bold');
    set(hx,'fontsize',18,'fontweight','bold');
    set(hy,'fontsize',18,'fontweight','bold');
    set(ht,'fontsize',20,'fontweight','bold');

    legend('GOM','GB','ENE','SNE','NYB1','NYB2','SS1','SS2','MAB','location','northwest');
    eval(['print -dpng -r300 ' figoutdir 'MAB_MonthlyTemp.png']);

    fmonthly = figure(2);
    set(gcf,'visible','on','unit','inches')
    set(gcf,'paperposition',[0 0 10 7]);
    hold on

    pmonth = datenum(0, ccc.month, 15);
    pregions = [9 7 8 1 2 3 4 5 6]; % north to south, then MAB
    for ri = 1:length(pregions)
        hp{ri} = plot(pmonth, ccc.saltmonthly(:,pregions(ri)));
        set(hp{ri},'color',[ri*0.1 0.3 1-ri*0.1],'linewidth',1);
    end %for
    set(gca,'box','on','xgrid','on','ygrid','on')
    datetick('x',5);
    set(hp{9},'color',[0 0 0],'linewidth',2)
    ht = title('Monthly Salinity Climatology')
    hx = xlabel('Month');
    hy = ylabel('Salinity');
    set(gca,'fontsize',18,'fontweight','bold');
    set(hx,'fontsize',18,'fontweight','bold');
    set(hy,'fontsize',18,'fontweight','bold');
    set(ht,'fontsize',20,'fontweight','bold');

    legend('GOM','GB','ENE','SNE','NYB1','NYB2','SS1','SS2','MAB','location','southwest');
    eval(['print -dpng -r300 ' figoutdir 'MAB_MonthlySalt.png']);

end %if
