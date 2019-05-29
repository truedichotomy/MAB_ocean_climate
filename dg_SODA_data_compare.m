% DG 2018-04-07
% this script compares the time series data for SODAcasts vs NEFSC CTD casts.

loadflag = 1
plotflag = 2

yind = [4:39];

if loadflag == 1
    sodaflag = 0
    dg_setup_MABclimate_dir
    load([workdir 'MABclimate2DaSODA_nobiasDG']);
    cccsoda = ccc;
    load([workdir 'MABclimate2Da_nobiasDG']);
end %if

if plotflag == 1
    for si = 1:3
        for ri = 1:9
            [si ri]
            gtempind = find(~isnan(ccc.temp_nobias(yind,si,ri)));
            gsaltind = find(~isnan(ccc.salt_nobias(yind,si,ri)));
            rtemp = corrcoef(ccc.temp_nobias(yind(gtempind),si,ri),cccsoda.temp_nobias(yind(gtempind),si,ri));
            rsalt = corrcoef(ccc.salt_nobias(yind(gsaltind),si,ri),cccsoda.salt_nobias(yind(gsaltind),si,ri));

            f1 = figure;
            set(gcf,'visible','off')

            subplot(2,1,1)
            hold on
            plot(ccc.yyyy,ccc.temp_nobias(:,si,ri),'k');
            plot(cccsoda.yyyy,cccsoda.temp_nobias(:,si,ri),'b-');
            hold off
            legend('CTD','SODA','location','southeast')
            ht1 = title([region(ri).label ' Temperature for ' season(si).label '  R=' num2str(round(rtemp(1,2),2))]);
            set(gca,'xgrid','on','ygrid','on')

            subplot(2,1,2)
            hold on
            plot(ccc.yyyy,ccc.salt_nobias(:,si,ri),'k');
            plot(cccsoda.yyyy,cccsoda.salt_nobias(:,si,ri),'b-');
            hold off
            legend('CTD','SODA','location','southeast')
            ht2 = title([region(ri).label ' Salinity for ' season(si).label '  R=' num2str(round(rsalt(1,2),2))]);
            set(gca,'xgrid','on','ygrid','on')

            eval(['print -depsc -r200 ' figoutdir 'tempsaltSODAcompare_' region(ri).label '_' season(si).name '.eps']);
            close(f1)
        end %for ri
    end %for si
elseif plotflag == 2
    for ri = 1:9
        rtemp = [];
        rsalt = [];

        for si = 1:3
            gtempind = find(~isnan(ccc.temp_nobias(yind,si,ri)));
            gsaltind = find(~isnan(ccc.salt_nobias(yind,si,ri)));
            rtemp{si} = corrcoef(ccc.temp_nobias(yind(gtempind),si,ri),cccsoda.temp_nobias(yind(gtempind),si,ri));
            rsalt{si} = corrcoef(ccc.salt_nobias(yind(gsaltind),si,ri),cccsoda.salt_nobias(yind(gsaltind),si,ri));
        end %for

        f1 = figure;
        set(gcf,'visible','off')

        subplot(3,2,1)
        hold on
        hp1 = plot(ccc.yyyy,ccc.temp_nobias(:,1,ri),'k');
        hp2 = plot(cccsoda.yyyy,cccsoda.temp_nobias(:,1,ri),'k--');
        set(hp1,'linewidth',1)
        set(hp2,'linewidth',1)
        hold off
        set(gca,'box','on')
        %legend('CTD','SODA','location','southeast')
        ht1 = title([region(ri).label ' Temperature for ' season(1).label '  R=' num2str(round(rtemp{1}(1,2),2))]);
        set(gca,'xgrid','on','ygrid','on')

        subplot(3,2,2)
        hold on
        hp1 = plot(ccc.yyyy,ccc.salt_nobias(:,1,ri),'k');
        hp2 = plot(cccsoda.yyyy,cccsoda.salt_nobias(:,1,ri),'k--');
        set(hp1,'linewidth',1)
        set(hp2,'linewidth',1)
        hold off
        set(gca,'box','on')
        legend('CTD','SODA','location','southeast')
        ht2 = title([region(ri).label ' Salinity for ' season(1).label '  R=' num2str(round(rsalt{1}(1,2),2))]);
        set(gca,'xgrid','on','ygrid','on')

        subplot(3,2,3)
        hold on
        hp1 = plot(ccc.yyyy,ccc.temp_nobias(:,2,ri),'k');
        hp2 = plot(cccsoda.yyyy,cccsoda.temp_nobias(:,2,ri),'k--');
        set(hp1,'linewidth',1)
        set(hp2,'linewidth',1)
        hold off
        set(gca,'box','on')
        %legend('CTD','SODA','location','southeast')
        ht1 = title([region(ri).label ' Temperature for ' season(2).label '  R=' num2str(round(rtemp{2}(1,2),2))]);
        set(gca,'xgrid','on','ygrid','on')

        subplot(3,2,4)
        hold on
        hp1 = plot(ccc.yyyy,ccc.salt_nobias(:,2,ri),'k');
        hp2 = plot(cccsoda.yyyy,cccsoda.salt_nobias(:,2,ri),'k--');
        set(hp1,'linewidth',1)
        set(hp2,'linewidth',1)
        hold off
        set(gca,'box','on')
        %legend('CTD','SODA','location','southeast')
        ht2 = title([region(ri).label ' Salinity for ' season(2).label '  R=' num2str(round(rsalt{2}(1,2),2))]);
        set(gca,'xgrid','on','ygrid','on')

        subplot(3,2,5)
        hold on
        hp1 = plot(ccc.yyyy,ccc.temp_nobias(:,3,ri),'k');
        hp2 = plot(cccsoda.yyyy,cccsoda.temp_nobias(:,3,ri),'k--');
        set(hp1,'linewidth',1)
        set(hp2,'linewidth',1)
        hold off
        set(gca,'box','on')
        %legend('CTD','SODA','location','southeast')
        ht1 = title([region(ri).label ' Temperature for ' season(3).label '  R=' num2str(round(rtemp{3}(1,2),2))]);
        set(gca,'xgrid','on','ygrid','on')

        subplot(3,2,6)
        hold on
        hp1 = plot(ccc.yyyy,ccc.salt_nobias(:,3,ri),'k');
        hp2 = plot(cccsoda.yyyy,cccsoda.salt_nobias(:,3,ri),'k--');
        set(hp1,'linewidth',1)
        set(hp2,'linewidth',1)
        hold off
        set(gca,'box','on')
        %legend('CTD','SODA','location','southeast')
        ht2 = title([region(ri).label ' Salinity for ' season(3).label '  R=' num2str(round(rsalt{3}(1,2),2))]);
        set(gca,'xgrid','on','ygrid','on')

        eval(['print -depsc -r300 ' figoutdir 'tempsaltSODAcompare6_' region(ri).label '.eps']);
        close(f1)
    end %for

end %if plotflag
