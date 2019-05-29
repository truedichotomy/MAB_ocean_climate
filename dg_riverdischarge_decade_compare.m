% dg_riverdischarge_decade_compare.m

%if ~exist('flow_chesapeake_monthly')
load('riverdischarge.mat');
%end %if

p1ind = find(yyyylist <= 1999);
p2ind = find(yyyylist > 1999);

flow_chesapeake_monthly_p1 = nanmean(flow_chesapeake_monthly(p1ind,:),1);
flow_chesapeake_monthly_p2 = nanmean(flow_chesapeake_monthly(p2ind,:),1);
flow_chesapeake_monthly_p1p2 = [flow_chesapeake_monthly_p1; flow_chesapeake_monthly_p2];

flow_hudson_monthly_p1 = nanmean(flow_hudson_monthly(p1ind,:),1);
flow_hudson_monthly_p2 = nanmean(flow_hudson_monthly(p2ind,:),1);
flow_hudson_monthly_p1p2 = [flow_hudson_monthly_p1; flow_hudson_monthly_p2];

%bar(flow_hudson_monthly_p1p2')
%bar(flow_chesapeake_monthly_p1p2')

[errp1_hudson] = dg_ci(flow_hudson_monthly(p1ind,:)');
[errp2_hudson] = dg_ci(flow_hudson_monthly(p2ind,:)');
[errp1_chesapeake] = dg_ci(flow_chesapeake_monthly(p1ind,:)');
[errp2_chesapeake] = dg_ci(flow_chesapeake_monthly(p2ind,:)');

figure(1)
h_hudson = barerrorbar({flow_hudson_monthly_p1p2'},{flow_hudson_monthly_p1p2',[errp1_hudson,errp2_hudson],'k.'});
title('Hudson River Monthly Discharge Comparison')
legend('1977-1999','2000-2016')
xlabel('Month')
ylabel('Discharge (m^{3}/s)')
eval(['print -depsc -r300 MonthlyHudsonDischargeDecadalCompare.eps'])

figure(2)
h_chesapeake = barerrorbar({flow_chesapeake_monthly_p1p2'},{flow_chesapeake_monthly_p1p2',[errp1_chesapeake,errp2_chesapeake],'k.'});
title('Chesapeake Monthly Discharge Comparison')
legend('1977-1999','2000-2016')
xlabel('Month')
ylabel('Discharge (m^{3}/s)')
eval(['print -depsc -r300 MonthlyChesapeakeDischargeDecadalCompare.eps'])

flow_chesapeake_all_p1 = nanmean(nanmean(flow_chesapeake_monthly(p1ind,:),2))
flow_chesapeake_all_p2 = nanmean(nanmean(flow_chesapeake_monthly(p2ind,:),2))

flow_hudson_all_p1 = nanmean(nanmean(flow_hudson_monthly(p1ind,:),2))
flow_hudson_all_p2 = nanmean(nanmean(flow_hudson_monthly(p2ind,:),2))



flow_chesapeake_JFMA_p1 = nanmean(flow_chesapeake_JFMA(p1ind))
flow_chesapeake_JFMA_p2 = nanmean(flow_chesapeake_JFMA(p2ind))

flow_hudson_JFMA_p1 = nanmean(flow_hudson_JFMA(p1ind))
flow_hudson_JFMA_p2 = nanmean(flow_hudson_JFMA(p2ind))
