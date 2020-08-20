datetick=[];
for month=2:24;
    	datetick=[datetick,datenum(2020,month,1)];
end
dstr=datestr(datetick,'mmm dd');
dstr(:,5)=[];

[N,statsd]=size(Stats);
    
clf
%set(gcf,'Position',[1800 680 750 225])
%hold on

for i=1:N
	Eout=Stats(i,end);
	X=Stats(i,end-1);
	k=Stats(i,end-2);
	M0=Stats(i,end-3);
	alp=Stats(i,end-4);
	gamma=Stats([1:N],end-5);
	beta=Stats(i,1);
	tbet=Stats(i,2:4);
	EndTime=Stats(i,5:7)+AddTime;
    %EndTime=EndTime+AddTime;
	KAPL=Stats(i,8);
	KAPS=Stats(i,9);
	LamMask=Stats(i,10);
	Zeta=Stats(i,11);
	LockRel=Stats(i,12);
end

alpha_plot=kapl_plot/kapl_max;
delta_plot=kaps_plot/kaps_max;
beta_plot=LAM_plot/lam_max;
gamma_plot = ((MeasureI./pop)./(TotI));


fig1 = figure(1);
set(gcf,'Position',[1800 680 750 225]);
ax1 = axes(fig1);
plot(Tth, alpha_plot,'b','LineWidth',2);
title("Lockdown Parameter \alpha(t)/\alpha-max")
%ax1.Xticks = datetick
%ax1.Xticklabel = dstr
%xticks(ax1, [datetick])
%xticklabels(ax1, {dstr})
ylim(ax1,[0,1.2])
h=gca;
set(h,'XTick',datetick);
xtickangle(40)
xtick=get(h,'XTick');
set(h,'XTickLabel',dstr);
ax1.FontWeight = 'bold';
box on
grid on

eval(['print -dpng -r',int2str(PlotRes),' ',lockdownstr]);



fig2 = figure(2);
set(gcf,'Position',[1800 680 750 225]);
ax2 = axes(fig2);
plot(Tth,delta_plot,'b','LineWidth',2);
title("Social Distancing Parameter \delta(t)/\delta-max")
%xticks(datetick)
%xticklabels(dstr)
ylim(ax2,[0,1.2])
h=gca;
set(h,'XTick',datetick);
xtickangle(40)
xtick=get(h,'XTick');
set(h,'XTickLabel',dstr);
ax2.FontWeight = 'bold';
box on
grid on

eval(['print -dpng -r',int2str(PlotRes),' ',distancestr]);



fig3 = figure(3);
set(gcf,'Position',[1800 680 750 225]);
ax3 = axes(fig3);
plot(Tth,beta_plot,'b','LineWidth',2);
title("Probability of Infection \beta/\beta-max")
%xticks(datetick)
%xticklabels(dstr)
ylim(ax3,[0,1.2])
h=gca;
set(h,'XTick',datetick);
xtickangle(40)
xtick=get(h,'XTick');
set(h,'XTickLabel',dstr);
ax3.FontWeight = 'bold';
box on
grid on

eval(['print -dpng -r',int2str(PlotRes),' ',betastr]);



fig4 = figure(4);
set(gcf,'Position',[1800 680 750 225]);
ax4 = axes(fig4);
plot(TS,gamma_plot(1:length(TS)),'b','LineWidth',2);
title("Ratio of Confirmed Infections \gamma(t)")
%xticks(datetick)
%xticklabels(dstr)
%ylim(ax4,[0,1.2])
h=gca;
set(h,'XTick',datetick);
xtickangle(40)
xtick=get(h,'XTick');
set(h,'XTickLabel',dstr);
ax4.FontWeight = 'bold';
box on
grid on

eval(['print -dpng -r',int2str(PlotRes),' ',gammastr]);
