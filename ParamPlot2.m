datetick=[];
for month=2:24;
    	datetick=[datetick,datenum(2020,month,1)];
end
dstr=datestr(datetick,'mmm dd');
dstr(:,5)=[];

[N,statsd]=size(Stats);
    
clf
set(gcf,'Position',[1800 680 750 450])
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


t = tiledlayout(2,2);
%xlabel=(t,"Days Since First Infection")
%xticks(datetick)
%xticklabels(dstr)

ax1 = nexttile;
plot(Tth, alpha_plot,'b')
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
box on
grid on


ax2 = nexttile;
plot(Tth,delta_plot,'b')
title("Social Distancing Parameter \delta(t)/\delta-max")
%xticks(datetick)
%xticklabels(dstr)
ylim(ax2,[0,1.2])
h=gca;
set(h,'XTick',datetick);
xtickangle(40)
xtick=get(h,'XTick');
set(h,'XTickLabel',dstr);
box on
grid on

ax3 = nexttile;
plot(Tth,beta_plot,'b')
title("Probability of Infection \beta/\beta-max")
%xticks(datetick)
%xticklabels(dstr)
ylim(ax3,[0,1.2])
h=gca;
set(h,'XTick',datetick);
xtickangle(40)
xtick=get(h,'XTick');
set(h,'XTickLabel',dstr);
box on
grid on

Gamma_time = length(TS) + 14; 
%Here we only plot gamma 2 weeks into the future from the data 
%Since gamma becomes less certain the farther we project into the future we
%limit it

ax4 = nexttile;
plot(TS,gamma_plot(1:length(TS)),'b')
title("Percentage of Confirmed Infections \gamma(t)")
%xticks(datetick)
%xticklabels(dstr)
%ylim(ax4,[0,1.2])
h=gca;
set(h,'XTick',datetick);
xtickangle(40)
xtick=get(h,'XTick');
set(h,'XTickLabel',dstr);
box on
grid on


%figh=plot(plot_time,alpha_plot,'b',plot_time,delta_plot,'k');
%legend(figh,'\alpha','\delta','Location','best')
%set(ax1,'LineWidth',3)
%set(ax2,'LineWidth',3)
%hold off
%axis([xmin,xmax,ymin,ymax])
%h=gca;
%set(h,'XTick',datetick);
%xtickangle(40)
%xtick=get(h,'XTick');
%set(h,'XTickLabel',dstr);
%h=gca;
%set(h,'XTick',datetick);
%xtickangle(40)
%xtick=get(h,'XTick');
%set(h,'XTickLabel',dstr);
%box on
%grid on
%tstr=['Model Parameters w/respect to Time'];
%title(tstr);
%set(gca,'FontSize',18)			
%eval(['print -depsc ',pstr])
eval(['print -dpng -r',int2str(PlotRes),' ',paramstr])


