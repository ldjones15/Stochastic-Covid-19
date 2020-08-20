datetick=[];
for month=2:24;
    	datetick=[datetick,datenum(2020,month,1)];
end
dstr=datestr(datetick,'mmm dd');
dstr(:,5)=[];

gamma_plot = ((MeasureI./pop)./(TotI));

clf
set(gcf,'Position',[1800 680 750 450])
hold on
xmin=min(Tth);xmax=max(Tth);ymin=min(Iplot)^lam;ymax=max(max(Iplot)^lam,max(IS)^lam);
if(usebounds==1)
	xmin=min(Tth);xmax=max(Tth);ymin=min(Imin)^lam;ymax=max(max(Imax)^lam,max(IS)^lam);
	figs=fill([Tth;flipud(Tth)],[Imin.^lam;flipud(Imax.^lam)],cgray);
	set(figs,'FaceColor',cgray,'EdgeColor',cgray)
end
figh=plot(TS,IS.^lam,'ro',Tth,Iplot.^lam,'k');
legend(figh,'Observed','Model','Location','best')
set(figh(2),'LineWidth',3)
set(figh(1),'MarkerSize',6,'LineWidth',2)
hold off
axis([xmin,xmax,ymin,ymax])
h=gca;
set(h,'XTick',datetick);
xtickangle(40)
xtick=get(h,'XTick');
set(h,'XTickLabel',dstr);
yabs=max(abs([ymin,ymax].^(1/lam)));
ytick=get(gca,'YTick');
if(yabs<=1000)
	yticklabel=round(ytick.^(1/lam));
	set(h,'YTick',ytick,'YTickLabel',yticklabel);	
	ystr='Infections';

else
	sigfig=ceil(log10(yabs));
	yticklabel=round(ytick.^(1/lam)/10^(sigfig-3));
	set(h,'YTick',ytick,'YTickLabel',yticklabel);
	ystr=['Infections (x ',int2str(10^(sigfig-3)),')'];
end
box on
grid on
tstr=['Confirmed Daily Infections'];
title(tstr);
ylabel(ystr)
set(gca,'FontSize',18)			
%eval(['print -depsc ',pstr])
eval(['print -dpng -r',int2str(PlotRes),' ',pstr])


