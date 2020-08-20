datetick=[];
for month=2:24;
    	datetick=[datetick,datenum(2020,month,1)];
end
dstr=datestr(datetick,'mmm dd');
dstr(:,5)=[];

%k=K0(idxnow,:);
%bet=BETmin(idxnow,1);
%tbet=TLOCK(idxnow,:);
%EndTime=ENDTIME(idxnow,:);
%gam=GAMmin(idxnow,:);
%alp=1-gam;
%M0=MFULLmin(idxnow,:);
%T=length(Sfit);
%X=X0;
%KAPL=KAPLmin(idxnow,:);
%KAPS=KAPSmin(idxnow,:);
%LamMask=LamMaskmin(idxnow,:);


[plot_time,kapl_plot,kaps_plot,LAM_plot]=ComputePlot(bet,tbet,gam,alp,M0,k,T,X,KAPL,KAPS,Zeta,EndTime,LamMask,LockRel);
%Tplot=3:length(Sfit);
Mi = min(kaps_plot);
Mx = max(kaps_plot);
Mxl = max(kapl_plot);
h=gca;
figh=plot(plot_time,kapl_plot/Mxl.^0.25,'k', plot_time, kaps_plot/Mx,'r');
legend(figh,'\alpha','\delta','Location','best')
xtickangle(40)
xtick=get(h,'XTick');
set(h,'XTickLabel',dstr);
title('Model Parameters w/ respect to time')
xlabel('Date')
%ylabel('\delta')
%axis([0 140 0.55,0.56]);

%ComputeE(bet,tbet,gam,alp,M0,k,Sfit,X,KAPL,KAPS,Zeta,EndTime,LamMask,LockRel)
