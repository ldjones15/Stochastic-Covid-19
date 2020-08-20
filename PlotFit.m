k=K0(idxk,:);
bet=BETmin(idxk,1);
tbet=TLOCK(idxk,:);
EndTime=ENDTIME(idxk,:);
gam=GAMmin(idxk,:);
alp=1-gam;
M0=MFULLmin(idxk,:);
T=length(Sfit);
X=X0;
KAPL=KAPLmin(idxk,:);
KAPS=KAPSmin(idxk,:);
LamMask=LamMaskmin(idxk,:);
Smooth=Smoothmin(idxk,:);


[Sth,Mth,Qth,Rth,Ith]=Compute(bet,tbet,gam,alp,M0,k,T,X,KAPL,KAPS,Smooth,Zeta,EndTime,LamMask,LockRel);
Tplot=2:length(Sfit);
Sth=Sth(k+2:end);
plot(Tplot,diff(Sfit).^0.25,'ro',Tplot,diff(Sth).^0.25,'k');

ComputeE(bet,tbet,gam,alp,M0,k,Sfit,X,KAPL,KAPS,Smooth,Zeta,EndTime,LamMask,LockRel)
