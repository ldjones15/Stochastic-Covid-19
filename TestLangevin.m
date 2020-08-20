clear

BET=0.7;
tbet=[14;15;16];
gam=0.05;
alp=1-gam;
M0=1;
k=8;
T=150;
X=9000000;
KAPL=0.1;
KAPS=0.1;
Zeta=2;
EndTime=[100;100;120];
LamMask=0.6;

[Sth,Mth,Qth,Rth,Ith]=Compute(BET,tbet,gam,alp,M0,k,T,X,KAPL,KAPS,Zeta,EndTime,LamMask);
ComputeE(BET,tbet,gam,alp,M0,k,Sth,X,KAPL,KAPS,Zeta,EndTime,LamMask)

lam=0.25;
figh=plot(diff(Sth).^lam);
%figh=plot(Sth.^lam);
%figh=plot(Ith);
ytick=get(gca,'YTick');
set(gca,'YTickLabel',round(ytick.^(1/lam)));
