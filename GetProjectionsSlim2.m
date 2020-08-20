function [Smin,Smax,Smean,Imin,Imax,Imean,Mmin,Mmax,Mmean,Rmin,Rmax,Rmean,ImaxDist,MDist,TotIDist,plot_time,kapl_plot,kaps_plot,LAM_plot,EndTime,gam]=GetProjectionsSlim(Stats,Tth,Tend,err,pop,tau,Itol,Strue,targidx,AddTime2)


idxcur=length(Strue);
[N,statsd]=size(Stats);
Smin=0*Tth;Imin=0*Tth;Mmin=0*Tth;Rmin=0*Tth;
Smax=0*Tth;Imax=0*Tth;Mmax=0*Tth;Rmax=0*Tth;
Smean=0*Tth;Imean=0*Tth;Mmean=0*Tth;Rmean=0*Tth;

Emin=min(Stats(:,end));
Etarg=(1+err)*Emin;
jpos=find(Stats(:,end)<=Etarg);
num=0;
ImaxDist=zeros(length(jpos),1);
MDist=zeros(length(jpos),2);
TotIDist=zeros(length(jpos),2);

for i=1:N
	Eout=Stats(i,end);
	X=Stats(i,end-1);
	k=Stats(i,end-2);
	M0=Stats(i,end-3);
	alp=Stats(i,end-4);
	gam=Stats(i,end-5);
	BET=Stats(i,1);
	tbet=Stats(i,2:4);
	EndTime=Stats(i,5:7);
    EndTime=EndTime+AddTime2;
	KAPL=Stats(i,8);
	KAPS=Stats(i,9);
	LamMask=Stats(i,10);
	Zeta=Stats(i,11);
	LockRel=Stats(i,12);
    Smooth = Stats(i,13);
    

	if(Eout<=Etarg);
		[Sth,Mth,Qth,Rth,Ith]=Compute(BET,tbet,gam,alp,M0,k,Tend,X,KAPL,KAPS,Smooth,Zeta,EndTime,LamMask,LockRel);
        [plot_time,kapl_plot,kaps_plot,LAM_plot]=ComputePlot(BET,tbet,gam,alp,M0,k,Tend,X,KAPL,KAPS,Smooth,Zeta,EndTime,LamMask,LockRel);
		Sth=Sth(k+2:end);Mth=Mth(k+2:end);Ith=[Sth(1);diff(Sth)];Rth=Rth(k+2:end);
		if(sum(Smax)<1)
			Smin=Sth;Imin=Ith;Mmin=Mth;Rmin=Rth;
			Smax=Sth;Imax=Ith;Mmax=Mth;Rmax=Rth;
		end
		Smin=(min([Smin';Sth']))';Imin=(min([Imin';Ith']))';Mmin=(min([Mmin';Mth']))';Rmin=(min([Rmin';Rth']))';
		Smax=(max([Smax';Sth']))';Imax=(max([Imax';Ith']))';Mmax=(max([Mmax';Mth']))';Rmax=(max([Rmax';Rth']))';
		num=num+1;
		Smean=Smean+(Sth-Smean)/num;
		Imean=Imean+(Ith-Imean)/num;
		Mmean=Mmean+(Mth-Mmean)/num;
		Rmean=Rmean+(Rth-Rmean)/num;
		[maxI,imax]=max(Ith);
		ImaxDist(num)=Tth(imax);
		MDist(num,1)=Mth(idxcur);
		TotIDist(num,1)=Mth(idxcur)+Sth(idxcur)+Rth(idxcur);
		MDist(num,2)=Mth(targidx);
		TotIDist(num,2)=Mth(targidx)+Sth(targidx)+Rth(targidx);
        
	end
end

