function XFULL=GetX(BETmin,TBETmin,GAMmin,ALPmin,MFULLmin,K0,Xmin,S,tp,pop,minfit);

XFULL=Xmin;
XTRY=[0.005,0.01:0.01:0.1,0.15:0.05:1]*pop;

%%%%%%%% Determine X %%%%%%%%%%%%%%%%%

for iK=1:length(K0)

	k=K0(iK);
	bet=BETmin(iK,:);
	tbet=TBETmin(iK,:);
	gam=GAMmin(iK,:);
	alp=ALPmin(iK,:);
	M0=MFULLmin(iK,:);
	for iX=1:length(XTRY)
		Xcur=XTRY(iX);
		Eout=ComputeE(bet,tbet,gam,alp,M0,k,S,Xcur);
		if(iX==1)
			Ecurmin=Eout;
			Xcurmin=Xcur;
		end
		if(Eout<Ecurmin)
			Ecurmin=Eout;
			Xcurmin=Xcur;
		end		
	end
	XFULL(iK)=Xcurmin;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

