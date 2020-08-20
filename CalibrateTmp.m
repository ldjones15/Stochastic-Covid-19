function [EMIN,BETmin,TBETmin,GAMmin,ALPmin,MFULLmin,Xmin]=Calibrate(BET,TBET,GAM,ALP,MFULL,K,X0,S,Etarg,sstr,NumTries,delbet)
%NumTries should be odd;

BETmin=zeros(length(K),length(TBET(1,:)));
TBETmin=BETmin;
GAMmin=zeros(length(K),1);
ALPmin=zeros(length(K),1);
MFULLmin=zeros(length(K),1);
Xmin=zeros(length(K),1);

Xtry=X0;

ngam=length(GAM);
nalp=1;
nK=length(K);
nX=length(Xtry(1,:));

EMIN=100000000000*mean(S.^2)*ones(nK,1);

for ik=1:nK
	k=K(ik);
	bet0=BET(ik,:);
	tbet=TBET(ik,:);tbet(end)=length(S);
	M0=MFULL(ik);
	nbet=NumTries^(length(tbet));
	niter=ngam*nalp*nX*nbet;
	iter=0;
	for igam=1:ngam
		if(niter>500000)	
%			[k,ik,nK,NumTries, iter/niter]
		end
	gam=GAM(igam);
	alp=1-gam;
	
	for iX=1:nX
	X=Xtry(ik,iX);

	for ibet=0:nbet-1
		iter=iter+1;
		A=ibet;
		for ic=1:length(tbet)
			Af=floor(A/NumTries);
			aic=(A-NumTries*Af);
			A=Af;
			bet(ic)=((1-delbet)+aic*2*delbet/(NumTries-1))*bet0(ic);
		end
		bet=min(bet,0.999);
		bet=max(bet,0.001);
		Ecur=ComputeE(bet,tbet,gam,alp,M0,k,S,X);
		if(Ecur<EMIN(ik))
			EMIN(ik)=Ecur;
			BETmin(ik,:)=bet;TBETmin(ik,:)=tbet;
			GAMmin(ik)=gam;ALPmin(ik)=alp;
			MFULLmin(ik)=M0;Xmin(ik)=X;
		end
end
end
end
end


