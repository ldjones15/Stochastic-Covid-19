function [BET,TBET,GAM,ALP,MFULL,ETARG,tp,ifit,TBETfull]=Initialize(Sraw,K,X,minfit,cmax,cmaxTP,niterC,pval,rthresh,tp_robust,threshneg,KAPL,KAPS,Smooth,Zeta,LamMask,GAMTRY,BETTRY,LockRel)
%S=cum infections from first infection;k=lag;minfit=min-length of changepoint;
%cmax=max number of phases;X=pop-at-risk;niterC=number of random to estimate
%significance;pval=confidence changepoint is non-random;rthresh=fraction of
%improvement random needs for changepoint to not be random.


%First identify the turning point.
Iraw=[Sraw(1);diff(Sraw)];
[tbet1,egrowth1]=ChangePoints(Iraw,cmaxTP,niterC,pval,rthresh,min(minfit+3,length(Iraw)));


threshpos=1.01;
ineg=find(egrowth1<=threshneg);
itmp=find(ineg>1);
ineg=ineg(itmp);

if(length(ineg)==0)
	ifit=length(Sraw);
else
	jfit=ineg(1);
	while(egrowth1(jfit)<=threshpos)
		ifit=tbet1(jfit);
		jfit=jfit+1;
		if(jfit>length(egrowth1))
			break;
		end
	end
end


tp=TurningPoint(Iraw(1:ifit),tp_robust);


Igrow=Iraw(1:tp);
Sgrow=Sraw(1:tp);

%Identify changepoints.

minfit=minfit+max(ceil((tp-38)/8),0);


[tbetOBS,egrowth]=ChangePoints(Sgrow,cmax,niterC,pval,rthresh,min(minfit,length(Sgrow)));
%tbetOBS is an array of the observed times where changepoints occur

if(length(Sraw)-ifit<=minfit)
	TBETfull=ones(length(K),1)*tbetOBS-K*ones(1,length(tbetOBS));
	TBETfull(:,end)=length(Sraw);
else
	[tbetOBS1,egrowth1]=ChangePoints(Sraw(ifit:end),cmax,niterC,pval,rthresh,minfit);
	tbetfull=tbetOBS;tbetfull(:,end)=ifit;
	tbetfull=[tbetfull,tbetOBS1+ifit-1];
	TBETfull=ones(length(K),1)*tbetfull-K*ones(1,length(tbetfull));
	TBETfull(:,end)=length(Sraw);
end

S=Sraw(1:ifit);

%length(tbetOBS) represents the number of changepoints
%length(K) represents the number of lags
%For the current data they both should be 1, so these values are 1x1 matrices
BET=zeros(length(K),length(tbetOBS)); 
TBET=zeros(length(K),length(tbetOBS));
GAM=zeros(length(K),1);
ALP=zeros(length(K),1);
MFULL=zeros(length(K),1);
ETARG=zeros(length(K),1);
XFULL=zeros(length(K),1);


for iK=1:length(K)
k=K(iK);

tbet=tbetOBS-k;     %changepoints in beta occur k steps before observation.
bethat=min(0.999,(egrowth-1)./(1-1./egrowth.^k));


%Ready to start the fitting. 

Sk=[zeros(k+1,1);S];
Tk=(1:length(Sk))';

A=zeros(length(tbet),2);
if(tbetOBS(1)-k-1>=minfit)
	x1=Tk(2*k+3:k+1+tbetOBS(1));
	y1=Sk(2*k+3:k+1+tbetOBS(1));
else
	x1=Tk(k+2:k+1+tbetOBS(1));
	y1=Sk(k+2:k+1+tbetOBS(1));
end

A(1,:)=exp(polyfit(x1,log(y1),1));
for cidx=2:length(tbetOBS)
	xidx=Tk(k+1+tbetOBS(cidx-1)+1:k+1+tbetOBS(cidx));
	yidx=Sk(k+1+tbetOBS(cidx-1)+1:k+1+tbetOBS(cidx));
	A(cidx,:)=exp(polyfit(xidx,log(yidx),1));
end


%Since we can only identify gam*M0 and (gam+alp), assume M0=ceil(S(1));
% Start by assuming alp=gam;
%Try a bunch of gammas.

M0=ceil(Sk(k+2));

tbet(end)=length(S);

tstart=tbet(1)*ones(3,1);
EndTime=(length(S)+k)*ones(3,1);



for i=1:length(GAMTRY)
	gam=GAMTRY(i);
	alp=1-gam;
	for ibet=1:length(BETTRY)
		bet=BETTRY(ibet);
		Eout=ComputeE(bet,tstart,gam,alp,M0,k,S,X,KAPL,KAPS,Smooth,Zeta,EndTime,LamMask,LockRel);
		if(i==1)
			gamopt=gam;
			Eopt=Eout;
			betopt=bet;
		elseif(Eout<Eopt)
			gamopt=gam;
			Eopt=Eout;
			betopt=bet;
		end
	end
end


gam=gamopt;
alp=1-gam;
bet=betopt;
Etarg=Eopt;

BET(iK,1)=bet;
TBET(iK,:)=tbet;
GAM(iK)=gam;
ALP(iK)=alp;
MFULL(iK)=M0;
ETARG(iK)=Etarg;

end
