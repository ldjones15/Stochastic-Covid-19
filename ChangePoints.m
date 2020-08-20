function [tbet,bethat]=ChangePoints(I,cmax,niter,pval,rthresh,taumin)
%I infection counts up to turning point.
%iter is the number of trials to determine statistical significance.
%cmax is the max number of beta regions (e.g. 10).
%tbet gives the changepoints in the observed. (includes lag)
%bethat is the exp-growth factor, so subtract 1 to get beta.

I=I+0.1;    		  %regularize for taking logs if I can be 0.

T=(1:length(I))';
jpos=find(I>0);
x=T(jpos);
y=log(I(jpos));
n=length(x);
cmaxposs=min(cmax,floor(n/taumin));

%%%%%First Fit optimal fits up to cmaxposs exponentials.%%%%%
[DP,iDP]=GetDPtaumin(x,y,cmaxposs,taumin);

%%%%%Compute optimal fits.%%%%%

IDP=zeros(cmaxposs,n);
yDP=zeros(cmaxposs,n);

for c=1:cmaxposs
	iend=n;
	for curc=c:-1:1
		ibeg=iDP(curc,iend);
		xj=x(ibeg:iend);
		yj=y(ibeg:iend);
		aj=polyfit(xj,yj,1);
		yhat=aj(1)*xj+aj(2);
		yDP(c,ibeg:iend)=yhat;
		IDP(c,ibeg:iend)=exp(yDP(c,ibeg:iend));
		iend=ibeg;		
	end
end


%%%%%Statistical test to see if going from c to c+1 gives improvement.%%%%%
copt=1;
for c=1:cmaxposs-1
	y0=yDP(c,:)';
	sig=sqrt(n/(n-2*c))*sqrt(mean((y-y0).^2));
	Egain=zeros(niter,1);
	for iter=1:niter
		ynoise=y0+sig*randn(size(y0));
		DPiter=GetDPtaumin(x,ynoise,c+1,taumin);
		Egain(iter)=DPiter(c,end)-DPiter(c+1,end);
	end
	if(quantile(Egain,pval)>=rthresh*(DP(c,end)-DP(c+1,end)))
		copt=c;
		break;
	else
		copt=c+1;
	end
end

%%%%%Obtain changepoint times and exponential growth parameter%%%%%
c=copt;
tbet=zeros(1,c);
bethat=zeros(1,c);
iend=n;
for curc=c:-1:1
	tbet(curc)=iend;
	ibeg=iDP(curc,iend);
	xj=x(ibeg:iend);
	yj=y(ibeg:iend);
	aj=polyfit(xj,yj,1);
	bethat(curc)=exp(aj(1));
	iend=ibeg;		
end


%{

%Merge changepoints that are too small (overfitted).
%Not needed if using GetDPtaumin

Dtbet=[tbet(1),diff(tbet)];
while(min(Dtbet)<taumin&length(Dtbet)>1)
	for i=[1,length(Dtbet),2:length(Dtbet)-1]
		if(Dtbet(i)<taumin&i==1)
			bethat(2)=0.5*(bethat(1)+bethat(2));
			tbet(1)=[];
			bethat(1)=[];
			break;
		end
		if(Dtbet(i)<taumin&i==length(Dtbet))
			bethat(end-1)=0.5*(bethat(end)+bethat(end-1));
			tbet(end-1)=tbet(end);
			tbet(end)=[];
			bethat(end)=[];
			break;
		end
		if(Dtbet(i)<taumin)
			if(abs(bethat(i)-bethat(i-1))<=abs(bethat(i)-bethat(i+1)))
				bethat(i-1)=0.5*(bethat(i-1)+bethat(i));
				tbet(i-1)=tbet(i);
				tbet(i)=[];
				bethat(i)=[];
			else
				bethat(i+1)=0.5*(bethat(i+1)+bethat(i));
				tbet(i)=[];
				bethat(i)=[];
			end
			break;
		end
	end
	Dtbet=[tbet(1),diff(tbet)];
end

%%%%%Re-Obtain changepoint times and exponential growth parameter%%%%%
c=length(tbet);
bethat=zeros(1,c);
for curc=1:length(tbet)
	if(curc==1)
		ibeg=1;
	else
		ibeg=tbet(curc-1);
	end
	iend=tbet(curc);
	xj=x(ibeg:iend);
	yj=y(ibeg:iend);
	aj=polyfit(xj,yj,1);
	bethat(curc)=exp(aj(1));
%	aj
%	plot(xj,yj,'ro',xj,aj(1)*xj+aj(2));pause
end

%}
