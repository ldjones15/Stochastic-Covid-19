function [DP,iDP]=GetDP(x,y,cmax,taumin);
%assumes x valuse are distinct and x(1)<x(2)<x(3)...
%at most cmax phases each of length at least taumin

  
n=length(x);

DP=zeros(cmax,n);
iDP=ones(cmax,n);

Y2=cumsum(y.^2);
X2=cumsum(x.^2);
X=cumsum(x);
Y=cumsum(y);
XY=cumsum(x.*y);

%c=1 means no changepoint, one phase.
%DP(c,t)=best error up to time t with c phases.
c=1;
for t=3:n
	sigy2=Y2(t)-Y(t)^2/t;
	sigx2=X2(t)-X(t)^2/t;
	sigxy=XY(t)-X(t)*Y(t)/t;
	DP(c,t)=sigy2-sigxy^2/sigx2;
end

for c=2:cmax
	iDP(c,2:c+1)=(2:c+1)-1;
	for t=c+2:n
		for j=2:t
			nj=t-j+1;
			if(nj<=2)
				Ej=0;
			else
				y2sum=Y2(t)-Y2(j-1);
				ysum=Y(t)-Y(j-1);
				x2sum=X2(t)-X2(j-1);
				xsum=X(t)-X(j-1);
				xysum=XY(t)-XY(j-1);
				sigy2=y2sum-ysum^2/nj;
				sigx2=x2sum-xsum^2/nj;
				sigxy=xysum-xsum*ysum/nj;
				Ej=sigy2-sigxy^2/sigx2;
			end
			if(j<=2)
				Eprefix=0;
			else
				Eprefix=DP(c-1,j);
			end
			Etot=Ej+Eprefix;
			if(j==2)
				DP(c,t)=Etot;
				iDP(c,t)=j;
			else
				if(Etot<DP(c,t))
					DP(c,t)=Etot;
					iDP(c,t)=j;
				end
			end
		end
	end
end


