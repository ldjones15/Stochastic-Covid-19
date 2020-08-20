(* ::Package:: *)

function [Sth,Mth,Qth,Rth,Ith]=Compute(BET,tbet,gam,alp,M0,k,T,X,KAPL,KAPS,Zeta,EndTime,LamMask,LockRel)
% tbet specifies the start of SD and LockDown.
%tbet and EndTime actions have effects k time steps later.

Sth=zeros(T+k+1,1);
Mth=zeros(T+k+1,1);
Qth=zeros(T+k+1,1);
Rth=zeros(T+k+1,1);
Ith=zeros(T+k+1,1);
Mth(2)=M0;
Qth(2)=M0;% New infections
Ith(2)=M0/X;

for t=3:length(Sth)
		kapL=KAPL;kapS=KAPS;bet=BET;
		if(t-k-1<=tbet(2))
			kapS=0;
		end
		if(t-k-1>EndTime(2))
			kapS=LockRel*kapS;
		end
		if(t-k-1<=tbet(1))
			kapL=0;
		end
		if(t-k-1>EndTime(1))
			kapL=LockRel*kapL;
		end
		if(t-k-1>tbet(3)&t-k-1<=EndTime(3))
			bet=LamMask*bet;
		end
		Ihat=Mth (t-1)/X;
		Qth(t)=X*(bet*Ihat*(1-Ith(t-1))-kapL*Ihat-kapS*(sqrt(Ihat))^(Zeta)*(1-2*Ihat));
%		Qth(t)=X*(bet*Ith(t-1)*(1-Ith(t-1))-kapL*Ith(t-1)-kapS*(sqrt(Ith(t-1)))^(Zeta)*(1-2*Ith(t-1)));
%		[t,kapS,kapL,bet,Qth(t),Ith(t-1)]
%		pause			
		if(Qth(t)<0.00000001)
			Qth(t)=0;
		end
		Ith(t)=Ith(t-1)+Qth (t)/X;
		%[t,Qth(t),Ith(t)]
		% pause
		if(t<=k+1)
			Mth(t)=Mth(t-1)+Qth(t);
		else
			Sth(t)=Sth(t-1)+gam*Qth(t-k);
			Rth(t)=Rth(t-1)+alp*Qth(t-k);
			Mth(t)=Mth(t-1)+Qth(t)-(gam+alp)*Qth(t-k);
		end		
	end

end

