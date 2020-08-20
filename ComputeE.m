function Eout=ComputeE(BET,tbet,gam,alp,M0,k,Strue,X,KAPL,KAPS,Smooth,Zeta,EndTime,LamMask,LockRel)
%Strue must be a column vector.

[S,M,Q,R,I]=Compute(BET,tbet,gam,alp,M0,k,length(Strue),X,KAPL,KAPS,Smooth,Zeta,EndTime,LamMask,LockRel);
S=(S(k+2:end));
Eout=norm(S-Strue)/100+norm((S-Strue)./max(Strue,1))/5;
