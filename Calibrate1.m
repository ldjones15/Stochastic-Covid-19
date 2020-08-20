function [EMIN,BETmin,GAMmin,ALPmin,MFULLmin,KAPLmin,KAPSmin,LamMaskmin,Smoothmin]=Calibrate1(BETTRY,GAMTRY,ALPTRY,MFULL,K,XCAL,KAPLTRY,KAPSTRY,LamMaskTRY,SmoothTRY,TLOCK,ENDTIME,S,Etarg,sstr,delE,Zeta,LockRel)
%NumTries should be odd;
BETmin=zeros(length(K),1);
GAMmin=zeros(length(K),1);
ALPmin=zeros(length(K),1);
MFULLmin=zeros(length(K),1);
KAPLmin=zeros(length(K),1);
KAPSmin=zeros(length(K),1);
LamMaskmin=zeros(length(K),1);
Smoothmin = zeros(length(K),1);

%[fout,msg]=fopen(sstr,'w+')
%assert(fout>=3,msg)

fout=fopen(sstr,'w+');

ngam=length(GAMTRY);
nalp=1;
nK=length(K);
nbet=length(BETTRY);
nKAPL=length(KAPLTRY);
nKAPS=length(KAPSTRY);
nLamMask=length(LamMaskTRY);
nSmooth=length(SmoothTRY);

niter=ngam*nalp*nK*nbet*nKAPL*nKAPS*nLamMask*nSmooth;


EMIN=100000000000*mean(S.^2)*ones(nK,1);

for ik=1:nK
	k=K(ik);
	tbet=TLOCK(ik,:);
	EndTime=ENDTIME(ik,:);
	M0=MFULL(ik);
	X=XCAL(ik);
	for igam=1:ngam
		gam=GAMTRY(igam);
		alp=1-gam;

		for ibet=1:nbet
			bet=BETTRY(ibet);

			for iKAPL=1:nKAPL
				KAPL=KAPLTRY(iKAPL);
				
				for iKAPS=1:nKAPS
					KAPS=KAPSTRY(iKAPS);

					for iLamMask=1:nLamMask
						LamMask=LamMaskTRY(iLamMask);
                        
                        for iSmooth=1:nSmooth
                            Smooth=SmoothTRY(iSmooth);

						Ecur=ComputeE(bet,tbet,gam,alp,M0,k,S,X,KAPL,KAPS,Smooth,Zeta,EndTime,LamMask,LockRel);
						if(Ecur<EMIN(ik))
							EMIN(ik)=Ecur;
							BETmin(ik,:)=bet;
							GAMmin(ik)=gam;ALPmin(ik)=alp;
							MFULLmin(ik)=M0;
							KAPLmin(ik)=KAPL;
							KAPSmin(ik)=KAPS;
                            Smoothmin(ik)=Smooth;
							LamMaskmin(ik)=LamMask;
							Etarg=min(Etarg,Ecur/(1-delE));
						end
						if(Ecur<=Etarg)
							fprintf(fout,'%3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d',[bet,tbet,EndTime,KAPL,KAPS,LamMask,Zeta,LockRel,Smooth]);
							fprintf(fout,' %3d %3d %3d %3d %3d %3d\n',[gam,alp,M0,k,X,Ecur]);
                        end
                        end%iSmooth
					end%iLamMask
				end%iKAPS
			end%iKAPL
		end%ibet
	end%igam
end%ik



fclose(fout);

