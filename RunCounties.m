clear
warning off
format short g
tic

lstr='C:\Users\Liam\Documents\MATLAB\COVID\DATA\JohnsHopkins\AllRegions.csv';
DATA=dlmread(lstr,',',1,2);
lstr='C:\Users\Liam\Documents\MATLAB\COVID\DATA\JohnsHopkins\RegionQuarItol.csv';
QPAR=dlmread(lstr,',',1,1);QPAR(3435,1)=0.5;%sweden
[nCounty,npar]=size(QPAR);
TotalInfThresh=500;
%LOCKlength=54; %Length of the lockdown
%39 days for NC, 54 for NY,
Tfull=fliplr(DATA(1,2:end));
delE=0.05;
bsmooth=2;bsmooth1=1;
fsmooth=2;fsmooth1=1;

%for county=[805] %Indiana All
%for county=[3306] %Canada
%for county=[3451] %United Kingdom
%for county=[2025] %North Carolina All
%for county=[3435] %Sweden (Never went full lockdown)
%for county=[1893]
%for county=[1924] %New York - All
%for county=[3357]
for county=[397] %Florida All
%for county=[557]
%for county=[3321] %Denmark
%for county=[1:nCounty]

%county

Quarantine=QPAR(county,1);


pop=DATA(county+1,1);
Sfull=fliplr(DATA(county+1,2:end));
jS=find(Sfull>0);
S=reshape(Sfull(jS),length(Sfull(jS)),1);

if(length(S)>0&S(end)>TotalInfThresh)

I=[S(1);diff(S)];

Scal=movemean1(cummax(S),[bsmooth,fsmooth]);
Scal1=movemean1(cummax(S),[bsmooth1,fsmooth1]);


%%%%Parameters for Initialize%%%%
cmaxTP=4;
cmax=3;
niterC=1000;
pval=0.99;
rthresh=1;
tp_robust=0.05;if(county==1921)tp_robust=0.2;end%Westchester
threshneg=0.9;if(county==3278)threshneg=0.965;end%IRAN
minfit=14;
if(county==1876)minfit=9;end%Duchess
if(county==3278)minfit=13;end%Iran
if(county==nCounty)minfit=20;end%World
X0=pop;
K0=(8:2:8)'; % This means the incubation period is 8 days (8:2:8)
KAPL=0.2;
%KAPL=0; %SEIR Model test
KAPS=0;
Smooth=0.33; %About 1 month
Zeta=2;
LockRel=0.5; 
%if(county==3435)LockRel=0.8;end %Willing Compliance of Swedish Population without Mandates
LamMask=0.6;
%if(county==3435)LamMask=0.85;end %Swedish very lax mask regulations, 
BETTRY=0.1:0.04:0.99;
GAMTRY=[0.005:0.005:0.015,0.02:0.04:0.3];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[BET,TBET,GAM0,ALP0,MFULL,ETARG,tp,ifit,TBETfull]=Initialize2(Scal,K0,X0,minfit,cmax,cmaxTP,niterC,pval,rthresh,tp_robust,threshneg,KAPL,KAPS,Smooth,Zeta,LamMask,GAMTRY,BETTRY,LockRel);
%Sfit=Scal(1:ifit);
Sfit=Scal1(1:ifit);
%New York
%LOCKlength=(ifit-(TBET(1))); %Lets code pick the length of Lockdown
LOCKlength=40;
EXT = 80; %extension of other measures 81+54=135
%I should use actual lockdown statistics
%39 days for NC, 54 for NY, Florida is 31

[county,Quarantine,toc]



%%%%%%Parameters for Calibrate%%%%
sstr='CalibrateTmp.dat';
ALPTRY=[0.005:0.005:0.015,0.02:0.04:0.3];
KAPLTRY=[0.00000001,0.001:0.001:0.005,0.01,0.05:0.05:0.3];
KAPSTRY=[0.00000001,0.001:0.001:0.005,0.01,0.05:0.05:0.2];
%KAPLTRY=[0]; %SEIR Model Test
%KAPSTRY=[0]; %SEIR Model Test
SmoothTRY=[0.01:0.01:1]; %Use if you wish to calculate per region
%SmoothTRY=[0.01:0.05:0.5]; %Faster run time
%SmoothTRY=[0.15:0.01:0.45]; %This way Smooth isn't determined by the machine
LamMaskTRY=[0.3:0.1:0.7];
%LamMaskTRY=[0.5:0.1:0.7];
%LamMaskTRY=[0.75:0.75:0.75];
%LamMaskTRY = 0.6;
%LamMaskTRY=1; %SEIR Model Test
Etarg=min(ETARG)/(1-delE);
XCAL=X0*ones(size(K0));
TLOCK=TBET(:,1)*ones(1,3); %TLOCK is time of lockdown
%if(county==3435)TLOCK(1)=0;end %For Sweden
ENDTIME=TLOCK+LOCKlength;
ENDTIME(2) = ENDTIME(2)+EXT;
ENDTIME(3) = ENDTIME(3)+EXT;
%if(county==3435)ENDTIME(2)=0;end %Sweden never locked down
%ENDTIME(:,end)=length(Sfit)+K0; %Extends Mask provisions past lockdown
%ENDTIME=TLOCK+100;ENDTIME(:,:)=length(Sfit)+K0;
%The Start time, Endtime, and length for a second implementation of 
%preventative measures can be found in Compute.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%Final round optimization -- fine%%%%%%%%

Etarg=min(ETARG)/(1-delE);
sstr=['LangevinAuto-20200618/CountyData/County',int2str(county),'.dat'];


[EMIN,BETmin,GAMmin,ALPmin,MFULLmin,KAPLmin,KAPSmin,LamMaskmin,Smoothmin]=Calibrate1(BETTRY,GAMTRY,ALPTRY,MFULL,K0,XCAL,KAPLTRY,KAPSTRY,LamMaskTRY,SmoothTRY,TLOCK,ENDTIME,Sfit,Etarg,sstr,delE,Zeta,LockRel);




Etarg=min(EMIN)/(1-delE);
BETTRY=0.95*BETmin:0.1*BETmin/10:1.05*BETmin;
GAMTRY=0.95*GAMmin:0.1*GAMmin/6:1.05*GAMmin;
KAPLTRY=0.95*KAPLmin:0.1*KAPLmin/6:1.05*KAPLmin;
KAPSTRY=0.95*KAPSmin:0.1*KAPSmin/6:1.05*KAPSmin;
LamMaskTRY=0.95*LamMaskmin:0.1*LamMaskmin/6:1.05*LamMaskmin;
SmoothTRY=0.95*Smoothmin:0.1*Smoothmin/6:1.05*Smoothmin;

[EMIN,BETmin,GAMmin,ALPmin,MFULLmin,KAPLmin,KAPSmin,LamMaskmin,Smoothmin]=Calibrate1(BETTRY,GAMTRY,ALPTRY,MFULL,K0,XCAL,KAPLTRY,KAPSTRY,LamMaskTRY,SmoothTRY,TLOCK,ENDTIME,Sfit,Etarg,sstr,delE,Zeta,LockRel);


Etarg=min(EMIN)/(1-delE);
BETTRY=0.95*BETmin:0.1*BETmin/10:1.05*BETmin;
GAMTRY=0.95*GAMmin:0.1*GAMmin/6:1.05*GAMmin;
KAPLTRY=0.95*KAPLmin:0.1*KAPLmin/6:1.05*KAPLmin;
KAPSTRY=0.95*KAPSmin:0.1*KAPSmin/6:1.05*KAPSmin;
LamMaskTRY=0.95*LamMaskmin:0.1*LamMaskmin/6:1.05*LamMaskmin;
SmoothTRY=0.95*Smoothmin:0.1*Smoothmin/6:1.05*Smoothmin;

[EMIN,BETmin,GAMmin,ALPmin,MFULLmin,KAPLmin,KAPSmin,LamMaskmin,Smoothmin]=Calibrate1(BETTRY,GAMTRY,ALPTRY,MFULL,K0,XCAL,KAPLTRY,KAPSTRY,LamMaskTRY,SmoothTRY,TLOCK,ENDTIME,Sfit,Etarg,sstr,delE,Zeta,LockRel);

sstr=['save C:\Users\Liam\Documents\MATLAB\COVID\LangevinAuto-20200618\FinalParams\CountyFinalParam',int2str(county),' EMIN BETmin GAMmin ALPmin MFULLmin K0 Etarg S Scal Sfit jS pop Quarantine tp ifit TBETfull TLOCK ENDTIME KAPLmin KAPSmin LamMaskmin Smoothmin'];

eval(sstr);
 

end
end


%%%%%%%%%%%%%%%Plot Fit

[minimumE,idxk]=min(EMIN);
%ParamPlot
PlotFit
pause(0.01)




