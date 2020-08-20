clear
warning off
format short g
tic

lstr='C:\Users\Liam\Documents\MATLAB\COVID\DATA\JohnsHopkins\AllRegions.csv';
DATA=dlmread(lstr,',',1,2);
lstr='C:\Users\Liam\Documents\MATLAB\COVID\DATA\JohnsHopkins\AllRegionsDeaths.csv';
DEATHS=dlmread(lstr,',',1,2);
lstr='C:\Users\Liam\Documents\MATLAB\COVID\DATA\JohnsHopkins\RegionQuarItol.csv';
QPAR=dlmread(lstr,',',1,1);QPAR(3435,1)=0.5;%swedenx
[nCounty,npar]=size(QPAR);
TotalInfThresh=500;
AddTime=[0,90,90]; %Extra time for each preventative measure [lockdown, social distancing, masks]
%Here we add AddTime2 so we can have the code run for two separate
%instances of preventative measure time extensions
AddTime2=[0,0,0]; %Used in GetProjectionsSlim2 
PlotRes=500;


lamdata=1;
lamdaily=0.35;
lamrel=0.35;
lamopt=0.35;

day1=datenum(2020,1,22);
daytarg=datenum(2020,07,31);
dmax=300; %The amount of days our model calculates since first infection
lagA=8;
ModelDir='C:\Users\Liam\Documents\MATLAB\COVID\LangevinAuto-20200618\';sstrdir='C:\Users\Liam\Documents\MATLAB\COVID\WarRoom\FiguresLangevin\';statsdir='C:\Users\Liam\Documents\MATLAB\COVID\WarRoom\StatsLangevin\';


err=0.05;

%for cidx=[805]; %Indiana All
%for cidx=[3306]; %Canada
%for cidx=[3451]; %United Kingdom
%for cidx=[1893];
%for cidx=[2016];
%for cidx=[2025]; %North Carolina - All
%for cidx=[1924]; %New York - All
%for cidx=[3435]; %Sweden
%for cidx=[557];
%for cidx=[3357];
for cidx=[397]; %Florida All
%for cidx=[3321]; %Denmark
%for cidx=[1:nCounty];


Itol=QPAR(cidx,2);
county=cidx;
[county,Itol,toc]


Sfull=fliplr(DATA(county+1,2:end));
Dfull=fliplr(DEATHS(county+1,2:end));
jS=find(Sfull>0);
S=reshape(Sfull(jS),length(Sfull(jS)),1);
D=reshape(Dfull(jS),length(Dfull(jS)),1);

sstr=[sstrdir,'County',int2str(cidx),'Plot'];
statssstr=[statsdir,'County',int2str(cidx),'Stats.csv'];

if(length(S)>0&S(end)>TotalInfThresh)
	load C:\Users\Liam\Documents\MATLAB\COVID\QuadHospitalModel\HospParamGen
	lstr=['load ',ModelDir,'FinalParams\CountyFinalParam',int2str(county)];
	eval(lstr)
  
	IS=[S(1);diff(S)];
	TS=day1+jS(1)-1+(0:length(S)-1)';
	Tend=length(S)+dmax;
	Tth=day1+jS(1)-1+(0:Tend-1)';
	idxnow=length(S);
	daynow= Tth(idxnow); %August 11
    dayend = Tth(idxnow+142)
	targidx=find(Tth==daytarg);
    

	lstr=[ModelDir,'CountyData\County',int2str(cidx),'.dat'];
    %lstr=['CountyData/County',int2str(cidx),'.dat'];
	Stats=load(lstr);
%% Stats has 18 cols:
%% 1   2   3  4   5     6     7    8  9  10  11   12    13     14  15 16 17 18  19
%% bet tSD tL tMk EndSD EndLk EndM KL KS LAM Zeta LRel Smooth  gam alp M0 k  X  E


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Low Quarantine %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	tau=1;

   % hold on %Plotting for the first case (AddTime)
    
	[Smin,Smax,Smean,Imin,Imax,Imean,Mmin,Mmax,Mmean,Rmin,Rmax,Rmean,ImaxDist,MDist,TotIDist,plot_time,kapl_plot,kaps_plot,LAM_plot,EndTime,gam]=GetProjectionsSlim(Stats,Tth,Tend,err,pop,tau,Itol,S,targidx,AddTime); 
    
    [N,statsd]=size(Stats);
   
    for i = 1:N
        lockstart = datestr(Tth(Stats(i,3)));
        lockend = datestr(Tth(Stats(i,6)));
    end
    
	TotInf1=Smean+Mmean+Rmean;
	CurConfInf=S(end);
	CurTotInf=TotInf1(length(D));
	CurFat=D(end);
	FatRate=CurFat/CurTotInf;
	TargTotInf=ceil(TotInf1(targidx));
	TargTotFat=CurFat+ceil(FatRate*(TotInf1(targidx)-TotInf1(length(D))));
	TargConfInf=CurConfInf+round(Smean(targidx)-Smean(idxnow));


	InfStats=[CurConfInf,CurTotInf,CurFat,FatRate,TargTotInf,TargTotFat,TargConfInf,daynow,daytarg];
	jpos=find(Stats(:,end)<=(1+err)*min(Stats(:,end)));
	ModelStats=Stats(jpos,[1,8,9,10,13,14]);
	if(length(jpos>1))
		ModelStats=mean(ModelStats);
	end

	dlmwrite(statssstr,[InfStats,ModelStats],'precision','%d');

	cbar=[0.7,0.7,1];
	cgray=[0.8,0.8,0.8];
	pstr=[sstr,'Daily'];
    paramstr=[sstr,'Parameters'];
    lockdownstr=[sstr, 'Alpha'];
    lockdownstrn=[sstr, 'Alpha_Param'];
    distancestr=[sstr,'Delta'];
    betastr=[sstr,'Beta'];
    gammastr=[sstr,'Gamma'];
	lam=lamdaily;
	usebounds=1;
	Iplot=abs(Imean);TotI=(Mmean+Smean+Rmean)/pop;
    TotalInf=(Mmean+Smean+Rmean);
    AugInf = TotalInf(idxnow) %Total Infections on Aug 11
    DecInf = TotalInf(end) %Total Infection on Dec 31
    MeasureI = abs(Smean);
    INow=Imean(idxnow);
	InfNow=Smean(idxnow)+Rmean(idxnow)+Mmean(idxnow);
	IDistPlot=TotIDist;
	PlotDaily;
    kapl_max=max(kapl_plot);
    kaps_max=max(kaps_plot);
    lam_max=max(LAM_plot);
    ParamPlot2 %Plots parameters as 4 plots on one figure through TiledLayout 
    ParamPlot3; %Plots parameters each as a separate figure (4 plots per county)
    
    
		%Needs pstr,Tth,Imin,Imax,TS,IS,Iplot,Hplot,Hospplot,Aplot
		%CDPHP,lam,hospfrac,cblue,cgray,Qdyn,usebounds
	pause(0.001)
%	PlotInfDist
%		%Needs IDistPlot and HDistPlot
	ModelStats
	InfStats
    
    %hold off 
    
   % hold on %Plotting for the second case (AddTime2)
    
    [Smin,Smax,Smean,Imin,Imax,Imean,Mmin,Mmax,Mmean,Rmin,Rmax,Rmean,ImaxDist,MDist,TotIDist,plot_time,kapl_plot,kaps_plot,LAM_plot,EndTime,gam]=GetProjectionsSlim2(Stats,Tth,Tend,err,pop,tau,Itol,S,targidx,AddTime2); 
    
    [N,statsd]=size(Stats);
    
	TotInf1=Smean+Mmean+Rmean;
	CurConfInf=S(end);
	CurTotInf=TotInf1(length(D));
	CurFat=D(end);
	FatRate=CurFat/CurTotInf;
	TargTotInf=ceil(TotInf1(targidx));
	TargTotFat=CurFat+ceil(FatRate*(TotInf1(targidx)-TotInf1(length(D))));
	TargConfInf=CurConfInf+round(Smean(targidx)-Smean(idxnow));


	InfStats2=[CurConfInf,CurTotInf,CurFat,FatRate,TargTotInf,TargTotFat,TargConfInf,daynow,daytarg];
	jpos=find(Stats(:,end)<=(1+err)*min(Stats(:,end)));
	ModelStats2=Stats(jpos,[1,8,9,10,13,14]);
	if(length(jpos>1))
		ModelStats2=mean(ModelStats2);
	end

	dlmwrite(statssstr,[InfStats2,ModelStats2],'precision','%d');

	cbar=[0.7,0.7,1];
	cgray=[0.8,0.8,0.8];
	pstr=[sstr,'Daily_2'];
    paramstr=[sstr,'Parameters_2'];
    lockdownstr=[sstr, 'Alpha_2'];
    lockdownstrn=[sstr, 'Alpha_Param'];
    distancestr=[sstr,'Delta_2'];
    betastr=[sstr,'Beta_2'];
    gammastr=[sstr,'Gamma_2'];
	lam=lamdaily;
	usebounds=1;
	Iplot=abs(Imean);TotI=(Mmean+Smean+Rmean)/pop;
    MeasureI = abs(Smean);
	InfNow=Smean(idxnow)+Rmean(idxnow)+Mmean(idxnow);
    TotalInf=(Mmean+Smean+Rmean);
    AugInf = TotalInf(idxnow) %Total Infections on Aug 11
    DecInf = TotalInf(end) %Total Infection on Dec 31
	IDistPlot=TotIDist;
	PlotDaily;
    kapl_max=max(kapl_plot);
    kaps_max=max(kaps_plot);
    lam_max=max(LAM_plot);
    ParamPlot2 %Plots parameters as 4 plots on one figure through TiledLayout 
    ParamPlot3; %Plots parameters each as a separate figure (4 plots per county)
    
    
		%Needs pstr,Tth,Imin,Imax,TS,IS,Iplot,Hplot,Hospplot,Aplot
		%CDPHP,lam,hospfrac,cblue,cgray,Qdyn,usebounds
	pause(0.001)
%	PlotInfDist
%		%Needs IDistPlot and HDistPlot
	ModelStats2
	InfStats2
    
    %hold off
    
      
end
end

