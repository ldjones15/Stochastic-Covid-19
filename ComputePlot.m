function [plot_time,kapl_plot,kaps_plot,LAM_plot]=ComputePlot(BET,tbet,gam,alp,M0,k,T,X,KAPL,KAPS,Smooth,Zeta,EndTime,LamMask,LockRel)
%tbet specifies the start of SD and LockDown.
%tbet and EndTime actions have effects k time steps later.

Gth=zeros(T,1); 
plot_time=1:T;
kapl_plot=zeros(T-k-1,1);
kaps_plot=zeros(T-k-1,1);
LAM_plot=zeros(T-k-1,1);
%gam_plot=zeros(T-K-1,1);


%Start time of second lockdown
tbet2 = EndTime + 0*ones(3,1); %Assuming it starts 90 days after the first lockdown
%Length of second lockdown
LOCKlength2 = 0*ones(3,1); %Assuming our new Lockdown is also 90 days long
%End Time of second Lockdown
EndTime2 = tbet2 + LOCKlength2;
%AddTime2 =[0,0,0] %Option to add extra time to any parameter



for t=3:length(Gth)
        a1=tbet(1);a2=tbet(2);a3=tbet(3); % first switchpoint
        b1=EndTime(1);b2=EndTime(2);b3=EndTime(3); %Second switchpoint
        c = Smooth;
        smooth_kapl_1 = (1+1*tanh(((t-k-1)-(a1))*(c)))/2;
        smooth_kapl_2 = (1+1*tanh(((t-k-1)-(b1))*(c)))/2;
        smooth_kaps_1 = (1+1*tanh(((t-k-1)-(a2))*(c)))/2;
        smooth_kaps_2 = (1+1*tanh(((t-k-1)-(b2))*(c)))/2;
        smooth_bet_1 = (1+1*tanh(((t-k-1)-(a3))*(c)))/2;
        smooth_bet_2 = (1+1*tanh(((t-k-1)-(b3))*(c)))/2;
        kapL=KAPL; 
        kapS=KAPS; 
        bet=BET; 
        %kapL=KAPL;kapS=KAPS;bet=BET;
		if(t-k-1<=tbet(2))
			kapS=0;
        end
        if(t-k-1>tbet(2) && t-k-1<=EndTime(2))
            %kapS = kapS*smooth_kaps_1;
            kapS = kapS*smooth_kaps_1;
            %kapS = kapS;
        end
		if(t-k-1>EndTime(2))
			kapS=(LockRel*kapS*smooth_kaps_2)+kapS*(1-smooth_kaps_2);
            %kapS = kapS + smooth_kaps_2*(LockRel*kapS - kapS);
            %kapS = LockRel * kapS;
		end
		if(t-k-1<=tbet(1))
			kapL=0;
        end
        if(t-k-1>tbet(1) && t-k-1<=EndTime(1))
            %kapL = kapL*smooth_kapl_1;
            kapL =  kapL*smooth_kapl_1;
            %kapL = kapL;
        end
		if(t-k-1>EndTime(1))
			kapL=(LockRel*kapL*smooth_kapl_2)+kapL*(1-smooth_kapl_2);
            %kapL = kapL + smooth_kapl_2*(LockRel*kapL - kapL);
            %kapL = LockRel*kapL;
        end
        %if(tbet(1)==EndTime(1))
            %kapL=0.0*kapL; %This is for sweden
        %end
        if(t-k-1<=tbet(3))
			bet=bet;
        end
		if(t-k-1>tbet(3)&&t-k-1<=EndTime(3))
			bet=(LamMask*bet*smooth_bet_1)+bet*(1-smooth_bet_1);
            %bet = bet + smooth_bet_1*(LamMask*bet-bet);
            %bet = LamMask*bet;
        end
        if(t-k-1>EndTime(3))
            bet=((LamMask + LockRel*LamMask)*bet*smooth_bet_2)+(LamMask*bet)*(1-smooth_bet_2);
            %bet = bet + smooth_bet_2*(bet - LamMask*bet);
        end
        
        %kapL = 0.5*kapL; 
        %kapS=0.5*kapS;
        %kapL = kapL*tanh(kapL*(t-k-1))+kapL*(1-tanh(kapL*(t-k-1)));
        %kapS = kapS*tanh(kapS*(t-k-1))+kapS*(1-tanh(kapS*(t-k-1)));
        %bet = bet*tanh(bet*(t-k-1))+bet*(1-tanh(bet*(t-k-1)));
        
        kapl_plot(t)=kapL;
        kaps_plot(t)=kapS;
        LAM_plot(t)=bet;
        %These are the plot parameters to be called in ParamPlot
	end

end