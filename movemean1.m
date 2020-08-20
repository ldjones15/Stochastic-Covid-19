function out=movmean1(S,smoothwindow);
bsmooth=smoothwindow(1);
fsmooth=smoothwindow(2);

%first deal with the noisy zero infections days that are likely non-reports.
for i=3:length(S)-2
	if(S(i)==0)
		xlocal=S(i-2:i+2);
		jpos=find(xlocal>0);
		S(i)=min(xlocal(jpos));
	end
end

out=movmean(S,[bsmooth,fsmooth]);
for i=1:bsmooth
	out(i)=mean(S(1:2*i-1));
end
for i=1:fsmooth
	out(end-i+1)=mean(S(end-2*i+2:end));
end
