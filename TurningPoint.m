function [tp]=TurningPoint(I,tp_robust)

DI=[0;diff(I)];

DI=tanh(DI/50);

jpos=find(DI>=0);
Spos=sum(DI(jpos));

Eopt=Spos;
E=Spos;
tp=0;
for i=1:length(DI)
	E=E-DI(i);
	if(E<=(1+tp_robust)*Eopt)
		tp=i;
	end
	if(E<=Eopt)
		Eopt=E;
	end
end
