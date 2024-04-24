function[peak2]=FDPK(CSY)
[peaka,peakb]=max(abs(CSY));
CSY(peakb)=0;
[a,b]=max(abs(CSY));
Pa=zeros(2,1);
Pb=zeros(2,1);
countp=1;
pkf=0.75;
while (abs(b-peakb)<4)||(abs(b-peakb)>7)
    if abs(b-peakb)==1
        Pa(countp)=a;
        Pb(countp)=b;
        countp=countp+1;
    end
    CSY(b)=0;
    [a,b]=max(abs(CSY));
    if a==0
        break
    end
end
if a>pkf*peaka
    peak2=b;
else
    peak2=peakb;
    pkc=0.98;
    if max(Pa)>pkc*peaka
        [~,b2]=max(Pa);
        peak2=Pb(b2);
    end
end
end