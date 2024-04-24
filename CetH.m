function[h]=CetH(vec)
Nt=32;
Nr=128;
L=4;
Ar=zeros(Nr,L);
At=zeros(Nt,L);
lambda=1e-3;
d=lambda/2;
aR_vector=zeros(Nr,1);
aT_vector=zeros(Nt,1);
AoA=vec(1:L);
AoD=vec(L+1:2*L);
            for A_j=1:L
                for A_k=1:Nr
                    aR_vector(A_k,1)=(1/sqrt(Nr))*exp((-1i)*2*pi/lambda*(A_k-1)*d*sin(AoA(A_j)));
                end
                Ar(:,A_j)=aR_vector;
                for A_x=1:Nt
                    aT_vector(A_x,1)=(1/sqrt(Nt))*exp((-1i)*2*pi/lambda*(A_x-1)*d*sin(AoD(A_j)));
                end
                At(:,A_j)=aT_vector;
            end
H=sqrt(Nr*Nt/L)*Ar*At';
h=H(:);