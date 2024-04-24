close all
clear all
%par.runld=1;
%rng(par.runld);
Tmax=15;
TMAX=10;
Ncan=128;
v_i=2;
v=[4,6,8];
sigma=[1,sqrt(0.1),0.1];
Nt=64;
Nr=64;
Nrf=4;
N=9;
M=5;
lambda=1e-3;
d=lambda/2;
L=4;
aru=zeros(L*N,1);
atu=zeros(L*N,1);
a=[1,1;0,1];
A=diag(diag(repmat(a,2*L,2*L)))+diag(diag(repmat(a,2*L,2*L),1),1);
X_r=zeros(2*L,1);
Xko=zeros(2*L,3);
rXko=zeros(2*L,3);
%_________________________________________________________________________
G=256;
Gw=13;
AoA=zeros(L,1);
deta_AoA=zeros(L,1);
AoD=zeros(L,1);
deta_AoD=zeros(L,1);
Acan_W_vector=zeros(Nr,1);
Acan_F_vector=zeros(Nt,1);
Ar=zeros(Nr,L);
At=zeros(Nt,L);
H=zeros(Nr,Nt);
W_ko=zeros(Nr,Nrf);
F_ko=zeros(Nt,Nrf);
W_ez=zeros(Nr,Nrf);
F_ez=zeros(Nt,Nrf);
a_W=[];
a_F=[];
Ya=zeros(L,1);
Xod=zeros(L,1);
dta=zeros(L-1,1);
dtd=zeros(L-1,1);
ARL1=zeros(Nrf,1);
ARL2=zeros(Nrf,1);
ARL3=zeros(Nrf,1);
ATL1=zeros(Nrf,1);
ATL2=zeros(Nrf,1);
ATL3=zeros(Nrf,1);
Acan_W=zeros(Nr,Ncan);
Acan_F=zeros(Nt,Ncan);
for Acan_j=1:Ncan
    psi=-pi/2+pi*(Acan_j-1)/Ncan;
    a_W(Acan_j)=psi;
    for Acan_k=1:Nr
        Acan_W_vector(Acan_k,1)=exp((-1i)*2*pi/lambda*(Acan_k-1)*d*sin(psi));
    end
    Acan_W(:,Acan_j)=Acan_W_vector;
end
for Acan_j=1:Ncan
    psi=-pi/2+pi*(Acan_j-1)/Ncan;
    a_F(Acan_j)=psi;
    for Acan_k=1:Nt
        Acan_F_vector(Acan_k,1)=exp((-1i)*2*pi/lambda*(Acan_k-1)*d*sin(psi));
    end
    Acan_F(:,Acan_j)=Acan_F_vector;
end
ARb=zeros(Nr,G);
ATb=zeros(Nt,G);
aR_v=zeros(Nr,1);
aT_v=zeros(Nt,1);
for g=1:G
    for A_k=1:Nr
        aR_v(A_k,1)=(1/sqrt(Nr))*exp((-1i)*pi*(2/G*(g-1)-1)*(A_k-1));
    end
    ARb(:,g)=aR_v;
    
    for A_x=1:Nt
        aT_v(A_x,1)=(1/sqrt(Nt))*exp((1i)*pi*(2/G*(g-1)-1)*(A_x-1));
    end
    ATb(:,g)=aT_v;
end
C=kron(ATb,ARb);
out=zeros(Tmax,3);
out1=zeros(TMAX,Tmax);
out2=zeros(TMAX,Tmax);
out3=zeros(TMAX,Tmax);
E1=zeros(Tmax,TMAX);
E2=zeros(Tmax,TMAX);
E3=zeros(Tmax,TMAX);
delta_A=v(v_i)*pi/180;
%G_wide=ceil(sin(delta_A)/2*G);
% Gw1=ceil(1.6*delta_A/pi*G);
Gw1=14;
Gw2=7;
for times=1:TMAX;
    T_AOA=zeros(Tmax,L);
    T_AOD=zeros(Tmax,L);
    P_AOA1=zeros(Tmax,L);
    P_AOA2=zeros(Tmax,L);
    P_AOA3=zeros(Tmax,L);
    P_AOD1=zeros(Tmax,L);
    P_AOD2=zeros(Tmax,L);
    P_AOD3=zeros(Tmax,L);
    D=1;%diag([0.8893,0.0953,0.0107,0.0047]);
    a=1;
    AoA(a)=0.85*pi*(rand-1/2);
    if AoA(a)>0
        deta_AoA(a)=-delta_A*rand;
    else
        deta_AoA(a)=delta_A*rand;
    end
    while a<4
        x=0.85*pi*(rand-1/2);
        if min(abs(x-AoA(1:a)))<pi/8
        else
            a=a+1;
            AoA(a)=x;
            if AoA(a)>0
                deta_AoA(a)=-delta_A*rand;
            else
                deta_AoA(a)=delta_A*rand;
            end
        end
    end
    a=1;
    AoD(a)=0.85*pi*(rand-1/2);
    if AoD(a)>0
        deta_AoD(a)=-delta_A*rand;
    else
        deta_AoD(a)=delta_A*rand;
    end
    while a<4
        x=0.85*pi*(rand-1/2);
        if min(abs(x-AoD(1:a)))<pi/8
        else
            a=a+1;
            AoD(a)=x;
            if AoD(a)>0
                deta_AoD(a)=-delta_A*rand;
            else
                deta_AoD(a)=delta_A*rand;
            end
        end
    end
    [AoA,aoao]=sort(AoA);
    AoD=AoD(aoao);
    deta_AoA=deta_AoA(aoao);
    deta_AoD=deta_AoD(aoao);
    for n=1:Tmax
        %Hn
        T_AOA(n,:)=AoA;
        T_AOD(n,:)=AoD;
        X_r(1:L)=AoA;
        X_r(L+1:2*L)=AoD;
        h_r=CetH(X_r);
        H=reshape(h_r,Nr,Nt);
        AoA=AoA+deta_AoA;
        AoD=AoD+deta_AoD;
        %---sigma1
        if n==1
            Nn_ko=sigma(1)/sqrt(2)*(randn(Nr,Ncan)+1i*randn(Nr,Ncan));
            bW_ko=Acan_W;
            bF_ko=Acan_F;
            xn=bW_ko'*H*bF_ko+bW_ko'*Nn_ko;
            y=xn(:);
            Q1=kron(bF_ko.',bW_ko')*C;
            CS=reshape(Q1'*y,G,G);
            for cs_i=1:L
                [arl,atl]=find(abs(CS)==max(max(abs(CS))));
                arl1=FDPK(CS(:,atl));
                ARL1(cs_i)=round((arl+arl1)/2);
                oA=(asin(2/G*(arl-1)-1)+asin(2/G*(arl1-1)-1))/2;
                P_AOA1(n,cs_i)=oA;
                [~,b_W]=min(abs(a_W-oA));
                W_ko(:,cs_i)=Acan_W(:,b_W);
                %F_omp
                atl1=FDPK(CS(ARL1(cs_i),:));
                ATL1(cs_i)=round((atl+atl1)/2);                
                oD=(asin(2/G*(atl-1)-1)+asin(2/G*(atl1-1)-1))/2;
                P_AOD1(n,cs_i)=oD;
                [~,b_F]=min(abs(a_F-oD));
                F_ko(:,cs_i)=Acan_F(:,b_F);
                if arl-Gw<1
                    CS(1:arl+Gw,:)=0;
                    CS(G+arl-Gw:G,:)=0;
                else
                    if arl+Gw>G
                        CS(arl-Gw:G,:)=0;
                        CS(1:arl-G+Gw,:)=0;
                    else
                        CS(arl-Gw:arl+Gw,:)=0;
                    end
                end
                if atl-Gw<1
                    CS(:,1:atl+Gw)=0;
                    CS(:,G-Gw+atl:G)=0;
                else
                    if atl+Gw>G
                        CS(:,atl-Gw:G)=0;
                        CS(:,1:atl+Gw-G)=0;
                    else
                        CS(:,atl-Gw:atl+Gw)=0;
                    end
                end
            end
            [P_AOA1(1,:),oo]=sort(P_AOA1(1,:));
            OD1=P_AOD1(1,:);
            OD1=OD1(oo);
            ARL1=ARL1(oo);
            P_AOD1(1,:)=OD1;
            ATL1=ATL1(oo);
        end
        if n==2
            Nn_ko=sigma(1)/sqrt(2)*(randn(Nr,Nrf*N)+1i*randn(Nr,Nrf*N));
            bW_ko=zeros(Nr,Nrf*N);
            bF_ko=zeros(Nt,Nrf*N);
            aro=zeros(Nrf*N,1);
            ato=zeros(Nrf*N,1);
            %W_omp
            for i=1:Nrf*N
                aro(i)=P_AOA1(n-1,ceil(i/N))+(rem(i,N)-(N-1)/2)*pi/Ncan;
                [~,b]=min(abs(a_W-aro(i)));
                bW_ko(:,i)=Acan_W(:,b);
                ato(i)=P_AOD1(n-1,ceil(i/N))+(rem(i,N)-(N-1)/2)*pi/Ncan;
                [~,b]=min(abs(a_F-ato(i)));
                bF_ko(:,i)=Acan_F(:,b);
            end
            xn=bW_ko'*H*bF_ko+bW_ko'*Nn_ko;
            y=xn(:);
            Q1=kron(bF_ko.',bW_ko')*C;
            CS=reshape(Q1'*y,G,G);
            CSS=CS;
            for cs_i=1:L
                CS=CSS;
                CS(ARL1(cs_i)+Gw1+1:G,:)=0;
                CS(1:ARL1(cs_i)-Gw1-1,:)=0;
                CS(:,ATL1(cs_i)+Gw1+1:G)=0;
                CS(:,1:ATL1(cs_i)-Gw1-1)=0;
                [arl,atl]=find(abs(CS)==max(max(abs(CS))));
                %W_omp
                arl1=FDPK(CS(:,atl));
                ARL1(cs_i)=round((arl+arl1)/2);
                oA=(asin(2/G*(arl-1)-1)+asin(2/G*(arl1-1)-1))/2;
                P_AOA1(n,cs_i)=oA;
                [~,b_W]=min(abs(a_W-oA));
                W_ko(:,cs_i)=Acan_W(:,b_W);
                %F_omp
                atl1=FDPK(CS(ARL1(cs_i),:));
                ATL1(cs_i)=round((atl+atl1)/2);
                oD=(asin(2/G*(atl-1)-1)+asin(2/G*(atl1-1)-1))/2;
                P_AOD1(n,cs_i)=oD;
                [~,b_F]=min(abs(a_F-oD));
                F_ko(:,cs_i)=Acan_F(:,b_F);
            end
            Rko1=zeros(2*L,2*L);
            for i=1:L
                Xko(i,1)=P_AOA1(n,i)-P_AOA1(n-1,i);
                Xko(i+L,1)=P_AOD1(n,i)-P_AOD1(n-1,i);
                Rko1(i)=(Xko(i,1)-deta_AoA(i))^2;
                Rko1(i+L)=(Xko(i+L,1)-deta_AoD(i))^2;
            end
        end
        if n>2
            Nn_ko=sigma(1)/sqrt(2)*(randn(Nr,Nrf*M)+1i*randn(Nr,Nrf*M));
            pXko=Xko(:,1);
            pRko=Rko1;
            Kt=Rko1*(pRko+(pi/Ncan)^2*eye(2*L,2*L))^-1;
            bW_ko=zeros(Nr,Nrf*M);
            bF_ko=zeros(Nt,Nrf*M);
            aro=zeros(Nrf*M,1);
            ato=zeros(Nrf*M,1);
            Yoa=P_AOA1(n-1,:)+(pXko(1:L))';
            Yod=P_AOD1(n-1,:)+(pXko(L+1:2*L))';
            for i=1:L
                ARL1(i)=round((sin(Yoa(i))+1)*G/2);
                ATL1(i)=round((sin(Yod(i))+1)*G/2);
            end
            [Ya,Yb]=sort(Yoa);
            sma=0;
            clos=0;
            for i=1:L-1
                dta(i)=abs(Ya(i)-Ya(i+1));
            end
            if min(dta)<5*pi/Ncan
                [~,bta]=min(dta);
                sma=[Yb(bta),Yb(bta+1)];
                if abs(Yod(Yb(bta))-Yod(Yb(bta+1)))<5*pi/Ncan
                    clos=sma;
                end
            end
            %W_omp
            for i=1:Nrf*M
                aro(i)=Yoa(ceil(i/M))+(rem(i,M)-(M-1)/2)*pi/Ncan;
                if abs(aro(i))>pi/2
                    aro(i)=0;
                end
                [~,b]=min(abs(a_W-aro(i)));
                bW_ko(:,i)=Acan_W(:,b);
                ato(i)=Yod(ceil(i/M))+(rem(i,M)-(M-1)/2)*pi/Ncan;
                if abs(ato(i))>pi/2
                    ato(i)=0;
                end
                [~,b]=min(abs(a_F-ato(i)));
                bF_ko(:,i)=Acan_F(:,b);
            end
            xn=bW_ko'*H*bF_ko+bW_ko'*Nn_ko;
            y=xn(:);
            Q1=kron(bF_ko.',bW_ko')*C;
            CS=reshape(Q1'*y,G,G);
            CSS=CS;
            for cs_i=1:L
                if  ismember(cs_i,clos)>0
                    P_AOA1(n,cs_i)=Yoa(cs_i);
                    [~,b_w]=min(abs(a_W-Yoa(cs_i)));
                    W_ko(:,cs_i)=Acan_W(:,b_w);
                    ARL1(cs_i)=ceil((sin(Yoa(cs_i))+1)*G/2);
                    P_AOD1(n,cs_i)=Yod(cs_i);
                    [~,b_F]=min(abs(a_F-Yod(cs_i)));
                    F_ko(:,cs_i)=Acan_F(:,b_F);
                    ATL1(cs_i)=ceil((sin(Yod(cs_i))+1)*G/2);
                else
                    CS=CSS;
                    CS(ARL1(cs_i)+Gw2+1:G,:)=0;
                    CS(1:ARL1(cs_i)-Gw2-1,:)=0;
                    CS(:,ATL1(cs_i)+Gw2+1:G)=0;
                    CS(:,1:ATL1(cs_i)-Gw2-1)=0;
                    CS(:,ATL1(1:cs_i-1))=0;
                    CS(ARL1(1:cs_i-1),:);
                    [arl,atl]=find(abs(CS)==max(max(abs(CS))));
                    arl1=FDPK(CS(:,atl));
                    ARL1(cs_i)=round((arl+arl1)/2);
                    oA=(asin(2/G*(arl-1)-1)+asin(2/G*(arl1-1)-1))/2;
                    P_AOA1(n,cs_i)=oA;
                    [~,b_W]=min(abs(a_W-oA));
                    W_ko(:,cs_i)=Acan_W(:,b_W);
                    atl1=FDPK(CS(ARL1(cs_i),:));
                    ATL1(cs_i)=round((atl+atl1)/2);
                    oD=(asin(2/G*(atl-1)-1)+asin(2/G*(atl1-1)-1))/2;
                    P_AOD1(n,cs_i)=oD;
                    [~,b_F]=min(abs(a_F-oD));
                    F_ko(:,cs_i)=Acan_F(:,b_F);
                end
            end
            rXko(1:L,1)=P_AOA1(n,:)-P_AOA1(n-1,:);
            rXko(L+1:2*L,1)=P_AOD1(n,:)-P_AOD1(n-1,:);
            Xko(:,1)=pXko+Kt*(rXko(:,1)-pXko);
            Rko1=pRko-Kt*pRko+(pi/(2*G))^2*eye(2*L,2*L);
        end
        bHeff=W_ko'*H*F_ko;
        [~,A_ko,~]=svd(bHeff);
        Rn=0;
        for i=1:Nrf
            Rn=Rn+log2(1+(A_ko(i,i)^2/(sigma(1)^2)));
        end
        out1(times,n)=real(Rn);
        %---sigma2
        if n==1
            Nn_ko=sigma(2)/sqrt(2)*(randn(Nr,Ncan)+1i*randn(Nr,Ncan));
            bW_ko=Acan_W;
            bF_ko=Acan_F;
            xn=bW_ko'*H*bF_ko+bW_ko'*Nn_ko;
            y=xn(:);
            Q1=kron(bF_ko.',bW_ko')*C;
            CS=reshape(Q1'*y,G,G);
            for cs_i=1:L
                [arl,atl]=find(abs(CS)==max(max(abs(CS))));
                arl1=FDPK(CS(:,atl));
                ARL2(cs_i)=round((arl+arl1)/2);
                oA=(asin(2/G*(arl-1)-1)+asin(2/G*(arl1-1)-1))/2;
                P_AOA2(n,cs_i)=oA;
                [~,b_W]=min(abs(a_W-oA));
                W_ko(:,cs_i)=Acan_W(:,b_W);
                %F_omp
                atl1=FDPK(CS(ARL2(cs_i),:));
                ATL2(cs_i)=round((atl+atl1)/2);
                oD=(asin(2/G*(atl-1)-1)+asin(2/G*(atl1-1)-1))/2;
                P_AOD2(n,cs_i)=oD;
                [~,b_F]=min(abs(a_F-oD));
                F_ko(:,cs_i)=Acan_F(:,b_F);
                if arl-Gw<1
                    CS(1:arl+Gw,:)=0;
                    CS(G+arl-Gw:G,:)=0;
                else
                    if arl+Gw>G
                        CS(arl-Gw:G,:)=0;
                        CS(1:arl-G+Gw,:)=0;
                    else
                        CS(arl-Gw:arl+Gw,:)=0;
                    end
                end
                if atl-Gw<1
                    CS(:,1:atl+Gw)=0;
                    CS(:,G-Gw+atl:G)=0;
                else
                    if atl+Gw>G
                        CS(:,atl-Gw:G)=0;
                        CS(:,1:atl+Gw-G)=0;
                    else
                        CS(:,atl-Gw:atl+Gw)=0;
                    end
                end
            end
            [P_AOA2(1,:),oo]=sort(P_AOA2(1,:));
            OD2=P_AOD2(1,:);
            OD2=OD2(oo);
            ARL2=ARL2(oo);
            P_AOD2(1,:)=OD2;
            ATL2=ATL2(oo);
        end
        if n==2
            Nn_ko=sigma(2)/sqrt(2)*(randn(Nr,Nrf*N)+1i*randn(Nr,Nrf*N));
            bW_ko=zeros(Nr,Nrf*N);
            bF_ko=zeros(Nt,Nrf*N);
            aro=zeros(Nrf*N,1);
            ato=zeros(Nrf*N,1);
            %W_omp
            for i=1:Nrf*N
                aro(i)=P_AOA2(n-1,ceil(i/N))+(rem(i,N)-(N-1)/2)*pi/Ncan;
                [~,b]=min(abs(a_W-aro(i)));
                bW_ko(:,i)=Acan_W(:,b);
                ato(i)=P_AOD2(n-1,ceil(i/N))+(rem(i,N)-(N-1)/2)*pi/Ncan;
                [~,b]=min(abs(a_F-ato(i)));
                bF_ko(:,i)=Acan_F(:,b);
            end
            xn=bW_ko'*H*bF_ko+bW_ko'*Nn_ko;
            y=xn(:);
            Q1=kron(bF_ko.',bW_ko')*C;
            CS=reshape(Q1'*y,G,G);
            CSS=CS;
            for cs_i=1:L
                CS=CSS;
                CS(ARL2(cs_i)+Gw1+1:G,:)=0;
                CS(1:ARL2(cs_i)-Gw1-1,:)=0;
                CS(:,ATL2(cs_i)+Gw1+1:G)=0;
                CS(:,1:ATL2(cs_i)-Gw1-1)=0;
                [arl,atl]=find(abs(CS)==max(max(abs(CS))));
                %W_omp
                arl1=FDPK(CS(:,atl));
                ARL2(cs_i)=round((arl+arl1)/2);
                oA=(asin(2/G*(arl-1)-1)+asin(2/G*(arl1-1)-1))/2;
                P_AOA2(n,cs_i)=oA;
                [~,b_W]=min(abs(a_W-oA));
                W_ko(:,cs_i)=Acan_W(:,b_W);
                %F_omp
                atl1=FDPK(CS(ARL2(cs_i),:));
                ATL2(cs_i)=round((atl+atl1)/2);
                oD=(asin(2/G*(atl-1)-1)+asin(2/G*(atl1-1)-1))/2;
                P_AOD2(n,cs_i)=oD;
                [~,b_F]=min(abs(a_F-oD));
                F_ko(:,cs_i)=Acan_F(:,b_F);
            end
            Rko2=zeros(2*L,2*L);
            for i=1:L
                Xko(i,2)=P_AOA2(n,i)-P_AOA2(n-1,i);
                Xko(i+L,2)=P_AOD2(n,i)-P_AOD2(n-1,i);
                Rko2(i)=(Xko(i,2)-deta_AoA(i))^2;
                Rko2(i+L)=(Xko(i+L,2)-deta_AoD(i))^2;
            end
        end
        if n>2
            Nn_ko=sigma(2)/sqrt(2)*(randn(Nr,Nrf*M)+1i*randn(Nr,Nrf*M));
            pXko=Xko(:,2);
            pRko=Rko2;
            Kt=Rko2*(pRko+(pi/Ncan)^2*eye(2*L,2*L))^-1;
            bW_ko=zeros(Nr,Nrf*M);
            bF_ko=zeros(Nt,Nrf*M);
            aro=zeros(Nrf*M,1);
            ato=zeros(Nrf*M,1);
            Yoa=P_AOA2(n-1,:)+(pXko(1:L))';
            Yod=P_AOD2(n-1,:)+(pXko(L+1:2*L))';
            for i=1:4
                ARL2(i)=round((sin(Yoa(i))+1)*G/2);
                ATL2(i)=round((sin(Yod(i))+1)*G/2);
            end
            [Ya,Yb]=sort(Yoa);
            sma=0;
            clos=0;
            for i=1:L-1
                dta(i)=abs(Ya(i)-Ya(i+1));
            end
            if min(dta)<5*pi/Ncan
                [~,bta]=min(dta);
                sma=[Yb(bta),Yb(bta+1)];
                if abs(Yod(Yb(bta))-Yod(Yb(bta+1)))<5*pi/Ncan
                    clos=sma;
                end
            end
            %W_omp
            for i=1:Nrf*M
                aro(i)=Yoa(ceil(i/M))+(rem(i,M)-(M-1)/2)*pi/Ncan;
                if abs(aro(i))>pi/2
                    aro(i)=0;
                end
                [~,b]=min(abs(a_W-aro(i)));
                bW_ko(:,i)=Acan_W(:,b);
                ato(i)=Yod(ceil(i/M))+(rem(i,M)-(M-1)/2)*pi/Ncan;
                if abs(ato(i))>pi/2
                    ato(i)=0;
                end
                [~,b]=min(abs(a_F-ato(i)));
                bF_ko(:,i)=Acan_F(:,b);
            end
            xn=bW_ko'*H*bF_ko+bW_ko'*Nn_ko;
            y=xn(:);
            Q1=kron(bF_ko.',bW_ko')*C;
            CS=reshape(Q1'*y,G,G);
            CSS=CS;
            for cs_i=1:L
                if  ismember(cs_i,clos)>0
                    P_AOA2(n,cs_i)=Yoa(cs_i);
                    [~,b_w]=min(abs(a_W-Yoa(cs_i)));
                    W_ko(:,cs_i)=Acan_W(:,b_w);
                    ARL2(cs_i)=ceil((sin(Yoa(cs_i))+1)*G/2+1);
                    P_AOD2(n,cs_i)=Yod(cs_i);
                    [~,b_F]=min(abs(a_F-Yod(cs_i)));
                    F_ko(:,cs_i)=Acan_F(:,b_F);
                    ATL2(cs_i)=ceil((sin(Yod(cs_i))+1)*G/2+1);
                else
                    CS=CSS;
                    CS(ARL2(cs_i)+Gw2+1:G,:)=0;
                    CS(1:ARL2(cs_i)-Gw2-1,:)=0;
                    CS(:,ATL2(cs_i)+Gw2+1:G)=0;
                    CS(:,1:ATL2(cs_i)-Gw2-1)=0;
                    CS(:,ATL2(1:cs_i-1))=0;
                    CS(ARL2(1:cs_i-1),:);
                    [arl,atl]=find(abs(CS)==max(max(abs(CS))));
                    arl1=FDPK(CS(:,atl));
                    ARL2(cs_i)=round((arl+arl1)/2);
                    oA=(asin(2/G*(arl-1)-1)+asin(2/G*(arl1-1)-1))/2;
                    P_AOA2(n,cs_i)=oA;
                    [~,b_W]=min(abs(a_W-oA));
                    W_ko(:,cs_i)=Acan_W(:,b_W);
                    atl1=FDPK(CS(ARL2(cs_i),:));
                    ATL2(cs_i)=round((atl+atl1)/2);
                    oD=(asin(2/G*(atl-1)-1)+asin(2/G*(atl1-1)-1))/2;
                    P_AOD2(n,cs_i)=oD;
                    [~,b_F]=min(abs(a_F-oD));
                    F_ko(:,cs_i)=Acan_F(:,b_F);
                end
            end
            rXko(1:L,2)=P_AOA2(n,:)-P_AOA2(n-1,:);
            rXko(L+1:2*L,2)=P_AOD2(n,:)-P_AOD2(n-1,:);
            Xko(:,2)=pXko+Kt*(rXko(:,2)-pXko);
            Rko2=pRko-Kt*pRko+(pi/(2*G))^2*eye(2*L,2*L);
        end
        bHeff=W_ko'*H*F_ko;
        [~,A_ko,~]=svd(bHeff);
        Rn=0;
        for i=1:Nrf
            Rn=Rn+log2(1+(A_ko(i,i)^2/(sigma(2)^2)));
        end
        out2(times,n)=real(Rn);
        %---sigma3
        if n==1
            Nn_ko=sigma(3)/sqrt(2)*(randn(Nr,Ncan)+1i*randn(Nr,Ncan));
            bW_ko=Acan_W;
            bF_ko=Acan_F;
            xn=bW_ko'*H*bF_ko+bW_ko'*Nn_ko;
            y=xn(:);
            Q1=kron(bF_ko.',bW_ko')*C;
            CS=reshape(Q1'*y,G,G);
            for cs_i=1:L
                [arl,atl]=find(abs(CS)==max(max(abs(CS))));
                arl1=FDPK(CS(:,atl));
                ARL3(cs_i)=round((arl+arl1)/2);
                oA=(asin(2/G*(arl-1)-1)+asin(2/G*(arl1-1)-1))/2;
                P_AOA3(n,cs_i)=oA;
                [~,b_W]=min(abs(a_W-oA));
                W_ko(:,cs_i)=Acan_W(:,b_W);
                %F_omp
                atl1=FDPK(CS(ARL3(cs_i),:));
                ATL3(cs_i)=round((atl+atl1)/2);
                oD=(asin(2/G*(atl-1)-1)+asin(2/G*(atl1-1)-1))/2;
                P_AOD3(n,cs_i)=oD;
                [~,b_F]=min(abs(a_F-oD));
                F_ko(:,cs_i)=Acan_F(:,b_F);
                if arl-Gw<1
                    CS(1:arl+Gw,:)=0;
                    CS(G+arl-Gw:G,:)=0;
                else
                    if arl+Gw>G
                        CS(arl-Gw:G,:)=0;
                        CS(1:arl-G+Gw,:)=0;
                    else
                        CS(arl-Gw:arl+Gw,:)=0;
                    end
                end
                if atl-Gw<1
                    CS(:,1:atl+Gw)=0;
                    CS(:,G-Gw+atl:G)=0;
                else
                    if atl+Gw>G
                        CS(:,atl-Gw:G)=0;
                        CS(:,1:atl+Gw-G)=0;
                    else
                        CS(:,atl-Gw:atl+Gw)=0;
                    end
                end
            end
            [P_AOA3(1,:),oo]=sort(P_AOA3(1,:));
            OD3=P_AOD3(1,:);
            OD3=OD3(oo);
            ARL3=ARL3(oo);
            P_AOD3(1,:)=OD3;
            ATL3=ATL3(oo);
        end
        if n==2
            Nn_ko=sigma(3)/sqrt(2)*(randn(Nr,Nrf*N)+1i*randn(Nr,Nrf*N));
            bW_ko=zeros(Nr,Nrf*N);
            bF_ko=zeros(Nt,Nrf*N);
            aro=zeros(Nrf*N,1);
            ato=zeros(Nrf*N,1);
            %W_omp
            for i=1:Nrf*N
                aro(i)=P_AOA3(n-1,ceil(i/N))+(rem(i,N)-(N-1)/2)*pi/Ncan;
                [~,b]=min(abs(a_W-aro(i)));
                bW_ko(:,i)=Acan_W(:,b);
                ato(i)=P_AOD3(n-1,ceil(i/N))+(rem(i,N)-(N-1)/2)*pi/Ncan;
                [~,b]=min(abs(a_F-ato(i)));
                bF_ko(:,i)=Acan_F(:,b);
            end
            xn=bW_ko'*H*bF_ko+bW_ko'*Nn_ko;
            y=xn(:);
            Q1=kron(bF_ko.',bW_ko')*C;
            CS=reshape(Q1'*y,G,G);
            CSS=CS;
            for cs_i=1:L
                CS=CSS;
                CS(ARL3(cs_i)+Gw1+1:G,:)=0;
                CS(1:ARL3(cs_i)-Gw1-1,:)=0;
                CS(:,ATL3(cs_i)+Gw1+1:G)=0;
                CS(:,1:ATL3(cs_i)-Gw1-1)=0;
                [arl,atl]=find(abs(CS)==max(max(abs(CS))));
                %W_omp
                arl1=FDPK(CS(:,atl));
                ARL3(cs_i)=round((arl+arl1)/2);
                oA=(asin(2/G*(arl-1)-1)+asin(2/G*(arl1-1)-1))/2;
                P_AOA3(n,cs_i)=oA;
                [~,b_W]=min(abs(a_W-oA));
                W_ko(:,cs_i)=Acan_W(:,b_W);
                %F_omp
                atl1=FDPK(CS(ARL3(cs_i),:));
                ATL3(cs_i)=round((atl+atl1)/2);
                oD=(asin(2/G*(atl-1)-1)+asin(2/G*(atl1-1)-1))/2;
                P_AOD3(n,cs_i)=oD;
                [~,b_F]=min(abs(a_F-oD));
                F_ko(:,cs_i)=Acan_F(:,b_F);
            end
            Rko3=zeros(2*L,2*L);
            for i=1:L
                Xko(i,3)=P_AOA3(n,i)-P_AOA3(n-1,i);
                Xko(i+L,3)=P_AOD3(n,i)-P_AOD3(n-1,i);
                Rko3(i)=(Xko(i,3)-deta_AoA(i))^2;
                Rko3(i+L)=(Xko(i+L,3)-deta_AoD(i))^2;
            end
        end
        if n>2
            Nn_ko=sigma(3)/sqrt(2)*(randn(Nr,Nrf*M)+1i*randn(Nr,Nrf*M));
            pXko=Xko(:,3);
            pRko=Rko3;
            Kt=Rko3*(pRko+(pi/Ncan)^2*eye(2*L,2*L))^-1;
            bW_ko=zeros(Nr,Nrf*M);
            bF_ko=zeros(Nt,Nrf*M);
            aro=zeros(Nrf*M,1);
            ato=zeros(Nrf*M,1);
            Yoa=P_AOA3(n-1,:)+(pXko(1:L))';
            Yod=P_AOD3(n-1,:)+(pXko(L+1:2*L))';
            for i=1:L
                ARL3(i)=round((sin(Yoa(i))+1)*G/2);
                ATL3(i)=round((sin(Yod(i))+1)*G/2);
            end
            [Ya,Yb]=sort(Yoa);
            sma=0;
            clos=0;
            for i=1:L-1
                dta(i)=abs(Ya(i)-Ya(i+1));
            end
            if min(dta)<5*pi/Ncan
                [~,bta]=min(dta);
                sma=[Yb(bta),Yb(bta+1)];
                if abs(Yod(Yb(bta))-Yod(Yb(bta+1)))<5*pi/Ncan
                    clos=sma;
                end
            end
            %W_omp
            for i=1:Nrf*M
                aro(i)=Yoa(ceil(i/M))+(rem(i,M)-(M-1)/2)*pi/Ncan;
                if abs(aro(i))>pi/2
                    aro(i)=0;
                end
                [~,b]=min(abs(a_W-aro(i)));
                bW_ko(:,i)=Acan_W(:,b);
                ato(i)=Yod(ceil(i/M))+(rem(i,M)-(M-1)/2)*pi/Ncan;
                if abs(ato(i))>pi/2
                    ato(i)=0;
                end
                [~,b]=min(abs(a_F-ato(i)));
                bF_ko(:,i)=Acan_F(:,b);
            end
            xn=bW_ko'*H*bF_ko+bW_ko'*Nn_ko;
            y=xn(:);
            Q1=kron(bF_ko.',bW_ko')*C;
            CS=reshape(Q1'*y,G,G);
            CSS=CS;
            for cs_i=1:L
                if ismember(cs_i,clos)>0
                    P_AOA3(n,cs_i)=Yoa(cs_i);
                    [~,b_w]=min(abs(a_W-Yoa(cs_i)));
                    W_ko(:,cs_i)=Acan_W(:,b_w);
                    ARL3(cs_i)=ceil((sin(Yoa(cs_i))+1)*G/2+1);
                    P_AOD3(n,cs_i)=Yod(cs_i);
                    [~,b_F]=min(abs(a_F-Yod(cs_i)));
                    F_ko(:,cs_i)=Acan_F(:,b_F);
                    ATL3(cs_i)=ceil((sin(Yod(cs_i))+1)*G/2+1);
                else
                    CS=CSS;
                    CS(ARL3(cs_i)+Gw2+1:G,:)=0;
                    CS(1:ARL3(cs_i)-Gw2-1,:)=0;
                    CS(:,ATL3(cs_i)+Gw2+1:G)=0;
                    CS(:,1:ATL3(cs_i)-Gw2-1)=0;
                    CS(:,ATL3(1:cs_i-1))=0;
                    CS(ARL3(1:cs_i-1),:)=0;
                    [arl,atl]=find(abs(CS)==max(max(abs(CS))));
                    arl1=FDPK(CS(:,atl));
                    ARL3(cs_i)=round((arl+arl1)/2);
                    oA=(asin(2/G*(arl-1)-1)+asin(2/G*(arl1-1)-1))/2;
                    P_AOA3(n,cs_i)=oA;
                    [~,b_W]=min(abs(a_W-oA));
                    W_ko(:,cs_i)=Acan_W(:,b_W);
                    atl1=FDPK(CS(ARL3(cs_i),:));
                    ATL3(cs_i)=round((atl+atl1)/2);
                    oD=(asin(2/G*(atl-1)-1)+asin(2/G*(atl1-1)-1))/2;
                    P_AOD3(n,cs_i)=oD;
                    [~,b_F]=min(abs(a_F-oD));
                    F_ko(:,cs_i)=Acan_F(:,b_F);
                end
            end
            rXko(1:L,3)=P_AOA3(n,:)-P_AOA3(n-1,:);
            rXko(L+1:2*L,3)=P_AOD3(n,:)-P_AOD3(n-1,:);
            Xko(:,3)=pXko+Kt*(rXko(:,3)-pXko);
            Rko3=pRko-Kt*pRko+(pi/(2*G))^2*eye(2*L,2*L);
        end
        bHeff=W_ko'*H*F_ko;
        [~,A_ko,~]=svd(bHeff);
        Rn=0;
        for i=1:Nrf
            Rn=Rn+log2(1+(A_ko(i,i)^2/(sigma(3)^2)));
        end
        out3(times,n)=real(Rn);
        %计算平均角度误差方差
        PA1=P_AOA1(n,:)-T_AOA(n,:);
        PA2=P_AOA2(n,:)-T_AOA(n,:);
        PA3=P_AOA3(n,:)-T_AOA(n,:);
        PD1=P_AOD1(n,:)-T_AOD(n,:);
        PD2=P_AOD2(n,:)-T_AOD(n,:);
        PD3=P_AOD3(n,:)-T_AOD(n,:);
        P1=PA1.^2+PD1.^2;
        P2=PA2.^2+PD2.^2;
        P3=PA3.^2+PD3.^2;
        E1(n,times)=sum(P1(:))/(2*L);
        E2(n,times)=sum(P2(:))/(2*L);
        E3(n,times)=sum(P3(:))/(2*L);
        if n==15
            if max(E1(1,:))==E1(1,times)
                EAT=T_AOA;
                EAP1=P_AOA1;
                EAP2=P_AOA2;
                EAP3=P_AOA3;
                EDT=T_AOD;
                EDP1=P_AOD1;
                EDP2=P_AOD2;
                EDP3=P_AOD3;
            end
        end
        disp([times,n])
    end
end
out(:,1)=sum(out1,1)/TMAX;
out(:,2)=sum(out2,1)/TMAX;
out(:,3)=sum(out3,1)/TMAX;
OUT=zeros(15,4);
OUT(:,1)=1:15;
OUT(:,2:4)=out;
px=(1:Tmax);
figure(1);
plot(px,OUT(px,2),'-ro');
hold on;
plot(px,OUT(px,3),'-g^');
plot(px,OUT(px,4),'-bo');
xlabel('Channel Block indices');
ylabel('Achievable Throughput(bits)');
legend('SNR=0dB','SNR=10dB','SNR=20dB');
MSE1=sum(E1,2)/TMAX;
MSE2=sum(E2,2)/TMAX;
MSE3=sum(E3,2)/TMAX;
MSE=zeros(15,4);
MSE(:,1)=1:15;
MSE(:,2)=MSE1;
MSE(:,3)=MSE2;
MSE(:,4)=MSE3;
figure(2)
plot(px,MSE(:,2),'-ro');
hold on
plot(px,MSE(:,3),'-g');
plot(px,MSE(:,4),'-bo');
legend('SNR=0dB','SNR=10dB','SNR=20dB');