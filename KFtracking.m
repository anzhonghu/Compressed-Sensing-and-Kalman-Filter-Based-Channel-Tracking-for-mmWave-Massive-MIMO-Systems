close all
clear all
%par.runld=1;
%rng(par.runld);
Tmax=15;%相干时间数
TMAX=200;%实验次数
v_i=1;
v=[4,6,8];
Nt=32;%发射天线数
Nr=128;%接收天线数
Nrf=4;%射频链数
Ncan=128;%可选波束数目
N=13;%n==2时扰动次数
M=5;%n>2时扰动次数
lambda=1e-3;%信号波长
d=lambda/2;%天线间距
L=4;%路径数量
aru=zeros(L*N,1);
atu=zeros(L*N,1);
%________UKF参数___________________________________________________________
Xuo=zeros(4*L,1);
a=[1,1;0,1];
A=diag(diag(repmat(a,2*L,2*L)))+diag(diag(repmat(a,2*L,2*L),1),1);
Eta_c=1e-4;
Eta_m=10;
Mu=2;
Kapa=0;
Lambda_m=Eta_m^2*(4*L+Kapa)-4*L;
Lambda_c=Eta_c^2*(4*L+Kapa)-4*L;
Omega_m=zeros(8*L+1,1);
Omega_m(1)=Lambda_m/(4*L+Lambda_m);
Omega_m(2:8*L+1)=1/(2*(4*L+Lambda_m));
Omega_c=zeros(8*L+1,1);
Omega_c(1)=Lambda_c/(4*L+Lambda_c)+(1-Eta_c^2+Mu);
Omega_c(2:8*L+1)=1/(2*(4*L+Lambda_c));
Chi=zeros(4*L,8*L+1);
Zeta=zeros(2*Nr*Nt,8*L+1);
bW_u=zeros(Nr,L*M);
bF_u=zeros(Nt,L*M);
W_u=zeros(Nr,Nrf);
F_u=zeros(Nt,Nrf);
X_r=zeros(1,2*L);
%KF参数
Xko=zeros(2*L,1);
rXko=zeros(2*L,1);
RXKO=zeros(Tmax,2*L);
XKO=zeros(Tmax,2*L);
%_________________________________________________________________________
Tb=ceil(Ncan/11);
G=256;
Gw=13;
%————————————————网格矩阵生成——————————————————
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
AoA=zeros(1,L);
deta_AoA=zeros(1,L);
AoD=zeros(1,L);
deta_AoD=zeros(1,L);
sigma=sqrt(0.1);%噪声方差
aR_vector=zeros(Nr,1);
aT_vector=zeros(Nt,1);
Acan_W_vector=zeros(Nr,1);
Acan_F_vector=zeros(Nt,1);
Ar=zeros(Nr,L);
At=zeros(Nt,L);
H=zeros(Nr,Nt);
W_ko=zeros(Nr,Nrf);
F_ko=zeros(Nt,Nrf);
W_ez=zeros(Nr,Nrf);
F_ez=zeros(Nt,Nrf);
AoA_W0=[];
AoD_F0=[];
a_W=[];
a_F=[];
Xoa=zeros(L,1);
Xod=zeros(L,1);
dta=zeros(L-1,1);
dtd=zeros(L-1,1);
ARL=zeros(Nrf,1);
ATL=zeros(Nrf,1);
W_real=zeros(Nr,Nrf);%理想接收端
F_real=zeros(Nt,Nrf);%理想发射端
%————————————————Acan的生成———————————————————
Acan_W=zeros(Nr,Ncan);
Acan_F=zeros(Nt,Ncan);
for Acan_j=1:Ncan
    psi=-pi/2+pi*(Acan_j-1)/Ncan; %角度φ随列序改变
    a_W(Acan_j)=psi;
    for Acan_k=1:Nr
        Acan_W_vector(Acan_k,1)=exp((-1i)*2*pi/lambda*(Acan_k-1)*d*sin(psi));
    end
    Acan_W(:,Acan_j)=Acan_W_vector;
end
for Acan_j=1:Ncan
    psi=-pi/2+pi*(Acan_j-1)/Ncan; %角度φ随列序改变
    a_F(Acan_j)=psi;
    for Acan_k=1:Nt
        Acan_F_vector(Acan_k,1)=exp((-1i)*2*pi/lambda*(Acan_k-1)*d*sin(psi));
    end
    Acan_F(:,Acan_j)=Acan_F_vector;
end
%输出结果矩阵
out1=zeros(TMAX,Tmax);
out2=zeros(TMAX,Tmax);
out3=zeros(TMAX,Tmax);
out4=zeros(TMAX,Tmax);
delta_A=v(v_i)*pi/180;%角度变化范围随速度改变
%G_wide=ceil(sin(delta_A)/2*G);
Gw1=ceil(1.6*delta_A/pi*G);
Gw2=7;
for times=1:TMAX;
    T_AOA=zeros(Tmax,L);
    T_AOD=zeros(Tmax,L);
    P_AOA=zeros(Tmax,L);
    P_AOA_ez=zeros(Tmax,L);
    P_AOD=zeros(Tmax,L);
    P_AOD_ez=zeros(Tmax,L);
    D=1;%diag([0.8893,0.0953,0.0107,0.0047]);%L条路径的增益
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
    Bw_ez=zeros(Nrf,Tmax);
    Bf_ez=zeros(Nrf,Tmax);
    [AoA,aoao]=sort(AoA);
    AoD=AoD(aoao);
    deta_AoA=deta_AoA(aoao);
    deta_AoD=deta_AoD(aoao);
    Bw_r=zeros(Nrf,Tmax);
    Bf_r=zeros(Nrf,Tmax);
    for n=1:Tmax
        W_rd=Acan_W(:,randi([1,Ncan],1,Nrf));
        F_rd=Acan_F(:,randi([1,Ncan],1,Nrf));
        %Hn
        T_AOA(n,:)=AoA;
        T_AOD(n,:)=AoD;
        X_r(1:L)=AoA;
        X_r(L+1:2*L)=AoD;
        h_r=CetH(X_r);
        H=reshape(h_r,Nr,Nt);
        AoA=AoA+deta_AoA;
        AoD=AoD+deta_AoD;
        %%%%理想情况
        for l=1:L
            [~,RWb]=min(abs(a_W-T_AOA(n,l)));
            W_real(:,l)=Acan_W(:,RWb);
            Bw_r(l,n)=RWb;
            [~,RFb]=min(abs(a_F-T_AOD(n,l)));
            F_real(:,l)=Acan_F(:,RFb);
            Bf_r(l,n)=RFb;
        end
        bHeff=W_real'*H*F_real;
        [~,A_R,~]=svd(bHeff);
        Rn=0;
        for i=1:Nrf
            Rn=Rn+log2(1+(A_R(i,i)^2/sigma^2));
        end
        out3(times,n)=real(Rn);
        %————————————接收端与发射端预编码更新————————————————
        if n==1
            Nn_ez=sigma/sqrt(2)*(randn(Nr,Ncan)+1i*randn(Nr,Ncan));          %噪声
            bW_ez=Acan_W;               %接收端使用所有波束
            bF_ez=Acan_F;
            Xn_ez=bW_ez'*H*bF_ez+bW_ez'*Nn_ez;
            Bw_ez=zeros(Nrf,Tmax);
            Bf_ez=zeros(Nrf,Tmax);
            for l=1:Nrf
                xn_ez=Xn_ez(:);
                [~,b]=max(abs(xn_ez));
                wb=rem(b,Ncan);
                if wb==0
                    wb=Ncan;
                end
                W_ez(:,l)=bW_ez(:,wb);                 %选取相关性最高的波束
                Bw_ez(l,n)=wb;
                P_AOA_ez(n,l)=a_W(wb);
                if wb>Ncan-Tb
                    Xn_ez(wb-Tb:Ncan,:)=0;
                    Xn_ez(1:wb+Tb-Ncan,:)=0;
                end
                if wb<Tb+1
                    Xn_ez(Ncan-Tb+wb:Ncan,:)=0;
                    Xn_ez(1:wb+Tb,:)=0;
                end
                if (Tb<wb)&&(wb<Ncan-Tb+1)
                    Xn_ez(wb-Tb:wb+Tb,:)=0;
                end
                fb=ceil(b/Ncan);
                F_ez(:,l)=bF_ez(:,fb);
                Bf_ez(l,n)=fb;
                P_AOD_ez(n,l)=a_F(fb);
                if fb>Ncan-Tb
                    Xn_ez(:,fb-Tb:Ncan)=0;
                    Xn_ez(:,1:fb+Tb-Ncan)=0;
                end
                if fb<Tb+1
                    Xn_ez(:,Ncan-Tb+fb:Ncan)=0;
                    Xn_ez(:,1:fb+Tb)=0;
                end
                if (Tb<fb)&&(fb<Ncan-Tb+1)
                    Xn_ez(:,fb-Tb:fb+Tb)=0;
                end
            end
        end
        if n==2
            %W
            Nn_Wez=sigma/sqrt(2)*(randn(Nr,1)+1i*randn(Nr,1));
            Nn_Fez=sigma/sqrt(2)*(randn(Nr,N)+1i*randn(Nr,N));
            bS=ones(Nrf,1);
            bW_ez=zeros(Nr,N);
            ar=zeros(N,1);
            at=zeros(N,1);
            for l=1:Nrf
                for i=1:N  %N次扰动
                    ar(i)=P_AOA_ez(n-1,l)+(i-(N+1)/2)*pi/Ncan;
                    if abs(ar(i))>pi/2
                        ar(i)=0;
                    end
                    [~,b]=min(abs(a_W-ar(i)));
                    bW_ez(:,i)=Acan_W(:,b);
                end
                xn_Wn=bW_ez'*H*F_ez*Vn*bS+bW_ez'*Nn_Wez;
                [~,b]=max(abs(xn_Wn));
                W_ez(:,l)=bW_ez(:,b);
                [~,B_w]=min(abs(a_W-ar(b)));
                Bw_ez(l,n)=B_w;
                P_AOA_ez(n,l)=ar(b);
            end
            %F
            bF_ez=zeros(Nt,N);
            Nn_F=sigma/sqrt(2)*(randn(Nt,N)+1i*randn(Nt,N));
            for l=1:Nrf
                for i=1:N
                    at(i)=P_AOD_ez(n-1,l)+(i-(N+1)/2)*pi/Ncan;
                    if abs(at(i))>pi/2
                        at(i)=0;
                    end
                    [~,F_b]=min(abs(a_F-at(i)));
                    bF_ez(:,i)=Acan_F(:,F_b);
                end
                xn_Fn=W_ez'*H*bF_ez+W_ez'*Nn_Fez;
                [~,b]=max(sum(abs(xn_Fn)));
                F_ez(:,l)=bF_ez(:,b);
                [~,B_f]=min(abs(a_F-at(b)));
                Bf_ez(l,n)=B_f;
                P_AOD_ez(n,l)=at(b);
            end
        end
        if n>2
            %W
            Nn_Wez=sigma/sqrt(2)*(randn(Nr,1)+1i*randn(Nr,1));
            Nn_Fez=sigma/sqrt(2)*(randn(Nr,M)+1i*randn(Nr,M));
            bS=ones(Nrf,1);
            bW_ez=zeros(Nr,M);
            ar=zeros(M,1);
            at=zeros(M,1);
            for l=1:Nrf
                for i=1:M  %N次扰动
                    ar(i)=P_AOA_ez(n-1,l)+(i-(M+1)/2)*pi/Ncan;
                    if abs(ar(i))>pi/2
                        ar(i)=0;
                    end
                    [~,b]=min(abs(a_W-ar(i)));
                    bW_ez(:,i)=Acan_W(:,b);
                end
                xn_Wn=bW_ez'*H*F_ez*Vn*bS+bW_ez'*Nn_Wez;
                [~,b]=max(abs(xn_Wn));
                W_ez(:,l)=bW_ez(:,b);
                [~,B_w]=min(abs(a_W-ar(b)));
                Bw_ez(l,n)=B_w;
                P_AOA_ez(n,l)=ar(b);
            end
            %F
            bF_ez=zeros(Nt,M);
            Nn_F=sigma/sqrt(2)*(randn(Nt,M)+1i*randn(Nt,M));
            for l=1:Nrf
                for i=1:M
                    at(i)=P_AOD_ez(n-1,l)+(i-(M+1)/2)*pi/Ncan;
                    if abs(at(i))>pi/2
                        at(i)=0;
                    end
                    [~,F_b]=min(abs(a_F-at(i)));
                    bF_ez(:,i)=Acan_F(:,F_b);
                end
                xn_Fn=W_ez'*H*bF_ez+W_ez'*Nn_Fez;
                [~,b]=max(sum(abs(xn_Fn)));
                F_ez(:,l)=bF_ez(:,b);
                [~,B_f]=min(abs(a_F-at(b)));
                Bf_ez(l,n)=B_f;
                P_AOD_ez(n,l)=at(b);
            end
        end
        %无预测波束训练法的信道容量
        bHeff=W_ez'*H*F_ez;
        [~,A_ez,Vn]=svd(bHeff);
        Rn=0;
        for i=1:Nrf
            Rn=Rn+log2(1+(A_ez(i,i)^2/sigma^2));
        end
        out1(times,n)=real(Rn);
        %-------------------------卡尔曼滤波预测--------------------------------
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
                ARL(cs_i)=round((arl+arl1)/2);
                oA=(asin(2/G*(arl-1)-1)+asin(2/G*(arl1-1)-1))/2;
                P_AOA(n,cs_i)=oA;
                [~,b_W]=min(abs(a_W-oA));
                W_ko(:,cs_i)=Acan_W(:,b_W);
                %F_omp
                atl1=FDPK(CS(ARL(cs_i),:));
                ATL(cs_i)=round((atl+atl1)/2);                
                oD=(asin(2/G*(atl-1)-1)+asin(2/G*(atl1-1)-1))/2;
                P_AOD(n,cs_i)=oD;
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
            [P_AOA(1,:),oo]=sort(P_AOA(1,:));
            OD1=P_AOD(1,:);
            OD1=OD1(oo);
            ARL=ARL(oo);
            P_AOD(1,:)=OD1;
            ATL=ATL(oo);
        end
        if n==2
            Nn_ko=sigma(1)/sqrt(2)*(randn(Nr,Nrf*N)+1i*randn(Nr,Nrf*N));
            bW_ko=zeros(Nr,Nrf*N);
            bF_ko=zeros(Nt,Nrf*N);
            aro=zeros(Nrf*N,1);
            ato=zeros(Nrf*N,1);
            %W_omp
            for i=1:Nrf*N
                aro(i)=P_AOA(n-1,ceil(i/N))+(rem(i,N)-(N-1)/2)*pi/Ncan;
                [~,b]=min(abs(a_W-aro(i)));
                bW_ko(:,i)=Acan_W(:,b);
                ato(i)=P_AOD(n-1,ceil(i/N))+(rem(i,N)-(N-1)/2)*pi/Ncan;
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
                CS(ARL(cs_i)+Gw1+1:G,:)=0;
                CS(1:ARL(cs_i)-Gw1-1,:)=0;
                CS(:,ATL(cs_i)+Gw1+1:G)=0;
                CS(:,1:ATL(cs_i)-Gw1-1)=0;
                [arl,atl]=find(abs(CS)==max(max(abs(CS))));
                %W_omp
                arl1=FDPK(CS(:,atl));
                ARL(cs_i)=round((arl+arl1)/2);
                oA=(asin(2/G*(arl-1)-1)+asin(2/G*(arl1-1)-1))/2;
                P_AOA(n,cs_i)=oA;
                [~,b_W]=min(abs(a_W-oA));
                W_ko(:,cs_i)=Acan_W(:,b_W);
                %F_omp
                atl1=FDPK(CS(ARL(cs_i),:));
                ATL(cs_i)=round((atl+atl1)/2);
                oD=(asin(2/G*(atl-1)-1)+asin(2/G*(atl1-1)-1))/2;
                P_AOD(n,cs_i)=oD;
                [~,b_F]=min(abs(a_F-oD));
                F_ko(:,cs_i)=Acan_F(:,b_F);
            end
            Rko1=zeros(2*L,2*L);
            for i=1:L
                Xko(i,1)=P_AOA(n,i)-P_AOA(n-1,i);
                Xko(i+L,1)=P_AOD(n,i)-P_AOD(n-1,i);
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
            Yoa=P_AOA(n-1,:)+(pXko(1:L))';
            Yod=P_AOD(n-1,:)+(pXko(L+1:2*L))';
            for i=1:4
                ARL(i)=round((sin(Yoa(i))+1)*G/2);
                ATL(i)=round((sin(Yod(i))+1)*G/2);
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
                    P_AOA(n,cs_i)=Yoa(cs_i);
                    [~,b_w]=min(abs(a_W-Yoa(cs_i)));
                    W_ko(:,cs_i)=Acan_W(:,b_w);
                    ARL(cs_i)=ceil((sin(Yoa(cs_i))+1)*G/2);
                    P_AOD(n,cs_i)=Yod(cs_i);
                    [~,b_F]=min(abs(a_F-Yod(cs_i)));
                    F_ko(:,cs_i)=Acan_F(:,b_F);
                    ATL(cs_i)=ceil((sin(Yod(cs_i))+1)*G/2);
                else
                    CS=CSS;
                    CS(ARL(cs_i)+Gw2+1:G,:)=0;
                    CS(1:ARL(cs_i)-Gw2-1,:)=0;
                    CS(:,ATL(cs_i)+Gw2+1:G)=0;
                    CS(:,1:ATL(cs_i)-Gw2-1)=0;
                    CS(:,ATL(1:cs_i-1))=0;
                    CS(ARL(1:cs_i-1),:);
                    [arl,atl]=find(abs(CS)==max(max(abs(CS))));
                    arl1=FDPK(CS(:,atl));
                    ARL(cs_i)=round((arl+arl1)/2);
                    oA=(asin(2/G*(arl-1)-1)+asin(2/G*(arl1-1)-1))/2;
                    P_AOA(n,cs_i)=oA;
                    [~,b_W]=min(abs(a_W-oA));
                    W_ko(:,cs_i)=Acan_W(:,b_W);
                    atl1=FDPK(CS(ARL(cs_i),:));
                    ATL(cs_i)=round((atl+atl1)/2);
                    oD=(asin(2/G*(atl-1)-1)+asin(2/G*(atl1-1)-1))/2;
                    P_AOD(n,cs_i)=oD;
                    [~,b_F]=min(abs(a_F-oD));
                    F_ko(:,cs_i)=Acan_F(:,b_F);
                end
            end
            rXko(1:L)=P_AOA(n,:)-P_AOA(n-1,:);
            rXko(L+1:2*L)=P_AOD(n,:)-P_AOD(n-1,:);
            Xko=pXko+Kt*(rXko-pXko);
            Rko1=pRko-Kt*pRko+(pi/(2*G))^2*eye(2*L,2*L);
        end
        %KF预测速度的信道容量
        bHeff=W_ko'*H*F_ko;
        [~,A_ko,~]=svd(bHeff);
        Rn=0;
        for i=1:Nrf
            Rn=Rn+log2(1+(A_ko(i,i)^2/sigma^2));
        end
        out2(times,n)=real(Rn);
        %以UKF的预测与更新-----------------------------
        if n==2
            for i=1:L
                Xuo(2*i-1)=P_AOA_ez(n,i);
                Xuo(2*i)=P_AOA_ez(n,i)-P_AOA_ez(n-1,i);
                Xuo(2*i-1+2*L)=P_AOD_ez(n,i);
                Xuo(2*i+2*L)=P_AOD_ez(n,i)-P_AOD_ez(n-1,i);
            end
            Rk=zeros(4*L,4*L);
            for i=1:L
                [~,bAoA]=min(abs(AoA-P_AOA_ez(n,i)));
                Rk(2*i-1,2*i-1)=(AoA(bAoA)-P_AOA_ez(n,i))^2;
                Rk(2*i,2*i)=(deta_AoA(bAoA)-Xuo(2*i))^2;
                [~,bAoD]=min(abs(AoD-P_AOD_ez(n,i)));
                Rk(2*i-1+2*L,2*i-1+2*L)=(AoD(bAoD)-P_AOD_ez(n,i))^2;
                Rk(2*i+2*L,2*i+2*L)=(deta_AoD(bAoD)-Xuo(2*i+2*L))^2;
            end
        end
        if n>2
            Xkk=A*Xuo;
            for i=1:L*M  %M次扰动
                aru(i)=Xkk(2*ceil(i/M)-1)+(rem(i,M)+1-(M+1)/2)*pi/Ncan;
                [~,b]=min(abs(a_W-aru(i)));
                bW_u(:,i)=Acan_W(:,b);
            end
            for i=1:L*M  %M次扰动
                atu(i)=Xkk(2*ceil(i/M)+2*L-1)+(rem(i,M)+1-(M+1)/2)*pi/Ncan;
                [~,b]=min(abs(a_F-atu(i)));
                bF_u(:,i)=Acan_F(:,b);
            end
            Gk=kron(bF_u.',bW_u');
            Rkk=Rk;%A*Rk*A.';
            bRk=sqrt((4*L+Lambda_c)*Rkk);
            Chi(:,1)=Xkk;
            for i=1:4*L
                Chi(:,i+1)=Xkk+bRk(:,i);
                Chi(:,i+4*L+1)=Xkk-bRk(:,i);
            end
            hk=0;
            Pik=0;
            for i=1:8*L+1
                Zeta(:,i)=[real(CetH(Chi(:,i)));imag(CetH(Chi(:,i)))];
                hk=hk+Omega_m(i)*Zeta(:,i);
            end
            for i=1:8*L+1
                Pik=Pik+Omega_c(i)*(Zeta(:,i)-hk)*(Zeta(:,i)-hk).';
            end
            rGk=[real(Gk) -imag(Gk);imag(Gk) real(Gk)];
            rh=[real(H(:));imag(H(:))];
            ryk=rGk*rh;
            byk=ryk-rGk*hk;
            Sk=rGk*Pik*rGk.'+1/2*diag(diag(ones(2*(M*Nrf)^2,2*(M*Nrf)^2)));
            %Sk=rGk*Pik*rGk.'+1/2*diag(diag(ones(2*Nrf^2,2*Nrf^2)));
            Puxi=0;
            for i=1:8*L+1
                Puxi=Puxi+Omega_c(i)*(Chi(:,i)-Xkk)*(Zeta(:,i)-hk).';
            end
            Xuo=Xkk+Puxi*rGk.'*Sk^-1*byk;
            %Rk=Rkk-Puxi*rGk.'*Sk^-1*rGk*Puxi.';
            for i=1:Nrf
                [~,bw]=min(abs(a_W-Xuo(2*i-1)));
                W_u(:,i)=Acan_W(:,bw);
                [~,bf]=min(abs(a_F-Xuo(2*L+2*i-1)));
                F_u(:,i)=Acan_F(:,bf);
            end
            %UKF预测的信道容量
            bHeff=W_u'*H*F_u;
            [~,A_ukf,~]=svd(bHeff);
            Rn=0;
            for i=1:Nrf
                Rn=Rn+log2(1+(A_ukf(i,i)^2/sigma^2));
            end
            out4(times,n)=real(Rn);
        end
        disp([times,n]) 
    end
    out4(times,1)=out1(times,1);
    out4(times,2)=out1(times,2);
end
out11=sum(out1,1)/TMAX;%波束训练
out22=sum(out2,1)/TMAX;%卡尔曼滤波算法
out33=sum(out3,1)/TMAX;%理想情况
out44=sum(out4,1)/TMAX;%波束UKF
OUT=zeros(Tmax,4);
OUT(:,1)=1:15;
OUT(:,2)=out11;
OUT(:,3)=out22;
OUT(:,4)=out33;
OUT(:,5)=out44;
px=(1:Tmax);
figure(1);
plot(px,OUT(:,2),'-ro');
hold on;
plot(px,OUT(:,3),'-g^');
plot(px,OUT(:,4),'-bo');
plot(px,OUT(:,5),'-mo');
xlabel('Channel Block indices');
ylabel('Achievable Throughput(bits)');
legend('Traditional estimation','KF','Ideal State','UKF');