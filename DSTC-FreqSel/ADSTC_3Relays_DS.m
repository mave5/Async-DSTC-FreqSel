% D-OFDM DSTC over Frequency Selective Channels for two relays
% refer to paper: M. R. Avendi and H. Jafarkhani, 
% "Differential Distributed Space-Time Coding with Imperfect 
% Synchronization in Frequency-Selective Channels," IEEE Transactions on 
% Wireless Communications, vol.14, no.4, pp.1811,1822, April 2015

% R=3,  three relays
% Double sampling

close all;
clear all;
clc;
addpath('functions')
%%

% Synch or ASynch
sync_type='Async';

% error threshold
err_th=1000;

% Flat-Fading or Frequency-Fading 
ch_type='freqsel';
if strcmp(ch_type,'flat')
    L=1;
else
    L=2;
end
% delay spread for frequency-selective channels
Tm=5; 

% number of relays
R=3; 

% synch errors
Ts=1; % symbol duaration
if strcmp(sync_type,'sync')
    tau1=0; % relay 1 delay
    tau2=0.0; % relay 2 delay 
    tau3=0.0; % relay 3 delay 
else
    tau1=0;   % relay 1 delay
    tau2=0.25; % relay 2 delay
    tau3=0.25; %  relay3 delay 
end
    
% Match-Filter outputs, sample at kTs
alfa10=raised_cosine(-tau1,0.9);
alfa11=raised_cosine(1-tau1,0.9);
alfa12=raised_cosine(-1-tau1,0.9);
alfa20=raised_cosine(-tau2,0.9);
alfa21=raised_cosine(1-tau2,.9);
alfa22=raised_cosine(-1-tau2,.9);
alfa30=raised_cosine(-tau3,0.9);
alfa31=raised_cosine(1-tau3,.9);
alfa32=raised_cosine(-1-tau3,.9);

% Match-Filter outputs, sample at kTs+Ts/2
beta10=raised_cosine(-tau1+Ts/2,0.9);
beta11=raised_cosine(1-tau1+Ts/2,0.9);
beta12=raised_cosine(-1-tau1+Ts/2,0.9);
beta20=raised_cosine(-tau2+Ts/2,0.9);
beta21=raised_cosine(1-tau2+Ts/2,.9);
beta22=raised_cosine(-1-tau2+Ts/2,.9);
beta30=raised_cosine(-tau3+Ts/2,0.9);
beta31=raised_cosine(1-tau3+Ts/2,.9);
beta32=raised_cosine(-1-tau3+Ts/2,.9);

% totla power
Ptot_dB=10:5:35;
N0=1;
Ptot=10.^(Ptot_dB/10)*N0;
sig_sr=1;% channel variances
% power allocation and amplification factor in relays
P0=Ptot./2;
Pr=Ptot./(2*R);
AF= Pr./(P0*sig_sr+N0);

% channels variation
ch_dis=1;
fdTs=1e-3;

M=2; % MPSK symbols

% code words
cw=cir_cw(M,R);
%%
for snr_ind=1:length(Ptot)
nerr=0;
nsyms=0;
clc
Ptot(snr_ind)


while nerr < err_th
nSim=0;

%number of OFDM sub-channels
N=16;

if Ptot_dB(snr_ind)<20
    Ns=N*floor(1E5/N);
else
    Ns=N*floor(1E5/N);
end

Ncp=10;% cyclic prefix length
Np=N+Ncp;

% generate channels
Ac1=1/L; % channel power profile
q11=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);
q21=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);
q31=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);
g11=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);
g21=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);
g31=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);

if strcmp(ch_type,'flat')
    q12=zeros(Ns,1);
    q22=zeros(Ns,1);
    q32=zeros(Ns,1);
    g12=zeros(Ns,1);
    g22=zeros(Ns,1);
    g32=zeros(Ns,1);
else
    q12=sqrt(1-Ac1)*flat_cos(Ns,fdTs,ch_dis);
    q22=sqrt(1-Ac1)*flat_cos(Ns,fdTs,ch_dis);
    q32=sqrt(1-Ac1)*flat_cos(Ns,fdTs,ch_dis);
    g12=sqrt(1-Ac1)*flat_cos(Ns,fdTs,ch_dis);
    g22=sqrt(1-Ac1)*flat_cos(Ns,fdTs,ch_dis);
    g32=sqrt(1-Ac1)*flat_cos(Ns,fdTs,ch_dis);
end

while  nSim<Ns 
nSim=nSim+1;

% generate random integer number between 1 to R*M
rcwn=randi(length(cw),1,N);

% space-time encoder
for k1=1:N
    V_in{k1}=cw{rcwn(k1)};
end

% differential encoder
if nSim==1
    %s_km1=[ones(1,N);ones(1,N);ones(1,N)]/sqrt(R);
    s_km1=[ones(1,N);zeros(1,N); zeros(1,N)];
    s_k=s_km1;
else
    s_k=diff_encoder_v(V_in,s_km1);
    s_km1=s_k;
end

% IDFT 
S1=sqrt(N)*ifft(s_k(1,:));
S2=sqrt(N)*ifft(s_k(2,:));
S3=sqrt(N)*ifft(s_k(3,:));

% Add cyclic prefix
S1_cp=[S1(end-Ncp+1:end),S1];
S2_cp=[S2(end-Ncp+1:end),S2];
S3_cp=[S3(end-Ncp+1:end),S3];

%% %%%%%%%%% RX signals at relays

% AWGN noise CN(0,N0)
n11=cxn(Np,N0);%noise at relay 1, TS1
n12=cxn(Np,N0);%noise at relay 1, TS2
n13=cxn(Np,N0);%noise at relay 1, TS3

n21=cxn(Np,N0);%noise at relay 2, TS1
n22=cxn(Np,N0);%noise at relay 2, TS2
n23=cxn(Np,N0);%noise at relay 2, TS3

n31=cxn(Np,N0);%noise at relay 3, TS1
n32=cxn(Np,N0);%noise at relay 3, TS2
n33=cxn(Np,N0);%noise at relay 3, TS3

% Relay 1
q1t=[q11(nSim),zeros(1,Tm),q12(nSim),zeros(1,N-2-Tm)];
temp11=conv(q1t,S1_cp);
R11_cp=sqrt(P0(snr_ind)*R)*temp11(1:Np)+n11;
temp12=conv(q1t,S2_cp);
R12_cp=sqrt(P0(snr_ind)*R)*temp12(1:Np)+n12;
temp13=conv(q1t,S3_cp);
R13_cp=sqrt(P0(snr_ind)*R)*temp13(1:Np)+n13;

% Relay 2
q2t=[q21(nSim),zeros(1,Tm),q22(nSim),zeros(1,N-2-Tm)];
temp21=conv(q2t,S1_cp);
R21_cp=sqrt(P0(snr_ind)*R)*temp21(1:Np)+n21;
temp22=conv(q2t,S2_cp);
R22_cp=sqrt(P0(snr_ind)*R)*temp22(1:Np)+n22;
temp23=conv(q2t,S3_cp);
R23_cp=sqrt(P0(snr_ind)*R)*temp23(1:Np)+n23;

% Relay 3
q3t=[q31(nSim),zeros(1,Tm),q32(nSim),zeros(1,N-2-Tm)];
temp31=conv(q3t,S1_cp);
R31_cp=sqrt(P0(snr_ind)*R)*temp31(1:Np)+n31;
temp32=conv(q3t,S2_cp);
R32_cp=sqrt(P0(snr_ind)*R)*temp32(1:Np)+n32;
temp33=conv(q3t,S3_cp);
R33_cp=sqrt(P0(snr_ind)*R)*temp33(1:Np)+n33;

%%% remove CP at relays
% Relay 1
R11=R11_cp(Ncp+1:Np);
R12=R12_cp(Ncp+1:Np);
R13=R13_cp(Ncp+1:Np);

% Relay 2
R21=R21_cp(Ncp+1:Np);
R22=R22_cp(Ncp+1:Np);
R23=R23_cp(Ncp+1:Np);
% circular time-reverse
cR21_tr=conj([R21(1) fliplr(R21(2:length(R21)))]);
cR22_tr=conj([R22(1) fliplr(R22(2:length(R22)))]);
cR23_tr=conj([R23(1) fliplr(R23(2:length(R23)))]);

% Relay 3
R31=R31_cp(Ncp+1:Np);
R32=R32_cp(Ncp+1:Np);
R33=R33_cp(Ncp+1:Np);
% circular time-reverse
cR31_tr=conj([R31(1) fliplr(R31(2:length(R31)))]);
cR32_tr=conj([R32(1) fliplr(R32(2:length(R32)))]);
cR33_tr=conj([R33(1) fliplr(R33(2:length(R33)))]);

%%% configuration at Relays
% relay 1
B1=eye(R);
X1=B1*[R11;R12;R13];

% Relay 2
B2=[0 1 0;0 0 1;1 0 0];
X2=B2*[R21;R22;R23];

% Relay 3
B3=[0 0 1;1 0 0;0 1 0];
X3=B3*[R31;R32;R33];

% Add Cyclic Prefix at Relays
% Relay 1
CP1=X1(:,end-Ncp+1:end);
X1_cp=sqrt(AF(snr_ind))*[CP1,X1];

% Relay 2
CP2=X2(:,end-Ncp+1:end);
X2_cp=sqrt(AF(snr_ind))*[CP2,X2];

% Relay 3
CP3=X3(:,end-Ncp+1:end);
X3_cp=sqrt(AF(snr_ind))*[CP3,X3];

%% RX signals at Destination

% channels at kTs
g1t=[alfa12*g11(nSim),alfa10*g11(nSim),alfa11*g11(nSim),zeros(1,Tm-3),alfa12*g12(nSim),alfa10*g12(nSim),alfa11*g12(nSim),zeros(1,N-3-Tm)];
g2t=[alfa22*g21(nSim),alfa20*g21(nSim),alfa21*g21(nSim),zeros(1,Tm-3),alfa22*g22(nSim),alfa20*g22(nSim),alfa21*g22(nSim),zeros(1,N-3-Tm)];
g3t=[alfa32*g31(nSim),alfa30*g31(nSim),alfa31*g31(nSim),zeros(1,Tm-3),alfa32*g32(nSim),alfa30*g32(nSim),alfa31*g32(nSim),zeros(1,N-3-Tm)];

% channels at kTs+Ts/2
g1dt=[beta12*g11(nSim),beta10*g11(nSim),beta11*g11(nSim),zeros(1,Tm-3),alfa12*g12(nSim),alfa10*g12(nSim),beta11*g12(nSim),zeros(1,N-3-Tm)];
g2dt=[beta22*g21(nSim),beta20*g21(nSim),beta21*g21(nSim),zeros(1,Tm-3),alfa22*g22(nSim),alfa20*g22(nSim),beta21*g22(nSim),zeros(1,N-3-Tm)];
g3dt=[beta32*g31(nSim),beta30*g31(nSim),beta31*g31(nSim),zeros(1,Tm-3),alfa32*g32(nSim),alfa30*g32(nSim),beta31*g32(nSim),zeros(1,N-3-Tm)];

% RX signals at kTs
Y1_ch=conv(g1t,X1_cp(1,:))+conv(g2t,X2_cp(1,:))+conv(g3t,X3_cp(1,:));
Y2_ch=conv(g1t,X1_cp(2,:))+conv(g2t,X2_cp(2,:))+conv(g3t,X3_cp(2,:));
Y3_ch=conv(g1t,X1_cp(3,:))+conv(g2t,X2_cp(3,:))+conv(g3t,X3_cp(3,:));

% discard the tail
Y1_ch=Y1_ch(1:Np);
Y2_ch=Y2_ch(1:Np);
Y3_ch=Y3_ch(1:Np);

% add AWGN
ly=length(Y1_ch);
w1=cxn(ly,N0);
w2=cxn(ly,N0);
w3=cxn(ly,N0);
Y1_cp=Y1_ch+w1;
Y2_cp=Y2_ch+w2;
Y3_cp=Y3_ch+w3;

% remove CP
Y1=Y1_cp(Ncp+1:end);
Y2=Y2_cp(Ncp+1:end);
Y3=Y3_cp(Ncp+1:end);

% RX signals sampled at kTs+Ts/2
Y1d_ch=conv(g1dt,X1_cp(1,:))+conv(g2dt,X2_cp(1,:))+conv(g3dt,X3_cp(1,:));
Y2d_ch=conv(g1dt,X1_cp(2,:))+conv(g2dt,X2_cp(2,:))+conv(g3dt,X3_cp(2,:));
Y3d_ch=conv(g1dt,X1_cp(3,:))+conv(g2dt,X2_cp(3,:))+conv(g3dt,X3_cp(3,:));

% discard the tail
Y1d_ch=Y1d_ch(1:Np);
Y2d_ch=Y2d_ch(1:Np);
Y3d_ch=Y3d_ch(1:Np);

% add AWGN
ly=length(Y1d_ch);
w1d=cxn(ly,N0);
w2d=cxn(ly,N0);
w3d=cxn(ly,N0);

Y1d_cp=Y1d_ch+w1d;
Y2d_cp=Y2d_ch+w2d;
Y3d_cp=Y3d_ch+w3d;

% remove CP
Y1d=Y1d_cp(Ncp+1:end);
Y2d=Y2d_cp(Ncp+1:end);
Y3d=Y3d_cp(Ncp+1:end);

% add two samples
Y1t=Y1+Y1d;
Y2t=Y2+Y2d;
Y3t=Y3+Y3d;

% take IDFT 
y1= sqrt(1/N)*fft(Y1t,N);
y2= sqrt(1/N)*fft(Y2t,N);
y3= sqrt(1/N)*fft(Y3t,N);
y_k=[y1;y2;y3];
 
 if nSim==1
     y_km1=y_k;
 else
    y=[y_km1,y_k]; 
    Vh=circulant_decoder(y,cw);
    y_km1=y_k;
    
    % count errors
    for k1=1:length(V_in)
        nerr=nerr+1-isequal(Vh{k1},V_in{k1});
    end
    nsyms=length(V_in)+nsyms;
 end
 
end
end
BER(snr_ind)=(nerr)/(nsyms);
end
%% PLOT 
Pb_th=1./(4*Ptot);
clr=['r-+'; 'k->'; 'y-*'; 'g-o'; 'r-<'];
figure
semilogy(Ptot_dB,BER,clr(1,:),'LineWidth',2,'MarkerSize',8); 
grid on
xlabel('Total Power');
ylabel('BER');
legend([ch_type,', \tau=',num2str(tau2) ]);

set(gca,'XTick',Ptot_dB(1):5:Ptot_dB(end),'FontSize',16,...
   'FontName','Times New Roman');
%%
%filename=[ch_type,'_M',num2str(M),'_R',num2str(R),'_tau0',num2str(100*tau2),'_DS'];
%save (filename)