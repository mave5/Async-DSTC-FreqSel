% D-OFDM DSTC over Frequency Selective Channels for two relays
% refer to paper: M. R. Avendi and H. Jafarkhani, 
% "Differential Distributed Space-Time Coding with Imperfect 
% Synchronization in Frequency-Selective Channels," IEEE Transactions on 
% Wireless Communications, vol.14, no.4, pp.1811,1822, April 2015

% R=4 Four relays
close all;
clear all;
clc;
%%

% Synch or ASynch
sync_type='Async';

% Flat-Fading or Frequency-Fading 
ch_type='freq sel';
if strcmp(ch_type,'flat')
    L=1;
else
    L=2;
end
% delay spread for frequency-selective channels
Tm=5; 

%number of OFDM sub-channels
N=16;
Ns=N*floor(2E5/N);
Ncp=7;% cyclic prefix length
Np=N+Ncp;

% number of relays
R=4; 

% synch errors
if strcmp(sync_type,'sync')
    tau1=0; % relay 1 delay
    tau2=0.0; % relay 2 delay 
    tau3=0.0; % relay 3 delay 
    tau4=0.0; % relay 4 delay   
else
    tau1=0;   % relay 1 delay
    tau2=0.3; % relay 2 delay
    tau3=0.3; %  relay3 delay 
    tau4=0.3; % relay 4 delay 
end
    
% Match-Filter outputs
alfa1=raised_cosine(tau1,0.9);
beta1=raised_cosine(1-tau1,0.9);
alfa2=raised_cosine(tau2,0.9);
beta2=raised_cosine(1-tau2,.9);
alfa3=raised_cosine(tau3,0.9);
beta3=raised_cosine(1-tau3,.9);
alfa4=raised_cosine(tau4,0.9);
beta4=raised_cosine(1-tau4,.9);

% totla power
Ptot_dB=0:5:25;
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
% rotation of second constellation
if M==2
    rot_ang=pi/2;
elseif M==4
    rot_ang=pi/4;
end

%%
for snr_ind=1:length(Ptot)
nerr1=0;
nerr2=0;
nerr3=0;
nerr4=0;
nbits=0;
clc
Ptot(snr_ind)


err_th=100;
while mean([nerr1,nerr2,nerr3,nerr4])<err_th
nSim=0;

% generate channels
Ac1=1/L; % channel power profile
q11=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);
q21=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);
q31=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);
q41=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);
g11=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);
g21=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);
g31=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);
g41=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);

if strcmp(ch_type,'flat')
    q12=zeros(Ns,1);
    q22=zeros(Ns,1);
    q32=zeros(Ns,1);
    q42=zeros(Ns,1);
    g12=zeros(Ns,1);
    g22=zeros(Ns,1);
    g32=zeros(Ns,1);
    g42=zeros(Ns,1);
else
    q12=sqrt(1-Ac1)*flat_cos(Ns,fdTs,ch_dis);
    q22=sqrt(1-Ac1)*flat_cos(Ns,fdTs,ch_dis);
    q32=sqrt(1-Ac1)*flat_cos(Ns,fdTs,ch_dis);
    q42=sqrt(1-Ac1)*flat_cos(Ns,fdTs,ch_dis);
    g12=sqrt(1-Ac1)*flat_cos(Ns,fdTs,ch_dis);
    g22=sqrt(1-Ac1)*flat_cos(Ns,fdTs,ch_dis);
    g32=sqrt(1-Ac1)*flat_cos(Ns,fdTs,ch_dis);
    g42=sqrt(1-Ac1)*flat_cos(Ns,fdTs,ch_dis);
end

while  nSim<Ns 
nSim=nSim+1;

% binary data
x_bin=bits(log2(M)*N*R);
b1_in=x_bin(1:4:end);
b2_in=x_bin(2:4:end);
b3_in=x_bin(3:4:end);
b4_in=x_bin(4:4:end);

% MPSK
v1=bin2mpsk(b1_in,M);
v2=bin2mpsk(b2_in,M);
v3=bin2mpsk(b3_in,M,rot_ang);
v4=bin2mpsk(b4_in,M,rot_ang);

% space-time encoder
V_in=QOSTC4TX_encoder(v1,v2,v3,v4);

% differential encoder
if nSim==1
    s_km1=[ones(1,N);ones(1,N);ones(1,N);ones(1,N)]/sqrt(R);
    s_k=s_km1;
else
    s_k=diff_encoder_v(V_in,s_km1);
    s_km1=s_k;
end

% IDFT 
S1=sqrt(N)*ifft(s_k(1,:));
S2=sqrt(N)*ifft(s_k(2,:));
S3=sqrt(N)*ifft(s_k(3,:));
S4=sqrt(N)*ifft(s_k(4,:));

% Add cyclic prefix
S1_cp=[S1(end-Ncp+1:end),S1];
S2_cp=[S2(end-Ncp+1:end),S2];
S3_cp=[S3(end-Ncp+1:end),S3];
S4_cp=[S4(end-Ncp+1:end),S4];


%% %%%%%%%%% RX signals at relays

% AWGN noise CN(0,N0)
n11=cxn(Np,N0);%noise at relay 1, TS1
n12=cxn(Np,N0);%noise at relay 1, TS2
n13=cxn(Np,N0);%noise at relay 1, TS3
n14=cxn(Np,N0);%noise at relay 1, TS4

n21=cxn(Np,N0);%noise at relay 2, TS1
n22=cxn(Np,N0);%noise at relay 2, TS2
n23=cxn(Np,N0);%noise at relay 2, TS3
n24=cxn(Np,N0);%noise at relay 2, TS4

n31=cxn(Np,N0);%noise at relay 3, TS1
n32=cxn(Np,N0);%noise at relay 3, TS2
n33=cxn(Np,N0);%noise at relay 3, TS3
n34=cxn(Np,N0);%noise at relay 3, TS4

n41=cxn(Np,N0);%noise at relay 4, TS1
n42=cxn(Np,N0);%noise at relay 4, TS2
n43=cxn(Np,N0);%noise at relay 4, TS3
n44=cxn(Np,N0);%noise at relay 4, TS4

% Relay 1
q1t=[q11(nSim),zeros(1,Tm),q12(nSim),zeros(1,N-2-Tm)];
temp11=conv(q1t,S1_cp);
R11_cp=sqrt(P0(snr_ind)*R)*temp11(1:Np)+n11;
temp12=conv(q1t,S2_cp);
R12_cp=sqrt(P0(snr_ind)*R)*temp12(1:Np)+n12;
temp13=conv(q1t,S3_cp);
R13_cp=sqrt(P0(snr_ind)*R)*temp13(1:Np)+n13;
temp14=conv(q1t,S4_cp);
R14_cp=sqrt(P0(snr_ind)*R)*temp14(1:Np)+n14;

% Relay 2
q2t=[q21(nSim),zeros(1,Tm),q22(nSim),zeros(1,N-2-Tm)];
temp21=conv(q2t,S1_cp);
R21_cp=sqrt(P0(snr_ind)*R)*temp21(1:Np)+n21;
temp22=conv(q2t,S2_cp);
R22_cp=sqrt(P0(snr_ind)*R)*temp22(1:Np)+n22;
temp23=conv(q2t,S3_cp);
R23_cp=sqrt(P0(snr_ind)*R)*temp23(1:Np)+n23;
temp24=conv(q2t,S4_cp);
R24_cp=sqrt(P0(snr_ind)*R)*temp24(1:Np)+n24;

% Relay 3
q3t=[q31(nSim),zeros(1,Tm),q32(nSim),zeros(1,N-2-Tm)];
temp31=conv(q3t,S1_cp);
R31_cp=sqrt(P0(snr_ind)*R)*temp31(1:Np)+n31;
temp32=conv(q3t,S2_cp);
R32_cp=sqrt(P0(snr_ind)*R)*temp32(1:Np)+n32;
temp33=conv(q3t,S3_cp);
R33_cp=sqrt(P0(snr_ind)*R)*temp33(1:Np)+n33;
temp34=conv(q3t,S4_cp);
R34_cp=sqrt(P0(snr_ind)*R)*temp34(1:Np)+n34;

% Relay 4
q4t=[q41(nSim),zeros(1,Tm),q42(nSim),zeros(1,N-2-Tm)];
temp41=conv(q4t,S1_cp);
R41_cp=sqrt(P0(snr_ind)*R)*temp41(1:Np)+n41;
temp42=conv(q4t,S2_cp);
R42_cp=sqrt(P0(snr_ind)*R)*temp42(1:Np)+n42;
temp43=conv(q4t,S3_cp);
R43_cp=sqrt(P0(snr_ind)*R)*temp43(1:Np)+n43;
temp44=conv(q4t,S4_cp);
R44_cp=sqrt(P0(snr_ind)*R)*temp44(1:Np)+n44;

%%% remove CP at relays
% Relay 1
R11=R11_cp(Ncp+1:Np);
R12=R12_cp(Ncp+1:Np);
R13=R13_cp(Ncp+1:Np);
R14=R14_cp(Ncp+1:Np);

% Relay 2
R21=R21_cp(Ncp+1:Np);
R22=R22_cp(Ncp+1:Np);
R23=R23_cp(Ncp+1:Np);
R24=R24_cp(Ncp+1:Np);
% circular time-reverse
cR21_tr=conj([R21(1) fliplr(R21(2:length(R21)))]);
cR22_tr=conj([R22(1) fliplr(R22(2:length(R22)))]);
cR23_tr=conj([R23(1) fliplr(R23(2:length(R23)))]);
cR24_tr=conj([R24(1) fliplr(R24(2:length(R24)))]);

% Relay 3
R31=R31_cp(Ncp+1:Np);
R32=R32_cp(Ncp+1:Np);
R33=R33_cp(Ncp+1:Np);
R34=R34_cp(Ncp+1:Np);
% circular time-reverse
cR31_tr=conj([R31(1) fliplr(R31(2:length(R31)))]);
cR32_tr=conj([R32(1) fliplr(R32(2:length(R32)))]);
cR33_tr=conj([R33(1) fliplr(R33(2:length(R33)))]);
cR34_tr=conj([R34(1) fliplr(R34(2:length(R34)))]);

% Relay 4
R41=R41_cp(Ncp+1:Np);
R42=R42_cp(Ncp+1:Np);
R43=R43_cp(Ncp+1:Np);
R44=R44_cp(Ncp+1:Np);

%%% configuration at Relays
% relay 1
B1=eye(R);
C1=zeros(R);
X1=B1*[R11;R12;R13;R14];

% Relay 2
B2=zeros(R);
C2=[0 -1 0 0;1 0 0 0;0 0 0 -1;0 0 1 0];
X2=C2*[cR21_tr;cR22_tr;cR23_tr;cR24_tr];

% Relay 3
B3=zeros(R);
C3=[0 0 -1 0;0 0 0 -1;1 0 0 0;0 1 0 0];
X3=C3*[cR31_tr;cR32_tr;cR33_tr;cR34_tr];

% Relay 4
B4=[0 0 0 1; 0 0 -1 0; 0 -1 0 0; 1 0 0 0];
C4=zeros(R);
X4=B4*[R41;R42;R43;R44];

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

% Relay 4
CP4=X4(:,end-Ncp+1:end);
X4_cp=sqrt(AF(snr_ind))*[CP4,X4];

%% RX signals at Destination

% channels 
g1t=[alfa1*g11(nSim),beta1*g11(nSim),zeros(1,Tm-1),alfa1*g12(nSim),beta1*g12(nSim),zeros(1,N-3-Tm)];
g2t=[alfa2*g21(nSim),beta2*g21(nSim),zeros(1,Tm-1),alfa2*g22(nSim),beta2*g22(nSim),zeros(1,N-3-Tm)];
g3t=[alfa3*g31(nSim),beta3*g31(nSim),zeros(1,Tm-1),alfa3*g32(nSim),beta3*g32(nSim),zeros(1,N-3-Tm)];
g4t=[alfa4*g41(nSim),beta4*g41(nSim),zeros(1,Tm-1),alfa4*g42(nSim),beta4*g42(nSim),zeros(1,N-3-Tm)];

Y1_ch=conv(g1t,X1_cp(1,:))+conv(g2t,X2_cp(1,:))+conv(g3t,X3_cp(1,:))+conv(g4t,X4_cp(1,:));
Y2_ch=conv(g1t,X1_cp(2,:))+conv(g2t,X2_cp(2,:))+conv(g3t,X3_cp(2,:))+conv(g4t,X4_cp(2,:));
Y3_ch=conv(g1t,X1_cp(3,:))+conv(g2t,X2_cp(3,:))+conv(g3t,X3_cp(3,:))+conv(g4t,X4_cp(3,:));
Y4_ch=conv(g1t,X1_cp(4,:))+conv(g2t,X2_cp(4,:))+conv(g3t,X3_cp(4,:))+conv(g4t,X4_cp(4,:));

% discard the tail
Y1_ch=Y1_ch(1:Np);
Y2_ch=Y2_ch(1:Np);
Y3_ch=Y3_ch(1:Np);
Y4_ch=Y4_ch(1:Np);

% add AWGN
ly=length(Y1_ch);
w1=cxn(ly,N0);
w2=cxn(ly,N0);
w3=cxn(ly,N0);
w4=cxn(ly,N0);
Y1_cp=Y1_ch+w1;
Y2_cp=Y2_ch+w2;
Y3_cp=Y3_ch+w3;
Y4_cp=Y4_ch+w4;

% remove CP
Y1=Y1_cp(Ncp+1:end);
Y2=Y2_cp(Ncp+1:end);
Y3=Y3_cp(Ncp+1:end);
Y4=Y4_cp(Ncp+1:end);

% take IDFT 
y1= sqrt(1/N)*fft(Y1,N);
y2= sqrt(1/N)*fft(Y2,N);
y3= sqrt(1/N)*fft(Y3,N);
y4= sqrt(1/N)*fft(Y4,N);
y_k=[y1;y2;y3;y4];
 
 if nSim==1
     y_km1=y_k;
 else
    [v1_h,v2_h,v3_h,v4_h]= QOSTC4TX_decoder2(y_k,y_km1,M,rot_ang);
    y_km1=y_k;
    
    % symbols to binary
    b1_hat=mpsk2bin(v1_h,M);
    b2_hat=mpsk2bin(v2_h,M);
    b3_hat=mpsk2bin(v3_h,M,rot_ang);
    b4_hat=mpsk2bin(v4_h,M,rot_ang);

    % count errors
    nerr1=nerr1+sum(abs(b1_in-b1_hat));
    nerr2=nerr2+sum(abs(b2_in-b2_hat));
    nerr3=nerr3+sum(abs(b3_in-b3_hat));    
    nerr4=nerr4+sum(abs(b4_in-b4_hat));    
    nbits=log2(M)*N+nbits;
 end
 
end
end
BER1(snr_ind)=(nerr1)/(nbits);
BER2(snr_ind)=(nerr2)/(nbits);
BER3(snr_ind)=(nerr3)/(nbits);
BER4(snr_ind)=(nerr4)/(nbits);
end
% average BER
BER=mean([BER1;BER2;BER3;BER4]);
%% PLOT 
Pb_th=1./(4*Ptot);
clr=['r-+'; 'k->'; 'y-*'; 'g-o'; 'r-<'];
figure
semilogy(Ptot_dB,BER1,clr(1,:),'LineWidth',2,'MarkerSize',8); 
hold on;
semilogy(Ptot_dB,BER2,clr(4,:),'LineWidth',2,'MarkerSize',8); 

semilogy(Ptot_dB,BER3,clr(2,:),'LineWidth',2,'MarkerSize',8); 
semilogy(Ptot_dB,BER4,clr(3,:),'LineWidth',2,'MarkerSize',8); 

semilogy(Ptot_dB,BER,clr(5,:),'LineWidth',2,'MarkerSize',8); 
semilogy(Ptot_dB,Pb_th,'k-.','LineWidth',2,'MarkerSize',8);
grid on
xlabel('Total Power');
ylabel('BER');
legend([ch_type,', \tau=',num2str(tau2) ]);

set(gca,'XTick',Ptot_dB(1):5:Ptot_dB(end),'FontSize',16,...
   'FontName','Times New Roman');
%save M2_R4_FF_tau04