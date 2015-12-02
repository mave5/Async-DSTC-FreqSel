% D-OFDM DSTC over Frequency Selective Channels for two relays
% refer to paper: M. R. Avendi and H. Jafarkhani, 
% "Differential Distributed Space-Time Coding with Imperfect 
% Synchronization in Frequency-Selective Channels," IEEE Transactions on 
% Wireless Communications, vol.14, no.4, pp.1811,1822, April 2015

% Double Sampling is added 
close all;
clear all;
clc;
addpath functions
%% 

%synch or ASynch 
sync_type='Async';

% flat or freq sel
ch_type='freq-sel';
if strcmp(ch_type,'flat')
    L=1;
else
    L=2;
end

% delay spread for frequency-selective channels
Tm=5; 

% totla power
Ptot_dB=0:5:30;
N0=1;
Ptot=10.^(Ptot_dB/10)*N0;

% power allocation
sig_sr=1;
P0=Ptot./2;
Pr=Ptot./4;
AF= Pr./(P0*sig_sr+N0);

% synch errors
Ts=1;
if strcmp(sync_type,'sync')
    tau1=0; 
    tau2=0;
else
    tau1=0.0; 
    tau2=.25;
end

% Matched Filter output, sample at kTs
alfa10=raised_cosine(-tau1,0.9);
alfa11=raised_cosine(1-tau1,0.9);
alfa12=raised_cosine(-1-tau1,0.9);
alfa20=raised_cosine(-tau2,0.9);
alfa21=raised_cosine(1-tau2,.9);
alfa22=raised_cosine(-1-tau2,.9);

% Matched Filter output, sample at Ts/2+kTs
s_d=Ts/2;% sampling delay
beta10=raised_cosine(s_d-tau1,0.9);
beta11=raised_cosine(s_d+1-tau1,0.9);
beta12=raised_cosine(s_d-1-tau1,0.9);
beta20=raised_cosine(s_d-tau2,0.9);
beta21=raised_cosine(s_d+1-tau2,.9);
beta22=raised_cosine(s_d-1-tau2,.9);

% number of relays
R=2; 

% channel variation
fdTs=1e-3;
ch_dis=1;

% MPSK symbols
M=2; 

% Alamouti type 
type=1;

for snr_ind=1:length(Ptot)
nerr1=0;
nerr2=0;
nbits=0;
clc
Ptot(snr_ind)

err_th=1000;
while mean([nerr1,nerr2])<err_th
nSim=0;


%number of OFDM sub-channels
N=16;
if Ptot_dB(snr_ind)<30
    Ns=N*floor(2E5/N);
else
    Ns=N*floor(3E5/N);
end
Ncp=10;% cyclic prefix length
Np=N+Ncp;

% generate channels
Ac1=1/L; % channel power profile
q11=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);
q12=sqrt(1-Ac1)*flat_cos(Ns,fdTs,ch_dis);
q21=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);
q22=sqrt(1-Ac1)*flat_cos(Ns,fdTs,ch_dis);
g11=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);
g12=sqrt(1-Ac1)*flat_cos(Ns,fdTs,ch_dis);
g21=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);
g22=sqrt(1-Ac1)*flat_cos(Ns,fdTs,ch_dis);

while  nSim<Ns 
nSim=nSim+1;

%input bits generation
b1_in=bits(log2(M)*N);
b2_in=bits(log2(M)*N);

%MPSK symbol
v1_in=bin2mpsk(b1_in,M); 
v2_in=bin2mpsk(b2_in,M); 
temp=[v1_in,v2_in];
v_in=reshape(temp.',2*N,1);

% space-time encoding
V_in=stc_alamouti(v_in,type);

% differential encoder
if nSim==1
    s_km1=[ones(1,N);ones(1,N)]/sqrt(2);
    s_k=s_km1;
else
    s_k=diff_encoder_v(V_in,s_km1);
    s_km1=s_k;
end

% IDFT 
S1=sqrt(N)*ifft(s_k(1,:));
S2=sqrt(N)*ifft(s_k(2,:));

% Add cyclic prefix
S1_cp=[S1(end-Ncp+1:end),S1];
S2_cp=[S2(end-Ncp+1:end),S2];

% AWGN noise CN(0,N0)
n11=cxn(Np,N0);%noise at relay 1, TS1
n12=cxn(Np,N0);%noise at relay 1, TS2
n21=cxn(Np,N0);%noise at relay 2, TS1
n22=cxn(Np,N0);%noise at relay 2, TS2

%% RX signals at relays

% Relay 1
q1t=[q11(nSim),zeros(1,Tm),q12(nSim),zeros(1,N-2-Tm)];
temp1=conv(q1t,S1_cp);
X11_cp=sqrt(P0(snr_ind)*R)*temp1(1:Np)+n11;
temp2=conv(q1t,S2_cp);
X12_cp=sqrt(P0(snr_ind)*R)*temp2(1:Np)+n12;

% Relay 2
q2t=[q21(nSim),zeros(1,Tm),q22(nSim),zeros(1,N-2-Tm)];
temp3=conv(q2t,S1_cp);
X21_cp=sqrt(P0(snr_ind)*R)*temp3(1:Np)+n21;
temp4=conv(q2t,S2_cp);
X22_cp=sqrt(P0(snr_ind)*R)*temp4(1:Np)+n22;

% remove CP at relays
X11=X11_cp(Ncp+1:Np);
X12=X12_cp(Ncp+1:Np);
X21=X21_cp(Ncp+1:Np);
X22=X22_cp(Ncp+1:Np);

% circular time-reverse
X21_tr=[X21(1) fliplr(X21(2:length(X21)))];
X22_tr=[X22(1) fliplr(X22(2:length(X22)))];

% Add Cyclic Prefrex
CP11=X11(end-Ncp+1:end);
CP12=X12(end-Ncp+1:end);
CP21=X21_tr(end-Ncp+1:end);
CP22=X22_tr(end-Ncp+1:end);
X11_cp=sqrt(AF(snr_ind))*[CP11,X11];
X12_cp=sqrt(AF(snr_ind))*[CP12,X12];
X21_cp=sqrt(AF(snr_ind))*[CP21,X21_tr];
X22_cp=sqrt(AF(snr_ind))*[CP22,X22_tr];

%% RX signals at Destination 
% channels at kTs
g1t=[alfa12*g11(nSim),alfa10*g11(nSim),alfa11*g11(nSim),zeros(1,Tm-3),alfa12*g12(nSim),alfa10*g12(nSim),alfa11*g12(nSim),zeros(1,N-3-Tm)];
g2t=[alfa22*g21(nSim),alfa20*g21(nSim),alfa21*g21(nSim),zeros(1,Tm-3),alfa22*g22(nSim),alfa20*g22(nSim),alfa21*g22(nSim),zeros(1,N-3-Tm)];

% channels at Ts/2+kTs
g1dt=[beta12*g11(nSim),beta10*g11(nSim),beta11*g11(nSim),zeros(1,Tm-3),alfa12*g12(nSim),alfa10*g12(nSim),beta11*g12(nSim),zeros(1,N-3-Tm)];
g2dt=[beta22*g21(nSim),beta20*g21(nSim),beta21*g21(nSim),zeros(1,Tm-3),alfa22*g22(nSim),alfa20*g22(nSim),beta21*g22(nSim),zeros(1,N-3-Tm)];

% sampled at kTs
Y1_ch=conv(g1t,X11_cp)-conv(g2t,conj(X22_cp));
Y2_ch=conv(g1t,(X12_cp))+conv(g2t,conj(X21_cp));
Y1_ch=Y1_ch(1:Np);
Y2_ch=Y2_ch(1:Np);
ly=length(Y1_ch);
w1=cxn(ly,N0);
w2=cxn(ly,N0);
Y1_cp=Y1_ch+w1;
Y2_cp=Y2_ch+w2;
% remove CP
Y1=Y1_cp(Ncp+1:end);
Y2=Y2_cp(Ncp+1:end);

% RX signals sampled at Ts/2+kTs
Y1d_ch=conv(g1dt,X11_cp)-conv(g2dt,conj(X22_cp));
Y2d_ch=conv(g1dt,(X12_cp))+conv(g2dt,conj(X21_cp));
Y1d_ch=Y1d_ch(1:Np);
Y2d_ch=Y2d_ch(1:Np);
ly=length(Y1d_ch);
w1d=cxn(ly,N0);
w2d=cxn(ly,N0);
Y1d_cp=Y1d_ch+w1d;
Y2d_cp=Y2d_ch+w2d;
% remove CP
Y1d=Y1d_cp(Ncp+1:end);
Y2d=Y2d_cp(Ncp+1:end);
 
% add samples together
Y1t=Y1+Y1d;
Y2t=Y2+Y2d;

% demodulation OFDM
y1= sqrt(1/N)*fft(Y1t,N);
y2= sqrt(1/N)*fft(Y2t,N);
y_k=[y1;y2];
 
 if nSim==1
     y_km1=y_k;
 else
    [v1_hat,v2_hat]= dstc_decoder(y_k,y_km1,type);
    y_km1=y_k;
    
    % MPSK demodulation
    b1_hat=mpsk2bin(v1_hat,M);
    b2_hat=(mpsk2bin(v2_hat,M));
    
    % count errors
    nerr1=nerr1+sum(abs(b1_in-b1_hat));
    nerr2=nerr2+sum(abs(b2_in-b2_hat));
    nbits=log2(M)*N+nbits;
 end
 
end
end
BER1(snr_ind)=(nerr1)/(nbits);
BER2(snr_ind)=(nerr2)/(nbits);

end

% average BER
BER=mean([BER1;BER2]);

%% PLOT 
Pb_th=1./(4*Ptot);
clr=['r-+'; 'k->'; 'y-*'; 'g-o'; 'r-<'];
figure
semilogy(Ptot_dB,BER1,clr(1,:),'LineWidth',2,'MarkerSize',8); 
hold on;
semilogy(Ptot_dB,BER2,clr(4,:),'LineWidth',2,'MarkerSize',8); 
semilogy(Ptot_dB,BER,clr(2,:),'LineWidth',2,'MarkerSize',8); 
semilogy(Ptot_dB,Pb_th,'k-.','LineWidth',2,'MarkerSize',8);
grid on
xlabel('Total Power');
ylabel('BER');
%legend(['Proposed, \tau=',num2str(tau2)]);
legend([ch_type,', \tau=',num2str(tau2) ]);

set(gca,'XTick',Ptot_dB(1):5:Ptot_dB(end),'FontSize',16,...
   'FontName','Times New Roman');
axis([0 Ptot_dB(end) 1e-5 1])
%axis([0 30 1e-4 1])