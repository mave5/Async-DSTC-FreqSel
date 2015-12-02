% D-OFDM DSTC over Frequency Selective Channels for two relays
% refer to paper: M. R. Avendi and H. Jafarkhani, 
% "Differential Distributed Space-Time Coding with Imperfect 
% Synchronization in Frequency-Selective Channels," IEEE Transactions on 
% Wireless Communications, vol.14, no.4, pp.1811,1822, April 2015

% two relays and coded scenario
% channels have three taps

close all;
clear all;
clc;
%% 

% OFDM length
N=16;

% code rate
code_rate=3;
d_intv=N*1000;    

%synch or ASynch 
sync_type='Async';

% flat or freq sel
ch_type='freq-sel';
if strcmp(ch_type,'flat')
    L=1;
else
    L=3;
end
% delay spread for frequency-selective channels
Tm=4; 

%number of OFDM sub-channels
Ns=N*floor(1E5/N);
Ncp=15;% cyclic prefix length
Np=N+Ncp;

% totla power
Ptot_dB=0:5:20;
N0=1;
Ptot=10.^(Ptot_dB/10)*N0;
% power allocation
sig_sr=1;
P0=Ptot./2;
Pr=Ptot./4;
AF= Pr./(P0*sig_sr+N0);

% synch errors
if strcmp(sync_type,'sync')
    tau1=0; 
    tau2=0;
else
    tau1=0.0; 
    tau2=.0;
end    
% Matched Filter output
alfa1=raised_cosine(tau1,0.9);
beta1=raised_cosine(1-tau1,0.9);
alfa2=raised_cosine(tau2,0.9);
beta2=raised_cosine(1-tau2,.9);
beta22=raised_cosine(-1-tau2,.9);

% number of relays
R=2; 

% channel variation
fdTs=1e-3;
ch_dis=1;

% MPSK symbols
M=2; 

% Alamouti type 
type=1;
%%
for snr_ind=1:length(Ptot)
nerr1=0;
nerr2=0;
nbits=0;
clc
Ptot(snr_ind)

err_th=100;
while mean([nerr1,nerr2])<err_th
nSim=1;

% generate channels
Ac1=1/L; % channel power profile
q11=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);
q12=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);
q13=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);

q21=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);
q22=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);
q23=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);

g11=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);
g12=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);
g13=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);

g21=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);
g22=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);
g23=sqrt(Ac1)*flat_cos(Ns,fdTs,ch_dis);

while  nSim<Ns-d_intv

if nSim==1

% reference symbols
s_km1=[ones(1,N);ones(1,N)]/sqrt(2);
q1=[q11(nSim),zeros(1,Tm-1),q12(nSim),0,q13(nSim)];
q2=[q21(nSim),zeros(1,Tm-1),q22(nSim),0,q23(nSim)];
g1=[0,alfa1*g11(nSim),beta1*g11(nSim),zeros(1,Tm-3),alfa1*g12(nSim),beta1*g12(nSim),alfa1*g13(nSim),beta1*g13(nSim)];
g2=[beta22*g21(nSim),alfa2*g21(nSim),beta2*g21(nSim),zeros(1,Tm-3),alfa2*g22(nSim),beta2*g22(nSim),alfa2*g23(nSim),beta2*g23(nSim)];
y_km1=OFDM_2R(P0(snr_ind),N0,AF(snr_ind),R,s_km1,N,Ncp,q1,q2,g1,g2);
nSim=2;
else

%input bits generation
b1_uncoded=bits(log2(M)*d_intv);
b2_uncoded=bits(log2(M)*d_intv);


% repetition encoder
b1_coded=rep_encoder(b1_uncoded',code_rate);
b2_coded=rep_encoder(b2_uncoded',code_rate);

% interleaving
b1_intv=reshape(b1_coded,code_rate,[]);
b2_intv=reshape(b2_coded,code_rate,[]);

% read interleaver by rows
b1_intv_row=reshape(b1_intv', 1,[]);
b2_intv_row=reshape(b2_intv', 1,[]);

v11_h=[];
v22_h=[];
for k1=1:N:length(b1_intv_row)

nSim=nSim+1;
%input bits generation
b1_in=b1_intv_row(k1:N+k1-1)';
b2_in=b2_intv_row(k1:N+k1-1)';

%MPSK symbol
v1_in=bin2mpsk(b1_in,M); 
v2_in=bin2mpsk(b2_in,M); 
temp=[v1_in,v2_in];
v_in=reshape(temp.',2*N,1);

% space-time encoding
V_in=stc_alamouti(v_in,type);

% differential encoder
s_k=diff_encoder_v(V_in,s_km1);
s_km1=s_k;

q1=[q11(nSim),zeros(1,Tm),q12(nSim),q13(nSim)];
q2=[q21(nSim),zeros(1,Tm),q22(nSim),q23(nSim)];
g1=[0,alfa1*g11(nSim),beta1*g11(nSim),zeros(1,Tm-2),alfa1*g12(nSim),beta1*g12(nSim),alfa1*g13(nSim),beta1*g13(nSim)];
g2=[beta22*g21(nSim),alfa2*g21(nSim),beta2*g21(nSim),zeros(1,Tm-2),alfa2*g22(nSim),beta2*g22(nSim),alfa2*g23(nSim),beta2*g23(nSim)];
y_k=OFDM_2R(P0(snr_ind),N0,AF(snr_ind),R,s_k,N,Ncp,q1,q2,g1,g2);

[v1_h,v2_h]= dstc_decoder(y_k,y_km1,type);
y_km1=y_k;
v11_h=[v11_h,v1_h];
v22_h=[v22_h,v2_h];

end

% applly MRC
v1r_h=reshape(v11_h,[],code_rate);
v1_hat=sum(v1r_h,2);

v2r_h=reshape(v22_h,[],code_rate);
v2_hat=sum(v2r_h,2);
    
    % MPSK demodulation
    b1_hat=mpsk2bin(v1_hat,M);
    b2_hat=(mpsk2bin(v2_hat,M));

    % count errors
    nerr1=nerr1+sum(abs(b1_uncoded-b1_hat));
    nerr2=nerr2+sum(abs(b2_uncoded-b2_hat));
    nbits=length(b1_uncoded)+nbits;
    
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
legend([ch_type,', code rate=',num2str(1/code_rate) ]);

set(gca,'XTick',Ptot_dB(1):5:Ptot_dB(end),'FontSize',16,...
   'FontName','Times New Roman');
axis([0 Ptot_dB(end) 1e-5 1])
