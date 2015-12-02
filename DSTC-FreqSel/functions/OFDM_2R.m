
function y_k=OFDM_2R(P0,N0,AF,R,s_k,N,Ncp,q1,q2,g1,g2)
% inputs
% P0    : transmit power at source
% N0    : noise power
% AF    : amplification factor
% R     : number of relays
% s_k   : input symbols
% N     : OFDM symbol length
% Ncp   : cyclic prefix length
% q1 and q2 : channels from source to relays
% g1 and g2 : channels from relays to destination

% output
% y_k : output symbols ready to decode


Np=N+Ncp;
% IDFT 
S1=sqrt(N)*ifft(s_k(1,:));
S2=sqrt(N)*ifft(s_k(2,:));

% Add cyclic prefix
S1_cp=[S1(end-Ncp+1:end),S1];
S2_cp=[S2(end-Ncp+1:end),S2];

% RX signals at relays

% AWGN noise CN(0,N0)
n11=cxn(Np,N0);%noise at relay 1, TS1
n12=cxn(Np,N0);%noise at relay 1, TS2
n21=cxn(Np,N0);%noise at relay 2, TS1
n22=cxn(Np,N0);%noise at relay 2, TS2

% Relay 1
q1t=[q1,zeros(1,N-length(q1))];
temp1=conv(q1t,S1_cp);
X11_cp=sqrt(P0*R)*temp1(1:Np)+n11;
temp2=conv(q1t,S2_cp);
X12_cp=sqrt(P0*R)*temp2(1:Np)+n12;

% Relay 2
q2t=[q2,zeros(1,N-length(q2))];
temp3=conv(q2t,S1_cp);
X21_cp=sqrt(P0*R)*temp3(1:Np)+n21;
temp4=conv(q2t,S2_cp);
X22_cp=sqrt(P0*R)*temp4(1:Np)+n22;

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
X11_cp=sqrt(AF)*[CP11,X11];
X12_cp=sqrt(AF)*[CP12,X12];
X21_cp=sqrt(AF)*[CP21,X21_tr];
X22_cp=sqrt(AF)*[CP22,X22_tr];

% channels 
g1t=[g1,zeros(1,N-length(g1))];
g2t=[g2,zeros(1,N-length(g2))];

% RX signals
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
 
  % demodulation OFDM
 y1= sqrt(1/N)*fft(Y1,N);
 y2= sqrt(1/N)*fft(Y2,N);
 y_k=[y1;y2];

end