function [ C ] = channel(N)
%CHANNEL Summary of this function goes here
%   Detailed explanation goes here
% number of frequency samples
%N=64;
v=15;
fc=10^9;
c=3*10^8;
fD=(fc/c)*v;
tu=[0 4];
%% Two Taps => tu=[0,4]
cn1=fd_ch(fD);
cn2=fd_ch(fD);

for k=1:N
     C(k)=cn1(1)+(cn2(1)*exp(-sqrt(-1)*2*pi*(k-1)*tu(2)/N));
 end



