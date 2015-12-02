% Received SNR as a function of n and tau
% 
clc
clear all
close all
%%

R=4;
Ptot_dB=25;% totla power
N0=1;
sig_sr=1;
Ptot=10.^(Ptot_dB/10)*N0;
P0=Ptot./2;
Pr=Ptot./(2*R);
AF= Pr./(P0*sig_sr+N0);

g1=1;
g2=1;

N=64;
vtau=[0 .2 .4 .5 .6 .8 1];
a=0.0;

for k=1:length(vtau)
    tau=vtau(k);
    %alfa=sinc(tau);
    %beta=sinc(1-tau);
    alfa=sinc(tau).*cos(pi*a*tau)./(1-4*a^2*tau^2);    
    tau=1-vtau(k);
    beta=sinc(tau).*cos(pi*a*tau)./(1-4*a^2*tau^2);
    
    n=0:N-1;
    c=abs(alfa+beta*exp(-1i*2*pi*n/N)).^2;
    sig2(k,:)=N0*(1+AF*(1+(R-1)*c));
    gama(k,:)=AF*P0*(1+(R-1)*c)./(1+AF*(1+(R-1)*c));
    
end
figure
plot(n,sig2,'o')
legend('\tau=0','\tau=0.2','\tau=.4')
xlabel('n')
ylabel('sig2')

a1=mean(sig2(1,:));
a2=mean(sig2(2,:));
a3=mean(sig2(3,:));

figure 
plot(n,10*log10(gama),'LineWidth',2);
legend('\tau=0','\tau=0.2','\tau=0.4','\tau=0.5')

grid on
xlabel('n');
ylabel('gama');
set(gca,'XTick',n(1):7:n(end),'FontSize',16,...
   'FontName','Times New Roman');
%axis([0 N-1  15 20])


