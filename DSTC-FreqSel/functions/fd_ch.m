% Project 1 filter method
function [cn_out] =fd_ch(fD,Ts)
%% parameters
Fs=1/Ts;
% filter window length
Ng=1000;
Nz=50000;
%% Complex Gaussian input
x=randn(2*Nz+1,1)/sqrt(2)+sqrt(-1)*randn(2*Nz+1,1)/sqrt(2);
%% compute g(t) and g_hat(t)
t=(-Ng:Ng)*Ts;
for k=1:2*Ng+1
    if(k==Ng+1)
        g(k)=((pi*fD)^(1/4))/gamma(5/4);
    else
        g(k)= besselj(1/4,2*pi*fD*abs(t(k)))/((abs(t(k)))^(1/4));
    end
end
% compute K such that ... 
g_hat=g./sqrt(sum(g.^2));
% verify K
sum_g_hat=sum(g_hat.^2);
% impulse response
cn=conv(x,g_hat)*Ts;
% reject transient
cn=cn(2*Ng+1:length(cn)-2*Ng);
cn=cn/sqrt(var(cn));
var_cn=var(cn);
cn_out=cn;
end
