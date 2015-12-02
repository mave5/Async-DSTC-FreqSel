%% Flat fading channel, improved Jakes model
% By M. R. Avendi

function h=flat_cos(N,fdTs,ch_dis)
% inputs
% N: numebr of samples
%fdTs: fD*Ts normalized doppler frequency
% ch_dis : channel uses distance
% output
% h : channel samples

% number of multipaths
L=8;

% generating uniform random variables
a=-pi;
b=pi;

t=0:N-1;
omega_d=2*pi*fdTs*ch_dis;

Zc=zeros(1,N);
Zs=zeros(1,N);
for n=1:L
    phi_n=a+(b-a).*rand;
    teta_c=a+(b-a).*rand;
    ac_n=(2*pi*n-pi+teta_c)./(4*L);  

    Zc=Zc+sqrt(2/L)*cos(omega_d*t.*cos(ac_n)+phi_n);
end

for n=1:L
    varphi=a+(b-a).*rand;
    teta_s=a+(b-a).*rand;
    as_n=(2*pi*n-pi+teta_s)./(4*L);  
    Zs=Zs+sqrt(2/L)*sin(omega_d*t.*sin(as_n)+varphi);
end    

h=(Zc+1i*Zs)/sqrt(2);
h=h./sqrt(var(h));
end



