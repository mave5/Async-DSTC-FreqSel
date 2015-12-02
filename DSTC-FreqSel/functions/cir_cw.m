% circulant codes- code words
function cw=cir_cw(M,R)
% M : MPSK modulation size
% R : number of relays  


if R==2
% rotation angle
teta1=0;
teta2=pi/2;

m=0:M-1;
u1=exp(1i*(2*pi*m/M+teta1));
cw{1}=u1(1)*eye(R);
cw{2}=u1(2)*eye(R);

u2=exp(1i*(2*pi*m/M+teta2));
cw{3}=u2(1)*[0 1;1 0];
cw{4}=u2(2)*[0 1;1 0];

end

if R==3
% rotation angle
teta1=0;
teta2=pi/9;
teta3=2*pi/9;

m=0:M-1;
u1=exp(1i*(2*pi*m/M+teta1));
cw{1}=u1(1)*eye(R);
cw{2}=u1(2)*eye(R);

u2=exp(1i*(2*pi*m/M+teta2));
cw{3}=u2(1)*[0 1 0;0 0 1;1 0 0];
cw{4}=u2(2)*[0 1 0;0 0 1;1 0 0];

u3=exp(1i*(2*pi*m/M+teta3));
cw{5}=u3(1)*[0 0 1;1 0 0;0 1 0];
cw{6}=u3(2)*[0 0 1;1 0 0;0 1 0];
end


end