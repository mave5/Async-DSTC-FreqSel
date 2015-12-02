% Quasi-Orthogonal Space-Time code encoder for 4-tx antennas
function V=QOSTC4TX_encoder(vu1,vu2,vu3,vu4)


V{length(vu1)}=[];% define an empty cell, memory of codewords

% creat the codewords
for k=1:length(vu1)
u1=vu1(k);
u2=vu2(k);
u3=vu3(k);
u4=vu4(k);

sumui=sqrt(abs(u1)^2+abs(u2)^2+abs(u3)^2+abs(u4)^2);

V1=[u1  -u2'    -u3'  u4];
V2=[u2   u1'    -u4'  -u3];
V3=[u3  -u4'    u1'   -u2];
V4=[u4   u3'    u2'   u1];


V{k}=[V1;V2;V3;V4]/sumui;

end

