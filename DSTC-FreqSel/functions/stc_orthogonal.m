
function out=stc_orthogonal(in)

u1=in(1:4:end);
u2=in(2:4:end);
u3=in(3:4:end);
u4=in(4:4:end);

V{length(u1)}=[];% define an empty cell, memory of codewords

% creat the codewords
for k=1:length(u1)

V1=[u1(k)  -u2(k)   -u3(k)  -u4(k)];
V2=[u2(k)   u1(k)    u4(k)  -u3(k)];
V3=[u3(k)  -u4(k)    u1(k)   u2(k)];
V4=[u4(k)   u3(k)   -u2(k)   u1(k)];

sqrtui=sqrt(u1(k)^2+u2(k)^2+u3(k)^2+u4(k)^2);

V{k}=[V1;V2;V3;V4]/sqrtui;

end
out=V;

end


