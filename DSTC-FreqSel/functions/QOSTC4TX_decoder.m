% decoder for Orthogonal STC with 4 antennas
function [v1_h,v2_h,v3_h,v4_h]= QOSTC4TX_decoder(y_k,y_km1,M,rot_ang)

[~,N]=size(y_k);

% codewords
m=0:M-1;
V1=exp(1i*2*pi*m/M);
V2=exp(1i*2*pi*m/M)*exp(1i*rot_ang);
V1c=allcomb(V1,V1,V2,V2);

for n=1:length(V1c)
    u1=V1c(n,1);
    u2=V1c(n,2);
    u3=V1c(n,3);
    u4=V1c(n,4);
    sumui=sqrt(abs(u1)^2+abs(u2)^2+abs(u3)^2+abs(u4)^2);
    
    V1=[u1  -u2'    -u3'  u4];
    V2=[u2   u1'    -u4'  -u3];
    V3=[u3  -u4'    u1'   -u2];
    V4=[u4   u3'    u2'   u1];
    
    U{n}=[V1;V2;V3;V4]/sumui;
end

for n=1:N      
    for l=1:length(V1c)
        temp(l)=norm(y_k(:,n)-U{l}*y_km1(:,n));
    end
    [t1,lh]=min(temp);
    v1_h(n)=V1c(lh,1);
    v2_h(n)=V1c(lh,2);
    v3_h(n)=V1c(lh,3);
    v4_h(n)=V1c(lh,4);
end
    
end
