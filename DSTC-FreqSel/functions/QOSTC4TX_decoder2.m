% decoder for Orthogonal STC with 4 antennas
function [v1_h,v2_h,v3_h,v4_h]= QOSTC4TX_decoder2(y_k,y_km1,M,rot_ang)

[~,N]=size(y_k);

for n=1:N      
    v1_h(n)=y_km1(1,n)*y_k(1,n)'+y_km1(2,n)'*y_k(2,n)+y_km1(3,n)'*y_k(3,n)+y_km1(4,n)*y_k(4,n)';
    v4_h(n)=-y_km1(4,n)*y_k(1,n)'+y_km1(3,n)'*y_k(2,n)+y_km1(2,n)'*y_k(3,n)-y_km1(1,n)*y_k(4,n)';
    v2_h(n)=-y_km1(2,n)*y_k(1,n)'+y_km1(1,n)'*y_k(2,n)-y_km1(4,n)'*y_k(3,n)+y_km1(3,n)*y_k(4,n)';
    v3_h(n)=-y_km1(3,n)*y_k(1,n)'-y_km1(4,n)'*y_k(2,n)+y_km1(1,n)'*y_k(3,n)+y_km1(2,n)*y_k(4,n)';
end
    
end
