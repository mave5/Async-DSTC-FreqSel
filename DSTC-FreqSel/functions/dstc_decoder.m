
function [v1,v2]= dstc_decoder(y_k,y_km1,type)
    [~,N]=size(y_k);
    if nargin==2
        type=1;
    end
    
    if type==1
    for k=1:N
        v1(k)=(conj(y_km1(1,k))*y_k(1,k)+y_km1(2,k)*conj(y_k(2,k)));
        v2(k)=(-y_km1(2,k)*conj(y_k(1,k))+conj(y_km1(1,k))*y_k(2,k));
    end
    
    else
    for k=1:N
        v1(k)=(conj(y_km1(1,k))*y_k(1,k)+y_km1(2,k)*conj(y_k(2,k)));
        v2(k)=(-y_km1(1,k)*conj(y_k(2,k))+conj(y_km1(2,k))*y_k(1,k));
    end
    end        
end

% the decoding formula for alamouti
% [ y1[k] ] = [ u1   -u2']* [y1[k-1]]+n1[k]
% [ y2[k] ] = [ u2   u1' ]* [y2[k-1]]+n2[k]

% u1= y1[k-1]'*y1[k]+y2[k-1]*y2[k]'
% u2= -y2[k-1]*y1[k]'+y1[k-1]'*y2[k]

