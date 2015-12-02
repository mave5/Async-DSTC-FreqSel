% decoder for circualnt decoder
function Vh= circulant_decoder(y,cw)

[~,N1]=size(y);
N=N1/2;

% two consecutive vector symbols
y_k=y(:,N+1:end);
y_km1=y(:,1:N);

for k=1:N

% compute |yk-V*y_km1|    
for k1=1:length(cw)
    temp(k1)=norm(y_k(:,k)-cw{k1}*y_km1(:,k));
end
% find closest symbol
[t2 k1_min]=min(temp);
Vh{k}=cw{k1_min};
    
end




